/* The copyright in this software is being made available under the BSD
 * License, included below. This software may be subject to other third party
 * and contributor rights, including patent rights, and no such rights are
 * granted under this license.
 *
 * Copyright (c) 2010-2021, ITU/ISO/IEC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *  * Neither the name of the ITU/ISO/IEC nor the names of its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/** \file     Prediction.cpp
    \brief    prediction class
*/

#include "IntraPrediction.h"

#include "Unit.h"
#include "UnitTools.h"
#include "Buffer.h"

#include "dtrace_next.h"
#include "Rom.h"

#include <memory.h>

#include "CommonLib/InterpolationFilter.h"

//! \ingroup CommonLib
//! \{

// ====================================================================================================================
// Tables
// ====================================================================================================================

const uint8_t IntraPrediction::m_aucIntraFilter[MAX_INTRA_FILTER_DEPTHS] =
{
  24, //   1xn
  24, //   2xn
  24, //   4xn
  14, //   8xn
  2,  //  16xn
  0,  //  32xn
  0,  //  64xn
  0   // 128xn
};

#if JVET_W0123_TIMD_FUSION
const uint8_t IntraPrediction::m_aucIntraFilterExt[MAX_INTRA_FILTER_DEPTHS] =
{
  48, //   1xn
  48, //   2xn
  48, //   4xn
  28, //   8xn
  4,  //  16xn
  0,  //  32xn
  0,  //  64xn
  0   // 128xn
};
#endif

// ====================================================================================================================
// Constructor / destructor / initialize
// ====================================================================================================================

IntraPrediction::IntraPrediction()
:
  m_currChromaFormat( NUM_CHROMA_FORMAT )
{
#if !MERGE_ENC_OPT
  for( uint32_t ch = 0; ch < MAX_NUM_COMPONENT; ch++ )
  {
    for( uint32_t buf = 0; buf < 4; buf++ )
    {
      m_yuvExt2[ch][buf] = nullptr;
    }
  }
#endif

#if JVET_W0123_TIMD_FUSION
  m_timdSatdCost = nullptr;
#endif
  m_piTemp = nullptr;
  m_pMdlmTemp = nullptr;
#if MMLM
  m_encPreRDRun = false;
#endif
}

IntraPrediction::~IntraPrediction()
{
  destroy();
}

void IntraPrediction::destroy()
{
#if !MERGE_ENC_OPT
  for( uint32_t ch = 0; ch < MAX_NUM_COMPONENT; ch++ )
  {
    for( uint32_t buf = 0; buf < 4; buf++ )
    {
      delete[] m_yuvExt2[ch][buf];
      m_yuvExt2[ch][buf] = nullptr;
    }
  }
#endif

#if JVET_W0123_TIMD_FUSION
  delete m_timdSatdCost;
#endif
  delete[] m_piTemp;
  m_piTemp = nullptr;
  delete[] m_pMdlmTemp;
  m_pMdlmTemp = nullptr;

  for( auto &buffer : m_tempBuffer )
  {
    buffer.destroy();
  }
  m_tempBuffer.clear();
}

void IntraPrediction::init(ChromaFormat chromaFormatIDC, const unsigned bitDepthY)
{
#if MERGE_ENC_OPT
  if (m_currChromaFormat != chromaFormatIDC)
  {
    destroy();
  }

  m_currChromaFormat = chromaFormatIDC;
#else
  if( m_yuvExt2[COMPONENT_Y][0] != nullptr && m_currChromaFormat != chromaFormatIDC )
  {
    destroy();
  }

  m_currChromaFormat = chromaFormatIDC;


  if( m_yuvExt2[COMPONENT_Y][0] == nullptr ) // check if first is null (in which case, nothing initialised yet)
  {
    m_yuvExtSize2 = ( MAX_CU_SIZE ) * ( MAX_CU_SIZE );

    for( uint32_t ch = 0; ch < MAX_NUM_COMPONENT; ch++ )
    {
      for( uint32_t buf = 0; buf < 4; buf++ )
      {
        m_yuvExt2[ch][buf] = new Pel[m_yuvExtSize2];
      }
    }
  }
#endif

#if LMS_LINEAR_MODEL
  int shift = bitDepthY + 4;
  for (int i = 32; i < 64; i++)
  {
    m_auShiftLM[i - 32] = ((1 << shift) + i / 2) / i;
  }
#endif

#if JVET_W0123_TIMD_FUSION
  if (m_timdSatdCost == nullptr)
  {
    m_timdSatdCost = new RdCost;
  }
#endif
  if (m_piTemp == nullptr)
  {
    m_piTemp = new Pel[(MAX_CU_SIZE + 1) * (MAX_CU_SIZE + 1)];
  }
  if (m_pMdlmTemp == nullptr)
  {
    m_pMdlmTemp = new Pel[(2 * MAX_CU_SIZE + 1)*(2 * MAX_CU_SIZE + 1)];//MDLM will use top-above and left-below samples.
  }


  for( auto &buffer : m_tempBuffer )
  {
    buffer.destroy();
  }

  // the number of total temporal buffers can be adjusted by chaning the number here
  m_tempBuffer.resize( 2 );

  for( auto &buffer : m_tempBuffer )
  {
    buffer.create( chromaFormatIDC, Area( 0, 0, MAX_CU_SIZE, MAX_CU_SIZE ) );
  }
}

#if JVET_W0123_TIMD_FUSION
void IntraPrediction::xIntraPredTimdAngPdpc(Pel* pDsty,const int dstStride,Pel* refSide,const int width,const int height, int xOffset, int yOffset, int scale,int invAngle)
{
  int xlim = std::min(3 << scale, width);
  for (int y = yOffset; y<height; y++)
  {
    int invAngleSum = 256;
    if (width < 4)
    {
      for (int x = xOffset; x < 2; x++)
      {
        invAngleSum += invAngle;
        int wL   = 32 >> (2 * x >> scale);
        Pel left = refSide[y + (invAngleSum >> 9) + 1];
        pDsty[x] = pDsty[x] + ((wL * (left - pDsty[x]) + 32) >> 6);
      }
    }
    else
    {
      for (int x = xOffset; x < xlim; x++)
      {
        invAngleSum += invAngle;
        int wL   = 32 >> (2 * x >> scale);
        Pel left = refSide[y + (invAngleSum >> 9) + 1];
        pDsty[x] = pDsty[x] + ((wL * (left - pDsty[x]) + 32) >> 6);
      }
    }
    pDsty += dstStride;
  }
}

#if GRAD_PDPC
void IntraPrediction::xIntraPredTimdAngGradPdpc(Pel* pDsty, const int dstStride, Pel* refMain, Pel* refSide, const int width, const int height, int xOffset, int yOffset, int scale, int deltaPos, int intraPredAngle, const ClpRng& clpRng)
{
  for (int y = yOffset; y<height; y++)
  {
    const int deltaInt   = deltaPos >> 6;
    const int deltaFract = deltaPos & 63;
    const Pel left = refSide[1 + y];
    const Pel topLeft = refMain[deltaInt] + ((deltaFract * (refMain[deltaInt + 1] - refMain[deltaInt]) + 32) >> 6);
    for (int x = xOffset; x < std::min(3 << scale, width); x++)
    {
      int wL = 32 >> (2 * (x - xOffset) >> scale);
      pDsty[x] = ClipPel(pDsty[x] + ((wL * (left - topLeft) + 32) >> 6), clpRng);
    }
    pDsty += dstStride;
    deltaPos += intraPredAngle;
  }
}
#endif

void IntraPrediction::xIntraPredTimdHorVerPdpc(Pel* pDsty,const int dstStride,Pel* refSide,const int width,const int height, int xOffset, int yOffset, int scale,const Pel* refMain, const ClpRng& clpRng)
{
  const Pel topLeft = refMain[0];

  for( int y = yOffset; y < height; y++ )
  {
    memcpy(pDsty,&refMain[1],width*sizeof(Pel));
    const Pel left    = refSide[1 + y];
    for (int x = xOffset; x < std::min(3 << scale, width); x++)
    {
      const int wL  = 32 >> (2 * x >> scale);
      const Pel val = pDsty[x];
      pDsty[x]      = ClipPel(val + ((wL * (left - topLeft) + 32) >> 6), clpRng);
    }
    pDsty += dstStride;
  }
}

void IntraPrediction::xIntraPredTimdPlanarDcPdpc(const CPelBuf &pSrc, Pel* pDst, int iDstStride, int width, int height, TEMPLATE_TYPE eTempType, int iTemplateWidth, int iTemplateHeight)
{
  if (eTempType == LEFT_ABOVE_NEIGHBOR)
  {
    int xOffset = 0;
    int yOffset = 0;
    // PDPC for above template
    {
      const int iWidth  = width;
      const int iHeight = iTemplateHeight;
      xOffset = iTemplateWidth;
      const int scale = ((floorLog2(width) - 2 + floorLog2(height) - 2 + 2) >> 2);
      for (int y = 0; y < iHeight; y++)
      {
        const int wT   = 32 >> std::min(31, ((y << 1) >> scale));
        const Pel left = pSrc.at(y + 1, 1);
        for (int x = xOffset; x < iWidth; x++)
        {
          const int wL    = 32 >> std::min(31, ((x << 1) >> scale));
          const Pel top   = pSrc.at(x + 1, 0);
          const Pel val   = pDst[y * iDstStride + x];
          pDst[y * iDstStride + x] = val + ((wL * (left - val) + wT * (top - val) + 32) >> 6);
        }
      }
    }

    // PDPC for left template
    {
      const int iWidth  = iTemplateWidth;
      const int iHeight = height;
      yOffset = iTemplateHeight;
      const int scale = ((floorLog2(width) - 2 + floorLog2(height) - 2 + 2) >> 2);
      for (int y = yOffset; y < iHeight; y++)
      {
        const int wT   = 32 >> std::min(31, ((y << 1) >> scale));
        const Pel left = pSrc.at(y + 1, 1);
        for (int x = 0; x < iWidth; x++)
        {
          const int wL    = 32 >> std::min(31, ((x << 1) >> scale));
          const Pel top   = pSrc.at(x + 1, 0);
          const Pel val   = pDst[y * iDstStride + x];
          pDst[y * iDstStride + x] = val + ((wL * (left - val) + wT * (top - val) + 32) >> 6);
        }
      }
    }
  }
  else if (eTempType == LEFT_NEIGHBOR)
  {
    const int iHeight = height;
    const int scale = ((floorLog2(width) - 2 + floorLog2(height) - 2 + 2) >> 2);
    for (int y = 0; y < iHeight; y++)
    {
      const int wT   = 32 >> std::min(31, ((y << 1) >> scale));
      const Pel left = pSrc.at(y + 1, 1);
      for (int x = 0; x < iTemplateWidth; x++)
      {
        const int wL    = 32 >> std::min(31, ((x << 1) >> scale));
        const Pel top   = pSrc.at(x + 1, 0);
        const Pel val   = pDst[y * iDstStride + x];
        pDst[y * iDstStride + x] = val + ((wL * (left - val) + wT * (top - val) + 32) >> 6);
      }
    }
  }
  else // eTempType == ABOVE_NEIGHBOR
  {
    const int iWidth  = width;
    const int scale = ((floorLog2(width) - 2 + floorLog2(height) - 2 + 2) >> 2);
    for (int y = 0; y < iTemplateHeight; y++)
    {
      const int wT   = 32 >> std::min(31, ((y << 1) >> scale));
      const Pel left = pSrc.at(y + 1, 1);
      for (int x = 0; x < iWidth; x++)
      {
        const int wL    = 32 >> std::min(31, ((x << 1) >> scale));
        const Pel top   = pSrc.at(x + 1, 0);
        const Pel val   = pDst[y * iDstStride + x];
        pDst[y * iDstStride + x] = val + ((wL * (left - val) + wT * (top - val) + 32) >> 6);
      }
    }
  }
}

void IntraPrediction::xIntraPredTimdAngLuma(Pel* pDstBuf, const ptrdiff_t dstStride, Pel* refMain, int width, int height, int deltaPos, int intraPredAngle, const ClpRng& clpRng, int xOffset, int yOffset)
{
  for (int y = yOffset; y<height; y++ )
  {
    const int deltaInt   = deltaPos >> 6;
    const int deltaFract = deltaPos & 63;
    const TFilterCoeff* const f = InterpolationFilter::getExtIntraCubicFilter(deltaFract);
    int refMainIndex = deltaInt + 1 + xOffset;
    for( int x = xOffset; x < width; x++, refMainIndex++ )
    {
      pDstBuf[y*dstStride + x] = (f[0] * refMain[refMainIndex - 1] + f[1] * refMain[refMainIndex] + f[2] * refMain[refMainIndex + 1] + f[3] * refMain[refMainIndex + 2] + 128) >> 8;
      pDstBuf[y*dstStride + x] = ClipPel( pDstBuf[y*dstStride + x], clpRng ); // always clip even though not always needed
    }
    deltaPos += intraPredAngle;
  }
}
#endif

// ====================================================================================================================
// Public member functions
// ====================================================================================================================

// Function for calculating DC value of the reference samples used in Intra prediction
//NOTE: Bit-Limit - 25-bit source
Pel IntraPrediction::xGetPredValDc( const CPelBuf &pSrc, const Size &dstSize )
{
  CHECK( dstSize.width == 0 || dstSize.height == 0, "Empty area provided" );

  int idx, sum = 0;
  Pel dcVal;
  const int width  = dstSize.width;
  const int height = dstSize.height;
  const auto denom     = (width == height) ? (width << 1) : std::max(width,height);
  const auto divShift  = floorLog2(denom);
  const auto divOffset = (denom >> 1);
  const Pel* src;

  if ( width >= height )
  {
    src = pSrc.bufAt(m_ipaParam.multiRefIndex + 1, 0);
    for (idx = 0; idx < width; idx++)
    {
      sum += src[idx];
    }
  }
  if ( width <= height )
  {
    src = pSrc.bufAt(m_ipaParam.multiRefIndex + 1, 1);
    for (idx = 0; idx < height; idx++)
    {
      sum += src[idx];
    }
  }

  dcVal = (sum + divOffset) >> divShift;
  return dcVal;
}

int IntraPrediction::getModifiedWideAngle( int width, int height, int predMode )
{
  //The function returns a 'modified' wide angle index, given that it is not necessary 
  //in this software implementation to reserve the values 0 and 1 for Planar and DC to generate the prediction signal.
  //It should only be used to obtain the intraPredAngle parameter.
  //To simply obtain the wide angle index, the function PU::getWideAngle should be used instead.
  if ( predMode > DC_IDX && predMode <= VDIA_IDX )
  {
    int modeShift[] = { 0, 6, 10, 12, 14, 15 };
    int deltaSize = abs(floorLog2(width) - floorLog2(height));
    if (width > height && predMode < 2 + modeShift[deltaSize])
    {
      predMode += (VDIA_IDX - 1);
    }
    else if (height > width && predMode > VDIA_IDX - modeShift[deltaSize])
    {
      predMode -= (VDIA_IDX - 1);
    }
  }
  return predMode;
}

#if JVET_W0123_TIMD_FUSION
int IntraPrediction::getWideAngleExt( int width, int height, int predMode )
{
  if ( predMode > DC_IDX && predMode <= EXT_VDIA_IDX )
  {
    int modeShift[] = { 0, 11, 19, 23, 27, 29 };
    int deltaSize = abs(floorLog2(width) - floorLog2(height));
    if (width > height && predMode < 2 + modeShift[deltaSize])
    {
      predMode += (EXT_VDIA_IDX - 1);
    }
    else if (height > width && predMode > EXT_VDIA_IDX - modeShift[deltaSize])
    {
      predMode -= (EXT_VDIA_IDX - 1);
    }
  }
  return predMode;
}
#endif

void IntraPrediction::setReferenceArrayLengths( const CompArea &area )
{
  // set Top and Left reference samples length
  const int  width    = area.width;
  const int  height   = area.height;

  m_leftRefLength     = (height << 1);
  m_topRefLength      = (width << 1);
}

void IntraPrediction::predIntraAng( const ComponentID compId, PelBuf &piPred, const PredictionUnit &pu)
{
#if CIIP_PDPC
  CHECK((pu.ciipFlag == false && pu.ciipPDPC == true),"ciip_PDPC can not be true for an non CIIP PU");
#endif
  const ComponentID    compID       = MAP_CHROMA( compId );
  const ChannelType    channelType  = toChannelType( compID );
  const int            iWidth       = piPred.width;
  const int            iHeight      = piPred.height;
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  CHECK(iWidth == 2, "Width of 2 is not supported");
#endif

  CHECK(PU::isMIP(pu, toChannelType(compId)), "We should not get here for MIP.");
  const uint32_t       uiDirMode    = isLuma( compId ) && pu.cu->bdpcmMode ? BDPCM_IDX : !isLuma(compId) && pu.cu->bdpcmModeChroma ? BDPCM_IDX : PU::getFinalIntraMode(pu, channelType);
#if JVET_W0123_TIMD_FUSION
  bool bExtIntraDir = pu.cu->timd && isLuma( compId );
#endif

  CHECK( floorLog2(iWidth) < 2 && pu.cs->pcv->noChroma2x2, "Size not allowed" );
  CHECK( floorLog2(iWidth) > 7, "Size not allowed" );

  const int srcStride  = m_refBufferStride[compID];
  const int srcHStride = 2;

  const CPelBuf & srcBuf = CPelBuf(getPredictorPtr(compID), srcStride, srcHStride);
  const ClpRng& clpRng(pu.cu->cs->slice->clpRng(compID));
#if CIIP_PDPC
  if (!pu.ciipPDPC)
  {
#endif
    switch (uiDirMode)
    {
    case(PLANAR_IDX): xPredIntraPlanar(srcBuf, piPred); break;
    case(DC_IDX):     xPredIntraDc(srcBuf, piPred, channelType, false); break;
    case(BDPCM_IDX):  xPredIntraBDPCM(srcBuf, piPred, isLuma(compID) ? pu.cu->bdpcmMode : pu.cu->bdpcmModeChroma, clpRng); break;
#if JVET_W0123_TIMD_FUSION
    default:          xPredIntraAng(srcBuf, piPred, channelType, clpRng, bExtIntraDir); break;
#else
    default:          xPredIntraAng(srcBuf, piPred, channelType, clpRng); break;
#endif
    }
#if CIIP_PDPC
  }
#endif

#if ENABLE_DIMD
  if (pu.cu->dimd && pu.cu->dimd_is_blend && isLuma(compID))
  {
    int width = piPred.width;
    int height = piPred.height;
    const UnitArea localUnitArea( pu.chromaFormat, Area( 0, 0, width, height ) );

    PelBuf planarBuffer = m_tempBuffer[0].getBuf( localUnitArea.Y() );
    PelBuf predAng = m_tempBuffer[1].getBuf( localUnitArea.Y() );

    xPredIntraPlanar( srcBuf, planarBuffer );

    const bool applyPdpc = m_ipaParam.applyPDPC;
#if JVET_V0087_DIMD_NO_ISP   // this is pure cleanup to make code easier to read. It generates identical resut to the else part
    PredictionUnit pu2 = pu;
    pu2.intraDir[0] = pu.cu->dimdBlendMode[0];
    initPredIntraParams(pu2, pu.Y(), *(pu.cs->sps));

#if JVET_W0123_TIMD_FUSION
    xPredIntraAng(srcBuf, predAng, channelType, clpRng, false);
#else
    xPredIntraAng(srcBuf, predAng, channelType, clpRng);
#endif
#else
    const bool   useISP = NOT_INTRA_SUBPARTITIONS != pu.cu->ispMode && isLuma( CHANNEL_TYPE_LUMA );//ok
    const Size   cuSize = Size( pu.cu->blocks[compId].width, pu.cu->blocks[compId].height ); //ok
    const Size   puSize = Size( piPred.width, piPred.height );
    const Size&  blockSize = useISP ? cuSize : puSize;
    int blendDir = pu.cu->dimdBlendMode[0];
    const int      dirMode = blendDir;
    const int     predMode = getModifiedWideAngle( blockSize.width, blockSize.height, dirMode ); // to check later
    m_ipaParam.isModeVer = predMode >= DIA_IDX;
    m_ipaParam.multiRefIndex = 0;
    m_ipaParam.refFilterFlag = false;
    m_ipaParam.interpolationFlag = false;
    m_ipaParam.applyPDPC = ( ( puSize.width >= MIN_TB_SIZEY && puSize.height >= MIN_TB_SIZEY ) || !isLuma( compId ) ) && m_ipaParam.multiRefIndex == 0;

    const int    intraPredAngleMode = ( m_ipaParam.isModeVer ) ? predMode - VER_IDX : -( predMode - HOR_IDX );//ok
    int absAng = 0;
    if( dirMode > DC_IDX && dirMode < NUM_LUMA_MODE ) // intraPredAngle for directional modes
    {
      static const int angTable[32] = { 0,    1,    2,    3,    4,    6,     8,   10,   12,   14,   16,   18,   20,   23,   26,   29,   32,   35,   39,  45,  51,  57,  64,  73,  86, 102, 128, 171, 256, 341, 512, 1024 };//ok
      static const int invAngTable[32] = { 0,   16384, 8192, 5461, 4096, 2731, 2048, 1638, 1365, 1170, 1024, 910, 819, 712, 630, 565, 512, 468,   420,  364,  321,  287,  256,  224,  191,  161,  128,  96,  64,  48,  32,  16 };   // (512 * 32) / Angle

      const int     absAngMode = abs( intraPredAngleMode );
      const int     signAng = intraPredAngleMode < 0 ? -1 : 1;
      absAng = angTable[absAngMode];

      m_ipaParam.absInvAngle = invAngTable[absAngMode];
      m_ipaParam.intraPredAngle = signAng * absAng;
      if( intraPredAngleMode < 0 )
      {
        m_ipaParam.applyPDPC = false;
      }
      else if( intraPredAngleMode > 0 )
      {
        const int sideSize = m_ipaParam.isModeVer ? puSize.height : puSize.width;
        const int maxScale = 2;
#if GRAD_PDPC
        m_ipaParam.useGradPDPC = false;
#endif

        m_ipaParam.angularScale = std::min( maxScale, floorLog2( sideSize ) - ( floorLog2( 3 * m_ipaParam.absInvAngle - 2 ) - 8 ) );
#if GRAD_PDPC
        if( ( m_ipaParam.angularScale < 0 ) && ( isLuma( compId ) ) )
        {
          m_ipaParam.angularScale = ( floorLog2( puSize.width ) + floorLog2( puSize.height ) - 2 ) >> 2;
          m_ipaParam.useGradPDPC = true;
        }
#endif
        m_ipaParam.applyPDPC &= m_ipaParam.angularScale >= 0;
      }
    }

    if( pu.cs->sps->getSpsRangeExtension().getIntraSmoothingDisabledFlag()
        || ( !isLuma( CHANNEL_TYPE_LUMA ) && pu.chromaFormat != CHROMA_444 )
        || useISP
        || m_ipaParam.multiRefIndex
        || DC_IDX == dirMode
        )
    {
      //do nothing
    }
    else if( !useISP )// HOR, VER and angular modes (MDIS)
    {
      bool filterFlag = false;
      const int diff = std::min<int>( abs( predMode - HOR_IDX ), abs( predMode - VER_IDX ) );
      const int log2Size = ( ( floorLog2( puSize.width ) + floorLog2( puSize.height ) ) >> 1 );
      CHECK( log2Size >= MAX_INTRA_FILTER_DEPTHS, "Size not supported" );
      filterFlag = ( diff > m_aucIntraFilter[log2Size] );


      if( filterFlag )
      {
        const bool isRefFilter = isIntegerSlope( absAng );
        CHECK( puSize.width * puSize.height <= 32, "DCT-IF interpolation filter is always used for 4x4, 4x8, and 8x4 luma CB" );
        m_ipaParam.refFilterFlag = isRefFilter;
        m_ipaParam.interpolationFlag = !isRefFilter;
      }
    }

#if JVET_W0123_TIMD_FUSION
    xPredIntraAng( srcBuf, predAng, channelType, clpRng, false );
#else
    xPredIntraAng( srcBuf, predAng, channelType, clpRng );
#endif
#endif
    m_ipaParam.applyPDPC = applyPdpc;

    // do blending
    const int log2WeightSum = 6;
    Pel *pelPred = piPred.buf;
    Pel *pelPlanar = planarBuffer.buf;
    Pel *pelPredAng = predAng.buf;
    int  w0 = pu.cu->dimdRelWeight[0], w1 = pu.cu->dimdRelWeight[1], w2 = pu.cu->dimdRelWeight[2];

    for( int y = 0; y < height; y++ )
    {
      for( int x = 0; x < width; x++ )
      {
        int blend = pelPred[x] * w0;
        blend += pelPlanar[x] * w1;
        blend += pelPredAng[x] * w2;
        pelPred[x] = (Pel)(blend >> log2WeightSum);
      }

      pelPred += piPred.stride;
      pelPlanar += planarBuffer.stride;
      pelPredAng += predAng.stride;
    }
  }
#endif

#if JVET_W0123_TIMD_FUSION
  if (pu.cu->timd && pu.cu->timdIsBlended && isLuma(compID))
  {
    int width = piPred.width;
    int height = piPred.height;
    const UnitArea localUnitArea( pu.chromaFormat, Area( 0, 0, width, height ) );

    PelBuf predFusion = m_tempBuffer[1].getBuf( localUnitArea.Y() );

    const bool applyPdpc = m_ipaParam.applyPDPC;
    PredictionUnit pu2 = pu;
    pu2.intraDir[0] = pu.cu->timdModeSecondary;
    initPredIntraParams(pu2, pu.Y(), *(pu.cs->sps));

    switch (pu.cu->timdModeSecondary)
    {
    case(PLANAR_IDX): xPredIntraPlanar(srcBuf, predFusion); break;
    case(DC_IDX):     xPredIntraDc(srcBuf, predFusion, channelType, false); break;
    default:          xPredIntraAng(srcBuf, predFusion, channelType, clpRng, bExtIntraDir); break;
    }

    m_ipaParam.applyPDPC = applyPdpc;

    // do blending
    const int log2WeightSum = 6;
    Pel *pelPred = piPred.buf;
    Pel *pelPredFusion = predFusion.buf;
    int  w0 = pu.cu->timdFusionWeight[0], w1 = pu.cu->timdFusionWeight[1];

    for( int y = 0; y < height; y++ )
    {
      for( int x = 0; x < width; x++ )
      {
        int blend = pelPred[x] * w0;
        blend += pelPredFusion[x] * w1;
        pelPred[x] = (Pel)(blend >> log2WeightSum);
      }

      pelPred += piPred.stride;
      pelPredFusion += predFusion.stride;
    }
  }
#endif

#if CIIP_PDPC
  if (m_ipaParam.applyPDPC || pu.ciipPDPC)
#else
  if (m_ipaParam.applyPDPC)
#endif
  {
    PelBuf dstBuf = piPred;
    const int scale = ((floorLog2(iWidth) - 2 + floorLog2(iHeight) - 2 + 2) >> 2);
    CHECK(scale < 0 || scale > 31, "PDPC: scale < 0 || scale > 31");

    if (uiDirMode == PLANAR_IDX || uiDirMode == DC_IDX)
    { 
      const Pel* srcLeft = srcBuf.bufAt(1, 1);
      Pel* dst = dstBuf.buf;
#if CIIP_PDPC
      if (pu.ciipPDPC)
      {
        for (int y = 0; y < iHeight; y++)
        {
          const int wT = 32 >> std::min(31, ((y << 1) >> scale));
          const Pel left = *srcLeft;
          const Pel* srcTop = srcBuf.buf + 1;
          for (int x = 0; x < iWidth; x++)
          {
            const int wL = 32 >> std::min(31, ((x << 1) >> scale));
            const Pel top = *srcTop;
            dst[x] = ((wL * left  + wT * top + 32) >> 6);

            srcTop++;
          }
          srcLeft++;
          dst += dstBuf.stride;
        }
      }
      else
#endif

      for (int y = 0; y < iHeight; y++)
      {
        const int wT   = 32 >> std::min(31, ((y << 1) >> scale));
        const Pel left = *srcLeft;
        const Pel* srcTop = srcBuf.buf + 1;

        for (int x = 0; x < iWidth; x++)
        {
          const int wL = 32 >> std::min(31, ((x << 1) >> scale));
          const Pel top = *srcTop;
          const Pel val = dst[x];

          dst[x] = val + ((wL * (left - val) + wT * (top - val) + 32) >> 6);
          srcTop++;
        }

        srcLeft++;
        dst += dstBuf.stride;
      }
    }
  }
}

void IntraPrediction::predIntraChromaLM(const ComponentID compID, PelBuf &piPred, const PredictionUnit &pu, const CompArea& chromaArea, int intraDir)
{
  int  iLumaStride = 0;
  PelBuf Temp;
#if MMLM
  if ((intraDir == MDLM_L_IDX) || (intraDir == MDLM_T_IDX) || (intraDir == MMLM_L_IDX) || (intraDir == MMLM_T_IDX) || (m_encPreRDRun && intraDir == MMLM_CHROMA_IDX))
#else
  if ((intraDir == MDLM_L_IDX) || (intraDir == MDLM_T_IDX))
#endif
  {
    iLumaStride = 2 * MAX_CU_SIZE + 1;
    Temp = PelBuf(m_pMdlmTemp + iLumaStride + 1, iLumaStride, Size(chromaArea));
  }
  else
  {
    iLumaStride = MAX_CU_SIZE + 1;
    Temp = PelBuf(m_piTemp + iLumaStride + 1, iLumaStride, Size(chromaArea));
  }
  int a, b, iShift;
#if MMLM
  int a2, b2, iShift2, yThres;
#if LMS_LINEAR_MODEL
  xGetLMParameters_LMS(pu, compID, chromaArea, a, b, iShift, a2, b2, iShift2, yThres);
#else
  xGetLMParameters(pu, compID, chromaArea, a, b, iShift, a2, b2, iShift2, yThres);
#endif
#else
#if LMS_LINEAR_MODEL
  xGetLMParameters_LMS(pu, compID, chromaArea, a, b, iShift);
#else
  xGetLMParameters(pu, compID, chromaArea, a, b, iShift);
#endif
#endif

  ////// final prediction
  piPred.copyFrom(Temp);
#if MMLM
  if (PU::isMultiModeLM(pu.intraDir[1]))
  {
    Pel*  pPred = piPred.bufAt(0, 0);
    Pel  *pLuma = Temp.bufAt(0, 0);
    int uiPredStride = piPred.stride;
    int uiCWidth = chromaArea.width;
    int uiCHeight = chromaArea.height;

    for (int i = 0; i < uiCHeight; i++)
    {
      for (int j = 0; j < uiCWidth; j++)
      {
        if (pLuma[j] <= yThres)
        {
          pPred[j] = (Pel)ClipPel(((a * pLuma[j]) >> iShift) + b, pu.cs->slice->clpRng(compID));
        }
        else
        {
          pPred[j] = (Pel)ClipPel(((a2 * pLuma[j]) >> iShift2) + b2, pu.cs->slice->clpRng(compID));
        }
      }
      pPred += uiPredStride;
      pLuma += iLumaStride;
    }
  }
  else
#endif
  piPred.linearTransform(a, iShift, b, true, pu.cs->slice->clpRng(compID));
}

/** Function for deriving planar intra prediction. This function derives the prediction samples for planar mode (intra coding).
 */

//NOTE: Bit-Limit - 24-bit source
void IntraPrediction::xPredIntraPlanar( const CPelBuf &pSrc, PelBuf &pDst )
{
  const uint32_t width  = pDst.width;
  const uint32_t height = pDst.height;

  const uint32_t log2W = floorLog2( width );
  const uint32_t log2H = floorLog2( height );

  int leftColumn[MAX_CU_SIZE + 1], topRow[MAX_CU_SIZE + 1], bottomRow[MAX_CU_SIZE], rightColumn[MAX_CU_SIZE];
  const uint32_t offset = 1 << (log2W + log2H);

  // Get left and above reference column and row
  const Pel* src = pSrc.bufAt(1, 0);

  for (int k = 0; k < width + 1; k++)
  {
    topRow[k] = src[k];
  }

  src = pSrc.bufAt(1, 1);

  for (int k = 0; k < height + 1; k++)
  {
    leftColumn[k] = src[k];
  }

  // Prepare intermediate variables used in interpolation
  int bottomLeft = leftColumn[height];
  int topRight = topRow[width];

  for( int k = 0; k < width; k++ )
  {
    bottomRow[k] = bottomLeft - topRow[k];
    topRow[k]    = topRow[k] << log2H;
  }

  for( int k = 0; k < height; k++ )
  {
    rightColumn[k] = topRight - leftColumn[k];
    leftColumn[k]  = leftColumn[k] << log2W;
  }

  const uint32_t finalShift = 1 + log2W + log2H;
  const uint32_t stride     = pDst.stride;
  Pel*       pred       = pDst.buf;
  for( int y = 0; y < height; y++, pred += stride )
  {
    int horPred = leftColumn[y];

    for( int x = 0; x < width; x++ )
    {
      horPred += rightColumn[y];
      topRow[x] += bottomRow[x];

      int vertPred = topRow[x];
      pred[x]      = ( ( horPred << log2H ) + ( vertPred << log2W ) + offset ) >> finalShift;
    }
  }
}

void IntraPrediction::xPredIntraDc( const CPelBuf &pSrc, PelBuf &pDst, const ChannelType channelType, const bool enableBoundaryFilter )
{
  const Pel dcval = xGetPredValDc( pSrc, pDst );
  pDst.fill( dcval );
}

// Function for initialization of intra prediction parameters
void IntraPrediction::initPredIntraParams(const PredictionUnit & pu, const CompArea area, const SPS& sps)
{
  const ComponentID compId = area.compID;
  const ChannelType chType = toChannelType(compId);
#if JVET_W0123_TIMD_FUSION
  bool bExtIntraDir = pu.cu->timd && isLuma( chType );
#endif

  const bool        useISP = NOT_INTRA_SUBPARTITIONS != pu.cu->ispMode && isLuma( chType );

  const Size   cuSize    = Size( pu.cu->blocks[compId].width, pu.cu->blocks[compId].height );
  const Size   puSize    = Size( area.width, area.height );
  const Size&  blockSize = useISP ? cuSize : puSize;
  const int      dirMode = PU::getFinalIntraMode(pu, chType);
#if JVET_W0123_TIMD_FUSION
  const int     predMode = bExtIntraDir ? getWideAngleExt( blockSize.width, blockSize.height, dirMode ) : getModifiedWideAngle( blockSize.width, blockSize.height, dirMode );
#else
  const int     predMode = getModifiedWideAngle( blockSize.width, blockSize.height, dirMode );
#endif

#if JVET_W0123_TIMD_FUSION
  m_ipaParam.isModeVer            = bExtIntraDir ? (predMode >= EXT_DIA_IDX) : (predMode >= DIA_IDX);
#else
  m_ipaParam.isModeVer            = predMode >= DIA_IDX;
#endif
  m_ipaParam.multiRefIndex        = isLuma (chType) ? pu.multiRefIdx : 0 ;
  m_ipaParam.refFilterFlag        = false;
  m_ipaParam.interpolationFlag    = false;
  m_ipaParam.applyPDPC            = (puSize.width >= MIN_TB_SIZEY && puSize.height >= MIN_TB_SIZEY) && m_ipaParam.multiRefIndex == 0;

#if JVET_W0123_TIMD_FUSION
  const int    intraPredAngleMode = (m_ipaParam.isModeVer) ? (predMode - (bExtIntraDir? EXT_VER_IDX : VER_IDX)) : (-(predMode - (bExtIntraDir ? EXT_HOR_IDX : HOR_IDX)));
#else
  const int    intraPredAngleMode = (m_ipaParam.isModeVer) ? predMode - VER_IDX : -(predMode - HOR_IDX);
#endif


  int absAng = 0;
#if JVET_W0123_TIMD_FUSION
  if (dirMode > DC_IDX && dirMode < (bExtIntraDir ? EXT_VDIA_IDX + 1 : NUM_LUMA_MODE)) // intraPredAngle for directional modes
#else
  if (dirMode > DC_IDX && dirMode < NUM_LUMA_MODE) // intraPredAngle for directional modes
#endif
  {
    static const int angTable[32]    = { 0,    1,    2,    3,    4,    6,     8,   10,   12,   14,   16,   18,   20,   23,   26,   29,   32,   35,   39,  45,  51,  57,  64,  73,  86, 102, 128, 171, 256, 341, 512, 1024 };
    static const int invAngTable[32] = {
      0,   16384, 8192, 5461, 4096, 2731, 2048, 1638, 1365, 1170, 1024, 910, 819, 712, 630, 565,
      512, 468,   420,  364,  321,  287,  256,  224,  191,  161,  128,  96,  64,  48,  32,  16
    };   // (512 * 32) / Angle
#if JVET_W0123_TIMD_FUSION
    static const int extAngTable[64]    = { 0, 1, 2, 3, 4, 5, 6,7, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 43, 46, 49, 52, 55, 58, 61, 64, 67, 70, 74, 78, 84, 90, 96, 102, 108, 114, 121, 128, 137, 146, 159, 172, 188, 204, 230, 256, 299, 342, 427, 512, 597, 682, 853, 1024, 1536, 2048, 3072 };
    static const int extInvAngTable[64] = {
        0, 32768, 16384, 10923, 8192, 6554, 5461, 4681, 4096, 3277, 2731, 2341, 2048, 1820, 1638, 1489, 1365, 1260, 1170, 1092, 1024, 964, 910, 862, 819, 762, 712, 669, 630, 596, 565, 537, 512, 489, 468, 443, 420, 390, 364, 341, 321, 303, 287, 271, 256, 239, 224, 206, 191, 174, 161, 142, 128, 110, 96, 77, 64, 55, 48, 38, 32, 21, 16, 11
    };   // (512 * 64) / Angle
#endif

    const int     absAngMode         = abs(intraPredAngleMode);
    const int     signAng            = intraPredAngleMode < 0 ? -1 : 1;
#if JVET_W0123_TIMD_FUSION
                  absAng             = bExtIntraDir ? extAngTable[absAngMode] : angTable[absAngMode];
    m_ipaParam.absInvAngle              = bExtIntraDir ? extInvAngTable[absAngMode] : invAngTable[absAngMode];
#else
                  absAng             = angTable  [absAngMode];

    m_ipaParam.absInvAngle           = invAngTable[absAngMode];
#endif
    m_ipaParam.intraPredAngle        = signAng * absAng;
    if (intraPredAngleMode < 0)
    {
      m_ipaParam.applyPDPC = false;
    }
    else if (intraPredAngleMode > 0)
    {
      const int sideSize = m_ipaParam.isModeVer ? puSize.height : puSize.width;
      const int maxScale = 2;
#if GRAD_PDPC
      m_ipaParam.useGradPDPC = false;
#endif

      m_ipaParam.angularScale = std::min(maxScale, floorLog2(sideSize) - (floorLog2(3 * m_ipaParam.absInvAngle - 2) - 8));
#if GRAD_PDPC
      if ((m_ipaParam.angularScale < 0) && (isLuma(compId)))
      {
        m_ipaParam.angularScale = (floorLog2(puSize.width) + floorLog2(puSize.height) - 2) >> 2;
        m_ipaParam.useGradPDPC = true;
      }
#endif
      m_ipaParam.applyPDPC &= m_ipaParam.angularScale >= 0;
    }
  }

  // high level conditions and DC intra prediction
  if(   sps.getSpsRangeExtension().getIntraSmoothingDisabledFlag()
    || !isLuma( chType )
    || useISP
#if JVET_V0130_INTRA_TMP
	  || PU::isTmp(pu, chType)
#endif
    || PU::isMIP( pu, chType )
    || m_ipaParam.multiRefIndex
    || DC_IDX == dirMode
    )
  {
  }
  else if ((isLuma(chType) && pu.cu->bdpcmMode) || (!isLuma(chType) && pu.cu->bdpcmModeChroma)) // BDPCM
  {
    m_ipaParam.refFilterFlag = false;
  }
  else if (dirMode == PLANAR_IDX) // Planar intra prediction
  {
    m_ipaParam.refFilterFlag = puSize.width * puSize.height > 32 ? true : false;
  }
  else if (!useISP)// HOR, VER and angular modes (MDIS)
  {
    bool filterFlag = false;
    {
#if JVET_W0123_TIMD_FUSION
      const int diff = std::min<int>( abs( predMode - (bExtIntraDir ? EXT_HOR_IDX : HOR_IDX) ), abs( predMode - (bExtIntraDir ? EXT_VER_IDX : VER_IDX) ) );
#else
      const int diff = std::min<int>( abs( predMode - HOR_IDX ), abs( predMode - VER_IDX ) );
#endif
      const int log2Size = ((floorLog2(puSize.width) + floorLog2(puSize.height)) >> 1);
      CHECK( log2Size >= MAX_INTRA_FILTER_DEPTHS, "Size not supported" );
#if JVET_W0123_TIMD_FUSION
      filterFlag = (diff > (bExtIntraDir ? m_aucIntraFilterExt[log2Size] : m_aucIntraFilter[log2Size]));
#else
      filterFlag = (diff > m_aucIntraFilter[log2Size]);
#endif
    }

    // Selelection of either ([1 2 1] / 4 ) refrence filter OR Gaussian 4-tap interpolation filter
    if (filterFlag)
    {
#if JVET_W0123_TIMD_FUSION
      const bool isRefFilter       =  bExtIntraDir ? isIntegerSlopeExt(absAng) : isIntegerSlope(absAng);
#else
      const bool isRefFilter       =  isIntegerSlope(absAng);
#endif
      CHECK( puSize.width * puSize.height <= 32, "DCT-IF interpolation filter is always used for 4x4, 4x8, and 8x4 luma CB" );
      m_ipaParam.refFilterFlag     =  isRefFilter;
      m_ipaParam.interpolationFlag = !isRefFilter;
    }
  }
}


/** Function for deriving the simplified angular intra predictions.
*
* This function derives the prediction samples for the angular mode based on the prediction direction indicated by
* the prediction mode index. The prediction direction is given by the displacement of the bottom row of the block and
* the reference row above the block in the case of vertical prediction or displacement of the rightmost column
* of the block and reference column left from the block in the case of the horizontal prediction. The displacement
* is signalled at 1/32 pixel accuracy. When projection of the predicted pixel falls inbetween reference samples,
* the predicted value for the pixel is linearly interpolated from the reference samples. All reference samples are taken
* from the extended main reference.
*/
//NOTE: Bit-Limit - 25-bit source

#if JVET_W0123_TIMD_FUSION
void IntraPrediction::xPredIntraAng( const CPelBuf &pSrc, PelBuf &pDst, const ChannelType channelType, const ClpRng& clpRng, const bool bExtIntraDir)
#else
void IntraPrediction::xPredIntraAng( const CPelBuf &pSrc, PelBuf &pDst, const ChannelType channelType, const ClpRng& clpRng)
#endif
{
  int width =int(pDst.width);
  int height=int(pDst.height);

  const bool bIsModeVer     = m_ipaParam.isModeVer;
  const int  multiRefIdx    = m_ipaParam.multiRefIndex;
  const int  intraPredAngle = m_ipaParam.intraPredAngle;
  const int  absInvAngle    = m_ipaParam.absInvAngle;

  Pel* refMain;
  Pel* refSide;

#if !INTRA_6TAP
  Pel  refAbove[2 * MAX_CU_SIZE + 3 + 33 * MAX_REF_LINE_IDX];
  Pel  refLeft [2 * MAX_CU_SIZE + 3 + 33 * MAX_REF_LINE_IDX];
#else
  // 2 pixels more for 6 tap filter.
  Pel  refAbove[2 * MAX_CU_SIZE + 5 + 33 * MAX_REF_LINE_IDX];
  Pel  refLeft[2 * MAX_CU_SIZE + 5 + 33 * MAX_REF_LINE_IDX];
  // initializing for safeguard.
  ::memset(refAbove, 0, sizeof(refAbove));
  ::memset(refLeft, 0, sizeof(refAbove));
#endif

  // Initialize the Main and Left reference array.
  if (intraPredAngle < 0)
  {
#if INTRA_6TAP
    // x, y range increase by 1 (right extend)
    const Pel *src = pSrc.buf;
    Pel *dst = refAbove + height + 1;
    ::memcpy(dst, src, sizeof(Pel) * (width + 2 + multiRefIdx + 1));

    src = pSrc.buf + pSrc.stride;
    dst = refLeft + width + 1;
    ::memcpy(dst, src, sizeof(Pel) * (height + 2 + multiRefIdx + 1));

    refMain = bIsModeVer ? refAbove + height + 1 : refLeft + width + 1;
    refSide = bIsModeVer ? refLeft + width + 1: refAbove + height + 1;

    // Extend the Main reference to the left.
    int sizeSide = bIsModeVer ? height : width;
    // left extend by 1
    for (int k = -(sizeSide + 1); k <= -1; k++)
    {
#if JVET_W0123_TIMD_FUSION
      int frac32precision = bExtIntraDir ? ((-k * absInvAngle + 16) >> 5) : ((-k * absInvAngle + 8) >> 4);
#else
      int frac32precision = (-k * absInvAngle + 8) >> 4;
#endif
      int intpel = frac32precision >> 5;
      int fracpel = frac32precision & 31;
      //std::cout << " fracPel: " << fracpel << std::endl;
      int left_minus1 = refSide[Clip3(0, sizeSide + 2 + multiRefIdx, intpel - 1)];
      int left        = refSide[Clip3(0, sizeSide + 2 + multiRefIdx, intpel)];
      int right       = refSide[Clip3(0, sizeSide + 2 + multiRefIdx, intpel + 1)];
      int right_plus1 = refSide[Clip3(0, sizeSide + 2 + multiRefIdx, intpel + 2)];

      const TFilterCoeff* f = InterpolationFilter::getWeak4TapFilterTable(fracpel);
      int val = ((int)f[0] * left_minus1 + (int)f[1] * left + (int)f[2] * right + f[3] * (int)right_plus1 + 32) >> 6;
      refMain[k] = (Pel)ClipPel(val, clpRng);
    }
#else
    for (int x = 0; x <= width + 1 + multiRefIdx; x++)
    {
      refAbove[x + height] = pSrc.at(x, 0);
    }
    for (int y = 0; y <= height + 1 + multiRefIdx; y++)
    {
      refLeft[y + width] = pSrc.at(y, 1);
    }
    refMain = bIsModeVer ? refAbove + height : refLeft + width;
    refSide = bIsModeVer ? refLeft + width : refAbove + height;

    // Extend the Main reference to the left.
    int sizeSide = bIsModeVer ? height : width;
    for (int k = -sizeSide; k <= -1; k++)
    {
      refMain[k] = refSide[std::min((-k * absInvAngle + 256) >> 9, sizeSide)];
    }
#endif
  }
  else
  {
#if INTRA_6TAP
    const Pel *src = pSrc.buf;
    Pel *dst = refAbove + 1;
    ::memcpy(dst, src, sizeof(Pel) * (m_topRefLength + multiRefIdx + 1));

    src = pSrc.buf + pSrc.stride;
    dst = refLeft + 1;
    ::memcpy(dst, src, sizeof(Pel) * (m_leftRefLength + multiRefIdx + 1));

    // left extended by 1
    refAbove[0] = refAbove[1];
    refLeft[0] = refLeft[1];
    refMain = bIsModeVer ? refAbove + 1 : refLeft + 1;
    refSide = bIsModeVer ? refLeft + 1 : refAbove + 1;

    // Extend main reference to right using replication
    const int log2Ratio = floorLog2(width) - floorLog2(height);
    const int s = std::max<int>(0, bIsModeVer ? log2Ratio : -log2Ratio);
#if JVET_W0123_TIMD_FUSION
    const int maxIndex  = (multiRefIdx << s) + 6;
#else
    const int maxIndex = (multiRefIdx << s) + 2;
#endif
    const int refLength = bIsModeVer ? m_topRefLength : m_leftRefLength;
    const Pel val = refMain[refLength + multiRefIdx];
    // right extended by 1 (z range)
    for (int z = 1; z <= (maxIndex + 1); z++)
    {
      refMain[refLength + multiRefIdx + z] = val;
    }
#else
    for (int x = 0; x <= m_topRefLength + multiRefIdx; x++)
    {
      refAbove[x] = pSrc.at(x, 0);
    }
    for (int y = 0; y <= m_leftRefLength + multiRefIdx; y++)
    {
      refLeft[y] = pSrc.at(y, 1);
    }

    refMain = bIsModeVer ? refAbove : refLeft;
    refSide = bIsModeVer ? refLeft : refAbove;

    // Extend main reference to right using replication
    const int log2Ratio = floorLog2(width) - floorLog2(height);
    const int s         = std::max<int>(0, bIsModeVer ? log2Ratio : -log2Ratio);
#if JVET_W0123_TIMD_FUSION
    const int maxIndex  = (multiRefIdx << s) + 6;
#else
    const int maxIndex  = (multiRefIdx << s) + 2;
#endif
    const int refLength = bIsModeVer ? m_topRefLength : m_leftRefLength;
    const Pel val       = refMain[refLength + multiRefIdx];
    for (int z = 1; z <= maxIndex; z++)
    {
      refMain[refLength + multiRefIdx + z] = val;
    }
#endif
  }

  // swap width/height if we are doing a horizontal mode:
  if (!bIsModeVer)
  {
    std::swap(width, height);
  }
  Pel       tempArray[MAX_CU_SIZE * MAX_CU_SIZE];
  const int dstStride = bIsModeVer ? pDst.stride : width;
  Pel *     pDstBuf   = bIsModeVer ? pDst.buf : tempArray;

  // compensate for line offset in reference line buffers
  refMain += multiRefIdx;
  refSide += multiRefIdx;

  Pel *pDsty = pDstBuf;

  if( intraPredAngle == 0 )  // pure vertical or pure horizontal
  {
    for( int y = 0; y < height; y++ )
    {
      ::memcpy(pDsty, refMain + 1, width * sizeof(Pel));

      if (m_ipaParam.applyPDPC)
      {
        const int scale   = (floorLog2(width) + floorLog2(height) - 2) >> 2;
        const Pel topLeft = refMain[0];
        const Pel left    = refSide[1 + y];
        for (int x = 0; x < std::min(3 << scale, width); x++)
        {
          const int wL  = 32 >> (2 * x >> scale);
          const Pel val = pDsty[x];
          pDsty[x]      = ClipPel(val + ((wL * (left - topLeft) + 32) >> 6), clpRng);
        }
      }

      pDsty += dstStride;
    }
  }
  else
  {
    for (int y = 0, deltaPos = intraPredAngle * (1 + multiRefIdx); y<height; y++, deltaPos += intraPredAngle, pDsty += dstStride)
    {
#if JVET_W0123_TIMD_FUSION
      const int deltaInt   = bExtIntraDir ? deltaPos >> 6 : deltaPos >> 5;
      const int deltaFract = bExtIntraDir ? deltaPos & 63 : deltaPos & 31;
#else
      const int deltaInt   = deltaPos >> 5;
      const int deltaFract = deltaPos & 31;
#endif

#if JVET_W0123_TIMD_FUSION
      bool bIntSlope = bExtIntraDir ? isIntegerSlopeExt( abs(intraPredAngle) ) : isIntegerSlope( abs(intraPredAngle) );
      if ( !bIntSlope )
#else
      if ( !isIntegerSlope( abs(intraPredAngle) ) )
#endif
      {
        if( isLuma(channelType) )
        {
          const bool useCubicFilter = !m_ipaParam.interpolationFlag;

#if INTRA_6TAP
          const TFilterCoeff        intraSmoothingFilter[6] = { TFilterCoeff(0), TFilterCoeff(64 - (deltaFract << 1)), TFilterCoeff(128 - (deltaFract << 1)), TFilterCoeff(64 + (deltaFract << 1)), TFilterCoeff(deltaFract << 1), TFilterCoeff(0) };
          const TFilterCoeff        intraSmoothingFilter2[6] = { TFilterCoeff(16 - (deltaFract >> 1)), TFilterCoeff(64 - 3*(deltaFract >> 1)), TFilterCoeff(96 - (deltaFract)), TFilterCoeff(64 + (deltaFract)),
            TFilterCoeff(16 + 3*(deltaFract >> 1)), TFilterCoeff((deltaFract >> 1)) };
#if JVET_W0123_TIMD_FUSION
          const TFilterCoeff        intraSmoothingFilterExt[6] = { TFilterCoeff(0), TFilterCoeff(64 - (deltaFract)), TFilterCoeff(128 - (deltaFract)), TFilterCoeff(64 + (deltaFract)), TFilterCoeff(deltaFract), TFilterCoeff(0) };
          const TFilterCoeff        intraSmoothingFilter2Ext[6] = { TFilterCoeff(16 - (deltaFract >> 2)), TFilterCoeff(64 - 3*(deltaFract >> 2)), TFilterCoeff(96 - (deltaFract >> 1)), TFilterCoeff(64 + (deltaFract >> 1)),
            TFilterCoeff(16 + 3*(deltaFract >> 2)), TFilterCoeff((deltaFract >> 2)) };
          const TFilterCoeff* const f = (useCubicFilter) ? ( bExtIntraDir ? InterpolationFilter::getIntraLumaFilterTableExt(deltaFract) : InterpolationFilter::getIntraLumaFilterTable(deltaFract)) : (width >=32 && height >=32)? (bExtIntraDir ? intraSmoothingFilter2Ext : intraSmoothingFilter2) : (bExtIntraDir ? intraSmoothingFilterExt : intraSmoothingFilter);
#else
          const TFilterCoeff* const f = (useCubicFilter) ? InterpolationFilter::getIntraLumaFilterTable(deltaFract) : (width >=32 && height >=32)? intraSmoothingFilter2 : intraSmoothingFilter;
#endif
#else
#if IF_12TAP
          const TFilterCoeff        intraSmoothingFilter[4] = { TFilterCoeff(64 - (deltaFract << 1)), TFilterCoeff(128 - (deltaFract << 1)), TFilterCoeff(64 + (deltaFract << 1)), TFilterCoeff(deltaFract << 1) };   
#if JVET_W0123_TIMD_FUSION
          const TFilterCoeff        intraSmoothingFilterExt[4] = { TFilterCoeff(64 - (deltaFract)), TFilterCoeff(128 - (deltaFract)), TFilterCoeff(64 + (deltaFract)), TFilterCoeff(deltaFract) };
#endif
#else
          const TFilterCoeff        intraSmoothingFilter[4] = {TFilterCoeff(16 - (deltaFract >> 1)), TFilterCoeff(32 - (deltaFract >> 1)), TFilterCoeff(16 + (deltaFract >> 1)), TFilterCoeff(deltaFract >> 1)};
#endif  

#if JVET_W0123_TIMD_FUSION
          const TFilterCoeff* const f                       = (useCubicFilter) ? (bExtIntraDir ? InterpolationFilter::getExtIntraCubicFilter(deltaFract) : InterpolationFilter::getChromaFilterTable(deltaFract)) : (bExtIntraDir ? InterpolationFilter::getExtIntraGaussFilter(deltaFract) : intraSmoothingFilter);
#else
          const TFilterCoeff* const f                       = (useCubicFilter) ? InterpolationFilter::getChromaFilterTable(deltaFract) : intraSmoothingFilter;
#endif
#endif

          for (int x = 0; x < width; x++)
          {
#if INTRA_6TAP
            Pel* p = refMain + deltaInt + x - 1;
            int val32 = ((int)f[0] * p[0] + (int)f[1] * p[1] + (int)f[2] * p[2] + (int)f[3] * p[3] + (int)f[4] * p[4] + (int)f[5] * p[5] + 128) >> 8;
            Pel val = (Pel)ClipPel(val32, clpRng);
#else
            Pel p[4];

            p[0] = refMain[deltaInt + x];
            p[1] = refMain[deltaInt + x + 1];
            p[2] = refMain[deltaInt + x + 2];
            p[3] = refMain[deltaInt + x + 3];

#if IF_12TAP
            Pel val = ( f[0] * p[0] + f[1] * p[1] + f[2] * p[2] + f[3] * p[3] + 128 ) >> 8;
#else
#if JVET_W0123_TIMD_FUSION
            int tOffset = 32;
            int tShift = 6;
            if (bExtIntraDir)
            {
              tOffset = 128;
              tShift = 8;
            }
            Pel val = (f[0] * p[0] + f[1] * p[1] + f[2] * p[2] + f[3] * p[3] + tOffset) >> tShift;
#else
            Pel val = (f[0] * p[0] + f[1] * p[1] + f[2] * p[2] + f[3] * p[3] + 32) >> 6;
#endif
#endif
#endif

            pDsty[x] = ClipPel(val, clpRng);   // always clip even though not always needed
          }
        }
        else
        {
          // Do linear filtering
          for (int x = 0; x < width; x++)
          {
            Pel* p = refMain + deltaInt + x + 1;
            pDsty[x] = p[0] + ((deltaFract * (p[1] - p[0]) + 16) >> 5);
          }
        }
      }
      else
      {
        // Just copy the integer samples
        ::memcpy(pDsty, refMain + deltaInt + 1, width * sizeof(Pel));
      }
#if GRAD_PDPC
      if (m_ipaParam.applyPDPC && m_ipaParam.useGradPDPC)
      {
        const int scale = m_ipaParam.angularScale;
        const Pel left = refSide[1 + y];
#if JVET_W0123_TIMD_FUSION
        int gradOffset = 16;
        int gradShift = 5;
        if (bExtIntraDir)
        {
          gradOffset = 32;
          gradShift = 6;
        }
        const Pel topLeft = refMain[deltaInt] + ((deltaFract * (refMain[deltaInt + 1] - refMain[deltaInt]) + gradOffset) >> gradShift);
#else
        const Pel topLeft = refMain[deltaInt] + ((deltaFract * (refMain[deltaInt + 1] - refMain[deltaInt]) + 16) >> 5);
#endif

        for (int x = 0; x < std::min(3 << scale, width); x++)
        {
          int wL = 32 >> (2 * x >> scale);
          const Pel val = pDsty[x];
          pDsty[x] = ClipPel(val + ((wL * (left - topLeft) + 32) >> 6), clpRng);
        }
      }
      else
#endif
      if (m_ipaParam.applyPDPC)
      {
        const int scale       = m_ipaParam.angularScale;
        int       invAngleSum = 256;

        for (int x = 0; x < std::min(3 << scale, width); x++)
        {
          invAngleSum += absInvAngle;

          int wL   = 32 >> (2 * x >> scale);
          Pel left = refSide[y + (invAngleSum >> 9) + 1];
          pDsty[x] = pDsty[x] + ((wL * (left - pDsty[x]) + 32) >> 6);
        }
      }
    }
  }

  // Flip the block if this is the horizontal mode
  if( !bIsModeVer )
  {
    for( int y = 0; y < height; y++ )
    {
      Pel *dst = pDst.buf + y;

      for( int x = 0; x < width; x++ )
      {
        *dst = pDstBuf[x];
        dst += pDst.stride;
      }
      pDstBuf += dstStride;
    }
  }
}

void IntraPrediction::xPredIntraBDPCM(const CPelBuf &pSrc, PelBuf &pDst, const uint32_t dirMode, const ClpRng& clpRng )
{
  const int wdt = pDst.width;
  const int hgt = pDst.height;

  const int strideP = pDst.stride;
  const int strideS = pSrc.stride;

  CHECK( !( dirMode == 1 || dirMode == 2 ), "Incorrect BDPCM mode parameter." );

  Pel* pred = &pDst.buf[0];
  if( dirMode == 1 )
  {
    Pel  val;
    for( int y = 0; y < hgt; y++ )
    {
      val = pSrc.buf[(y + 1) + strideS];
      for( int x = 0; x < wdt; x++ )
      {
        pred[x] = val;
      }
      pred += strideP;
    }
  }
  else
  {
    for( int y = 0; y < hgt; y++ )
    {
      for( int x = 0; x < wdt; x++ )
      {
        pred[x] = pSrc.buf[x + 1];
      }
      pred += strideP;
    }
  }
}

// Explicit instanciation since the implementation is in cpp file
template void IntraPrediction::geneWeightedPred<false>( const ComponentID compId, AreaBuf<Pel>& pred, const PredictionUnit &pu, const AreaBuf<Pel> &interPred, const AreaBuf<Pel> &intraPred, const Pel* pLUT );
template void IntraPrediction::geneWeightedPred<true>( const ComponentID compId, AreaBuf<Pel>& pred, const PredictionUnit &pu, const AreaBuf<Pel> &interPred, const AreaBuf<Pel> &intraPred, const Pel* pLUT );

template<bool lmcs>
void IntraPrediction::geneWeightedPred( const ComponentID compId, PelBuf& pred, const PredictionUnit &pu, const PelBuf& interPred, const PelBuf& intraPred, const Pel* pLUT )
{
  const int            width = pred.width;
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  CHECK(width == 2, "Width of 2 is not supported");
#endif
  const int            height = pred.height;
  Pel*                 interPredBuf = interPred.buf;
  const int            interPredStride = interPred.stride;
  Pel*                 intraPredBuf = intraPred.buf;
  const int            intraPredStride = intraPred.stride;
  const int            dstStride = pred.stride;

  Pel*                 dstBuf = pred.buf;
#if CIIP_PDPC
  if (!pu.ciipPDPC)
  {
#endif
  int wIntra, wMerge;

  const Position posBL = pu.Y().bottomLeft();
  const Position posTR = pu.Y().topRight();
  const PredictionUnit *neigh0 = pu.cs->getPURestricted(posBL.offset(-1, 0), pu, CHANNEL_TYPE_LUMA);
  const PredictionUnit *neigh1 = pu.cs->getPURestricted(posTR.offset(0, -1), pu, CHANNEL_TYPE_LUMA);
  bool isNeigh0Intra = neigh0 && (CU::isIntra(*neigh0->cu));
  bool isNeigh1Intra = neigh1 && (CU::isIntra(*neigh1->cu));

  if (isNeigh0Intra && isNeigh1Intra)
  {
    wIntra = 3; wMerge = 1;
  }
  else
  {
    if (!isNeigh0Intra && !isNeigh1Intra)
    {
      wIntra = 1; wMerge = 3;
    }
    else
    {
      wIntra = 2; wMerge = 2;
    }
  }
  for (int y = 0; y < height; y++)
  {
    for( int x = 0; x < width; x++ )
    {
      if( lmcs )
      {
        dstBuf[x] = ( wMerge * pLUT[interPredBuf[x]] + wIntra * intraPredBuf[x] + 2 ) >> 2;
      }
      else
      {
        dstBuf[x] = ( wMerge * interPredBuf[x] + wIntra * intraPredBuf[x] + 2 ) >> 2;
      }
    }

    dstBuf += dstStride;
    interPredBuf += interPredStride;
    intraPredBuf += intraPredStride;
  }
#if CIIP_PDPC
  }
  else
  {
    const int scale = ((floorLog2(width) - 2 + floorLog2(height) - 2 + 2) >> 2);
    for (int y = 0; y < height; y++)
    {
      const int wT = 32 >> std::min(31, ((y << 1) >> scale));
      for (int x = 0; x < width; x++)
      {
        const int wL = 32 >> std::min(31, ((x << 1) >> scale));

        if( lmcs )
        {
          dstBuf[x] = ( Pel ) ClipPel( ( ( int( intraPredBuf[x] ) << 6 ) + ( 64 - wT - wL ) * int( pLUT[interPredBuf[x]] ) + 32 ) >> 6, pu.cs->slice->clpRng( compId ) );
        }
        else
        {
          dstBuf[x] = ( Pel ) ClipPel( ( ( int( intraPredBuf[x] ) << 6 ) + ( 64 - wT - wL ) * int( interPredBuf[x] ) + 32 ) >> 6, pu.cs->slice->clpRng( compId ) );
        }
      }

      dstBuf += dstStride;
      interPredBuf += interPredStride;
      intraPredBuf += intraPredStride;
    }
  }
#endif
}

void IntraPrediction::geneIntrainterPred(const CodingUnit &cu, PelStorage& pred)
{
  if (!cu.firstPU->ciipFlag)
  {
    return;
  }

  const PredictionUnit* pu = cu.firstPU;

  initIntraPatternChType(cu, pu->Y());

  const UnitArea localUnitArea(pu->cs->area.chromaFormat, Area(0, 0, pu->Y().width, pu->Y().height));
  PelBuf ciipBuff = pred.getBuf(localUnitArea.Y());
  predIntraAng(COMPONENT_Y, ciipBuff, *pu);

  if (isChromaEnabled(pu->chromaFormat))
  {
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  if (pu->chromaSize().width > 2)
  {
#endif
    initIntraPatternChType(cu, pu->Cb());
    PelBuf ciipBuffCb = pred.getBuf(localUnitArea.Cb());
    predIntraAng(COMPONENT_Cb, ciipBuffCb, *pu);

    initIntraPatternChType(cu, pu->Cr());
    PelBuf ciipBuffCr = pred.getBuf(localUnitArea.Cr());
    predIntraAng(COMPONENT_Cr, ciipBuffCr, *pu);
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  }
#endif
  }
}

#if !MERGE_ENC_OPT
void IntraPrediction::geneWeightedPred( const ComponentID compId, PelBuf &pred, const PredictionUnit &pu, Pel *srcBuf )
{
  const int            width = pred.width;
  CHECK( width == 2, "Width of 2 is not supported" );
  const int            height = pred.height;
  const int            srcStride = width;
  const int            dstStride = pred.stride;

  Pel*                 dstBuf = pred.buf;
  int wIntra, wMerge;

  const Position posBL = pu.Y().bottomLeft();
  const Position posTR = pu.Y().topRight();
  const PredictionUnit *neigh0 = pu.cs->getPURestricted( posBL.offset( -1, 0 ), pu, CHANNEL_TYPE_LUMA );
  const PredictionUnit *neigh1 = pu.cs->getPURestricted( posTR.offset( 0, -1 ), pu, CHANNEL_TYPE_LUMA );
  bool isNeigh0Intra = neigh0 && ( CU::isIntra( *neigh0->cu ) );
  bool isNeigh1Intra = neigh1 && ( CU::isIntra( *neigh1->cu ) );

  if( isNeigh0Intra && isNeigh1Intra )
  {
    wIntra = 3; wMerge = 1;
  }
  else
  {
    if( !isNeigh0Intra && !isNeigh1Intra )
    {
      wIntra = 1; wMerge = 3;
    }
    else
    {
      wIntra = 2; wMerge = 2;
    }
  }
  for( int y = 0; y < height; y++ )
  {
    for( int x = 0; x < width; x++ )
    {
      dstBuf[y*dstStride + x] = ( wMerge * dstBuf[y*dstStride + x] + wIntra * srcBuf[y*srcStride + x] + 2 ) >> 2;
    }
  }
}

void IntraPrediction::switchBuffer( const PredictionUnit &pu, ComponentID compID, PelBuf srcBuff, Pel *dst )
{
  Pel  *src = srcBuff.bufAt( 0, 0 );
  int compWidth = compID == COMPONENT_Y ? pu.Y().width : pu.Cb().width;
  int compHeight = compID == COMPONENT_Y ? pu.Y().height : pu.Cb().height;
  for( int i = 0; i < compHeight; i++ )
  {
    memcpy( dst, src, compWidth * sizeof( Pel ) );
    src += srcBuff.stride;
    dst += compWidth;
  }
}
#endif

inline bool isAboveLeftAvailable  ( const CodingUnit &cu, const ChannelType &chType, const Position &posLT );
inline int  isAboveAvailable      ( const CodingUnit &cu, const ChannelType &chType, const Position &posLT, const uint32_t uiNumUnitsInPU, const uint32_t unitWidth, bool *validFlags );
inline int  isLeftAvailable       ( const CodingUnit &cu, const ChannelType &chType, const Position &posLT, const uint32_t uiNumUnitsInPU, const uint32_t unitWidth, bool *validFlags );
inline int  isAboveRightAvailable ( const CodingUnit &cu, const ChannelType &chType, const Position &posRT, const uint32_t uiNumUnitsInPU, const uint32_t unitHeight, bool *validFlags );
inline int  isBelowLeftAvailable  ( const CodingUnit &cu, const ChannelType &chType, const Position &posLB, const uint32_t uiNumUnitsInPU, const uint32_t unitHeight, bool *validFlags );

void IntraPrediction::initIntraPatternChType(const CodingUnit &cu, const CompArea &area, const bool forceRefFilterFlag)
{
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  CHECK(area.width == 2, "Width of 2 is not supported");
#endif
  const CodingStructure& cs   = *cu.cs;

  if (!forceRefFilterFlag)
  {
    initPredIntraParams(*cu.firstPU, area, *cs.sps);
  }

  Pel *refBufUnfiltered = m_refBuffer[area.compID][PRED_BUF_UNFILTERED];
  Pel *refBufFiltered   = m_refBuffer[area.compID][PRED_BUF_FILTERED];

  setReferenceArrayLengths( area );

  // ----- Step 1: unfiltered reference samples -----
  xFillReferenceSamples( cs.picture->getRecoBuf( area ), refBufUnfiltered, area, cu );
  // ----- Step 2: filtered reference samples -----
  if( m_ipaParam.refFilterFlag || forceRefFilterFlag )
  {
    xFilterReferenceSamples( refBufUnfiltered, refBufFiltered, area, *cs.sps, cu.firstPU->multiRefIdx );
  }
}

void IntraPrediction::initIntraPatternChTypeISP(const CodingUnit& cu, const CompArea& area, PelBuf& recBuf, const bool forceRefFilterFlag)
{
  const CodingStructure& cs = *cu.cs;

  if (!forceRefFilterFlag)
  {
    initPredIntraParams(*cu.firstPU, area, *cs.sps);
  }

  const Position posLT = area;
  bool           isLeftAvail  = (cs.getCURestricted(posLT.offset(-1, 0), cu, CHANNEL_TYPE_LUMA) != NULL) && cs.isDecomp(posLT.offset(-1, 0), CHANNEL_TYPE_LUMA);
  bool           isAboveAvail = (cs.getCURestricted(posLT.offset(0, -1), cu, CHANNEL_TYPE_LUMA) != NULL) && cs.isDecomp(posLT.offset(0, -1), CHANNEL_TYPE_LUMA);
  // ----- Step 1: unfiltered reference samples -----
  if (cu.blocks[area.compID].x == area.x && cu.blocks[area.compID].y == area.y)
  {
    Pel *refBufUnfiltered = m_refBuffer[area.compID][PRED_BUF_UNFILTERED];
    // With the first subpartition all the CU reference samples are fetched at once in a single call to xFillReferenceSamples
    if (cu.ispMode == HOR_INTRA_SUBPARTITIONS)
    {
      m_leftRefLength = cu.Y().height << 1;
      m_topRefLength = cu.Y().width + area.width;
    }
    else //if (cu.ispMode == VER_INTRA_SUBPARTITIONS)
    {
      m_leftRefLength = cu.Y().height + area.height;
      m_topRefLength = cu.Y().width << 1;
    }


    xFillReferenceSamples(cs.picture->getRecoBuf(cu.Y()), refBufUnfiltered, cu.Y(), cu);

    // After having retrieved all the CU reference samples, the number of reference samples is now adjusted for the current subpartition
    m_topRefLength = cu.blocks[area.compID].width + area.width;
    m_leftRefLength = cu.blocks[area.compID].height + area.height;
  }
  else
  {
    m_topRefLength = cu.blocks[area.compID].width + area.width;
    m_leftRefLength = cu.blocks[area.compID].height + area.height;

    const int predSizeHor = m_topRefLength;
    const int predSizeVer = m_leftRefLength;
    if (cu.ispMode == HOR_INTRA_SUBPARTITIONS)
    {
      Pel* src = recBuf.bufAt(0, -1);
      Pel *ref = m_refBuffer[area.compID][PRED_BUF_UNFILTERED] + m_refBufferStride[area.compID];
      if (isLeftAvail)
      {
        for (int i = 0; i <= 2 * cu.blocks[area.compID].height - area.height; i++)
        {
          ref[i] = ref[i + area.height];
        }
      }
      else
      {
        for (int i = 0; i <= predSizeVer; i++)
        {
          ref[i] = src[0];
        }
      }
      Pel *dst = m_refBuffer[area.compID][PRED_BUF_UNFILTERED] + 1;
      dst[-1]  = ref[0];
      for (int i = 0; i < area.width; i++)
      {
        dst[i] = src[i];
      }
      Pel sample = src[area.width - 1];
      dst += area.width;
      for (int i = 0; i < predSizeHor - area.width; i++)
      {
        dst[i] = sample;
      }
    }
    else
    {
      Pel* src = recBuf.bufAt(-1, 0);
      Pel *ref = m_refBuffer[area.compID][PRED_BUF_UNFILTERED];
      if (isAboveAvail)
      {
        for (int i = 0; i <= 2 * cu.blocks[area.compID].width - area.width; i++)
        {
          ref[i] = ref[i + area.width];
        }
      }
      else
      {
        for (int i = 0; i <= predSizeHor; i++)
        {
          ref[i] = src[0];
        }
      }
      Pel *dst = m_refBuffer[area.compID][PRED_BUF_UNFILTERED] + m_refBufferStride[area.compID] + 1;
      dst[-1]  = ref[0];
      for (int i = 0; i < area.height; i++)
      {
        *dst = *src;
        src += recBuf.stride;
        dst++;
      }
      Pel sample = src[-recBuf.stride];
      for (int i = 0; i < predSizeVer - area.height; i++)
      {
        *dst = sample;
        dst++;
      }
    }
  }
  // ----- Step 2: filtered reference samples -----
  if (m_ipaParam.refFilterFlag || forceRefFilterFlag)
  {
    Pel *refBufUnfiltered = m_refBuffer[area.compID][PRED_BUF_UNFILTERED];
    Pel *refBufFiltered   = m_refBuffer[area.compID][PRED_BUF_FILTERED];
    xFilterReferenceSamples(refBufUnfiltered, refBufFiltered, area, *cs.sps, cu.firstPU->multiRefIdx);
  }
}

#if JVET_V0130_INTRA_TMP
#if JVET_W0069_TMP_BOUNDARY
RefTemplateType IntraPrediction::GetRefTemplateType(CodingUnit& cu, CompArea& area)
#else
bool IntraPrediction::isRefTemplateAvailable(CodingUnit& cu, CompArea& area)
#endif
{
	const ChannelType      chType = toChannelType(area.compID);
	const CodingStructure& cs = *cu.cs;
	const SPS& sps = *cs.sps;
	const PreCalcValues& pcv = *cs.pcv;

	const int  tuWidth = area.width;
	const int  tuHeight = area.height;
	const int  predSize = m_topRefLength;
	const int  predHSize = m_leftRefLength;
	//const int predStride = predSize;

	const int  unitWidth = pcv.minCUWidth >> getComponentScaleX(area.compID, sps.getChromaFormatIdc());
	const int  unitHeight = pcv.minCUHeight >> getComponentScaleY(area.compID, sps.getChromaFormatIdc());

	const int  totalAboveUnits = (predSize + (unitWidth - 1)) / unitWidth;
	const int  totalLeftUnits = (predHSize + (unitHeight - 1)) / unitHeight;
	const int  totalUnits = totalAboveUnits + totalLeftUnits + 1; //+1 for top-left
	const int  numAboveUnits = std::max<int>(tuWidth / unitWidth, 1);
	const int  numLeftUnits = std::max<int>(tuHeight / unitHeight, 1);
	const int  numAboveRightUnits = totalAboveUnits - numAboveUnits;
	const int  numLeftBelowUnits = totalLeftUnits - numLeftUnits;

  if( numAboveUnits <= 0 || numLeftUnits <= 0 || numAboveRightUnits <= 0 || numLeftBelowUnits <= 0 )
  {
#if JVET_W0069_TMP_BOUNDARY
	  return NO_TEMPLATE;
#else
    return false;
#endif
  }

	// ----- Step 1: analyze neighborhood -----
	const Position posLT = area;
	//const Position posRT = area.topRight();
	//const Position posLB = area.bottomLeft();

	bool  neighborFlags[4 * MAX_NUM_PART_IDXS_IN_CTU_WIDTH + 1];
	//int   numIntraNeighbor = 0;

	memset(neighborFlags, 0, totalUnits);

	//bool retVal = 1;

#if JVET_W0069_TMP_BOUNDARY
	if (isAboveLeftAvailable(cu, chType, posLT) && isAboveAvailable(cu, chType, posLT, numAboveUnits, unitWidth, (neighborFlags + totalLeftUnits + 1)) && isLeftAvailable(cu, chType, posLT, numLeftUnits, unitHeight, (neighborFlags + totalLeftUnits - 1)))
		return L_SHAPE_TEMPLATE;
	else if (isAboveLeftAvailable(cu, chType, posLT))
		return LEFT_TEMPLATE;
	else if (isAboveAvailable(cu, chType, posLT, numAboveUnits, unitWidth, (neighborFlags + totalLeftUnits + 1)))
		return ABOVE_TEMPLATE;
	else
		return NO_TEMPLATE;
	CHECK(1, "un defined template type");
#else
	return isAboveLeftAvailable(cu, chType, posLT) && isAboveAvailable(cu, chType, posLT, numAboveUnits, unitWidth, (neighborFlags + totalLeftUnits + 1)) && isLeftAvailable(cu, chType, posLT, numLeftUnits, unitHeight, (neighborFlags + totalLeftUnits - 1));
#endif

	//return retVal;
}
#endif
void IntraPrediction::xFillReferenceSamples( const CPelBuf &recoBuf, Pel* refBufUnfiltered, const CompArea &area, const CodingUnit &cu )
{
  const ChannelType      chType = toChannelType( area.compID );
  const CodingStructure &cs     = *cu.cs;
  const SPS             &sps    = *cs.sps;
  const PreCalcValues   &pcv    = *cs.pcv;

  const int multiRefIdx         = (area.compID == COMPONENT_Y) ? cu.firstPU->multiRefIdx : 0;

  const int  tuWidth            = area.width;
  const int  tuHeight           = area.height;
  const int  predSize           = m_topRefLength;
  const int  predHSize          = m_leftRefLength;
  const int predStride = predSize + 1 + multiRefIdx;
  m_refBufferStride[area.compID] = predStride;

  const bool noShift            = pcv.noChroma2x2 && area.width == 4; // don't shift on the lowest level (chroma not-split)
  const int  unitWidth          = tuWidth  <= 2 && cu.ispMode && isLuma(area.compID) ? tuWidth  : pcv.minCUWidth  >> (noShift ? 0 : getComponentScaleX(area.compID, sps.getChromaFormatIdc()));
  const int  unitHeight         = tuHeight <= 2 && cu.ispMode && isLuma(area.compID) ? tuHeight : pcv.minCUHeight >> (noShift ? 0 : getComponentScaleY(area.compID, sps.getChromaFormatIdc()));

  const int  totalAboveUnits    = (predSize + (unitWidth - 1)) / unitWidth;
  const int  totalLeftUnits     = (predHSize + (unitHeight - 1)) / unitHeight;
  const int  totalUnits         = totalAboveUnits + totalLeftUnits + 1; //+1 for top-left
  const int  numAboveUnits      = std::max<int>( tuWidth / unitWidth, 1 );
  const int  numLeftUnits       = std::max<int>( tuHeight / unitHeight, 1 );
  const int  numAboveRightUnits = totalAboveUnits - numAboveUnits;
  const int  numLeftBelowUnits  = totalLeftUnits - numLeftUnits;

  CHECK( numAboveUnits <= 0 || numLeftUnits <= 0 || numAboveRightUnits <= 0 || numLeftBelowUnits <= 0, "Size not supported" );

  // ----- Step 1: analyze neighborhood -----
  const Position posLT          = area;
  const Position posRT          = area.topRight();
  const Position posLB          = area.bottomLeft();

  bool  neighborFlags[4 * MAX_NUM_PART_IDXS_IN_CTU_WIDTH + 1];
  int   numIntraNeighbor = 0;

  memset( neighborFlags, 0, totalUnits );

  neighborFlags[totalLeftUnits] = isAboveLeftAvailable( cu, chType, posLT );
  numIntraNeighbor += neighborFlags[totalLeftUnits] ? 1 : 0;
  numIntraNeighbor += isAboveAvailable     ( cu, chType, posLT, numAboveUnits,      unitWidth,  (neighborFlags + totalLeftUnits + 1) );
  numIntraNeighbor += isAboveRightAvailable( cu, chType, posRT, numAboveRightUnits, unitWidth,  (neighborFlags + totalLeftUnits + 1 + numAboveUnits) );
  numIntraNeighbor += isLeftAvailable      ( cu, chType, posLT, numLeftUnits,       unitHeight, (neighborFlags + totalLeftUnits - 1) );
  numIntraNeighbor += isBelowLeftAvailable ( cu, chType, posLB, numLeftBelowUnits,  unitHeight, (neighborFlags + totalLeftUnits - 1 - numLeftUnits) );

  // ----- Step 2: fill reference samples (depending on neighborhood) -----

  const Pel*  srcBuf    = recoBuf.buf;
  const int   srcStride = recoBuf.stride;
        Pel*  ptrDst    = refBufUnfiltered;
  const Pel*  ptrSrc;
  const Pel   valueDC   = 1 << (sps.getBitDepth( chType ) - 1);

  if( numIntraNeighbor == 0 )
  {
    // Fill border with DC value
    for (int j = 0; j <= predSize + multiRefIdx; j++)
    {
      ptrDst[j] = valueDC;
    }
    for (int i = 0; i <= predHSize + multiRefIdx; i++)
    {
      ptrDst[i + predStride] = valueDC;
    }
  }
  else if( numIntraNeighbor == totalUnits )
  {
    // Fill top-left border and top and top right with rec. samples
    ptrSrc = srcBuf - (1 + multiRefIdx) * srcStride - (1 + multiRefIdx);

    std::memcpy( ptrDst, ptrSrc, sizeof( *ptrDst ) * ( predSize + multiRefIdx + 1 ) );

    for (int i = 0; i <= predHSize + multiRefIdx; i++)
    {
      ptrDst[i + predStride] = ptrSrc[i * srcStride];
    }
  }
  else // reference samples are partially available
  {
    // Fill top-left sample(s) if available
    ptrSrc = srcBuf - (1 + multiRefIdx) * srcStride - (1 + multiRefIdx);
    ptrDst = refBufUnfiltered;
    if (neighborFlags[totalLeftUnits])
    {
      ptrDst[0] = ptrSrc[0];
      ptrDst[predStride] = ptrSrc[0];
      for (int i = 1; i <= multiRefIdx; i++)
      {
        ptrDst[i] = ptrSrc[i];
        ptrDst[i + predStride] = ptrSrc[i * srcStride];
      }
    }

    // Fill left & below-left samples if available (downwards)
    ptrSrc += (1 + multiRefIdx) * srcStride;
    ptrDst += (1 + multiRefIdx) + predStride;
    for (int unitIdx = totalLeftUnits - 1; unitIdx > 0; unitIdx--)
    {
      if (neighborFlags[unitIdx])
      {
        for (int i = 0; i < unitHeight; i++)
        {
          ptrDst[i] = ptrSrc[i * srcStride];
        }
      }
      ptrSrc += unitHeight * srcStride;
      ptrDst += unitHeight;
    }
    // Fill last below-left sample(s)
    if (neighborFlags[0])
    {
      int lastSample = (predHSize % unitHeight == 0) ? unitHeight : predHSize % unitHeight;
      for (int i = 0; i < lastSample; i++)
      {
        ptrDst[i] = ptrSrc[i * srcStride];
      }
    }

    // Fill above & above-right samples if available (left-to-right)
    ptrSrc = srcBuf - srcStride * (1 + multiRefIdx);
    ptrDst = refBufUnfiltered + 1 + multiRefIdx;
    for (int unitIdx = totalLeftUnits + 1; unitIdx < totalUnits - 1; unitIdx++)
    {
      if (neighborFlags[unitIdx])
      {
        for (int j = 0; j < unitWidth; j++)
        {
          ptrDst[j] = ptrSrc[j];
        }
      }
      ptrSrc += unitWidth;
      ptrDst += unitWidth;
    }
    // Fill last above-right sample(s)
    if (neighborFlags[totalUnits - 1])
    {
      int lastSample = (predSize % unitWidth == 0) ? unitWidth : predSize % unitWidth;
      for (int j = 0; j < lastSample; j++)
      {
        ptrDst[j] = ptrSrc[j];
      }
    }

    // pad from first available down to the last below-left
    ptrDst = refBufUnfiltered;
    int lastAvailUnit = 0;
    if (!neighborFlags[0])
    {
      int firstAvailUnit = 1;
      while (firstAvailUnit < totalUnits && !neighborFlags[firstAvailUnit])
      {
        firstAvailUnit++;
      }

      // first available sample
      int firstAvailRow = -1;
      int firstAvailCol = 0;
      if (firstAvailUnit < totalLeftUnits)
      {
        firstAvailRow = (totalLeftUnits - firstAvailUnit) * unitHeight + multiRefIdx;
      }
      else if (firstAvailUnit == totalLeftUnits)
      {
        firstAvailRow = multiRefIdx;
      }
      else
      {
        firstAvailCol = (firstAvailUnit - totalLeftUnits - 1) * unitWidth + 1 + multiRefIdx;
      }
      const Pel firstAvailSample = ptrDst[firstAvailRow < 0 ? firstAvailCol : firstAvailRow + predStride];

      // last sample below-left (n.a.)
      int lastRow = predHSize + multiRefIdx;

      // fill left column
      for (int i = lastRow; i > firstAvailRow; i--)
      {
        ptrDst[i + predStride] = firstAvailSample;
      }
      // fill top row
      if (firstAvailCol > 0)
      {
        for (int j = 0; j < firstAvailCol; j++)
        {
          ptrDst[j] = firstAvailSample;
        }
      }
      lastAvailUnit = firstAvailUnit;
    }

    // pad all other reference samples.
    int currUnit = lastAvailUnit + 1;
    while (currUnit < totalUnits)
    {
      if (!neighborFlags[currUnit]) // samples not available
      {
        // last available sample
        int lastAvailRow = -1;
        int lastAvailCol = 0;
        if (lastAvailUnit < totalLeftUnits)
        {
          lastAvailRow = (totalLeftUnits - lastAvailUnit - 1) * unitHeight + multiRefIdx + 1;
        }
        else if (lastAvailUnit == totalLeftUnits)
        {
          lastAvailCol = multiRefIdx;
        }
        else
        {
          lastAvailCol = (lastAvailUnit - totalLeftUnits) * unitWidth + multiRefIdx;
        }
        const Pel lastAvailSample = ptrDst[lastAvailRow < 0 ? lastAvailCol : lastAvailRow + predStride];

        // fill current unit with last available sample
        if (currUnit < totalLeftUnits)
        {
          for (int i = lastAvailRow - 1; i >= lastAvailRow - unitHeight; i--)
          {
            ptrDst[i + predStride] = lastAvailSample;
          }
        }
        else if (currUnit == totalLeftUnits)
        {
          for (int i = 0; i < multiRefIdx + 1; i++)
          {
            ptrDst[i + predStride] = lastAvailSample;
          }
          for (int j = 0; j < multiRefIdx + 1; j++)
          {
            ptrDst[j] = lastAvailSample;
          }
        }
        else
        {
          int numSamplesInUnit = (currUnit == totalUnits - 1) ? ((predSize % unitWidth == 0) ? unitWidth : predSize % unitWidth) : unitWidth;
          for (int j = lastAvailCol + 1; j <= lastAvailCol + numSamplesInUnit; j++)
          {
            ptrDst[j] = lastAvailSample;
          }
        }
      }
      lastAvailUnit = currUnit;
      currUnit++;
    }
  }
}

void IntraPrediction::xFilterReferenceSamples(const Pel *refBufUnfiltered, Pel *refBufFiltered, const CompArea &area,
                                              const SPS &sps, int multiRefIdx)
{
  if (area.compID != COMPONENT_Y)
  {
    multiRefIdx = 0;
  }
  const int predSize = m_topRefLength + multiRefIdx;
  const int predHSize = m_leftRefLength + multiRefIdx;
  const size_t predStride = m_refBufferStride[area.compID];

  const Pel topLeft =
    (refBufUnfiltered[0] + refBufUnfiltered[1] + refBufUnfiltered[predStride] + refBufUnfiltered[predStride + 1] + 2)
    >> 2;

  refBufFiltered[0] = topLeft;

  for (int i = 1; i < predSize; i++)
  {
    refBufFiltered[i] = (refBufUnfiltered[i - 1] + 2 * refBufUnfiltered[i] + refBufUnfiltered[i + 1] + 2) >> 2;
  }
  refBufFiltered[predSize] = refBufUnfiltered[predSize];

  refBufFiltered += predStride;
  refBufUnfiltered += predStride;

  refBufFiltered[0] = topLeft;

  for (int i = 1; i < predHSize; i++)
  {
    refBufFiltered[i] = (refBufUnfiltered[i - 1] + 2 * refBufUnfiltered[i] + refBufUnfiltered[i + 1] + 2) >> 2;
  }
  refBufFiltered[predHSize] = refBufUnfiltered[predHSize];
}

#if JVET_W0123_TIMD_FUSION
Pel IntraPrediction::xGetPredTimdValDc( const CPelBuf &pSrc, const Size &dstSize, TEMPLATE_TYPE eTempType, int iTempHeight, int iTempWidth )
{
  int idx, sum = 0;
  Pel dcVal;
  const int width  = dstSize.width;
  const int height = dstSize.height;
  auto denom     = (width == height) ? (width << 1) : std::max(width,height);
  auto divShift  = floorLog2(denom);
  auto divOffset = (denom >> 1);

  if (eTempType == LEFT_NEIGHBOR)
  {
    denom = height;
    divShift = floorLog2(denom);
    divOffset = (denom >> 1);
    for(idx = 0; idx < height; idx++)
      sum += pSrc.at(1 + idx, 1);
    dcVal = (sum + divOffset) >> divShift;
    return dcVal;
  }
  else if (eTempType == ABOVE_NEIGHBOR)
  {
    denom = width;
    divShift = floorLog2(denom);
    divOffset = (denom >> 1);
    for(idx = 0; idx < width; idx++)
      sum += pSrc.at( 1 + idx, 0);
    dcVal = (sum + divOffset) >> divShift;
    return dcVal;
  }

  if ( width >= height )
  {
    for( idx = 0; idx < width; idx++ )
    {
      sum += pSrc.at(iTempWidth + 1 + idx, 0);
    }
  }
  if ( width <= height )
  {
    for( idx = 0; idx < height; idx++ )
    {
      sum += pSrc.at(iTempHeight + 1 + idx, 1);
    }
  }
  dcVal = (sum + divOffset) >> divShift;
  return dcVal;
}

void IntraPrediction::predTimdIntraAng( const ComponentID compId, const PredictionUnit &pu, uint32_t uiDirMode, Pel* pPred, uint32_t uiStride, uint32_t iWidth, uint32_t iHeight, TEMPLATE_TYPE eTempType, int32_t iTemplateWidth, int32_t iTemplateHeight)
{
  const ComponentID compID       = MAP_CHROMA( compId );

  const int srcStride  = m_refBufferStride[compID];
  const int srcHStride = 2;

  const CPelBuf & srcBuf = CPelBuf(getPredictorPtr(compID), srcStride, srcHStride);
  const ClpRng& clpRng(pu.cu->cs->slice->clpRng(compID));

  switch (uiDirMode)
  {
    case(PLANAR_IDX): xPredTimdIntraPlanar(srcBuf, pPred, uiStride, iWidth, iHeight, eTempType, iTemplateWidth, iTemplateHeight); break;
    case(DC_IDX):     xPredTimdIntraDc(pu, srcBuf, pPred, uiStride, iWidth, iHeight, eTempType, iTemplateWidth, iTemplateHeight); break;
    default:          xPredTimdIntraAng(srcBuf, clpRng, pPred, uiStride, iWidth, iHeight, eTempType, iTemplateWidth, iTemplateHeight, uiDirMode); break;
  }

  if (m_ipaParam.applyPDPC && (uiDirMode == PLANAR_IDX || uiDirMode == DC_IDX))
  {
    xIntraPredTimdPlanarDcPdpc(srcBuf, pPred, uiStride, iWidth, iHeight, eTempType, iTemplateWidth, iTemplateHeight);
  }
}

void IntraPrediction::xPredTimdIntraPlanar( const CPelBuf &pSrc, Pel* rpDst, int iDstStride, int width, int height, TEMPLATE_TYPE eTempType, int iTemplateWidth, int iTemplateHeight )
{
  static int leftColumn[MAX_CU_SIZE+DIMD_MAX_TEMP_SIZE+1] = {0}, topRow[MAX_CU_SIZE+DIMD_MAX_TEMP_SIZE+1] ={0}, bottomRow[MAX_CU_SIZE+DIMD_MAX_TEMP_SIZE] = {0}, rightColumn[MAX_CU_SIZE+DIMD_MAX_TEMP_SIZE]={0};
  if(eTempType == LEFT_ABOVE_NEIGHBOR)
  {
    //predict above template
    {
      uint32_t w = width - iTemplateWidth;
      const uint32_t log2W = floorLog2( w );
      const uint32_t log2H = floorLog2( iTemplateHeight );
      const uint32_t offset = 1 << (log2W + log2H);
      for(int k = 0; k < w + 1; k++)
      {
        topRow[k] = pSrc.at( k + iTemplateWidth + 1, 0 );
      }
      for (int k=0; k < iTemplateHeight + 1; k++)
      {
        leftColumn[k] = pSrc.at( k + 1, 1 );
      }

      int bottomLeft = leftColumn[iTemplateHeight];
      int topRight = topRow[w];
      for(int k = 0; k < w; k++)
      {
        bottomRow[k]  = bottomLeft - topRow[k];
        topRow[k]     = topRow[k] << log2H;
      }
      for(int k = 0; k < iTemplateHeight; k++)
      {
        rightColumn[k]  = topRight - leftColumn[k];
        leftColumn[k]   = leftColumn[k] << log2W;
      }

      const uint32_t finalShift = 1 + log2W + log2H;
      for (int y = 0; y < iTemplateHeight; y++)
      {
        int horPred = leftColumn[y];
        for (int x = 0; x < w; x++)
        {
          horPred   += rightColumn[y];
          topRow[x] += bottomRow[x];
          int vertPred = topRow[x];
          rpDst[y*iDstStride+x + iTemplateWidth] = ( ( horPred << log2H ) + ( vertPred << log2W ) + offset ) >> finalShift;
        }
      }
    }

    //predict left template
    {
      uint32_t h = height - iTemplateHeight;
      const uint32_t log2W = floorLog2( iTemplateWidth );
      const uint32_t log2H = floorLog2( h );
      const uint32_t offset = 1 << (log2W + log2H);
      for (int k = 0; k < h + 1; k++)
      {
        leftColumn[k] = pSrc.at( k + iTemplateHeight + 1, 1 );
      }
      for(int k = 0; k < iTemplateWidth + 1; k++)
      {
        topRow[k] = pSrc.at( k + 1, 0 );
      }
      int bottomLeft = leftColumn[h];
      int topRight = topRow[iTemplateWidth];
      for(int k = 0; k < iTemplateWidth; k++)
      {
        bottomRow[k]  = bottomLeft - topRow[k];
        topRow[k]     = topRow[k] << log2H;
      }
      for(int k = 0; k < h; k++)
      {
        rightColumn[k]  = topRight - leftColumn[k];
        leftColumn[k]   = leftColumn[k] << log2W;
      }
      const uint32_t finalShift = 1 + log2W + log2H;
      for (int y = 0; y < height; y++)
      {
        int horPred = leftColumn[y];
        for (int x = 0; x < iTemplateWidth; x++)
        {
          horPred   += rightColumn[y];
          topRow[x] += bottomRow[x];
          int vertPred = topRow[x];
          rpDst[(y + iTemplateHeight)*iDstStride+x] = ( ( horPred << log2H ) + ( vertPred << log2W ) + offset ) >> finalShift;
        }
      }
    }
  }
  else if(eTempType == LEFT_NEIGHBOR)
  {
    const uint32_t log2W = floorLog2( iTemplateWidth );
    const uint32_t log2H = floorLog2( height );
    const uint32_t offset = 1 << (log2W + log2H);
    for (int k = 0; k < height + 1; k++)
    {
      leftColumn[k] = pSrc.at( k + iTemplateHeight + 1, 1 );
    }
    for(int k = 0; k < iTemplateWidth + 1; k++)
    {
      topRow[k] = pSrc.at( k + 1, 0 );
    }

    int bottomLeft = leftColumn[height];
    int topRight = topRow[iTemplateWidth];
    for(int k = 0; k < iTemplateWidth; k++)
    {
      bottomRow[k]  = bottomLeft - topRow[k];
      topRow[k]     = topRow[k] << log2H;
    }
    for(int k = 0; k < height; k++)
    {
      rightColumn[k]  = topRight - leftColumn[k];
      leftColumn[k]   = leftColumn[k] << log2W;
    }

    const uint32_t finalShift = 1 + log2W + log2H;
    for (int y = 0; y < height; y++)
    {
      int horPred = leftColumn[y];
      for (int x = 0; x < iTemplateWidth; x++)
      {
        horPred   += rightColumn[y];
        topRow[x] += bottomRow[x];
        int vertPred = topRow[x];
        rpDst[y*iDstStride+x] = ( ( horPred << log2H ) + ( vertPred << log2W ) + offset ) >> finalShift;
      }
    }
  }
  else if(eTempType == ABOVE_NEIGHBOR)
  {
    const uint32_t log2W = floorLog2( width );
    const uint32_t log2H = floorLog2( iTemplateHeight );
    const uint32_t offset = 1 << (log2W + log2H);
    for(int k = 0; k < width + 1; k++)
    {
      topRow[k] = pSrc.at( k + iTemplateWidth + 1, 0 );
    }
    for (int k=0; k < iTemplateHeight + 1; k++)
    {
      leftColumn[k] = pSrc.at( k + 1, 1 );
    }

    int bottomLeft = leftColumn[iTemplateHeight];
    int topRight = topRow[width];
    for(int k=0;k<width;k++)
    {
      bottomRow[k]  = bottomLeft - topRow[k];
      topRow[k]     = topRow[k] << log2H;
    }
    for(int k = 0; k < iTemplateHeight; k++)
    {
      rightColumn[k]  = topRight - leftColumn[k];
      leftColumn[k]   = leftColumn[k] << log2W;
    }

    const uint32_t finalShift = 1 + log2W + log2H;
    for (int y = 0; y < iTemplateHeight; y++)
    {
      int horPred = leftColumn[y];
      for (int x = 0; x < width; x++)
      {
        horPred   += rightColumn[y];
        topRow[x] += bottomRow[x];
        int vertPred = topRow[x];
        rpDst[y*iDstStride+x] = ( ( horPred << log2H ) + ( vertPred << log2W ) + offset ) >> finalShift;
      }
    }
  }
  else
  {
    assert(0);
  }
}

void IntraPrediction::xPredTimdIntraDc( const PredictionUnit &pu, const CPelBuf &pSrc, Pel* pDst, int iDstStride, int iWidth, int iHeight, TEMPLATE_TYPE eTempType, int iTemplateWidth, int iTemplateHeight )
{
  const Size &dstSize = Size(pu.lwidth(), pu.lheight());
  const Pel dcval = xGetPredTimdValDc( pSrc, dstSize, eTempType, iTemplateHeight, iTemplateWidth );
  if(eTempType == LEFT_ABOVE_NEIGHBOR)
  {
    for (int y = 0; y < iHeight; y++,pDst += iDstStride)
    {
      if(y < iTemplateHeight)
      {
        for (int x = iTemplateWidth; x < iWidth; x++)
        {
          pDst[x] = dcval;
        }
      }
      else
      {
        for (int x = 0; x < iTemplateWidth; x++)
        {
          pDst[x] = dcval;
        }
      }
    }
  }
  else if(eTempType == LEFT_NEIGHBOR)
  {
    for (int y = 0; y < iHeight; y++, pDst += iDstStride)
    {
      for (int x = 0; x < iTemplateWidth; x++)
      {
        pDst[x] = dcval;
      }
    }
  }
  else if(eTempType == ABOVE_NEIGHBOR)
  {
    for (int y = 0; y < iTemplateHeight; y++, pDst+=iDstStride)
    {
      for (int x = 0; x < iWidth; x++)
      {
        pDst[x] = dcval;
      }
    }
  }
  else
  {
    assert(0);
  }
}

void IntraPrediction::initPredTimdIntraParams(const PredictionUnit & pu, const CompArea area, int dirMode)
{
  const Size   puSize    = Size( area.width, area.height );
  const Size&  blockSize = puSize;
  const int     predMode = getWideAngleExt( blockSize.width, blockSize.height, dirMode );

  m_ipaParam.isModeVer            = predMode >= EXT_DIA_IDX;
  m_ipaParam.refFilterFlag        = false;
  m_ipaParam.interpolationFlag    = false;
  m_ipaParam.applyPDPC            = puSize.width >= MIN_TB_SIZEY && puSize.height >= MIN_TB_SIZEY;
  const int    intraPredAngleMode = (m_ipaParam.isModeVer) ? predMode - EXT_VER_IDX : -(predMode - EXT_HOR_IDX);

  int absAng = 0;
  static const int extAngTable[64]    = { 0, 1, 2, 3, 4, 5, 6,7, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 43, 46, 49, 52, 55, 58, 61, 64, 67, 70, 74, 78, 84, 90, 96, 102, 108, 114, 121, 128, 137, 146, 159, 172, 188, 204, 230, 256, 299, 342, 427, 512, 597, 682, 853, 1024, 1536, 2048, 3072 };
  static const int extInvAngTable[64] = { 0, 32768, 16384, 10923, 8192, 6554, 5461, 4681, 4096, 3277, 2731, 2341, 2048, 1820, 1638, 1489, 1365, 1260, 1170, 1092, 1024, 964, 910, 862, 819, 762, 712, 669, 630, 596, 565, 537, 512, 489, 468, 443, 420, 390, 364, 341, 321, 303, 287, 271, 256, 239, 224, 206, 191, 174, 161, 142, 128, 110, 96, 77, 64, 55, 48, 38, 32, 21, 16, 11 }; // (512 * 64) / Angle

  const int     absAngMode         = abs(intraPredAngleMode);
  const int     signAng            = intraPredAngleMode < 0 ? -1 : 1;
                absAng             = extAngTable  [absAngMode];

  m_ipaParam.absInvAngle              = extInvAngTable[absAngMode];
  m_ipaParam.intraPredAngle        = signAng * absAng;

  if (dirMode > 1)
  {
    if (intraPredAngleMode < 0)
    {
      m_ipaParam.applyPDPC = false;
    }
    else if (intraPredAngleMode > 0)
    {
      const int sideSize = m_ipaParam.isModeVer ? puSize.height : puSize.width;
      const int maxScale = 2;
#if GRAD_PDPC
      m_ipaParam.useGradPDPC = false;
#endif
      m_ipaParam.angularScale = std::min(maxScale, floorLog2(sideSize) - (floorLog2(3 * m_ipaParam.absInvAngle - 2) - 8));
#if GRAD_PDPC
      if (m_ipaParam.angularScale < 0)
      {
        m_ipaParam.angularScale = (floorLog2(puSize.width) + floorLog2(puSize.height) - 2) >> 2;
        m_ipaParam.useGradPDPC = true;
      }
#endif
      m_ipaParam.applyPDPC &= m_ipaParam.angularScale >= 0;
    }
  }
}

void IntraPrediction::xPredTimdIntraAng( const CPelBuf &pSrc, const ClpRng& clpRng, Pel* pTrueDst, int iDstStride, int iWidth, int iHeight, TEMPLATE_TYPE eTempType, int iTemplateWidth , int iTemplateHeight, uint32_t dirMode)
{
  int width = iWidth;
  int height = iHeight;
  const bool bIsModeVer     = m_ipaParam.isModeVer;
  const int  intraPredAngle = m_ipaParam.intraPredAngle;
  const int  invAngle       = m_ipaParam.absInvAngle;
  Pel* refMain;
  Pel* refSide;
  static Pel  refAbove[2 * MAX_CU_SIZE + 5 + 33 * MAX_REF_LINE_IDX];
  static Pel  refLeft[2 * MAX_CU_SIZE + 5 + 33 * MAX_REF_LINE_IDX];

  // Initialize the Main and Left reference array.
  if (intraPredAngle < 0)
  {
    for (int x = 0; x <= width + 1; x++)
    {
      refAbove[x + height] = pSrc.at(x, 0);
    }
    for (int y = 0; y <= height + 1; y++)
    {
      refLeft[y + width] = pSrc.at(y, 1);
    }
    refMain = bIsModeVer ? refAbove + height : refLeft + width;
    refSide = bIsModeVer ? refLeft + width : refAbove + height;
    // Extend the Main reference to the left.
    int sizeSide = bIsModeVer ? height : width;
    for (int k = -sizeSide; k <= -1; k++)
    {
      refMain[k] = refSide[std::min((-k * invAngle + 256) >> 9, sizeSide)];
    }
  }
  else
  {
    for (int x = 0; x <= m_topRefLength; x++)
    {
      refAbove[x] = pSrc.at(x, 0);
    }
    for (int y = 0; y <= m_leftRefLength; y++)
    {
      refLeft[y] = pSrc.at(y, 1);
    }
    refMain = bIsModeVer ? refAbove : refLeft;
    refSide = bIsModeVer ? refLeft : refAbove;
    // Extend main reference to right using replication
    const int log2Ratio = floorLog2(width - iTemplateWidth) - floorLog2(height - iTemplateHeight);
    const int s         = std::max<int>(0, bIsModeVer ? log2Ratio : -log2Ratio);
    const int maxIndex  = (std::max(iTemplateWidth, iTemplateHeight) << s) + 2 + std::max(iTemplateWidth, iTemplateHeight);
    const int refLength = bIsModeVer ? m_topRefLength : m_leftRefLength;
    const Pel val       = refMain[refLength];
    for (int z = 1; z <= maxIndex; z++)
    {
      refMain[refLength + z] = val;
    }
  }

  // swap width/height if we are doing a horizontal mode:
  static Pel tempArray[(MAX_CU_SIZE+DIMD_MAX_TEMP_SIZE)*(MAX_CU_SIZE+DIMD_MAX_TEMP_SIZE)];  ///< buffer size may not be big enough
  const int dstStride = bIsModeVer ? iDstStride : (MAX_CU_SIZE+DIMD_MAX_TEMP_SIZE);
  Pel *pDst = bIsModeVer ? pTrueDst : tempArray;
  if (!bIsModeVer)
  {
    std::swap(width, height);
    std::swap(iTemplateWidth, iTemplateHeight);
  }

  if( intraPredAngle == 0 )  // pure vertical or pure horizontal
  {
    if(eTempType == LEFT_ABOVE_NEIGHBOR)
    {
      if (m_ipaParam.applyPDPC)
      {
        int scale = (floorLog2(width) + floorLog2(height) - 2) >> 2;
        xIntraPredTimdHorVerPdpc(pDst, dstStride, refSide, width, iTemplateHeight, iTemplateWidth, 0, scale, refMain, clpRng);
        xIntraPredTimdHorVerPdpc(pDst+iTemplateHeight*iDstStride, dstStride, refSide, iTemplateWidth, height, 0, iTemplateHeight, scale, refMain, clpRng);
      }
      else
      {
        for (int y = 0; y < iTemplateHeight; y++)
        {
          memcpy(pDst + y * dstStride + iTemplateWidth, &refMain[iTemplateWidth + 1], (width - iTemplateWidth) * sizeof(Pel));
        }
        for (int y = iTemplateHeight; y < height; y++)
        {
          memcpy(pDst + y * dstStride, &refMain[1], iTemplateWidth * sizeof(Pel));
        }
      }
    }
    else if(eTempType == LEFT_NEIGHBOR || eTempType == ABOVE_NEIGHBOR)
    {
      if((eTempType == LEFT_NEIGHBOR && bIsModeVer)||(eTempType == ABOVE_NEIGHBOR && !bIsModeVer))
      {
        if (m_ipaParam.applyPDPC)
        {
          const int scale   = (floorLog2(width) + floorLog2(height) - 2) >> 2;
          xIntraPredTimdHorVerPdpc(pDst, dstStride, refSide, iTemplateWidth, height, 0, 0, scale, refMain, clpRng);
        }
        else
        {
          for (int y = 0; y < height; y++)
          {
            for (int x = 0; x < iTemplateWidth; x++)
            {
              pDst[y * dstStride+x] = refMain[x + 1];
            }
          }
        }
      }
      else
      {
        if (m_ipaParam.applyPDPC)
        {
          const int scale   = (floorLog2(width) + floorLog2(height) - 2) >> 2;
          xIntraPredTimdHorVerPdpc(pDst, dstStride, refSide, width, iTemplateHeight, 0, 0, scale, refMain, clpRng);
        }
        else
        {
          for (int y = 0; y < iTemplateHeight; y++)
          {
            memcpy(pDst + y * dstStride, &refMain[1], width * sizeof(Pel));
          }
        }
      }
    }
    else
    {
      assert(0);
    }
  }
  else
  {
    Pel *pDsty=pDst;
    if ( !isIntegerSlopeExt( abs(intraPredAngle) ) )
    {
      int deltaPos = intraPredAngle;
      if (eTempType == LEFT_ABOVE_NEIGHBOR)
      {
        Pel *pDsty=pDst;
        // Above template
        xIntraPredTimdAngLuma(pDsty, dstStride, refMain, width, iTemplateHeight, deltaPos, intraPredAngle, clpRng, iTemplateWidth, 0);
        // Left template
        for (int y = 0; y < iTemplateHeight; y++)
          deltaPos += intraPredAngle;
        xIntraPredTimdAngLuma(pDsty, dstStride, refMain, iTemplateWidth, height, deltaPos, intraPredAngle, clpRng, 0, iTemplateHeight);
#if GRAD_PDPC
        if (m_ipaParam.applyPDPC && m_ipaParam.useGradPDPC)
        {
          int deltaPos2 = intraPredAngle;
          const int scale = m_ipaParam.angularScale;
          xIntraPredTimdAngGradPdpc(pDst, dstStride, refMain, refSide, width, iTemplateHeight, iTemplateWidth, 0, scale, deltaPos2, intraPredAngle, clpRng);
          for (int y = 0; y < iTemplateHeight; y++)
            deltaPos2 += intraPredAngle;
          xIntraPredTimdAngGradPdpc(pDst+iTemplateHeight*dstStride, dstStride, refMain, refSide, iTemplateWidth, height, 0, iTemplateHeight, scale, deltaPos2, intraPredAngle, clpRng);
        }
        else
#endif
        if (m_ipaParam.applyPDPC)
        {
          const int scale = m_ipaParam.angularScale;
          xIntraPredTimdAngPdpc(pDst, dstStride, refSide, width, iTemplateHeight, iTemplateWidth, 0, scale, invAngle);
          xIntraPredTimdAngPdpc(pDst+iTemplateHeight*dstStride, dstStride, refSide, iTemplateWidth, height, 0, iTemplateHeight, scale, invAngle);
        }
      }
      else if (eTempType == LEFT_NEIGHBOR || eTempType == ABOVE_NEIGHBOR)
      {
        int iRegionWidth, iRegionHeight;
        if((eTempType == LEFT_NEIGHBOR && bIsModeVer)||(eTempType == ABOVE_NEIGHBOR && !bIsModeVer))
        {
          iRegionWidth  = iTemplateWidth;
          iRegionHeight = height;
        }
        else
        {
          iRegionWidth  = width;
          iRegionHeight = iTemplateHeight;
        }
        xIntraPredTimdAngLuma(pDsty, dstStride, refMain, iRegionWidth, iRegionHeight, deltaPos, intraPredAngle, clpRng, 0, 0);
#if GRAD_PDPC
        if (m_ipaParam.applyPDPC && m_ipaParam.useGradPDPC)
        {
          int deltaPos2 = intraPredAngle;
          const int scale = m_ipaParam.angularScale;
          xIntraPredTimdAngGradPdpc(pDst, dstStride, refMain, refSide, iRegionWidth, iRegionHeight, 0, 0, scale, deltaPos2, intraPredAngle, clpRng);
        }
        else
#endif
        if (m_ipaParam.applyPDPC)
        {
          const int scale = m_ipaParam.angularScale;
          xIntraPredTimdAngPdpc(pDst, dstStride, refSide, iRegionWidth, iRegionHeight, 0, 0, scale, invAngle);
        }
      }
    }
    else
    {
      if(eTempType == LEFT_ABOVE_NEIGHBOR)
      {
        Pel *pDsty=pDst;
        for (int y = 0, deltaPos = intraPredAngle; y<height; y++, deltaPos += intraPredAngle, pDsty += dstStride)
        {
          const int deltaInt   = deltaPos >> 6;
          int iStartIdx, iEndIdx;
          if(y < iTemplateHeight)
          {
            iStartIdx = iTemplateWidth;
            iEndIdx   = width - 1;
          }
          else
          {
            iStartIdx = 0;
            iEndIdx   = iTemplateWidth - 1;
          }
          memcpy(pDsty + iStartIdx, &refMain[iStartIdx + deltaInt + 1], (iEndIdx - iStartIdx + 1) * sizeof(Pel));
        }
#if GRAD_PDPC
        if (m_ipaParam.applyPDPC && m_ipaParam.useGradPDPC)
        {
          int deltaPos2 = intraPredAngle;
          const int scale = m_ipaParam.angularScale;
          xIntraPredTimdAngGradPdpc(pDst, dstStride, refMain, refSide, width, iTemplateHeight, iTemplateWidth, 0, scale, deltaPos2, intraPredAngle, clpRng);
          for (int y = 0; y < iTemplateHeight; y++)
            deltaPos2 += intraPredAngle;
          xIntraPredTimdAngGradPdpc(pDst+iTemplateHeight*dstStride, dstStride, refMain, refSide, iTemplateWidth, height, 0, iTemplateHeight, scale, deltaPos2, intraPredAngle, clpRng);
        }
        else
#endif
        if (m_ipaParam.applyPDPC)
        {
          const int scale = m_ipaParam.angularScale;
          xIntraPredTimdAngPdpc(pDst, dstStride, refSide, width, iTemplateHeight, iTemplateWidth, 0, scale, invAngle);
          xIntraPredTimdAngPdpc(pDst+iTemplateHeight*dstStride, dstStride, refSide, iTemplateWidth, height, 0, iTemplateHeight, scale, invAngle);
        }
      }
      else // if (eTempType == LEFT_NEIGHBOR || eTempType == ABOVE_NEIGHBOR)
      {
        Pel *pDsty=pDst;
        assert(eTempType == LEFT_NEIGHBOR || eTempType == ABOVE_NEIGHBOR);
        int iRegionWidth, iRegionHeight;
        if((eTempType == LEFT_NEIGHBOR && bIsModeVer)||(eTempType == ABOVE_NEIGHBOR && !bIsModeVer))
        {
          iRegionWidth  = iTemplateWidth;
          iRegionHeight = height;
        }
        else
        {
          iRegionWidth  = width;
          iRegionHeight = iTemplateHeight;
        }
        for (int y = 0, deltaPos = intraPredAngle; y<iRegionHeight; y++, deltaPos += intraPredAngle, pDsty += dstStride)
        {
          const int deltaInt   = deltaPos >> 6;
          memcpy(pDsty, &refMain[deltaInt + 1], iRegionWidth * sizeof(Pel));
        }
#if GRAD_PDPC
        if (m_ipaParam.applyPDPC && m_ipaParam.useGradPDPC)
        {
          int deltaPos2 = intraPredAngle;
          const int scale = m_ipaParam.angularScale;
          xIntraPredTimdAngGradPdpc(pDst, dstStride, refMain, refSide, iRegionWidth, iRegionHeight, 0, 0, scale, deltaPos2, intraPredAngle, clpRng);
        }
        else
#endif
        if (m_ipaParam.applyPDPC)
        {
          const int scale   = m_ipaParam.angularScale;
          xIntraPredTimdAngPdpc(pDst, dstStride, refSide, iRegionWidth, iRegionHeight, 0, 0, scale, invAngle);
        }
      }
    }
  }

  // Flip the block if this is the horizontal mode
  if (!bIsModeVer)
  {
    if(eTempType == LEFT_ABOVE_NEIGHBOR)
    {
      for (int y = 0; y < height; y++)
      {
        int iStartIdx, iEndIdx;
        if(y < iTemplateHeight)
        {
          iStartIdx = iTemplateWidth;
          iEndIdx   = width - 1;
        }
        else
        {
          iStartIdx = 0;
          iEndIdx   = iTemplateWidth - 1;
        }
        for (int x = iStartIdx; x <= iEndIdx; x++)
        {
          pTrueDst[x*iDstStride+y] = pDst[y*dstStride+x];
        }
      }
    }
    else if(eTempType == LEFT_NEIGHBOR)
    {
      for (int y = 0; y < iTemplateHeight; y++)
      {
        for (int x = 0; x < width; x++)
        {
          pTrueDst[x*iDstStride+y] = pDst[y*dstStride+x];
        }
      }
    }
    else if(eTempType == ABOVE_NEIGHBOR)
    {
      for (int y = 0; y < height; y++)
      {
        for (int x = 0; x < iTemplateWidth; x++)
        {
          pTrueDst[x*iDstStride+y] = pDst[y*dstStride+x];
        }
      }
    }
    else
    {
      assert(0);
    }
  }
}

void IntraPrediction::initTimdIntraPatternLuma(const CodingUnit &cu, const CompArea &area, int iTemplateWidth, int iTemplateHeight, uint32_t uiRefWidth, uint32_t uiRefHeight)
{
  const CodingStructure& cs   = *cu.cs;
  Pel *refBufUnfiltered = m_refBuffer[area.compID][PRED_BUF_UNFILTERED];
  bool bLeftAbove = iTemplateHeight > 0 && iTemplateWidth > 0;
  m_leftRefLength     = bLeftAbove ? (uiRefHeight << 1) : ((uiRefHeight + iTemplateHeight) << 1);
  m_topRefLength      = bLeftAbove ? (uiRefWidth << 1) : ((uiRefWidth + iTemplateWidth) << 1);
  xFillTimdReferenceSamples(cs.picture->getRecoBuf(area), refBufUnfiltered, area, cu, iTemplateWidth, iTemplateHeight);
}

void IntraPrediction::xFillTimdReferenceSamples(const CPelBuf &recoBuf, Pel* refBufUnfiltered, const CompArea &area, const CodingUnit &cu, int iTemplateWidth, int iTemplateHeight)
{
  const ChannelType      chType = toChannelType( area.compID );
  const CodingStructure &cs     = *cu.cs;
  const SPS             &sps    = *cs.sps;
  const PreCalcValues   &pcv    = *cs.pcv;

  const int  tuWidth            = area.width;
  const int  tuHeight           = area.height;
  const int  predSize           = m_topRefLength;
  const int  predHSize          = m_leftRefLength;
  const int predStride = predSize + 1;
  m_refBufferStride[area.compID] = predStride;

  const bool noShift            = pcv.noChroma2x2 && area.width == 4; // don't shift on the lowest level (chroma not-split)
  const int  unitWidth          = tuWidth  <= 2 && cu.ispMode && isLuma(area.compID) ? tuWidth  : pcv.minCUWidth  >> (noShift ? 0 : getComponentScaleX(area.compID, sps.getChromaFormatIdc()));
  const int  unitHeight         = tuHeight <= 2 && cu.ispMode && isLuma(area.compID) ? tuHeight : pcv.minCUHeight >> (noShift ? 0 : getComponentScaleY(area.compID, sps.getChromaFormatIdc()));
  int leftTempUnitNum = 0;
  int aboveTempUnitNum = 0;
  if (iTemplateHeight >= 4)
  {
    leftTempUnitNum = iTemplateHeight / unitHeight;
  }
  if (iTemplateWidth >= 4)
  {
    aboveTempUnitNum = iTemplateWidth / unitWidth;
  }

  const int  totalAboveUnits = (predSize + (unitWidth - 1)) / unitWidth - aboveTempUnitNum;
  const int  totalLeftUnits = (predHSize + (unitHeight - 1)) / unitHeight - leftTempUnitNum;
  const int  totalUnits = totalAboveUnits + totalLeftUnits + 1 + aboveTempUnitNum + leftTempUnitNum; //+1 for top-left
  const int  numAboveUnits      = std::max<int>( tuWidth / unitWidth, 1 );
  const int  numLeftUnits       = std::max<int>( tuHeight / unitHeight, 1 );
  const int  numAboveRightUnits = totalAboveUnits - numAboveUnits;
  const int  numLeftBelowUnits  = totalLeftUnits - numLeftUnits;

  // ----- Step 1: analyze neighborhood -----
  const Position posLT          = area;
  const Position posRT          = area.topRight();
  const Position posLB          = area.bottomLeft();

  bool  neighborFlags[4 * MAX_NUM_PART_IDXS_IN_CTU_WIDTH + 1];
  int   numIntraNeighbor = 0;

  memset( neighborFlags, 0, totalUnits );

  neighborFlags[totalLeftUnits] = isAboveLeftAvailable( cu, chType, posLT.offset(-iTemplateWidth, -iTemplateHeight) );
  neighborFlags[totalLeftUnits + leftTempUnitNum] = neighborFlags[totalLeftUnits];
  neighborFlags[totalLeftUnits + leftTempUnitNum + aboveTempUnitNum] = neighborFlags[totalLeftUnits];
  numIntraNeighbor += neighborFlags[totalLeftUnits] ? 1 : 0;
  numIntraNeighbor += leftTempUnitNum > 0 && neighborFlags[totalLeftUnits] ? 1 : 0;
  numIntraNeighbor += aboveTempUnitNum > 0 && neighborFlags[totalLeftUnits] ? 1 : 0;
  numIntraNeighbor += isAboveAvailable     ( cu, chType, posLT.offset(0, -iTemplateHeight), numAboveUnits,      unitWidth,  (neighborFlags + totalLeftUnits + 1 + leftTempUnitNum + aboveTempUnitNum) );
  numIntraNeighbor += isAboveRightAvailable( cu, chType, posRT.offset(0, -iTemplateHeight), numAboveRightUnits, unitWidth,  (neighborFlags + totalLeftUnits + 1 + leftTempUnitNum + aboveTempUnitNum + numAboveUnits) );
  numIntraNeighbor += isLeftAvailable      ( cu, chType, posLT.offset(-iTemplateWidth, 0), numLeftUnits,       unitHeight, (neighborFlags + totalLeftUnits - 1) );
  numIntraNeighbor += isBelowLeftAvailable ( cu, chType, posLB.offset(-iTemplateWidth, 0), numLeftBelowUnits,  unitHeight, (neighborFlags + totalLeftUnits - 1 - numLeftUnits) );

  // ----- Step 2: fill reference samples (depending on neighborhood) -----

  const Pel*  srcBuf    = recoBuf.buf;
  const int   srcStride = recoBuf.stride;
        Pel*  ptrDst    = refBufUnfiltered;
  const Pel*  ptrSrc;
  const Pel   valueDC   = 1 << (sps.getBitDepth( chType ) - 1);


  if( numIntraNeighbor == 0 )
  {
    // Fill border with DC value
    for (int j = 0; j <= predSize; j++) { ptrDst[j] = valueDC; }
    for (int i = 0; i <= predHSize; i++)
    {
      ptrDst[i + predStride] = valueDC;
    }
  }
  else if( numIntraNeighbor == totalUnits )
  {
    // Fill top-left border and top and top right with rec. samples
    ptrSrc = srcBuf - (1 + iTemplateHeight) * srcStride - (1 + iTemplateWidth);
    for (int j = 0; j <= predSize; j++)
    {
      ptrDst[j] = ptrSrc[j];
    }
    for (int i = 0; i <= predHSize; i++)
    {
      ptrDst[i + predStride] = ptrSrc[i * srcStride];
    }
  }
  else // reference samples are partially available
  {
    // Fill top-left sample(s) if available
    ptrSrc = srcBuf - (1 + iTemplateHeight) * srcStride - (1 + iTemplateWidth);
    ptrDst = refBufUnfiltered;
    if (neighborFlags[totalLeftUnits])
    {
      for (int i = 0; i <= iTemplateWidth; i++)
        ptrDst[i] = ptrSrc[i];
      for (int i = 0; i <= iTemplateHeight; i++)
        ptrDst[i + predStride] = ptrSrc[i * srcStride];
    }

    // Fill left & below-left samples if available (downwards)
    ptrSrc += (1 + iTemplateHeight) * srcStride;
    ptrDst += (1 + iTemplateHeight) + predStride;
    for (int unitIdx = totalLeftUnits - 1; unitIdx > 0; unitIdx--)
    {
      if (neighborFlags[unitIdx])
      {
        for (int i = 0; i < unitHeight; i++)
        {
          ptrDst[i] = ptrSrc[i * srcStride];
        }
      }
      ptrSrc += unitHeight * srcStride;
      ptrDst += unitHeight;
    }
    // Fill last below-left sample(s)
    if (neighborFlags[0])
    {
      int lastSample = ((predHSize - iTemplateHeight) % unitHeight == 0) ? unitHeight : (predHSize - iTemplateHeight) % unitHeight;
      for (int i = 0; i < lastSample; i++)
      {
        ptrDst[i] = ptrSrc[i * srcStride];
      }
    }

    // Fill above & above-right samples if available (left-to-right)
    ptrSrc = srcBuf - srcStride * (1 + iTemplateHeight);
    ptrDst = refBufUnfiltered + 1 + iTemplateWidth;
    for (int unitIdx = totalLeftUnits + 1 + leftTempUnitNum + aboveTempUnitNum; unitIdx < totalUnits - 1; unitIdx++)
    {
      if (neighborFlags[unitIdx])
      {
        for (int j = 0; j < unitWidth; j++)
        {
          ptrDst[j] = ptrSrc[j];
        }
      }
      ptrSrc += unitWidth;
      ptrDst += unitWidth;
    }
    // Fill last above-right sample(s)
    if (neighborFlags[totalUnits - 1])
    {
      int lastSample = ((predSize - iTemplateWidth) % unitWidth == 0) ? unitWidth : (predSize - iTemplateWidth) % unitWidth;
      for (int j = 0; j < lastSample; j++)
      {
        ptrDst[j] = ptrSrc[j];
      }
    }

    // pad from first available down to the last below-left
    ptrDst = refBufUnfiltered;
    int lastAvailUnit = 0;
    if (!neighborFlags[0])
    {
      int firstAvailUnit = 1;
      while (firstAvailUnit < totalUnits && !neighborFlags[firstAvailUnit])
      {
        firstAvailUnit++;
      }

      // first available sample
      int firstAvailRow = -1;
      int firstAvailCol = 0;
      if (firstAvailUnit < totalLeftUnits)
      {
        firstAvailRow = (totalLeftUnits - firstAvailUnit) * unitHeight + iTemplateHeight;
      }
      else if (firstAvailUnit == totalLeftUnits)
      {
        firstAvailRow = iTemplateHeight;
      }
      else
      {
        firstAvailCol = (firstAvailUnit - (totalLeftUnits + leftTempUnitNum + aboveTempUnitNum) - 1) * unitWidth + 1 + iTemplateWidth;
      }
      const Pel firstAvailSample = ptrDst[firstAvailRow < 0 ? firstAvailCol : firstAvailRow + predStride];

      // last sample below-left (n.a.)
      int lastRow = predHSize;

      // fill left column
      for (int i = lastRow; i > firstAvailRow; i--)
      {
        ptrDst[i + predStride] = firstAvailSample;
      }
      // fill top row
      if (firstAvailCol > 0)
      {
        for (int j = 0; j < firstAvailCol; j++)
        {
          ptrDst[j] = firstAvailSample;
        }
      }
      lastAvailUnit = firstAvailUnit;
    }

    // pad all other reference samples.
    int currUnit = lastAvailUnit + 1;
    while (currUnit < totalUnits)
    {
      if (!neighborFlags[currUnit]) // samples not available
      {
        // last available sample
        int lastAvailRow = -1;
        int lastAvailCol = 0;
        if (lastAvailUnit < totalLeftUnits)
        {
          lastAvailRow = (totalLeftUnits - lastAvailUnit - 1) * unitHeight + iTemplateHeight + 1;
        }
        else if (lastAvailUnit == totalLeftUnits)
        {
          lastAvailCol = iTemplateWidth;
        }
        else
        {
          lastAvailCol = (lastAvailUnit - (totalLeftUnits + leftTempUnitNum + aboveTempUnitNum)) * unitWidth + iTemplateWidth;
        }
        const Pel lastAvailSample = ptrDst[lastAvailRow < 0 ? lastAvailCol : lastAvailRow + predStride];

        // fill current unit with last available sample
        if (currUnit < totalLeftUnits)
        {
          for (int i = lastAvailRow - 1; i >= lastAvailRow - unitHeight; i--)
          {
            ptrDst[i + predStride] = lastAvailSample;
          }
        }
        else if (currUnit == totalLeftUnits)
        {
          for (int i = 0; i < iTemplateHeight + 1; i++)
          {
            ptrDst[i + predStride] = lastAvailSample;
          }
          for (int j = 0; j < iTemplateWidth + 1; j++)
          {
            ptrDst[j] = lastAvailSample;
          }
        }
        else
        {
          int numSamplesInUnit = (currUnit == totalUnits - 1) ? (((predSize - iTemplateWidth) % unitWidth == 0) ? unitWidth : (predSize - iTemplateWidth) % unitWidth) : unitWidth;
          for (int j = lastAvailCol + 1; j <= lastAvailCol + numSamplesInUnit; j++)
          {
            ptrDst[j] = lastAvailSample;
          }
        }
      }
      lastAvailUnit = currUnit;
      currUnit++;
    }
  }
}

int IntraPrediction::deriveTimdMode( const CPelBuf &recoBuf, const CompArea &area, CodingUnit &cu )
{
  int channelBitDepth = cu.slice->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA);
  SizeType uiWidth = cu.lwidth();
  SizeType uiHeight = cu.lheight();

  static Pel PredLuma[(MAX_CU_SIZE + DIMD_MAX_TEMP_SIZE) * (MAX_CU_SIZE + DIMD_MAX_TEMP_SIZE)];
  memset(PredLuma, 0, (MAX_CU_SIZE + DIMD_MAX_TEMP_SIZE) * (MAX_CU_SIZE + DIMD_MAX_TEMP_SIZE) * sizeof(Pel));
  Pel* piPred = PredLuma;
  uint32_t uiPredStride = MAX_CU_SIZE + DIMD_MAX_TEMP_SIZE;

  int  iCurX  = cu.lx();
  int  iCurY  = cu.ly();
  int  iRefX  = -1, iRefY = -1;
  uint32_t uiRefWidth = 0, uiRefHeight = 0;

  int iTempWidth = 4, iTempHeight = 4;
  if(uiWidth <= 8)
  {
    iTempWidth = 2;
  }
  if(uiHeight <= 8)
  {
    iTempHeight = 2;
  }

  TEMPLATE_TYPE eTempType = CU::deriveTimdRefType(iCurX, iCurY, uiWidth, uiHeight, iTempWidth, iTempHeight, iRefX, iRefY, uiRefWidth, uiRefHeight);

  if (eTempType != NO_NEIGHBOR)
  {
    const CodingStructure& cs   = *cu.cs;
    m_ipaParam.multiRefIndex        = iTempWidth;
    Pel* piOrg = cs.picture->getRecoBuf( area ).buf;
    int iOrgStride = cs.picture->getRecoBuf( area ).stride;
    piOrg += (iRefY - iCurY) * iOrgStride + (iRefX - iCurX);
    DistParam distParamSad[2]; // above, left
    distParamSad[0].applyWeight = false;
    distParamSad[0].useMR = false;
    distParamSad[1].applyWeight = false;
    distParamSad[1].useMR = false;
    if(eTempType == LEFT_ABOVE_NEIGHBOR)
    {
      m_timdSatdCost->setTimdDistParam(distParamSad[0], piOrg + iTempWidth, piPred + iTempWidth, iOrgStride, uiPredStride, channelBitDepth, COMPONENT_Y, uiWidth, iTempHeight, 0, 1, true); // Use HAD (SATD) cost
      m_timdSatdCost->setTimdDistParam(distParamSad[1], piOrg + iTempHeight * iOrgStride, piPred + iTempHeight * uiPredStride, iOrgStride, uiPredStride, channelBitDepth, COMPONENT_Y, iTempWidth, uiHeight, 0, 1, true); // Use HAD (SATD) cost
    }
    else if(eTempType == LEFT_NEIGHBOR)
    {
      m_timdSatdCost->setTimdDistParam(distParamSad[1], piOrg, piPred, iOrgStride, uiPredStride, channelBitDepth, COMPONENT_Y, iTempWidth, uiHeight, 0, 1, true);
    }
    else if(eTempType == ABOVE_NEIGHBOR)
    {
      m_timdSatdCost->setTimdDistParam(distParamSad[0], piOrg, piPred, iOrgStride, uiPredStride, channelBitDepth, COMPONENT_Y, uiWidth, iTempHeight, 0, 1, true);
    }
    initTimdIntraPatternLuma(cu, area, eTempType != ABOVE_NEIGHBOR ? iTempWidth : 0, eTempType != LEFT_NEIGHBOR ? iTempHeight : 0, uiRefWidth, uiRefHeight);

    uint32_t uiIntraDirNeighbor[5] = {0}, modeIdx = 0;
    bool includedMode[EXT_VDIA_IDX + 1];
    memset(includedMode, false, (EXT_VDIA_IDX + 1) * sizeof(bool));
    auto &pu = *cu.firstPU;
    uint32_t uiRealW = uiRefWidth + (eTempType == LEFT_NEIGHBOR? iTempWidth : 0);
    uint32_t uiRealH = uiRefHeight + (eTempType == ABOVE_NEIGHBOR? iTempHeight : 0);
    uint64_t maxCost = (uint64_t)(iTempWidth * cu.lheight() + iTempHeight * cu.lwidth());

    uint64_t uiBestCost = MAX_UINT64;
    int iBestMode = PLANAR_IDX;
    uint64_t uiSecondaryCost = MAX_UINT64;
    int iSecondaryMode = PLANAR_IDX;

    const Position posLTx = pu.Y().topLeft();
    const Position posRTx = pu.Y().topRight();
    const Position posLBx = pu.Y().bottomLeft();

    // left
    const PredictionUnit *puLeftx = pu.cs->getPURestricted(posLBx.offset(-1, 0), pu, pu.chType);
    if (puLeftx && CU::isIntra(*puLeftx->cu))
    {
      uiIntraDirNeighbor[modeIdx] = PU::getIntraDirLuma( *puLeftx );
      if (!puLeftx->cu->timd)
      {
        uiIntraDirNeighbor[modeIdx] = MAP67TO131(uiIntraDirNeighbor[modeIdx]);
      }
      if( !includedMode[uiIntraDirNeighbor[modeIdx]] )
      {
        includedMode[uiIntraDirNeighbor[modeIdx]] = true;
        modeIdx++;
      }
    }
    // above
    const PredictionUnit *puAbovex = pu.cs->getPURestricted(posRTx.offset(0, -1), pu, pu.chType);
    if (puAbovex && CU::isIntra(*puAbovex->cu) && CU::isSameCtu(*pu.cu, *puAbovex->cu))
    {
      uiIntraDirNeighbor[modeIdx] =PU::getIntraDirLuma( *puAbovex );
      if (!puAbovex->cu->timd)
      {
        uiIntraDirNeighbor[modeIdx] = MAP67TO131(uiIntraDirNeighbor[modeIdx]);
      }
      if( !includedMode[uiIntraDirNeighbor[modeIdx]] )
      {
        includedMode[uiIntraDirNeighbor[modeIdx]] = true;
        modeIdx++;
      }
    }
    // below left
    const PredictionUnit *puLeftBottomx = cs.getPURestricted( posLBx.offset( -1, 1 ), pu, pu.chType );
    if (puLeftBottomx && CU::isIntra(*puLeftBottomx->cu))
    {
      uiIntraDirNeighbor[modeIdx] = PU::getIntraDirLuma( *puLeftBottomx );
      if (!puLeftBottomx->cu->timd)
      {
        uiIntraDirNeighbor[modeIdx] = MAP67TO131(uiIntraDirNeighbor[modeIdx]);
      }
      if( !includedMode[uiIntraDirNeighbor[modeIdx]] )
      {
        includedMode[uiIntraDirNeighbor[modeIdx]] = true;
        modeIdx++;
      }
    }
    // above right
    const PredictionUnit *puAboveRightx = cs.getPURestricted( posRTx.offset( 1, -1 ), pu, pu.chType );
    if (puAboveRightx && CU::isIntra(*puAboveRightx->cu))
    {
      uiIntraDirNeighbor[modeIdx] = PU::getIntraDirLuma( *puAboveRightx );
      if (!puAboveRightx->cu->timd)
      {
        uiIntraDirNeighbor[modeIdx] = MAP67TO131(uiIntraDirNeighbor[modeIdx]);
      }
      if( !includedMode[uiIntraDirNeighbor[modeIdx]] )
      {
        includedMode[uiIntraDirNeighbor[modeIdx]] = true;
        modeIdx++;
      }
    }
    //above left
    const PredictionUnit *puAboveLeftx = cs.getPURestricted( posLTx.offset( -1, -1 ), pu, pu.chType );
    if( puAboveLeftx && CU::isIntra(*puAboveLeftx->cu) )
    {
      uiIntraDirNeighbor[modeIdx] = PU::getIntraDirLuma( *puAboveLeftx );
      if (!puAboveLeftx->cu->timd)
      {
        uiIntraDirNeighbor[modeIdx] = MAP67TO131(uiIntraDirNeighbor[modeIdx]);
      }
      if( !includedMode[uiIntraDirNeighbor[modeIdx]] )
      {
        includedMode[uiIntraDirNeighbor[modeIdx]] = true;
        modeIdx++;
      }
    }
    bool bNoAngular = false;
    if(modeIdx >= 2)
    {
      bNoAngular = true;
      for(uint32_t i = 0; i < modeIdx; i++)
      {
        if(uiIntraDirNeighbor[i] > DC_IDX)
        {
          bNoAngular = false;
          break;
        }
      }
    }

    if (bNoAngular)
    {
      for(int iMode = 0; iMode <= 1; iMode ++)
      {
        uint64_t uiCost = 0;
        initPredTimdIntraParams(pu, area, iMode);
        predTimdIntraAng(COMPONENT_Y, pu, iMode, piPred, uiPredStride, uiRealW, uiRealH, eTempType, (eTempType == ABOVE_NEIGHBOR)? 0: iTempWidth, (eTempType == LEFT_NEIGHBOR)? 0: iTempHeight);
        if(eTempType == LEFT_ABOVE_NEIGHBOR)
        {
          uiCost += distParamSad[0].distFunc(distParamSad[0]);
          uiCost += distParamSad[1].distFunc(distParamSad[1]);
        }
        else if(eTempType == LEFT_NEIGHBOR)
        {
          uiCost = distParamSad[1].distFunc(distParamSad[1]);
        }
        else if(eTempType == ABOVE_NEIGHBOR)
        {
          uiCost += distParamSad[0].distFunc(distParamSad[0]);
        }
        else
        {
          assert(0);
        }

        if(uiCost < uiBestCost)
        {
          uiBestCost = uiCost;
          iBestMode = iMode;
        }
        if(uiBestCost <= maxCost)
        {
          break;
        }
      }
      cu.timdMode = iBestMode;
      cu.timdIsBlended = false;

      return iBestMode;
    }
#if SECONDARY_MPM
    uint8_t mpmList[NUM_MOST_PROBABLE_MODES];
    uint8_t intraNonMPM[NUM_NON_MPM_MODES];
    PU::getIntraMPMs(pu, mpmList, intraNonMPM);
#else
    unsigned mpmList[NUM_MOST_PROBABLE_MODES];
    PU::getIntraMPMs(pu, mpmList);
#endif
    unsigned mpmExtraList[NUM_MOST_PROBABLE_MODES + 3]; // +DC/VER/HOR
    int maxModeNum = NUM_MOST_PROBABLE_MODES;
    unsigned modeCandList[3] = {DC_IDX, HOR_IDX, VER_IDX};
    bool bNotExist[3] = {true, true, true};
    for (int i = 0; i < NUM_MOST_PROBABLE_MODES; i++)
    {
      mpmExtraList[i] = mpmList[i];
      if (bNotExist[0] && mpmList[i] == DC_IDX)
      {
        bNotExist[0] = false;
      }
      if (bNotExist[1] && mpmList[i] == HOR_IDX)
      {
        bNotExist[1] = false;
      }
      if (bNotExist[2] && mpmList[i] == VER_IDX)
      {
        bNotExist[2] = false;
      }
    }
    for (int i = 0; i < 3; i++)
    {
      if (bNotExist[i])
      {
        mpmExtraList[maxModeNum++] = modeCandList[i];
      }
    }
    for(int i = 0; i < maxModeNum; i ++)
    {
      uint64_t uiCost = 0;
      int iMode = mpmExtraList[i];
      if (iMode > DC_IDX)
      {
        iMode = MAP67TO131(iMode);
      }
      initPredTimdIntraParams(pu, area, iMode);
      predTimdIntraAng(COMPONENT_Y, pu, iMode, piPred, uiPredStride, uiRealW, uiRealH, eTempType, (eTempType == ABOVE_NEIGHBOR)? 0: iTempWidth, (eTempType == LEFT_NEIGHBOR)? 0: iTempHeight);
      if(eTempType == LEFT_ABOVE_NEIGHBOR)
      {
        uiCost += distParamSad[0].distFunc(distParamSad[0]);
        uiCost += distParamSad[1].distFunc(distParamSad[1]);
      }
      else if(eTempType == LEFT_NEIGHBOR)
      {
        uiCost = distParamSad[1].distFunc(distParamSad[1]);
      }
      else if(eTempType == ABOVE_NEIGHBOR)
      {
        uiCost += distParamSad[0].distFunc(distParamSad[0]);
      }
      else
      {
        assert(0);
      }

      if( uiCost < uiBestCost )
      {
        uiSecondaryCost = uiBestCost;
        iSecondaryMode  = iBestMode;
        uiBestCost  = uiCost;
        iBestMode = iMode;
      }
      else if (uiCost < uiSecondaryCost)
      {
        uiSecondaryCost = uiCost;
        iSecondaryMode  = iMode;
      }
      if (uiSecondaryCost <= maxCost)
      {
        break;
      }
    }

    int midMode = iBestMode;
    if (midMode > DC_IDX && uiBestCost > maxCost)
    {
      for (int i = -1; i <= 1; i+=2)
      {
        int iMode = midMode + i;
        if (iMode <= DC_IDX || iMode > EXT_VDIA_IDX)
        {
          continue;
        }
        initPredTimdIntraParams(pu, area, iMode);
        predTimdIntraAng(COMPONENT_Y, pu, iMode, piPred, uiPredStride, uiRealW, uiRealH, eTempType, (eTempType == ABOVE_NEIGHBOR)? 0: iTempWidth, (eTempType == LEFT_NEIGHBOR)? 0: iTempHeight);
        uint64_t uiCost = 0;
        if(eTempType == LEFT_ABOVE_NEIGHBOR)
        {
          uiCost += distParamSad[0].distFunc(distParamSad[0]);
          uiCost += distParamSad[1].distFunc(distParamSad[1]);
        }
        else if(eTempType == LEFT_NEIGHBOR)
        {
          uiCost = distParamSad[1].distFunc(distParamSad[1]);
        }
        else if(eTempType == ABOVE_NEIGHBOR)
        {
          uiCost += distParamSad[0].distFunc(distParamSad[0]);
        }
        else
        {
          assert(0);
        }

        if(uiCost < uiBestCost)
        {
          uiBestCost  = uiCost;
          iBestMode = iMode;
        }
        if(uiBestCost <= maxCost)
        {
          break;
        }
      }
    }

    midMode = iSecondaryMode;
    if (midMode > DC_IDX && uiSecondaryCost > maxCost)
    {
      for (int i = -1; i <= 1; i+=2)
      {
        int iMode = midMode + i;
        if (iMode <= DC_IDX || iMode > EXT_VDIA_IDX)
        {
          continue;
        }
        initPredTimdIntraParams(pu, area, iMode);
        predTimdIntraAng(COMPONENT_Y, pu, iMode, piPred, uiPredStride, uiRealW, uiRealH, eTempType, (eTempType == ABOVE_NEIGHBOR)? 0: iTempWidth, (eTempType == LEFT_NEIGHBOR)? 0: iTempHeight);
        uint64_t uiCost = 0;
        if(eTempType == LEFT_ABOVE_NEIGHBOR)
        {
          uiCost += distParamSad[0].distFunc(distParamSad[0]);
          uiCost += distParamSad[1].distFunc(distParamSad[1]);
        }
        else if(eTempType == LEFT_NEIGHBOR)
        {
          uiCost = distParamSad[1].distFunc(distParamSad[1]);
        }
        else if(eTempType == ABOVE_NEIGHBOR)
        {
          uiCost += distParamSad[0].distFunc(distParamSad[0]);
        }
        else
        {
          assert(0);
        }

        if(uiCost < uiSecondaryCost)
        {
          uiSecondaryCost  = uiCost;
          iSecondaryMode = iMode;
        }
        if(uiSecondaryCost <= maxCost)
        {
          break;
        }
      }
    }

    if ((uiSecondaryCost - uiBestCost) < uiBestCost)
  {
    cu.timdMode         = iBestMode;
    cu.timdIsBlended    = true;
    cu.timdModeSecondary = iSecondaryMode;

    const int blend_sum_weight = 6;
    int       sum_weight       = 1 << blend_sum_weight;

    double dRatio       = 0.0;
    dRatio              = (double) uiSecondaryCost / (double) (uiBestCost + uiSecondaryCost);
    int iRatio          = static_cast<int>(dRatio * sum_weight + 0.5);
    cu.timdFusionWeight[0] = iRatio;
    cu.timdFusionWeight[1] = sum_weight - iRatio;
  }
  else
  {
    cu.timdMode      = iBestMode;
    cu.timdIsBlended = false;
  }

    return iBestMode;
  }
  else
  {
    cu.timdMode = PLANAR_IDX;
    cu.timdIsBlended = false;

    return PLANAR_IDX;
  }
}
#endif
#if ENABLE_DIMD
void IntraPrediction::deriveDimdMode(const CPelBuf &recoBuf, const CompArea &area, CodingUnit &cu)
{
  if( !cu.slice->getSPS()->getUseDimd() )
  {
    return;
  }

  int sigcnt = 0;
  const CodingStructure  &cs = *cu.cs;
  const SPS             &sps = *cs.sps;
  const PreCalcValues   &pcv = *cs.pcv;
  const ChannelType   chType = toChannelType(area.compID);

  const Pel *pReco = recoBuf.buf;
  const uint32_t uiWidth = area.width;
  const uint32_t uiHeight = area.height;
  const int iStride = recoBuf.stride;
  const int predSize = (uiWidth << 1);
  const int predHSize = (uiHeight << 1);

  const bool noShift = pcv.noChroma2x2 && uiWidth == 4; // don't shift on the lowest level (chroma not-split)
  const int  unitWidth = pcv.minCUWidth >> (noShift ? 0 : getComponentScaleX(area.compID, sps.getChromaFormatIdc()));
  const int  unitHeight = pcv.minCUHeight >> (noShift ? 0 : getComponentScaleY(area.compID, sps.getChromaFormatIdc()));

  const int  totalAboveUnits = (predSize + (unitWidth - 1)) / unitWidth;
  const int  totalLeftUnits = (predHSize + (unitHeight - 1)) / unitHeight;
  const int  totalUnits = totalAboveUnits + totalLeftUnits + 1; //+1 for top-left
  const int  numAboveUnits = std::max<int>(uiWidth / unitWidth, 1);
  const int  numLeftUnits = std::max<int>(uiHeight / unitHeight, 1);
  const int  numAboveRightUnits = totalAboveUnits - numAboveUnits;
  const int  numLeftBelowUnits = totalLeftUnits - numLeftUnits;

  CHECK(numAboveUnits <= 0 || numLeftUnits <= 0 || numAboveRightUnits <= 0 || numLeftBelowUnits <= 0, "Size not supported");

  // ----- Step 1: analyze neighborhood -----
  const Position posLT = area;

  bool  neighborFlags[4 * MAX_NUM_PART_IDXS_IN_CTU_WIDTH + 1];
  memset(neighborFlags, 0, totalUnits);

  int numIntraAbove = isAboveAvailable(cu, chType, posLT, numAboveUnits, unitWidth, (neighborFlags + totalLeftUnits + 1));
  int numIntraLeft = isLeftAvailable(cu, chType, posLT, numLeftUnits, unitHeight, (neighborFlags + totalLeftUnits - 1));

  // ----- Step 2: build histogram of gradients -----
  int piHistogram_clean[NUM_LUMA_MODE] = { 0 };
  
  if (numIntraLeft)
  {
    uint32_t uiHeightLeft = numIntraLeft * unitHeight - 1 - (!numIntraAbove ? 1 : 0);
    const Pel *pRecoLeft = pReco - 2 + iStride * (!numIntraAbove ? 1 : 0);
    sigcnt += buildHistogram(pRecoLeft, iStride, uiHeightLeft, 1, piHistogram_clean, 1, uiWidth, uiHeight);
  }

  if (numIntraAbove)
  {
    uint32_t uiWidthAbove = numIntraAbove * unitWidth - 1 - (!numIntraLeft ? 1 : 0);
    const Pel *pRecoAbove = pReco - iStride * 2 + (!numIntraLeft ? 1 : 0);
    sigcnt += buildHistogram(pRecoAbove, iStride, 1, uiWidthAbove, piHistogram_clean, 2, uiWidth, uiHeight);
  }

  if (numIntraLeft && numIntraAbove)
  {
    const Pel *pRecoAboveLeft = pReco - 2 - iStride * 2;
    sigcnt += buildHistogram(pRecoAboveLeft, iStride, 2, 2, piHistogram_clean, 3, uiWidth, uiHeight);
  }

  int first_amp = 0, second_amp = 0, cur_amp = 0;
  int first_mode = 0, second_mode = 0, cur_mode = 0;
  for (int i = 0; i < NUM_LUMA_MODE; i++)
  {
    cur_amp = piHistogram_clean[i];
    cur_mode = i;
    if (cur_amp > first_amp)
    {
      second_amp = first_amp;
      second_mode = first_mode;
      first_amp = cur_amp;
      first_mode = cur_mode;
    }
    else
    {
      if (cur_amp > second_amp)
      {
        second_amp = cur_amp;
        second_mode = cur_mode;
      }
    }
  }

  // ----- Step 3: derive best mode from histogram of gradients -----
  cu.dimdMode = first_mode;

  cu.dimd_is_blend = true;
  cu.dimd_is_blend &= second_amp > 0;
  cu.dimd_is_blend &= second_mode > DC_IDX;
  cu.dimd_is_blend &= first_mode > DC_IDX;

  if( cu.dimd_is_blend )
  {
    cu.dimdBlendMode[0] = second_mode;
  }

  const int blend_sum_weight = 6;
  int sum_weight = 1 << blend_sum_weight;
  if (cu.dimd_is_blend)
  {
    double dRatio = 0.0;
    sum_weight -= static_cast<int>((double)sum_weight / 3); // ~ 1/3 of the weight to be reserved for planar
    dRatio = (double)first_amp / (double)(first_amp + second_amp);
    int iRatio = static_cast<int>(dRatio * sum_weight);
    cu.dimdRelWeight[0] = iRatio;
    cu.dimdRelWeight[2] = sum_weight - iRatio;
    cu.dimdRelWeight[1] = (1 << blend_sum_weight) - sum_weight;
  }
  else
  {
    cu.dimdRelWeight[0] = sum_weight;
    cu.dimdRelWeight[1] = 0;
    cu.dimdRelWeight[2] = 0;
  }

}

int buildHistogram(const Pel *pReco, int iStride, uint32_t uiHeight, uint32_t uiWidth, int* piHistogram, int direction, int bw, int bh)
{
  int w_step = 1, h_step = 1;
  int angTable[17] = { 0, 2048, 4096, 6144, 8192, 12288, 16384, 20480, 24576, 28672, 32768, 36864, 40960, 47104, 53248, 59392, 65536 };
  int offsets[4] = { HOR_IDX, HOR_IDX, VER_IDX, VER_IDX };
  int dirs[4] = { -1, 1, -1, 1 };
  int map_x_gr_y_1[2][2] = { { 1, 0 },{ 0, 1 } };
  int map_x_gr_y_0[2][2] = { { 2, 3 },{ 3, 2 } };

  for (uint32_t y = 0; y < uiHeight; y += h_step)
  {
    for (uint32_t x = 0; x < uiWidth; x += w_step)
    {
      if ((direction == 3) && x == (uiWidth - 1) && y == (uiHeight - 1))
        continue;

      const Pel *pRec = pReco + y * iStride + x;

      int iDy = pRec[-iStride - 1] + 2 * pRec[-1] + pRec[iStride - 1] - pRec[-iStride + 1] - 2 * pRec[+1] - pRec[iStride + 1];
      int iDx = pRec[iStride - 1] + 2 * pRec[iStride] + pRec[iStride + 1] - pRec[-iStride - 1] - 2 * pRec[-iStride] - pRec[-iStride + 1];

      if (iDy == 0 && iDx == 0)
        continue;

      int iAmp = (int)(abs(iDx) + abs(iDy));
      int iAng_uneven = -1;
      // for determining region
      if (iDx != 0 && iDy != 0) // pure angles are not concerned
      {
        // get the region
        int signx = iDx < 0 ? 1 : 0;
        int signy = iDy < 0 ? 1 : 0;
        int absx = iDx < 0 ? -iDx : iDx;
        int absy = iDy < 0 ? -iDy : iDy;
        int x_gr_y = absx > absy ? 1 : 0;
        int region = x_gr_y ? map_x_gr_y_1[signy][signx] : map_x_gr_y_0[signy][signx];
        //region = (region == 1 ? 2 : (region == 2 ? 1 : (region == 3 ? 4 : 3)));
        float fRatio = x_gr_y ? static_cast<float>(absy) / static_cast<float>(absx) : static_cast<float>(absx) / static_cast<float>(absy);
        float fRatio_scaled = fRatio * (1 << 16);
        int iRatio = static_cast<int>(fRatio_scaled);
        // get ang_idx
        int idx = -1;
        for( int i = 0; i < 17; i++ )
        {
          if( iRatio < angTable[i] )
          {
            idx = iRatio - angTable[i - 1] < angTable[i] - iRatio ? i - 1 : i;
            break;
          }
        }

        iAng_uneven = offsets[region] + dirs[region] * idx;
        //iAng_uneven = offsets[region - 1] + dirs[region - 1] * idx;
      }
      else
      {
        iAng_uneven = iDx == 0 ? VER_IDX : HOR_IDX;
      }
      piHistogram[iAng_uneven] += iAmp;
    }
  }
  return 0;
}
#endif
bool isAboveLeftAvailable(const CodingUnit &cu, const ChannelType &chType, const Position &posLT)
{
  const CodingStructure& cs = *cu.cs;
  const Position refPos = posLT.offset(-1, -1);

  if (!cs.isDecomp(refPos, chType))
  {
    return false;
  }

  return (cs.getCURestricted(refPos, cu, chType) != NULL);
}

int isAboveAvailable(const CodingUnit &cu, const ChannelType &chType, const Position &posLT, const uint32_t uiNumUnitsInPU, const uint32_t unitWidth, bool *bValidFlags)
{
  const CodingStructure& cs = *cu.cs;

  bool *    validFlags = bValidFlags;
  int       numIntra   = 0;
  const int maxDx      = uiNumUnitsInPU * unitWidth;

  for (int dx = 0; dx < maxDx; dx += unitWidth)
  {
    const Position refPos = posLT.offset(dx, -1);

    if (!cs.isDecomp(refPos, chType))
    {
      break;
    }

    const bool valid = (cs.getCURestricted(refPos, cu, chType) != NULL);
    numIntra += valid ? 1 : 0;
    *validFlags = valid;

    validFlags++;
  }

  return numIntra;
}

int isLeftAvailable(const CodingUnit &cu, const ChannelType &chType, const Position &posLT, const uint32_t uiNumUnitsInPU, const uint32_t unitHeight, bool *bValidFlags)
{
  const CodingStructure& cs = *cu.cs;

  bool *    validFlags = bValidFlags;
  int       numIntra   = 0;
  const int maxDy      = uiNumUnitsInPU * unitHeight;

  for (int dy = 0; dy < maxDy; dy += unitHeight)
  {
    const Position refPos = posLT.offset(-1, dy);

    if (!cs.isDecomp(refPos, chType))
    {
      break;
    }

    const bool valid = (cs.getCURestricted(refPos, cu, chType) != NULL);
    numIntra += valid ? 1 : 0;
    *validFlags = valid;

    validFlags--;
  }

  return numIntra;
}

int isAboveRightAvailable(const CodingUnit &cu, const ChannelType &chType, const Position &posRT, const uint32_t uiNumUnitsInPU, const uint32_t unitWidth, bool *bValidFlags )
{
  const CodingStructure& cs = *cu.cs;

  bool *    validFlags = bValidFlags;
  int       numIntra   = 0;
  const int maxDx      = uiNumUnitsInPU * unitWidth;

  for (int dx = 0; dx < maxDx; dx += unitWidth)
  {
    const Position refPos = posRT.offset(unitWidth + dx, -1);

    if (!cs.isDecomp(refPos, chType))
    {
      break;
    }

    const bool valid = (cs.getCURestricted(refPos, cu, chType) != NULL);
    numIntra += valid ? 1 : 0;
    *validFlags = valid;

    validFlags++;
  }

  return numIntra;
}

int isBelowLeftAvailable(const CodingUnit &cu, const ChannelType &chType, const Position &posLB, const uint32_t uiNumUnitsInPU, const uint32_t unitHeight, bool *bValidFlags )
{
  const CodingStructure& cs = *cu.cs;

  bool *    validFlags = bValidFlags;
  int       numIntra   = 0;
  const int maxDy      = uiNumUnitsInPU * unitHeight;

  for (int dy = 0; dy < maxDy; dy += unitHeight)
  {
    const Position refPos = posLB.offset(-1, unitHeight + dy);

    if (!cs.isDecomp(refPos, chType))
    {
      break;
    }

    const bool valid = (cs.getCURestricted(refPos, cu, chType) != NULL);
    numIntra += valid ? 1 : 0;
    *validFlags = valid;

    validFlags--;
  }

  return numIntra;
}

// LumaRecPixels
void IntraPrediction::xGetLumaRecPixels(const PredictionUnit &pu, CompArea chromaArea)
{
  int iDstStride = 0;
  Pel* pDst0 = 0;
  int curChromaMode = pu.intraDir[1];
#if MMLM
  if ((curChromaMode == MDLM_L_IDX) || (curChromaMode == MDLM_T_IDX) || (curChromaMode == MMLM_L_IDX) || (curChromaMode == MMLM_T_IDX))
#else
  if ((curChromaMode == MDLM_L_IDX) || (curChromaMode == MDLM_T_IDX))
#endif
  {
    iDstStride = 2 * MAX_CU_SIZE + 1;
    pDst0 = m_pMdlmTemp + iDstStride + 1;
  }
  else
  {
    iDstStride = MAX_CU_SIZE + 1;
    pDst0 = m_piTemp + iDstStride + 1; //MMLM_SAMPLE_NEIGHBOR_LINES;
  }
  //assert 420 chroma subsampling
  CompArea lumaArea = CompArea( COMPONENT_Y, pu.chromaFormat, chromaArea.lumaPos(), recalcSize( pu.chromaFormat, CHANNEL_TYPE_CHROMA, CHANNEL_TYPE_LUMA, chromaArea.size() ) );//needed for correct pos/size (4x4 Tus)

  CHECK(lumaArea.width == chromaArea.width && CHROMA_444 != pu.chromaFormat, "");
  CHECK(lumaArea.height == chromaArea.height && CHROMA_444 != pu.chromaFormat && CHROMA_422 != pu.chromaFormat, "");

  const SizeType uiCWidth = chromaArea.width;
  const SizeType uiCHeight = chromaArea.height;

  const CPelBuf Src = pu.cs->picture->getRecoBuf( lumaArea );
  Pel const* pRecSrc0   = Src.bufAt( 0, 0 );
  int iRecStride        = Src.stride;
  int logSubWidthC  = getChannelTypeScaleX(CHANNEL_TYPE_CHROMA, pu.chromaFormat);
  int logSubHeightC = getChannelTypeScaleY(CHANNEL_TYPE_CHROMA, pu.chromaFormat);

  int iRecStride2       = iRecStride << logSubHeightC;

  const CodingUnit& lumaCU = isChroma( pu.chType ) ? *pu.cs->picture->cs->getCU( lumaArea.pos(), CH_L ) : *pu.cu;
  const CodingUnit&     cu = *pu.cu;

  const CompArea& area = isChroma( pu.chType ) ? chromaArea : lumaArea;

  const uint32_t uiTuWidth  = area.width;
  const uint32_t uiTuHeight = area.height;

  int iBaseUnitSize = ( 1 << MIN_CU_LOG2 );

  const int  iUnitWidth       = iBaseUnitSize >> getComponentScaleX( area.compID, area.chromaFormat );
  const int  iUnitHeight = iBaseUnitSize >> getComponentScaleY(area.compID, area.chromaFormat);

  const int  iTUWidthInUnits = uiTuWidth / iUnitWidth;
  const int  iTUHeightInUnits = uiTuHeight / iUnitHeight;
  const int  iAboveUnits      = iTUWidthInUnits;
  const int  iLeftUnits       = iTUHeightInUnits;
  const int  chromaUnitWidth = iBaseUnitSize >> getComponentScaleX(COMPONENT_Cb, area.chromaFormat);
  const int  chromaUnitHeight = iBaseUnitSize >> getComponentScaleY(COMPONENT_Cb, area.chromaFormat);
  const int  topTemplateSampNum = 2 * uiCWidth; // for MDLM, the number of template samples is 2W or 2H.
  const int  leftTemplateSampNum = 2 * uiCHeight;
  assert(m_topRefLength >= topTemplateSampNum);
  assert(m_leftRefLength >= leftTemplateSampNum);
  const int  totalAboveUnits = (topTemplateSampNum + (chromaUnitWidth - 1)) / chromaUnitWidth;
  const int  totalLeftUnits = (leftTemplateSampNum + (chromaUnitHeight - 1)) / chromaUnitHeight;
  const int  totalUnits = totalLeftUnits + totalAboveUnits + 1;
  const int  aboveRightUnits = totalAboveUnits - iAboveUnits;
  const int  leftBelowUnits = totalLeftUnits - iLeftUnits;

  int avaiAboveRightUnits = 0;
  int avaiLeftBelowUnits = 0;
  bool  bNeighborFlags[4 * MAX_NUM_PART_IDXS_IN_CTU_WIDTH + 1];
  memset(bNeighborFlags, 0, totalUnits);
  bool aboveIsAvailable, leftIsAvailable;

  int availlableUnit = isLeftAvailable(isChroma(pu.chType) ? cu : lumaCU, toChannelType(area.compID), area.pos(),
                                       iLeftUnits, iUnitHeight, (bNeighborFlags + iLeftUnits + leftBelowUnits - 1));

  leftIsAvailable = availlableUnit == iTUHeightInUnits;

  availlableUnit = isAboveAvailable(isChroma(pu.chType) ? cu : lumaCU, toChannelType(area.compID), area.pos(),
                                    iAboveUnits, iUnitWidth, (bNeighborFlags + iLeftUnits + leftBelowUnits + 1));

  aboveIsAvailable = availlableUnit == iTUWidthInUnits;

  if (leftIsAvailable)   // if left is not available, then the below left is not available
  {
    avaiLeftBelowUnits = isBelowLeftAvailable(isChroma(pu.chType) ? cu : lumaCU, toChannelType(area.compID), area.bottomLeftComp(area.compID), leftBelowUnits, iUnitHeight, (bNeighborFlags + leftBelowUnits - 1));
  }

  if (aboveIsAvailable)   // if above is not available, then  the above right is not available.
  {
    avaiAboveRightUnits = isAboveRightAvailable(isChroma(pu.chType) ? cu : lumaCU, toChannelType(area.compID), area.topRightComp(area.compID), aboveRightUnits, iUnitWidth, (bNeighborFlags + iLeftUnits + leftBelowUnits + iAboveUnits + 1));
  }

  Pel*       pDst  = nullptr;
  Pel const* piSrc = nullptr;

  bool isFirstRowOfCtu = (lumaArea.y & ((pu.cs->sps)->getCTUSize() - 1)) == 0;

  if (aboveIsAvailable)
  {
    pDst  = pDst0    - iDstStride;
    int addedAboveRight = 0;
#if MMLM
    if ((curChromaMode == MDLM_L_IDX) || (curChromaMode == MDLM_T_IDX) || (curChromaMode == MMLM_L_IDX) || (curChromaMode == MMLM_T_IDX))
#else
    if ((curChromaMode == MDLM_L_IDX) || (curChromaMode == MDLM_T_IDX))
#endif
    {
      addedAboveRight = avaiAboveRightUnits*chromaUnitWidth;
    }
    for (int i = 0; i < uiCWidth + addedAboveRight; i++)
    {
      const bool leftPadding = i == 0 && !leftIsAvailable;
      if (pu.chromaFormat == CHROMA_444)
      {
        piSrc = pRecSrc0 - iRecStride;
        pDst[i] = piSrc[i];
      }
      else if (isFirstRowOfCtu)
      {
        piSrc   = pRecSrc0 - iRecStride;
        pDst[i] = (piSrc[2 * i] * 2 + piSrc[2 * i - (leftPadding ? 0 : 1)] + piSrc[2 * i + 1] + 2) >> 2;
      }
      else if (pu.chromaFormat == CHROMA_422)
      {
        piSrc = pRecSrc0 - iRecStride2;

        int s = 2;
        s += piSrc[2 * i] * 2;
        s += piSrc[2 * i - (leftPadding ? 0 : 1)];
        s += piSrc[2 * i + 1];
        pDst[i] = s >> 2;
      }
      else if (pu.cs->sps->getCclmCollocatedChromaFlag())
      {
        piSrc = pRecSrc0 - iRecStride2;

        int s = 4;
        s += piSrc[2 * i - iRecStride];
        s += piSrc[2 * i] * 4;
        s += piSrc[2 * i - (leftPadding ? 0 : 1)];
        s += piSrc[2 * i + 1];
        s += piSrc[2 * i + iRecStride];
        pDst[i] = s >> 3;
      }
      else
      {
        piSrc = pRecSrc0 - iRecStride2;
        int s = 4;
        s += piSrc[2 * i] * 2;
        s += piSrc[2 * i + 1];
        s += piSrc[2 * i - (leftPadding ? 0 : 1)];
        s += piSrc[2 * i + iRecStride] * 2;
        s += piSrc[2 * i + 1 + iRecStride];
        s += piSrc[2 * i + iRecStride - (leftPadding ? 0 : 1)];
        pDst[i] = s >> 3;
      }
    }
  }

  if (leftIsAvailable)
  {
    pDst  = pDst0    - 1;
    piSrc = pRecSrc0 - 1 - logSubWidthC;

    int addedLeftBelow = 0;
#if MMLM
    if ((curChromaMode == MDLM_L_IDX) || (curChromaMode == MDLM_T_IDX) || (curChromaMode == MMLM_L_IDX) || (curChromaMode == MMLM_T_IDX))
#else
    if ((curChromaMode == MDLM_L_IDX) || (curChromaMode == MDLM_T_IDX))
#endif
    {
      addedLeftBelow = avaiLeftBelowUnits*chromaUnitHeight;
    }

    for (int j = 0; j < uiCHeight + addedLeftBelow; j++)
    {
      if (pu.chromaFormat == CHROMA_444)
      {
        pDst[0] = piSrc[0];
      }
      else if (pu.chromaFormat == CHROMA_422)
      {
        int s = 2;
        s += piSrc[0] * 2;
        s += piSrc[-1];
        s += piSrc[1];
        pDst[0] = s >> 2;
      }
      else if (pu.cs->sps->getCclmCollocatedChromaFlag())
      {
        const bool abovePadding = j == 0 && !aboveIsAvailable;

        int s = 4;
        s += piSrc[-(abovePadding ? 0 : iRecStride)];
        s += piSrc[0] * 4;
        s += piSrc[-1];
        s += piSrc[1];
        s += piSrc[iRecStride];
        pDst[0] = s >> 3;
      }
      else
      {
        int s = 4;
        s += piSrc[0] * 2;
        s += piSrc[1];
        s += piSrc[-1];
        s += piSrc[iRecStride] * 2;
        s += piSrc[iRecStride + 1];
        s += piSrc[iRecStride - 1];
        pDst[0] = s >> 3;
      }

      piSrc += iRecStride2;
      pDst  += iDstStride;
    }
  }

  // inner part from reconstructed picture buffer
  for( int j = 0; j < uiCHeight; j++ )
  {
    for( int i = 0; i < uiCWidth; i++ )
    {
      if (pu.chromaFormat == CHROMA_444)
      {
        pDst0[i] = pRecSrc0[i];
      }
      else if (pu.chromaFormat == CHROMA_422)
      {
        const bool leftPadding  = i == 0 && !leftIsAvailable;

        int s = 2;
        s += pRecSrc0[2 * i] * 2;
        s += pRecSrc0[2 * i - (leftPadding ? 0 : 1)];
        s += pRecSrc0[2 * i + 1];
        pDst0[i] = s >> 2;
      }
      else if (pu.cs->sps->getCclmCollocatedChromaFlag())
      {
        const bool leftPadding  = i == 0 && !leftIsAvailable;
        const bool abovePadding = j == 0 && !aboveIsAvailable;

        int s = 4;
        s += pRecSrc0[2 * i - (abovePadding ? 0 : iRecStride)];
        s += pRecSrc0[2 * i] * 4;
        s += pRecSrc0[2 * i - (leftPadding ? 0 : 1)];
        s += pRecSrc0[2 * i + 1];
        s += pRecSrc0[2 * i + iRecStride];
        pDst0[i] = s >> 3;
      }
      else
      {
        CHECK(pu.chromaFormat != CHROMA_420, "Chroma format must be 4:2:0 for vertical filtering");
        const bool leftPadding = i == 0 && !leftIsAvailable;

        int s = 4;
        s += pRecSrc0[2 * i] * 2;
        s += pRecSrc0[2 * i + 1];
        s += pRecSrc0[2 * i - (leftPadding ? 0 : 1)];
        s += pRecSrc0[2 * i + iRecStride] * 2;
        s += pRecSrc0[2 * i + 1 + iRecStride];
        s += pRecSrc0[2 * i + iRecStride - (leftPadding ? 0 : 1)];
        pDst0[i] = s >> 3;
      }
    }

    pDst0    += iDstStride;
    pRecSrc0 += iRecStride2;
  }
}

#if !LMS_LINEAR_MODEL
void IntraPrediction::xGetLMParameters(const PredictionUnit &pu, const ComponentID compID,
                                              const CompArea &chromaArea,
#if MMLM
                                              int &a, int &b, int &iShift,
                                              int &a2, int &b2, int &iShift2, int &yThres)
#else
                                              int &a, int &b, int &iShift)
#endif
{
  CHECK(compID == COMPONENT_Y, "");

  const SizeType cWidth  = chromaArea.width;
  const SizeType cHeight = chromaArea.height;

  const Position posLT = chromaArea;

  CodingStructure & cs = *(pu.cs);
  const CodingUnit &cu = *(pu.cu);

  const SPS &        sps           = *cs.sps;
  const uint32_t     tuWidth     = chromaArea.width;
  const uint32_t     tuHeight    = chromaArea.height;
  const ChromaFormat nChromaFormat = sps.getChromaFormatIdc();

  const int baseUnitSize = 1 << MIN_CU_LOG2;
  const int unitWidth    = baseUnitSize >> getComponentScaleX(chromaArea.compID, nChromaFormat);
  const int unitHeight   = baseUnitSize >> getComponentScaleY(chromaArea.compID, nChromaFormat);

  const int tuWidthInUnits  = tuWidth / unitWidth;
  const int tuHeightInUnits = tuHeight / unitHeight;
  const int aboveUnits      = tuWidthInUnits;
  const int leftUnits       = tuHeightInUnits;
  int topTemplateSampNum = 2 * cWidth; // for MDLM, the template sample number is 2W or 2H;
  int leftTemplateSampNum = 2 * cHeight;
  assert(m_topRefLength >= topTemplateSampNum);
  assert(m_leftRefLength >= leftTemplateSampNum);
  int totalAboveUnits = (topTemplateSampNum + (unitWidth - 1)) / unitWidth;
  int totalLeftUnits = (leftTemplateSampNum + (unitHeight - 1)) / unitHeight;
  int totalUnits = totalLeftUnits + totalAboveUnits + 1;
  int aboveRightUnits = totalAboveUnits - aboveUnits;
  int leftBelowUnits = totalLeftUnits - leftUnits;
  int avaiAboveRightUnits = 0;
  int avaiLeftBelowUnits = 0;
  int avaiAboveUnits = 0;
  int avaiLeftUnits = 0;

  int curChromaMode = pu.intraDir[1];
  bool neighborFlags[4 * MAX_NUM_PART_IDXS_IN_CTU_WIDTH + 1];
  memset(neighborFlags, 0, totalUnits);

  bool aboveAvailable, leftAvailable;

  int availableUnit =
    isAboveAvailable(cu, CHANNEL_TYPE_CHROMA, posLT, aboveUnits, unitWidth,
    (neighborFlags + leftUnits + leftBelowUnits + 1));
  aboveAvailable = availableUnit == tuWidthInUnits;

  availableUnit =
    isLeftAvailable(cu, CHANNEL_TYPE_CHROMA, posLT, leftUnits, unitHeight,
    (neighborFlags + leftUnits + leftBelowUnits - 1));
  leftAvailable = availableUnit == tuHeightInUnits;
  if (leftAvailable) // if left is not available, then the below left is not available
  {
    avaiLeftUnits = tuHeightInUnits;
    avaiLeftBelowUnits = isBelowLeftAvailable(cu, CHANNEL_TYPE_CHROMA, chromaArea.bottomLeftComp(chromaArea.compID), leftBelowUnits, unitHeight, (neighborFlags + leftBelowUnits - 1));
  }
  if (aboveAvailable) // if above is not available, then  the above right is not available.
  {
    avaiAboveUnits = tuWidthInUnits;
    avaiAboveRightUnits = isAboveRightAvailable(cu, CHANNEL_TYPE_CHROMA, chromaArea.topRightComp(chromaArea.compID), aboveRightUnits, unitWidth, (neighborFlags + leftUnits + leftBelowUnits + aboveUnits + 1));
  }
  Pel *srcColor0, *curChroma0;
  int srcStride;

  PelBuf temp;
#if MMLM
  if ((curChromaMode == MDLM_L_IDX) || (curChromaMode == MDLM_T_IDX) || (curChromaMode == MMLM_L_IDX) || (curChromaMode == MMLM_T_IDX)
    || (m_encPreRDRun && curChromaMode == MMLM_CHROMA_IDX))
#else
  if ((curChromaMode == MDLM_L_IDX) || (curChromaMode == MDLM_T_IDX))
#endif
  {
    srcStride = 2 * MAX_CU_SIZE + 1;
    temp = PelBuf(m_pMdlmTemp + srcStride + 1, srcStride, Size(chromaArea));
  }
  else
  {
    srcStride = MAX_CU_SIZE + 1;
    temp        = PelBuf(m_piTemp + srcStride + 1, srcStride, Size(chromaArea));
  }
  srcColor0 = temp.bufAt(0, 0);
  curChroma0 = getPredictorPtr(compID);

  unsigned internalBitDepth = sps.getBitDepth(CHANNEL_TYPE_CHROMA);

  int minLuma[2] = {  MAX_INT, 0 };
  int maxLuma[2] = { -MAX_INT, 0 };
#if MMLM
  int minLuma2[2][2] = { { MAX_INT, 0 },{ MAX_INT, 0 } };
  int maxLuma2[2][2] = { { -MAX_INT, 0 },{ -MAX_INT, 0 } };
  int minDim = 1;
#endif
  Pel *src = srcColor0 - srcStride;
  int actualTopTemplateSampNum = 0;
  int actualLeftTemplateSampNum = 0;
#if MMLM
  if (curChromaMode == MDLM_T_IDX || curChromaMode == MMLM_T_IDX)
#else
  if (curChromaMode == MDLM_T_IDX)
#endif
  {
    leftAvailable = 0;
    avaiAboveRightUnits = avaiAboveRightUnits > (cHeight/unitWidth) ?  cHeight/unitWidth : avaiAboveRightUnits;
    actualTopTemplateSampNum = unitWidth*(avaiAboveUnits + avaiAboveRightUnits);
#if MMLM
    minDim = actualTopTemplateSampNum;
#endif
  }
#if MMLM
  else if (curChromaMode == MDLM_L_IDX || curChromaMode == MMLM_L_IDX)
#else
  else if (curChromaMode == MDLM_L_IDX)
#endif
  {
    aboveAvailable = 0;
    avaiLeftBelowUnits = avaiLeftBelowUnits > (cWidth/unitHeight) ? cWidth/unitHeight : avaiLeftBelowUnits;
    actualLeftTemplateSampNum = unitHeight*(avaiLeftUnits + avaiLeftBelowUnits);
#if MMLM
    minDim = actualLeftTemplateSampNum;
#endif
  }
#if MMLM
  else if (curChromaMode == LM_CHROMA_IDX || curChromaMode == MMLM_CHROMA_IDX)
#else
  else if (curChromaMode == LM_CHROMA_IDX)
#endif
  {
    actualTopTemplateSampNum = cWidth;
    actualLeftTemplateSampNum = cHeight;
#if MMLM
    minDim = leftAvailable && aboveAvailable ? 1 << g_aucPrevLog2[std::min(actualLeftTemplateSampNum, actualTopTemplateSampNum)]
      : 1 << g_aucPrevLog2[leftAvailable ? actualLeftTemplateSampNum : actualTopTemplateSampNum];
#endif
  }
#if MMLM
  int numSteps = minDim;
  int yAvg = 0;
  int avgCnt = 0;
#endif
  int startPos[2]; //0:Above, 1: Left
  int pickStep[2];

  int aboveIs4 = leftAvailable  ? 0 : 1;
  int leftIs4 =  aboveAvailable ? 0 : 1;

  startPos[0] = actualTopTemplateSampNum >> (2 + aboveIs4);
  pickStep[0] = std::max(1, actualTopTemplateSampNum >> (1 + aboveIs4));

  startPos[1] = actualLeftTemplateSampNum >> (2 + leftIs4);
  pickStep[1] = std::max(1, actualLeftTemplateSampNum >> (1 + leftIs4));

  Pel selectLumaPix[4] = { 0, 0, 0, 0 };
  Pel selectChromaPix[4] = { 0, 0, 0, 0 };

  int cntT, cntL;
  cntT = cntL = 0;
  int cnt = 0;
  if (aboveAvailable)
  {
    cntT = std::min(actualTopTemplateSampNum, (1 + aboveIs4) << 1);
    src = srcColor0 - srcStride;
    const Pel *cur = curChroma0 + 1;
    for (int pos = startPos[0]; cnt < cntT; pos += pickStep[0], cnt++)
    {
      selectLumaPix[cnt] = src[pos];
      selectChromaPix[cnt] = cur[pos];
    }
#if MMLM
    for (int j = 0; j < numSteps; j++)
    {
      int idx = (j * actualTopTemplateSampNum) / minDim;

      if (minLuma2[0][0] > src[idx])
      {
        minLuma2[0][0] = src[idx];
        minLuma2[0][1] = cur[idx];
      }
      if (maxLuma2[1][0] < src[idx])
      {
        maxLuma2[1][0] = src[idx];
        maxLuma2[1][1] = cur[idx];
      }

      yAvg += src[idx];
      avgCnt++;
    }
#endif
  }

  if (leftAvailable)
  {
    cntL = std::min(actualLeftTemplateSampNum, ( 1 + leftIs4 ) << 1 );
    src = srcColor0 - 1;
    const Pel *cur = curChroma0 + m_refBufferStride[compID] + 1;
    for (int pos = startPos[1], cnt = 0; cnt < cntL; pos += pickStep[1], cnt++)
    {
      selectLumaPix[cnt + cntT] = src[pos * srcStride];
      selectChromaPix[cnt + cntT] = cur[pos];
    }
#if MMLM
    for (int i = 0; i < numSteps; i++)
    {
      int idx = (i * actualLeftTemplateSampNum) / minDim;


      if (minLuma2[0][0] > src[srcStride * idx])
      {
        minLuma2[0][0] = src[srcStride * idx];
        minLuma2[0][1] = cur[idx];
      }
      if (maxLuma2[1][0] < src[srcStride * idx])
      {
        maxLuma2[1][0] = src[srcStride * idx];
        maxLuma2[1][1] = cur[idx];
      }

      yAvg += src[srcStride * idx];
      avgCnt++;
    }
#endif
  }
  cnt = cntL + cntT;

  if (cnt == 2)
  {
    selectLumaPix[3] = selectLumaPix[0]; selectChromaPix[3] = selectChromaPix[0];
    selectLumaPix[2] = selectLumaPix[1]; selectChromaPix[2] = selectChromaPix[1];
    selectLumaPix[0] = selectLumaPix[1]; selectChromaPix[0] = selectChromaPix[1];
    selectLumaPix[1] = selectLumaPix[3]; selectChromaPix[1] = selectChromaPix[3];
  }

  int minGrpIdx[2] = { 0, 2 };
  int maxGrpIdx[2] = { 1, 3 };
  int *tmpMinGrp = minGrpIdx;
  int *tmpMaxGrp = maxGrpIdx;
  if (selectLumaPix[tmpMinGrp[0]] > selectLumaPix[tmpMinGrp[1]])
  {
    std::swap(tmpMinGrp[0], tmpMinGrp[1]);
  }
  if (selectLumaPix[tmpMaxGrp[0]] > selectLumaPix[tmpMaxGrp[1]])
  {
    std::swap(tmpMaxGrp[0], tmpMaxGrp[1]);
  }
  if (selectLumaPix[tmpMinGrp[0]] > selectLumaPix[tmpMaxGrp[1]])
  {
    std::swap(tmpMinGrp, tmpMaxGrp);
  }
  if (selectLumaPix[tmpMinGrp[1]] > selectLumaPix[tmpMaxGrp[0]])
  {
    std::swap(tmpMinGrp[1], tmpMaxGrp[0]);
  }

  minLuma[0] = (selectLumaPix[tmpMinGrp[0]] + selectLumaPix[tmpMinGrp[1]] + 1 )>>1;
  minLuma[1] = (selectChromaPix[tmpMinGrp[0]] + selectChromaPix[tmpMinGrp[1]] + 1) >> 1;
  maxLuma[0] = (selectLumaPix[tmpMaxGrp[0]] + selectLumaPix[tmpMaxGrp[1]] + 1 )>>1;
  maxLuma[1] = (selectChromaPix[tmpMaxGrp[0]] + selectChromaPix[tmpMaxGrp[1]] + 1) >> 1;
#if MMLM
  if (avgCnt)
  {
    int x = floorLog2(avgCnt);
    static const uint8_t DivSigTable[1 << 4] = {
      // 4bit significands - 8 ( MSB is omitted )
      0,  7,  6,  5,  5,  4,  4,  3,  3,  2,  2,  1,  1,  1,  1,  0
    };
    int normDiff = (avgCnt << 4 >> x) & 15;
    int v = DivSigTable[normDiff] | 8;
    x += normDiff != 0;

    yAvg = (yAvg * v) >> (x + 3);

  }
  int cntMMLM[2] = { 0,0 };

  //minLuma2[0][0] = minLuma[0]; minLuma2[0][1] = minLuma[1];
  //maxLuma2[1][0] = maxLuma[0]; maxLuma2[1][1] = maxLuma[1];
  src = srcColor0 - srcStride;

  const Pel *curMMLM = curChroma0 + 1;
  if (aboveAvailable)
  {
    for (int j = 0; j < numSteps; j++)
    {
      int idx = (j * actualTopTemplateSampNum) / minDim;
      if (src[idx] >= yAvg)
      {
        if (minLuma2[1][0] > src[idx])
        {
          minLuma2[1][0] = src[idx];
          minLuma2[1][1] = curMMLM[idx];
        }
        cntMMLM[1]++;
      }
      else
      {
        if (maxLuma2[0][0] < src[idx])
        {
          maxLuma2[0][0] = src[idx];
          maxLuma2[0][1] = curMMLM[idx];
        }
        cntMMLM[0]++;
      }
    }
  }

  if (leftAvailable)
  {
    src = srcColor0 - 1;
    //curMMLM = curChroma0 - 1; // should check here -1 or not
    const Pel *curMMLM = curChroma0 + m_refBufferStride[compID] + 1;
    for (int i = 0; i < numSteps; i++)
    {
      int idx = (i * actualLeftTemplateSampNum) / minDim;

      if (src[srcStride * idx] >= yAvg)
      {
        if (minLuma2[1][0] > src[srcStride * idx])
        {
          minLuma2[1][0] = src[srcStride * idx];
          minLuma2[1][1] = curMMLM[idx];
        }
        cntMMLM[1]++;
      }
      else
      {
        if (maxLuma2[0][0] < src[srcStride * idx])
        {
          maxLuma2[0][0] = src[srcStride * idx];
          maxLuma2[0][1] = curMMLM[idx];
        }
        cntMMLM[0]++;
      }
    }
  }

  if (PU::isMultiModeLM(curChromaMode))
  {
    CHECK(cntMMLM[0] && minLuma2[0][0] > maxLuma2[0][0], "Invalid class");
    CHECK(cntMMLM[1] && minLuma2[1][0] > maxLuma2[1][0], "Invalid class");
    CHECK(cntMMLM[0] && cntMMLM[1] && maxLuma2[0][0] > minLuma2[1][0], "Invalid boundary");

    for (int i = 0; i < 2; i++)
    {
      int ax = 0, bx = 0, iShiftx = 0;
      if (cntMMLM[i])
      {
        int diff = maxLuma2[i][0] - minLuma2[i][0];
        if (diff > 0)
        {
          int diffC = maxLuma2[i][1] - minLuma2[i][1];
          int x = floorLog2(diff);
          static const uint8_t DivSigTable[1 << 4] = {
            // 4bit significands - 8 ( MSB is omitted )
            0,  7,  6,  5,  5,  4,  4,  3,  3,  2,  2,  1,  1,  1,  1,  0
          };
          int normDiff = (diff << 4 >> x) & 15;
          int v = DivSigTable[normDiff] | 8;
          x += normDiff != 0;

          int y = floorLog2(abs(diffC)) + 1;
          int add = 1 << y >> 1;
          ax = (diffC * v + add) >> y;
          iShiftx = 3 + x - y;
          if (iShiftx < 1) {
            iShiftx = 1;
            ax = ((ax == 0) ? 0 : (ax < 0) ? -15 : 15);   // a=Sign(a)*15
          }
          bx = minLuma2[i][1] - ((ax * minLuma2[i][0]) >> iShiftx);
        }
        else
        {
          ax = 0;
          bx = minLuma2[i][1];
          iShiftx = 0;
        }
      }
      else
      {
        ax = 0; bx = 1 << (internalBitDepth - 1); iShiftx = 0;
      }
      if (i == 0) { a = ax; b = bx; iShift = iShiftx; }
      else { a2 = ax; b2 = bx; iShift2 = iShiftx; }
    }
    yThres = yAvg;
  }
  else
  {
#endif // Non MMLM mode
  if (leftAvailable || aboveAvailable)
  {
    int diff = maxLuma[0] - minLuma[0];
    if (diff > 0)
    {
      int diffC = maxLuma[1] - minLuma[1];
      int x = floorLog2( diff );
      static const uint8_t DivSigTable[1 << 4] = {
        // 4bit significands - 8 ( MSB is omitted )
        0,  7,  6,  5,  5,  4,  4,  3,  3,  2,  2,  1,  1,  1,  1,  0
      };
      int normDiff = (diff << 4 >> x) & 15;
      int v = DivSigTable[normDiff] | 8;
      x += normDiff != 0;

      int y = floorLog2( abs( diffC ) ) + 1;
      int add = 1 << y >> 1;
      a = (diffC * v + add) >> y;
      iShift = 3 + x - y;
      if ( iShift < 1 )
      {
        iShift = 1;
        a = ( (a == 0)? 0: (a < 0)? -15 : 15 );   // a=Sign(a)*15
      }
      b = minLuma[1] - ((a * minLuma[0]) >> iShift);
    }
    else
    {
      a = 0;
      b = minLuma[1];
      iShift = 0;
    }
  }
  else
  {
    a = 0;

    b = 1 << (internalBitDepth - 1);

    iShift = 0;
  }
#if MMLM
  }
#endif
}
#endif

void IntraPrediction::initIntraMip( const PredictionUnit &pu, const CompArea &area )
{
  CHECK( area.width > MIP_MAX_WIDTH || area.height > MIP_MAX_HEIGHT, "Error: block size not supported for MIP" );

  // prepare input (boundary) data for prediction
  CHECK( m_ipaParam.refFilterFlag, "ERROR: unfiltered refs expected for MIP" );
  Pel       *ptrSrc     = getPredictorPtr(area.compID);
  const int  srcStride  = m_refBufferStride[area.compID];
  const int  srcHStride = 2;

  m_matrixIntraPred.prepareInputForPred(CPelBuf(ptrSrc, srcStride, srcHStride), area,
                                        pu.cu->slice->getSPS()->getBitDepth(toChannelType(area.compID)), area.compID);
}

void IntraPrediction::predIntraMip( const ComponentID compId, PelBuf &piPred, const PredictionUnit &pu )
{
  CHECK( piPred.width > MIP_MAX_WIDTH || piPred.height > MIP_MAX_HEIGHT, "Error: block size not supported for MIP" );
  CHECK( piPred.width != (1 << floorLog2(piPred.width)) || piPred.height != (1 << floorLog2(piPred.height)), "Error: expecting blocks of size 2^M x 2^N" );

  // generate mode-specific prediction
  uint32_t modeIdx       = MAX_NUM_MIP_MODE;
  bool     transposeFlag = false;
  if (compId == COMPONENT_Y)
  {
    modeIdx       = pu.intraDir[CHANNEL_TYPE_LUMA];
    transposeFlag = pu.mipTransposedFlag;
  }
  else
  {
    const PredictionUnit &coLocatedLumaPU = PU::getCoLocatedLumaPU(pu);

    CHECK(pu.intraDir[CHANNEL_TYPE_CHROMA] != DM_CHROMA_IDX, "Error: MIP is only supported for chroma with DM_CHROMA.");
    CHECK(!coLocatedLumaPU.cu->mipFlag, "Error: Co-located luma CU should use MIP.");

    modeIdx       = coLocatedLumaPU.intraDir[CHANNEL_TYPE_LUMA];
    transposeFlag = coLocatedLumaPU.mipTransposedFlag;
  }
  const int bitDepth = pu.cu->slice->getSPS()->getBitDepth(toChannelType(compId));

  CHECK(modeIdx >= getNumModesMip(piPred), "Error: Wrong MIP mode index");

  static_vector<int, MIP_MAX_WIDTH* MIP_MAX_HEIGHT> predMip( piPred.width * piPred.height );
  m_matrixIntraPred.predBlock(predMip.data(), modeIdx, transposeFlag, bitDepth, compId);

  Pel *pred = piPred.buf;
  int idx = 0;

  for( int y = 0; y < piPred.height; y++ )
  {
    for( int x = 0; x < piPred.width; x++ )
    {
      pred[x] = Pel(predMip[idx++]);
    }

    pred += piPred.stride;
  }
}

void IntraPrediction::reorderPLT(CodingStructure& cs, Partitioner& partitioner, ComponentID compBegin, uint32_t numComp)
{
  CodingUnit &cu = *cs.getCU(partitioner.chType);

  uint8_t        reusePLTSizetmp = 0;
  uint8_t        pltSizetmp = 0;
  Pel            curPLTtmp[MAX_NUM_COMPONENT][MAXPLTSIZE];
  bool           curPLTpred[MAXPLTPREDSIZE];

  for (int idx = 0; idx < MAXPLTPREDSIZE; idx++)
  {
    curPLTpred[idx] = false;
    cu.reuseflag[compBegin][idx] = false;
  }
  for (int idx = 0; idx < MAXPLTSIZE; idx++)
  {
    curPLTpred[idx] = false;
  }

  for (int predidx = 0; predidx < cs.prevPLT.curPLTSize[compBegin]; predidx++)
  {
    bool match = false;
    int curidx = 0;

    for (curidx = 0; curidx < cu.curPLTSize[compBegin]; curidx++)
    {
      if( curPLTpred[curidx] )
      {
        continue;
      }
      bool matchTmp = true;
      for (int comp = compBegin; comp < (compBegin + numComp); comp++)
      {
        matchTmp = matchTmp && (cu.curPLT[comp][curidx] == cs.prevPLT.curPLT[comp][predidx]);
      }
      if (matchTmp)
      {
        match = true;
        break;
      }
    }

    if (match)
    {
      cu.reuseflag[compBegin][predidx] = true;
      curPLTpred[curidx] = true;
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
      if( cu.isLocalSepTree() )
#else
      if (CS::isDualITree(*cu.cs))
#endif
      {
        cu.reuseflag[COMPONENT_Y][predidx] = true;
        for( int comp = COMPONENT_Y; comp < MAX_NUM_COMPONENT; comp++ )
        {
          curPLTtmp[comp][reusePLTSizetmp] = cs.prevPLT.curPLT[comp][predidx];
        }
      }
      else
      {
        for (int comp = compBegin; comp < (compBegin + numComp); comp++)
        {
          curPLTtmp[comp][reusePLTSizetmp] = cs.prevPLT.curPLT[comp][predidx];
        }
      }
      reusePLTSizetmp++;
      pltSizetmp++;
    }
  }
  cu.reusePLTSize[compBegin] = reusePLTSizetmp;
  for (int curidx = 0; curidx < cu.curPLTSize[compBegin]; curidx++)
  {
    if (!curPLTpred[curidx])
    {
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
      if( cu.isLocalSepTree() )
      {
        for( int comp = compBegin; comp < (compBegin + numComp); comp++ )
        {
          curPLTtmp[comp][pltSizetmp] = cu.curPLT[comp][curidx];
        }
        if( isLuma(partitioner.chType) )
        {
          curPLTtmp[COMPONENT_Cb][pltSizetmp] = 1 << (cs.sps->getBitDepth(CHANNEL_TYPE_CHROMA) - 1);
          curPLTtmp[COMPONENT_Cr][pltSizetmp] = 1 << (cs.sps->getBitDepth(CHANNEL_TYPE_CHROMA) - 1);
        }
        else
        {
          curPLTtmp[COMPONENT_Y][pltSizetmp] = 1 << (cs.sps->getBitDepth(CHANNEL_TYPE_LUMA) - 1);
        }
      }
      else
      {
#endif
        for (int comp = compBegin; comp < (compBegin + numComp); comp++)
        {
          curPLTtmp[comp][pltSizetmp] = cu.curPLT[comp][curidx];
        }
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
      }
#endif
      pltSizetmp++;
    }
  }
  assert(pltSizetmp == cu.curPLTSize[compBegin]);
  for (int curidx = 0; curidx < cu.curPLTSize[compBegin]; curidx++)
  {
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
    if( cu.isLocalSepTree() )
#else
    if (CS::isDualITree(*cu.cs))
#endif
    {
      for( int comp = COMPONENT_Y; comp < MAX_NUM_COMPONENT; comp++ )
      {
        cu.curPLT[comp][curidx] = curPLTtmp[comp][curidx];
      }
    }
    else
    {
    for (int comp = compBegin; comp < (compBegin + numComp); comp++)
    {
      cu.curPLT[comp][curidx] = curPLTtmp[comp][curidx];
    }
    }
  }
}

#if MMLM && LMS_LINEAR_MODEL
int IntraPrediction::xCalcLMParametersGeneralized(int x, int y, int xx, int xy, int count, int bitDepth, int &a, int &b, int &iShift)
{

  uint32_t uiInternalBitDepth = bitDepth;
  if (count == 0)
  {
    a = 0;
    b = 1 << (uiInternalBitDepth - 1);
    iShift = 0;
    return -1;
  }
  CHECK(count > 512, "");


  int iCountShift = g_aucLog2[count];

  int iTempShift = uiInternalBitDepth + iCountShift - 15;

  if (iTempShift > 0)
  {
    x = (x + (1 << (iTempShift - 1))) >> iTempShift;
    y = (y + (1 << (iTempShift - 1))) >> iTempShift;
    xx = (xx + (1 << (iTempShift - 1))) >> iTempShift;
    xy = (xy + (1 << (iTempShift - 1))) >> iTempShift;
    iCountShift -= iTempShift;
  }
  /////// xCalcLMParameters

  int avgX = x >> iCountShift;
  int avgY = y >> iCountShift;

  int RErrX = x & ((1 << iCountShift) - 1);
  int RErrY = y & ((1 << iCountShift) - 1);

  int iB = 7;
  iShift = 13 - iB;

  if (iCountShift == 0)
  {
    a = 0;
    b = 1 << (uiInternalBitDepth - 1);
    iShift = 0;
  }
  else
  {
    int a1 = xy - (avgX * avgY << iCountShift) - avgX * RErrY - avgY * RErrX;
    int a2 = xx - (avgX * avgX << iCountShift) - 2 * avgX * RErrX;
    const int iShiftA1 = uiInternalBitDepth - 2;
    const int iShiftA2 = 5;
    const int iAccuracyShift = uiInternalBitDepth + 4;

    int iScaleShiftA2 = 0;
    int iScaleShiftA1 = 0;
    int a1s = a1;
    int a2s = a2;

    iScaleShiftA1 = a1 == 0 ? 0 : floorLog2(abs(a1)) - iShiftA1;
    iScaleShiftA2 = a2 == 0 ? 0 : floorLog2(abs(a2)) - iShiftA2;

    if (iScaleShiftA1 < 0)
    {
      iScaleShiftA1 = 0;
    }

    if (iScaleShiftA2 < 0)
    {
      iScaleShiftA2 = 0;
    }

    int iScaleShiftA = iScaleShiftA2 + iAccuracyShift - iShift - iScaleShiftA1;

    a2s = a2 >> iScaleShiftA2;

    a1s = a1 >> iScaleShiftA1;

    if (a2s >= 32)
    {
      uint32_t a2t = m_auShiftLM[a2s - 32];
      a = a1s * a2t;
    }
    else
    {
      a = 0;
    }

    if (iScaleShiftA < 0)
    {
      a = a << -iScaleShiftA;
    }
    else
    {
      a = a >> iScaleShiftA;
    }
    a = Clip3(-(1 << (15 - iB)), (1 << (15 - iB)) - 1, a);
    a = a << iB;

    int16_t n = 0;
    if (a != 0)
    {
      n = floorLog2(abs(a) + ((a < 0 ? -1 : 1) - 1) / 2) - 5;
    }

    iShift = (iShift + iB) - n;
    a = a >> n;

    b = avgY - ((a * avgX) >> iShift);
  }
  return 0;
}

int IntraPrediction::xLMSampleClassifiedTraining(int count, int mean, int meanC, int LumaSamples[], int ChrmSamples[], int bitDepth, MMLM_parameter parameters[])
{

  //Initialize

  for (int i = 0; i < 2; i++)
  {
    parameters[i].a = 0;
    parameters[i].b = 1 << (bitDepth - 1);
    parameters[i].shift = 0;
  }

  if (count < 4)//
  {
    return -1;
  }
  int GroupCount[2] = { 0, 0 };

  CHECK(count > 512, "");

  int meanDiff = meanC - mean;
  mean = std::max(1, mean);

  int LumaPower2[2][128];
  int ChromaPower2[2][128];
  //int GroupCount[2] = { 0, 0 };
  for (int i = 0; i < count; i++)
  {
    if (LumaSamples[i] <= mean)
    {
      LumaPower2[0][GroupCount[0]] = LumaSamples[i];
      ChromaPower2[0][GroupCount[0]] = ChrmSamples[i];
      GroupCount[0]++;
    }
    else
    {
      LumaPower2[1][GroupCount[1]] = LumaSamples[i];
      ChromaPower2[1][GroupCount[1]] = ChrmSamples[i];
      GroupCount[1]++;
    }
  }

  // Take power of two
  for (int group = 0; group < 2; group++)
  {
    int existSampNum = GroupCount[group];
    if (existSampNum < 2)
    {
      continue;
    }

    int upperPower2 = 1 << (g_aucLog2[existSampNum - 1] + 1);
    int lowerPower2 = 1 << (g_aucLog2[existSampNum]);

    if (upperPower2 != lowerPower2)
    {
      int numPaddedSamples = std::min(existSampNum, upperPower2 - existSampNum);
      GroupCount[group] = upperPower2;
      int step = (int)(existSampNum / numPaddedSamples);
      for (int i = 0; i < numPaddedSamples; i++)
      {
        LumaPower2[group][existSampNum + i] = LumaPower2[group][i * step];
        ChromaPower2[group][existSampNum + i] = ChromaPower2[group][i * step];

      }
    }
  }

  int x[2], y[2], xy[2], xx[2];
  for (int group = 0; group < 2; group++)
  {
    x[group] = y[group] = xy[group] = xx[group] = 0;
  }

  for (int group = 0; group < 2; group++)
  {

    for (int i = 0; i < GroupCount[group]; i++)
    {
      x[group] += LumaPower2[group][i];
      y[group] += ChromaPower2[group][i];
      xx[group] += LumaPower2[group][i] * LumaPower2[group][i];
      xy[group] += LumaPower2[group][i] * ChromaPower2[group][i];
    }
  }
  for (int group = 0; group < 2; group++)
  {
    int a, b, iShift;
    if (GroupCount[group] > 1)
    {
      xCalcLMParametersGeneralized(x[group], y[group], xx[group], xy[group], GroupCount[group], bitDepth, a, b, iShift);

      parameters[group].a = a;
      parameters[group].b = b;
      parameters[group].shift = iShift;
    }
    else
    {
      parameters[group].a = 0;
      parameters[group].b = meanDiff;
      parameters[group].shift = 0;
    }
  }
  return 0;
}
#endif
#if LMS_LINEAR_MODEL
void IntraPrediction::xPadMdlmTemplateSample(Pel*pSrc, Pel*pCur, int cWidth, int cHeight, int existSampNum, int targetSampNum)
{
  int sampNumToBeAdd = targetSampNum - existSampNum;
  Pel*pTempSrc = pSrc + existSampNum;
  Pel*pTempCur = pCur + existSampNum;

  int step = (int)(existSampNum / sampNumToBeAdd);

  for (int i = 0; i < sampNumToBeAdd; i++)
  {
    pTempSrc[i] = pSrc[i * step];
    pTempCur[i] = pCur[i * step];
  }
}
void IntraPrediction::xGetLMParameters_LMS(const PredictionUnit &pu, const ComponentID compID, const CompArea& chromaArea,
#if MMLM
  int &a, int &b, int &iShift,
  int &a2, int &b2, int &iShift2, int &yThres)
#else
int &a, int &b, int &iShift)
#endif
{
  CHECK(compID == COMPONENT_Y, "");
  const SizeType cWidth = chromaArea.width;
  const SizeType cHeight = chromaArea.height;

  const Position posLT = chromaArea;

  CodingStructure & cs = *(pu.cs);
  const CodingUnit &cu = *(pu.cu);

  const SPS &        sps = *cs.sps;
  const uint32_t     tuWidth = chromaArea.width;
  const uint32_t     tuHeight = chromaArea.height;
  const ChromaFormat nChromaFormat = sps.getChromaFormatIdc();

  const int baseUnitSize = 1 << MIN_CU_LOG2;
  const int unitWidth = baseUnitSize >> getComponentScaleX(chromaArea.compID, nChromaFormat);
  const int unitHeight = baseUnitSize >> getComponentScaleY(chromaArea.compID, nChromaFormat);

  const int tuWidthInUnits = tuWidth / unitWidth;
  const int tuHeightInUnits = tuHeight / unitHeight;
  const int aboveUnits = tuWidthInUnits;
  const int leftUnits = tuHeightInUnits;
  int topTemplateSampNum = 2 * cWidth; // for MDLM, the template sample number is 2W or 2H;
  int leftTemplateSampNum = 2 * cHeight;
  assert(m_topRefLength >= topTemplateSampNum);
  assert(m_leftRefLength >= leftTemplateSampNum);
  int totalAboveUnits = (topTemplateSampNum + (unitWidth - 1)) / unitWidth;
  int totalLeftUnits = (leftTemplateSampNum + (unitHeight - 1)) / unitHeight;
  int totalUnits = totalLeftUnits + totalAboveUnits + 1;
  int aboveRightUnits = totalAboveUnits - aboveUnits;
  int leftBelowUnits = totalLeftUnits - leftUnits;
  int avaiAboveRightUnits = 0;
  int avaiLeftBelowUnits = 0;
  int avaiAboveUnits = 0;
  int avaiLeftUnits = 0;

  int curChromaMode = pu.intraDir[1];
  bool neighborFlags[4 * MAX_NUM_PART_IDXS_IN_CTU_WIDTH + 1];
  memset(neighborFlags, 0, totalUnits);

  bool aboveAvailable, leftAvailable;

  int availableUnit =
    isAboveAvailable(cu, CHANNEL_TYPE_CHROMA, posLT, aboveUnits, unitWidth,
    (neighborFlags + leftUnits + leftBelowUnits + 1));
  aboveAvailable = availableUnit == tuWidthInUnits;

  availableUnit =
    isLeftAvailable(cu, CHANNEL_TYPE_CHROMA, posLT, leftUnits, unitHeight,
    (neighborFlags + leftUnits + leftBelowUnits - 1));
  leftAvailable = availableUnit == tuHeightInUnits;
  if (leftAvailable) // if left is not available, then the below left is not available
  {
    avaiLeftUnits = tuHeightInUnits;
    avaiLeftBelowUnits = isBelowLeftAvailable(cu, CHANNEL_TYPE_CHROMA, chromaArea.bottomLeftComp(chromaArea.compID), leftBelowUnits, unitHeight, (neighborFlags + leftBelowUnits - 1));
  }
  if (aboveAvailable) // if above is not available, then  the above right is not available.
  {
    avaiAboveUnits = tuWidthInUnits;
    avaiAboveRightUnits = isAboveRightAvailable(cu, CHANNEL_TYPE_CHROMA, chromaArea.topRightComp(chromaArea.compID), aboveRightUnits, unitWidth, (neighborFlags + leftUnits + leftBelowUnits + aboveUnits + 1));
  }
  Pel *srcColor0, *curChroma0;
  int srcStride;

  PelBuf temp;
#if MMLM
  if ((curChromaMode == MDLM_L_IDX) || (curChromaMode == MDLM_T_IDX) || (curChromaMode == MMLM_L_IDX) || (curChromaMode == MMLM_T_IDX)
    || (m_encPreRDRun && curChromaMode == MMLM_CHROMA_IDX))
#else
  if ((curChromaMode == MDLM_L_IDX) || (curChromaMode == MDLM_T_IDX))
#endif
  {
    srcStride = 2 * MAX_CU_SIZE + 1;
    temp = PelBuf(m_pMdlmTemp + srcStride + 1, srcStride, Size(chromaArea));
  }
  else
  {
    srcStride = MAX_CU_SIZE + 1;
    temp = PelBuf(m_piTemp + srcStride + 1, srcStride, Size(chromaArea));
  }
  srcColor0 = temp.bufAt(0, 0);
  curChroma0 = getPredictorPtr(compID);

  int x = 0, y = 0, xx = 0, xy = 0;
  int iCountShift = 0;
  unsigned uiInternalBitDepth = sps.getBitDepth(CHANNEL_TYPE_CHROMA);

  Pel *src = srcColor0 - srcStride;
  int actualTopTemplateSampNum = 0;
  int actualLeftTemplateSampNum = 0;

  //get the temp buffer to store the downsampled luma and chroma
  Pel pTempBufferSrc[2 * MAX_CU_SIZE]; // for MDLM, use tempalte size 2W or 2H,
  Pel pTempBufferCur[2 * MAX_CU_SIZE];
  int minDim = 1;
  int cntT = 0; int cntL = 0;

#if MMLM
  if (curChromaMode == MDLM_T_IDX || curChromaMode == MMLM_T_IDX)
#else
  if (curChromaMode == MDLM_T_IDX)
#endif
  {
    leftAvailable = 0;
#if !LMS_LINEAR_MODEL
    avaiAboveRightUnits = avaiAboveRightUnits > (cHeight / unitWidth) ? cHeight / unitWidth : avaiAboveRightUnits;
#endif
    actualTopTemplateSampNum = unitWidth * (avaiAboveUnits + avaiAboveRightUnits);
    minDim = actualTopTemplateSampNum;
  }
#if MMLM
  else if (curChromaMode == MDLM_L_IDX || curChromaMode == MMLM_L_IDX)
#else
  else if (curChromaMode == MDLM_L_IDX)
#endif
  {
    aboveAvailable = 0;
#if !LMS_LINEAR_MODEL
    avaiLeftBelowUnits = avaiLeftBelowUnits > (cWidth / unitHeight) ? cWidth / unitHeight : avaiLeftBelowUnits;
#endif
    actualLeftTemplateSampNum = unitHeight * (avaiLeftUnits + avaiLeftBelowUnits);
    minDim = actualLeftTemplateSampNum;
  }
#if MMLM
  else if (curChromaMode == LM_CHROMA_IDX || curChromaMode == MMLM_CHROMA_IDX)
#else
  else if (curChromaMode == LM_CHROMA_IDX)
#endif
  {
    actualTopTemplateSampNum = cWidth;
    actualLeftTemplateSampNum = cHeight;
    minDim = leftAvailable && aboveAvailable ? 1 << g_aucPrevLog2[std::min(actualLeftTemplateSampNum, actualTopTemplateSampNum)]
      : 1 << g_aucPrevLog2[leftAvailable ? actualLeftTemplateSampNum : actualTopTemplateSampNum];
  }
  int numSteps = minDim;

  if (aboveAvailable)
  {
    cntT = numSteps;
    src = srcColor0 - srcStride;
    const Pel *cur = curChroma0 + 1;

    for (int j = 0; j < numSteps; j++)
    {
      int idx = (j * actualTopTemplateSampNum) / minDim;

      pTempBufferSrc[j] = src[idx];
      pTempBufferCur[j] = cur[idx];
    }

  }

  if (leftAvailable)
  {
    cntL = numSteps;
    src = srcColor0 - 1;
    const Pel *cur = curChroma0 + m_refBufferStride[compID] + 1;

    for (int i = 0; i < numSteps; i++)
    {
      int idx = (i * actualLeftTemplateSampNum) / minDim;

      pTempBufferSrc[i + cntT] = src[srcStride * idx];
      pTempBufferCur[i + cntT] = cur[idx];

    }
  }

  if ((curChromaMode == MDLM_L_IDX) || (curChromaMode == MDLM_T_IDX))
  {
    //pad the temple sample to targetSampNum.
    int orgNumSample = (curChromaMode == MDLM_T_IDX) ? (avaiAboveUnits*unitWidth) : (avaiLeftUnits*unitHeight);
    int existSampNum = (curChromaMode == MDLM_T_IDX) ? actualTopTemplateSampNum : actualLeftTemplateSampNum;

    if( !orgNumSample || !existSampNum )
    {
      a = 0;
      b = 1 << (uiInternalBitDepth - 1);
      iShift = 0;
      return;
    }

    int targetSampNum = 1 << ( floorLog2( existSampNum - 1 ) + 1 );

    if (targetSampNum != existSampNum)//if existSampNum not a value of power of 2
    {
      xPadMdlmTemplateSample(pTempBufferSrc, pTempBufferCur, cWidth, cHeight, existSampNum, targetSampNum);
    }
    for (int j = 0; j < targetSampNum; j++)
    {
      x += pTempBufferSrc[j];
      y += pTempBufferCur[j];
      xx += pTempBufferSrc[j] * pTempBufferSrc[j];
      xy += pTempBufferSrc[j] * pTempBufferCur[j];
    }
    iCountShift = g_aucLog2[targetSampNum];

  }
  else if (curChromaMode == LM_CHROMA_IDX)
  {
    int       minStep = 1;
    //int       numSteps = minDim;

    if( aboveAvailable )
    {
      iCountShift = g_aucLog2[minDim / minStep];
    }

    if( leftAvailable )
    {
      iCountShift += aboveAvailable ? 1 : g_aucLog2[minDim / minStep];
    }

    for (int i = 0; i < (cntT + cntL); i++)
    {
      x += pTempBufferSrc[i];
      y += pTempBufferCur[i];
      xx += pTempBufferSrc[i] * pTempBufferSrc[i];
      xy += pTempBufferSrc[i] * pTempBufferCur[i];
    }
  }
#if MMLM
  if (PU::isMultiModeLM(pu.intraDir[1]))
  {
    // Classify and training
    MMLM_parameter parameters[2];
    int LumaSamples[512];
    int ChrmSamples[512];
    int meanC = 0; int mean = 0;
    int avgCnt = cntT + cntL;

    for (int i = 0; i < avgCnt; i++)
    {
      mean += pTempBufferSrc[i];
      meanC += pTempBufferCur[i];
      LumaSamples[i] = pTempBufferSrc[i];
      ChrmSamples[i] = pTempBufferCur[i];
    }

    if (avgCnt)
    {
       int x = floorLog2(avgCnt);
       static const uint8_t DivSigTable[1 << 4] = {
         // 4bit significands - 8 ( MSB is omitted )
        0,  7,  6,  5,  5,  4,  4,  3,  3,  2,  2,  1,  1,  1,  1,  0
       };
       int normDiff = (avgCnt << 4 >> x) & 15;
       int v = DivSigTable[normDiff] | 8;
       x += normDiff != 0;

       mean = (mean * v) >> (x + 3);
       meanC = (meanC * v) >> (x + 3);
    }

    xLMSampleClassifiedTraining(avgCnt, mean, meanC, LumaSamples, ChrmSamples, uiInternalBitDepth, parameters);
    a = parameters[0].a; b = parameters[0].b; iShift = parameters[0].shift;
    a2 = parameters[1].a; b2 = parameters[1].b; iShift2 = parameters[1].shift;
    yThres = mean;
    return;
  }
#endif

  if ((curChromaMode == MDLM_L_IDX) || (curChromaMode == MDLM_T_IDX))
  {
    if ((curChromaMode == MDLM_L_IDX) ? (!leftAvailable) : (!aboveAvailable))
    {
      a = 0;
      b = 1 << (uiInternalBitDepth - 1);
      iShift = 0;
      return;
    }
  }
  else
  {
    if (!leftAvailable && !aboveAvailable)
    {
      a = 0;
      b = 1 << (uiInternalBitDepth - 1);
      iShift = 0;
      return;
    }
  }

  int iTempShift = uiInternalBitDepth + iCountShift - 15;

  if (iTempShift > 0)
  {
    x = (x + (1 << (iTempShift - 1))) >> iTempShift;
    y = (y + (1 << (iTempShift - 1))) >> iTempShift;
    xx = (xx + (1 << (iTempShift - 1))) >> iTempShift;
    xy = (xy + (1 << (iTempShift - 1))) >> iTempShift;
    iCountShift -= iTempShift;
  }

  /////// xCalcLMParameters

  int avgX = x >> iCountShift;
  int avgY = y >> iCountShift;

  int RErrX = x & ((1 << iCountShift) - 1);
  int RErrY = y & ((1 << iCountShift) - 1);

  int iB = 7;
  iShift = 13 - iB;

  if (iCountShift == 0)
  {
    a = 0;
    b = 1 << (uiInternalBitDepth - 1);
    iShift = 0;
  }
  else
  {
    int a1 = xy - (avgX * avgY << iCountShift) - avgX * RErrY - avgY * RErrX;
    int a2 = xx - (avgX * avgX << iCountShift) - 2 * avgX * RErrX;
    const int iShiftA1 = uiInternalBitDepth - 2;
    const int iShiftA2 = 5;
    const int iAccuracyShift = uiInternalBitDepth + 4;

    int iScaleShiftA2 = 0;
    int iScaleShiftA1 = 0;
    int a1s = a1;
    int a2s = a2;

    iScaleShiftA1 = a1 == 0 ? 0 : floorLog2(abs(a1)) - iShiftA1;
    iScaleShiftA2 = a2 == 0 ? 0 : floorLog2(abs(a2)) - iShiftA2;

    if (iScaleShiftA1 < 0)
    {
      iScaleShiftA1 = 0;
    }

    if (iScaleShiftA2 < 0)
    {
      iScaleShiftA2 = 0;
    }

    int iScaleShiftA = iScaleShiftA2 + iAccuracyShift - iShift - iScaleShiftA1;

    a2s = a2 >> iScaleShiftA2;

    a1s = a1 >> iScaleShiftA1;

    if (a2s >= 32)
    {
      uint32_t a2t = m_auShiftLM[a2s - 32];
      a = a1s * a2t;
    }
    else
    {
      a = 0;
    }

    if (iScaleShiftA < 0)
    {
      a = a << -iScaleShiftA;
    }
    else
    {
      a = a >> iScaleShiftA;
    }
    a = Clip3(-(1 << (15 - iB)), (1 << (15 - iB)) - 1, a);
    a = a << iB;

    int16_t n = 0;
    if (a != 0)
    {
      n = floorLog2(abs(a) + ((a < 0 ? -1 : 1) - 1) / 2) - 5;
    }

    iShift = (iShift + iB) - n;
    a = a >> n;

    b = avgY - ((a * avgX) >> iShift);
  }
}
#endif

//! \}
