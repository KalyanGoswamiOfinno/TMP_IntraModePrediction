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

/** \file     TrQuant.cpp
    \brief    transform and quantization class
*/

#include "TrQuant.h"
#include "TrQuant_EMT.h"

#include "UnitTools.h"
#include "ContextModelling.h"
#include "CodingStructure.h"

#include "dtrace_buffer.h"

#include <stdlib.h>
#include <limits>
#include <memory.h>

#include "QuantRDOQ.h"
#include "DepQuant.h"

#if RExt__DECODER_DEBUG_TOOL_STATISTICS
#include "CommonLib/CodingStatistics.h"
#endif

#if ENABLE_SIMD_TMP
#include "CommonDefX86.h"
#endif

#if JVET_V0130_INTRA_TMP
unsigned int g_uiDepth2Width[5] = { 4, 8, 16, 32, 64 };
#endif

struct coeffGroupRDStats
{
  int    iNNZbeforePos0;
  double d64CodedLevelandDist; // distortion and level cost only
  double d64UncodedDist;    // all zero coded block distortion
  double d64SigCost;
  double d64SigCost_0;
};

//! \ingroup CommonLib
//! \{

static inline int64_t square( const int d ) { return d * (int64_t)d; }

template<int signedMode> std::pair<int64_t,int64_t> fwdTransformCbCr( const PelBuf &resCb, const PelBuf &resCr, PelBuf& resC1, PelBuf& resC2 )
{
  const Pel*  cb  = resCb.buf;
  const Pel*  cr  = resCr.buf;
  Pel*        c1  = resC1.buf;
  Pel*        c2  = resC2.buf;
  int64_t     d1  = 0;
  int64_t     d2  = 0;
  for( SizeType y = 0; y < resCb.height; y++, cb += resCb.stride, cr += resCr.stride, c1 += resC1.stride, c2 += resC2.stride )
  {
    for( SizeType x = 0; x < resCb.width; x++ )
    {
      int cbx = cb[x], crx = cr[x];
      if      ( signedMode ==  1 )
      {
        c1[x] = Pel( ( 4*cbx + 2*crx ) / 5 );
        d1   += square( cbx - c1[x] ) + square( crx - (c1[x]>>1) );
      }
      else if ( signedMode == -1 )
      {
        c1[x] = Pel( ( 4*cbx - 2*crx ) / 5 );
        d1   += square( cbx - c1[x] ) + square( crx - (-c1[x]>>1) );
      }
      else if ( signedMode ==  2 )
      {
        c1[x] = Pel( ( cbx + crx ) / 2 );
        d1   += square( cbx - c1[x] ) + square( crx - c1[x] );
      }
      else if ( signedMode == -2 )
      {
        c1[x] = Pel( ( cbx - crx ) / 2 );
        d1   += square( cbx - c1[x] ) + square( crx + c1[x] );
      }
      else if ( signedMode ==  3 )
      {
        c2[x] = Pel( ( 4*crx + 2*cbx ) / 5 );
        d1   += square( cbx - (c2[x]>>1) ) + square( crx - c2[x] );
      }
      else if ( signedMode == -3 )
      {
        c2[x] = Pel( ( 4*crx - 2*cbx ) / 5 );
        d1   += square( cbx - (-c2[x]>>1) ) + square( crx - c2[x] );
      }
      else
      {
        d1   += square( cbx );
        d2   += square( crx );
      }
    }
  }
  return std::make_pair(d1,d2);
}

template<int signedMode> void invTransformCbCr( PelBuf &resCb, PelBuf &resCr )
{
  Pel*  cb  = resCb.buf;
  Pel*  cr  = resCr.buf;
  for( SizeType y = 0; y < resCb.height; y++, cb += resCb.stride, cr += resCr.stride )
  {
    for( SizeType x = 0; x < resCb.width; x++ )
    {
      if (signedMode == 1)
      {
        cr[x] = cb[x] >> 1;
      }
      else if (signedMode == -1)
      {
        cr[x] = -cb[x] >> 1;
      }
      else if (signedMode == 2)
      {
        cr[x] = cb[x];
      }
      else if (signedMode == -2)
      {
        // non-normative clipping to prevent 16-bit overflow
        cr[x] = (cb[x] == -32768 && sizeof(Pel) == 2) ? 32767 : -cb[x];
      }
      else if (signedMode == 3)
      {
        cb[x] = cr[x] >> 1;
      }
      else if (signedMode == -3)
      {
        cb[x] = -cr[x] >> 1;
      }
    }
  }
}

// ====================================================================================================================
// TrQuant class member functions
// ====================================================================================================================
#if TRANSFORM_SIMD_OPT
std::array<std::array<const TMatrixCoeff*, g_numTransformMatrixSizes>, NUM_TRANS_TYPE> TrQuant::m_forwardTransformKernels;
std::array<std::array<const TMatrixCoeff*, g_numTransformMatrixSizes>, NUM_TRANS_TYPE> TrQuant::m_inverseTransformKernels;
#endif

TrQuant::TrQuant() : m_quant( nullptr )
{
  // allocate temporary buffers
  {
    m_invICT      = m_invICTMem + maxAbsIctMode;
    m_invICT[ 0]  = invTransformCbCr< 0>;
    m_invICT[ 1]  = invTransformCbCr< 1>;
    m_invICT[-1]  = invTransformCbCr<-1>;
    m_invICT[ 2]  = invTransformCbCr< 2>;
    m_invICT[-2]  = invTransformCbCr<-2>;
    m_invICT[ 3]  = invTransformCbCr< 3>;
    m_invICT[-3]  = invTransformCbCr<-3>;
    m_fwdICT      = m_fwdICTMem + maxAbsIctMode;
    m_fwdICT[ 0]  = fwdTransformCbCr< 0>;
    m_fwdICT[ 1]  = fwdTransformCbCr< 1>;
    m_fwdICT[-1]  = fwdTransformCbCr<-1>;
    m_fwdICT[ 2]  = fwdTransformCbCr< 2>;
    m_fwdICT[-2]  = fwdTransformCbCr<-2>;
    m_fwdICT[ 3]  = fwdTransformCbCr< 3>;
    m_fwdICT[-3]  = fwdTransformCbCr<-3>;
#if JVET_V0130_INTRA_TMP
	  m_pppTarPatch = NULL;
#endif
  }
}

TrQuant::~TrQuant()
{
  if( m_quant )
  {
    delete m_quant;
    m_quant = nullptr;
  }

#if JVET_V0130_INTRA_TMP
  if (m_pppTarPatch != NULL)
  {
	  for (unsigned int uiDepth = 0; uiDepth < USE_MORE_BLOCKSIZE_DEPTH_MAX; uiDepth++)
	  {
		  unsigned int blkSize = g_uiDepth2Width[uiDepth];

		  unsigned int patchSize = blkSize + TMP_TEMPLATE_SIZE;
		  for (unsigned int uiRow = 0; uiRow < patchSize; uiRow++)
		  {
			  if (m_pppTarPatch[uiDepth][uiRow] != NULL)
			  {
				  delete[]m_pppTarPatch[uiDepth][uiRow]; m_pppTarPatch[uiDepth][uiRow] = NULL;
			  }
		  }
		  if (m_pppTarPatch[uiDepth] != NULL)
		  {
			  delete[]m_pppTarPatch[uiDepth]; m_pppTarPatch[uiDepth] = NULL;
		  }
	  }
	  delete[] m_pppTarPatch;
	  m_pppTarPatch = NULL;
  }
#endif
}

#if ENABLE_SPLIT_PARALLELISM
void TrQuant::copyState( const TrQuant& other )
{
  m_quant->copyState( *other.m_quant );
}
#endif

void TrQuant::xDeQuant(const TransformUnit &tu,
                             CoeffBuf      &dstCoeff,
                       const ComponentID   &compID,
                       const QpParam       &cQP)
{
  m_quant->dequant( tu, dstCoeff, compID, cQP );
}

void TrQuant::init( const Quant* otherQuant,
                    const uint32_t uiMaxTrSize,
                    const bool bUseRDOQ,
                    const bool bUseRDOQTS,
#if T0196_SELECTIVE_RDOQ
                    const bool useSelectiveRDOQ,
#endif
                    const bool bEnc
)
{
  delete m_quant;
  m_quant = nullptr;

  m_quant = new DepQuant(otherQuant, bEnc);

  if( m_quant )
  {
    m_quant->init( uiMaxTrSize, bUseRDOQ, bUseRDOQTS, useSelectiveRDOQ );
  }


#if JVET_V0130_INTRA_TMP
  unsigned int blkSize;
  
  if (m_pppTarPatch == NULL)
  {
	  m_pppTarPatch = new Pel * *[USE_MORE_BLOCKSIZE_DEPTH_MAX];
	  for (unsigned int uiDepth = 0; uiDepth < USE_MORE_BLOCKSIZE_DEPTH_MAX; uiDepth++)
	  {
		  blkSize = g_uiDepth2Width[uiDepth];

		  unsigned int patchSize = blkSize + TMP_TEMPLATE_SIZE;
		  m_pppTarPatch[uiDepth] = new Pel * [patchSize];
		  for (unsigned int uiRow = 0; uiRow < patchSize; uiRow++)
		  {
			  m_pppTarPatch[uiDepth][uiRow] = new Pel[patchSize];
		  }
	  }
}
#endif

#if TU_256
  fastFwdTrans =
  { {
    { fastForwardDCT2_B2, fastForwardDCT2_B4, fastForwardDCT2_B8, fastForwardDCT2_B16, fastForwardDCT2_B32, fastForwardDCT2_B64, fastForwardDCT2_B128, fastForwardDCT2_B256 },
    { nullptr,            fastForwardDCT8_B4, fastForwardDCT8_B8, fastForwardDCT8_B16, fastForwardDCT8_B32, fastForwardDCT8_B64, fastForwardDCT8_B128, fastForwardDCT8_B256 },
    { nullptr,            fastForwardDST7_B4, fastForwardDST7_B8, fastForwardDST7_B16, fastForwardDST7_B32, fastForwardDST7_B64, fastForwardDST7_B128, fastForwardDST7_B256 },
#if JVET_W0103_INTRA_MTS
    { nullptr,            fastForwardDCT5_B4, fastForwardDCT5_B8, fastForwardDCT5_B16, fastForwardDCT5_B32, fastForwardDCT5_B64, fastForwardDCT5_B128, fastForwardDCT5_B256 },
    { nullptr,            fastForwardDST4_B4, fastForwardDST4_B8, fastForwardDST4_B16, fastForwardDST4_B32, fastForwardDST4_B64, fastForwardDST4_B128, fastForwardDST4_B256 },
    { nullptr,            fastForwardDST1_B4, fastForwardDST1_B8, fastForwardDST1_B16, fastForwardDST1_B32, fastForwardDST1_B64, fastForwardDST1_B128, fastForwardDST1_B256 },
    { nullptr,            fastForwardIDTR_B4, fastForwardIDTR_B8, fastForwardIDTR_B16, fastForwardIDTR_B32, fastForwardIDTR_B64, fastForwardIDTR_B128, fastForwardIDTR_B256 },
#endif
  } };

  fastInvTrans =
  { {
    { fastInverseDCT2_B2, fastInverseDCT2_B4, fastInverseDCT2_B8, fastInverseDCT2_B16, fastInverseDCT2_B32, fastInverseDCT2_B64, fastInverseDCT2_B128, fastInverseDCT2_B256 },
    { nullptr,            fastInverseDCT8_B4, fastInverseDCT8_B8, fastInverseDCT8_B16, fastInverseDCT8_B32, fastInverseDCT8_B64, fastInverseDCT8_B128, fastInverseDCT8_B256 },
    { nullptr,            fastInverseDST7_B4, fastInverseDST7_B8, fastInverseDST7_B16, fastInverseDST7_B32, fastInverseDST7_B64, fastInverseDST7_B128, fastInverseDST7_B256 },
#if JVET_W0103_INTRA_MTS
    { nullptr,            fastInverseDCT5_B4, fastInverseDCT5_B8, fastInverseDCT5_B16, fastInverseDCT5_B32, fastInverseDCT5_B64, fastInverseDCT5_B128, fastInverseDCT5_B256 },
    { nullptr,            fastInverseDST4_B4, fastInverseDST4_B8, fastInverseDST4_B16, fastInverseDST4_B32, fastInverseDST4_B64, fastInverseDST4_B128, fastInverseDST4_B256 },
    { nullptr,            fastInverseDST1_B4, fastInverseDST1_B8, fastInverseDST1_B16, fastInverseDST1_B32, fastInverseDST1_B64, fastInverseDST1_B128, fastInverseDST1_B256 },
    { nullptr,            fastInverseIDTR_B4, fastInverseIDTR_B8, fastInverseIDTR_B16, fastInverseIDTR_B32, fastInverseIDTR_B64, fastInverseIDTR_B128, fastInverseIDTR_B256 },
#endif
  } };
#else
  fastFwdTrans =
  { {
    { fastForwardDCT2_B2, fastForwardDCT2_B4, fastForwardDCT2_B8, fastForwardDCT2_B16, fastForwardDCT2_B32, fastForwardDCT2_B64 },
    { nullptr,            fastForwardDCT8_B4, fastForwardDCT8_B8, fastForwardDCT8_B16, fastForwardDCT8_B32, nullptr },
    { nullptr,            fastForwardDST7_B4, fastForwardDST7_B8, fastForwardDST7_B16, fastForwardDST7_B32, nullptr },
#if JVET_W0103_INTRA_MTS
    { nullptr,            fastForwardDCT5_B4, fastForwardDCT5_B8, fastForwardDCT5_B16, fastForwardDCT5_B32, nullptr },
    { nullptr,            fastForwardDST4_B4, fastForwardDST4_B8, fastForwardDST4_B16, fastForwardDST4_B32, nullptr },
    { nullptr,            fastForwardDST1_B4, fastForwardDST1_B8, fastForwardDST1_B16, fastForwardDST1_B32, nullptr },
    { nullptr,            fastForwardIDTR_B4, fastForwardIDTR_B8, fastForwardIDTR_B16, fastForwardIDTR_B32, nullptr },
#endif
  } };

  fastInvTrans =
  { {
    { fastInverseDCT2_B2, fastInverseDCT2_B4, fastInverseDCT2_B8, fastInverseDCT2_B16, fastInverseDCT2_B32, fastInverseDCT2_B64 },
    { nullptr,            fastInverseDCT8_B4, fastInverseDCT8_B8, fastInverseDCT8_B16, fastInverseDCT8_B32, nullptr },
    { nullptr,            fastInverseDST7_B4, fastInverseDST7_B8, fastInverseDST7_B16, fastInverseDST7_B32, nullptr },
#if JVET_W0103_INTRA_MTS
    { nullptr,            fastInverseDCT5_B4, fastInverseDCT5_B8, fastInverseDCT5_B16, fastInverseDCT5_B32, nullptr },
    { nullptr,            fastInverseDST4_B4, fastInverseDST4_B8, fastInverseDST4_B16, fastInverseDST4_B32, nullptr },
    { nullptr,            fastInverseDST1_B4, fastInverseDST1_B8, fastInverseDST1_B16, fastInverseDST1_B32, nullptr },
    { nullptr,            fastInverseIDTR_B4, fastInverseIDTR_B8, fastInverseIDTR_B16, fastInverseIDTR_B32, nullptr },
#endif
  } };
#endif

#if ENABLE_SIMD_SIGN_PREDICTION
  m_computeSAD = xComputeSAD;
#endif
#if JVET_V0130_INTRA_TMP
  m_calcTemplateDiff = calcTemplateDiff;
#endif

#if ENABLE_SIMD_SIGN_PREDICTION || TRANSFORM_SIMD_OPT
#ifdef TARGET_SIMD_X86
  initTrQuantX86();
#endif
#endif
}

#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
void TrQuant::fwdLfnstNxN( TCoeff* src, TCoeff* dst, const uint32_t mode, const uint32_t index, const uint32_t size, int zeroOutSize )
#else
void TrQuant::fwdLfnstNxN( int* src, int* dst, const uint32_t mode, const uint32_t index, const uint32_t size, int zeroOutSize )
#endif
{
#if JVET_W0119_LFNST_EXTENSION
  const int8_t* trMat  = ( size > 8 ) ? g_lfnst16x16[ mode ][ index ][ 0 ] : ( ( size > 4 ) ? g_lfnst8x8[ mode ][ index ][ 0 ] : g_lfnst4x4[ mode ][ index ][ 0 ] );
  const int     trSize = ( size > 8 ) ? L16W_ZO : ( ( size > 4 ) ? L8W_ZO : 16 );
#else
  const int8_t* trMat  = ( size > 4 ) ? g_lfnst8x8[ mode ][ index ][ 0 ] : g_lfnst4x4[ mode ][ index ][ 0 ];
#if EXTENDED_LFNST
  const int     trSize = ( size > 4 ) ? 64 : 16;
#else
  const int     trSize = ( size > 4 ) ? 48 : 16;
#endif
#endif
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  TCoeff           coef;
  TCoeff*          out    = dst;
#else
  int           coef;
  int*          out    = dst;
#endif  
#if EXTENDED_LFNST || JVET_W0119_LFNST_EXTENSION
  assert( index < 4 );
#else
  assert( index < 3 );
#endif

  for( int j = 0; j < zeroOutSize; j++ )
  {
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
    TCoeff*          srcPtr   = src;
#else
    int*          srcPtr   = src;
#endif
    const int8_t* trMatTmp = trMat;
    coef = 0;
    for( int i = 0; i < trSize; i++ )
    {
      coef += *srcPtr++ * *trMatTmp++;
    }
    *out++ = ( coef + 64 ) >> 7;
    trMat += trSize;
  }

  ::memset( out, 0, ( trSize - zeroOutSize ) * sizeof( int ) );
}

#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
void TrQuant::invLfnstNxN( TCoeff* src, TCoeff* dst, const uint32_t mode, const uint32_t index, const uint32_t size, int zeroOutSize, const int maxLog2TrDynamicRange )
{
#else
void TrQuant::invLfnstNxN( int* src, int* dst, const uint32_t mode, const uint32_t index, const uint32_t size, int zeroOutSize )
{
  int             maxLog2TrDynamicRange =  15;
#endif
  const TCoeff    outputMinimum         = -( 1 << maxLog2TrDynamicRange );
  const TCoeff    outputMaximum         =  ( 1 << maxLog2TrDynamicRange ) - 1;

#if JVET_W0119_LFNST_EXTENSION
  const int8_t*   trMat  = ( size > 8 ) ? g_lfnst16x16[ mode ][ index ][ 0 ] : ( ( size > 4 ) ? g_lfnst8x8[ mode ][ index ][ 0 ] : g_lfnst4x4[ mode ][ index ][ 0 ] );
  const int       trSize = ( size > 8 ) ? L16W_ZO : ( ( size > 4 ) ? L8W_ZO : 16 );
#else
  const int8_t*   trMat                 =  ( size > 4 ) ? g_lfnst8x8[ mode ][ index ][ 0 ] : g_lfnst4x4[ mode ][ index ][ 0 ];
#if EXTENDED_LFNST
  const int       trSize                =  ( size > 4 ) ? 64 : 16;
#else
  const int       trSize                =  ( size > 4 ) ? 48 : 16;
#endif
#endif
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  TCoeff          resi;
  TCoeff*         out                   =  dst;
#else
  int             resi;
  int*            out                   =  dst;
#endif  

#if EXTENDED_LFNST || JVET_W0119_LFNST_EXTENSION
  assert( index < 4 );
#else
  assert( index < 3 );
#endif

  for( int j = 0; j < trSize; j++ )
  {
    resi = 0;
    const int8_t* trMatTmp = trMat;
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
    TCoeff*       srcPtr   = src;
#else
    int*          srcPtr   = src;
#endif
    for( int i = 0; i < zeroOutSize; i++ )
    {
      resi += *srcPtr++ * *trMatTmp;
      trMatTmp += trSize;
    }
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
    *out++ = Clip3<TCoeff>( outputMinimum, outputMaximum, ( resi + 64 ) >> 7 );
#else
    *out++ = Clip3( outputMinimum, outputMaximum, ( int ) ( resi + 64 ) >> 7 );
#endif
    trMat++;
  }
}

#if JVET_V0130_INTRA_TMP
void insertNode(int diff, int& iXOffset, int& iYOffset, int& pDiff, int& pX, int& pY, short& pId, unsigned int& setId)
{
	pDiff = diff;
	pX = iXOffset;
	pY = iYOffset;
	pId = setId;
}

void clipMvIntraConstraint(CodingUnit* pcCU, int regionId, int& iHorMin, int& iHorMax, int& iVerMin, int& iVerMax, unsigned int uiTemplateSize, unsigned int uiBlkWidth, unsigned int uiBlkHeight, int iCurrY, int iCurrX, int offsetLCUY, int offsetLCUX)
{
	int searchRangeWidth  = TMP_SEARCH_RANGE_MULT_FACTOR * uiBlkWidth;
	int searchRangeHeight = TMP_SEARCH_RANGE_MULT_FACTOR * uiBlkHeight;
	int iMvShift = 0;
	int iTemplateSize = uiTemplateSize;
	int iBlkWidth = uiBlkWidth;
	int iBlkHeight = uiBlkHeight;
	if (regionId == 0) //above outside LCU
	{
		iHorMax = std::min((iCurrX + searchRangeWidth) << iMvShift, (int)((pcCU->cs->sps->getMaxPicWidthInLumaSamples() - iBlkWidth) << iMvShift));
		iHorMin = std::max((iTemplateSize) << iMvShift, (iCurrX - searchRangeWidth) << iMvShift);

		iVerMax = (iCurrY - iBlkHeight - offsetLCUY) << iMvShift;
		iVerMin = std::max(((iTemplateSize) << iMvShift), ((iCurrY - searchRangeHeight) << iMvShift));

		iHorMin = iHorMin - iCurrX;
		iHorMax = iHorMax - iCurrX;
		iVerMax = iVerMax - iCurrY;
		iVerMin = iVerMin - iCurrY;
	}
	else if (regionId == 1) //left outside LCU
	{
		iHorMax = (iCurrX - offsetLCUX - iBlkWidth) << iMvShift;
		iHorMin = std::max((iTemplateSize) << iMvShift, (iCurrX - searchRangeWidth) << iMvShift);

		iVerMin = std::max((iTemplateSize) << iMvShift, (iCurrY - iBlkHeight - offsetLCUY) << iMvShift);
		iVerMax = (iCurrY) << iMvShift;

		iHorMin = iHorMin - iCurrX;
		iHorMax = iHorMax - iCurrX;
		iVerMax = iVerMax - iCurrY;
		iVerMin = iVerMin - iCurrY;
	}
	else if (regionId == 2) //left outside LCU (can reach the bottom row of LCU)
	{
		iHorMin = std::max((iTemplateSize) << iMvShift, (iCurrX - searchRangeWidth) << iMvShift);
		iHorMax = (iCurrX - offsetLCUX - iBlkWidth) << iMvShift;
		iVerMin = (iCurrY + 1) << iMvShift;
		iVerMax = std::min(pcCU->cs->sps->getMaxPicHeightInLumaSamples() - iBlkHeight, (iCurrY - offsetLCUY + pcCU->cs->sps->getCTUSize() - iBlkHeight) << iMvShift);

		iHorMin = iHorMin - iCurrX;
		iHorMax = iHorMax - iCurrX;
		iVerMax = iVerMax - iCurrY;
		iVerMin = iVerMin - iCurrY;
	}
}

TempLibFast::TempLibFast()
{
}

TempLibFast::~TempLibFast()
{
}

void TempLibFast::initTemplateDiff(unsigned int uiPatchWidth, unsigned int uiPatchHeight, unsigned int uiBlkWidth, unsigned int uiBlkHeight, int bitDepth)
{
	int maxValue = ((1 << bitDepth) >> (INIT_THRESHOULD_SHIFTBITS)) * (uiPatchHeight * uiPatchWidth - uiBlkHeight * uiBlkWidth);
	m_diffMax = maxValue;
	{
		m_pDiff = maxValue;
	}
}

#if JVET_W0069_TMP_BOUNDARY
void TrQuant::getTargetTemplate(CodingUnit* pcCU, unsigned int uiBlkWidth, unsigned int uiBlkHeight, RefTemplateType TempType)
#else
void TrQuant::getTargetTemplate(CodingUnit* pcCU, unsigned int uiBlkWidth, unsigned int uiBlkHeight)
#endif
{
	const ComponentID compID = COMPONENT_Y;
	unsigned int uiPatchWidth = uiBlkWidth + TMP_TEMPLATE_SIZE;
	unsigned int uiPatchHeight = uiBlkHeight + TMP_TEMPLATE_SIZE;
	unsigned int uiTarDepth = floorLog2(std::max(uiBlkHeight, uiBlkWidth)) - 2;
	Pel** tarPatch = m_pppTarPatch[uiTarDepth];
	CompArea area = pcCU->blocks[compID];
	Pel* pCurrStart = pcCU->cs->picture->getRecoBuf(area).buf;
	unsigned int  uiPicStride = pcCU->cs->picture->getRecoBuf(compID).stride;
	unsigned int uiY, uiX;

	//fill template
	//up-left & up 
	Pel* tarTemp;
#if JVET_W0069_TMP_BOUNDARY
	if (TempType == L_SHAPE_TEMPLATE)
	{
#endif
	Pel* pCurrTemp = pCurrStart - TMP_TEMPLATE_SIZE * uiPicStride - TMP_TEMPLATE_SIZE;
	for (uiY = 0; uiY < TMP_TEMPLATE_SIZE; uiY++)
	{
		tarTemp = tarPatch[uiY]; 
		for (uiX = 0; uiX < uiPatchWidth; uiX++)
		{
			tarTemp[uiX] = pCurrTemp[uiX];
		}
		pCurrTemp += uiPicStride;
	}
	//left
	for (uiY = TMP_TEMPLATE_SIZE; uiY < uiPatchHeight; uiY++)
	{
		tarTemp = tarPatch[uiY];
		for (uiX = 0; uiX < TMP_TEMPLATE_SIZE; uiX++)
		{
			tarTemp[uiX] = pCurrTemp[uiX];
		}
		pCurrTemp += uiPicStride;
	}
#if JVET_W0069_TMP_BOUNDARY
	}
  else if (TempType == ABOVE_TEMPLATE)
  {
    Pel* pCurrTemp = pCurrStart - TMP_TEMPLATE_SIZE * uiPicStride;
    for (uiY = 0; uiY < TMP_TEMPLATE_SIZE; uiY++)
    {
      tarTemp = tarPatch[uiY];
      for (uiX = 0; uiX < uiBlkWidth; uiX++)
      {
        tarTemp[uiX] = pCurrTemp[uiX];
      }
      pCurrTemp += uiPicStride;
    }
  }
	else if (TempType == LEFT_TEMPLATE)
	{
		Pel* pCurrTemp = pCurrStart - TMP_TEMPLATE_SIZE;
		for (uiY = TMP_TEMPLATE_SIZE; uiY < uiPatchHeight; uiY++)
		{
			tarTemp = tarPatch[uiY];
			for (uiX = 0; uiX < TMP_TEMPLATE_SIZE; uiX++)
			{
				tarTemp[uiX] = pCurrTemp[uiX];
			}
			pCurrTemp += uiPicStride;
		}
	}
#endif
}

#if JVET_W0069_TMP_BOUNDARY
void TrQuant::candidateSearchIntra(CodingUnit* pcCU, unsigned int uiBlkWidth, unsigned int uiBlkHeight, RefTemplateType TempType)
#else
void TrQuant::candidateSearchIntra(CodingUnit* pcCU, unsigned int uiBlkWidth, unsigned int uiBlkHeight)
#endif
{
	const ComponentID compID = COMPONENT_Y;
	const int channelBitDepth = pcCU->cs->sps->getBitDepth(toChannelType(compID));
	unsigned int uiPatchWidth = uiBlkWidth + TMP_TEMPLATE_SIZE;
	unsigned int uiPatchHeight = uiBlkHeight + TMP_TEMPLATE_SIZE;
	unsigned int uiTarDepth = floorLog2(std::max(uiBlkWidth, uiBlkHeight)) - 2;
	Pel** tarPatch = getTargetPatch(uiTarDepth);
	//Initialize the library for saving the best candidates
	m_tempLibFast.initTemplateDiff(uiPatchWidth, uiPatchHeight, uiBlkWidth, uiBlkHeight, channelBitDepth);
	short setId = 0; //record the reference picture.
#if JVET_W0069_TMP_BOUNDARY
	searchCandidateFromOnePicIntra(pcCU, tarPatch, uiPatchWidth, uiPatchHeight, setId, TempType);
#else
	searchCandidateFromOnePicIntra(pcCU, tarPatch, uiPatchWidth, uiPatchHeight, setId);
#endif
	//count collected candidate number
	int pDiff = m_tempLibFast.getDiff();
	int maxDiff = m_tempLibFast.getDiffMax();
	

  if( pDiff < maxDiff )
  {
    m_uiVaildCandiNum = 1;
  }
  else
  {
    m_uiVaildCandiNum = 0;
  }
}

#if JVET_W0069_TMP_BOUNDARY
void  TrQuant::searchCandidateFromOnePicIntra(CodingUnit* pcCU, Pel** tarPatch, unsigned int uiPatchWidth, unsigned int uiPatchHeight, unsigned int setId, RefTemplateType TempType)
#else
void  TrQuant::searchCandidateFromOnePicIntra(CodingUnit* pcCU, Pel** tarPatch, unsigned int uiPatchWidth, unsigned int uiPatchHeight, unsigned int setId)
#endif
{
	const ComponentID compID = COMPONENT_Y;
	unsigned int uiBlkWidth = uiPatchWidth - TMP_TEMPLATE_SIZE;
	unsigned int uiBlkHeight = uiPatchHeight - TMP_TEMPLATE_SIZE;

	int pX = m_tempLibFast.getX();
	int pY = m_tempLibFast.getY();
	int pDiff = m_tempLibFast.getDiff();
	short pId = m_tempLibFast.getId();
	CompArea area = pcCU->blocks[compID];
	int  refStride = pcCU->cs->picture->getRecoBuf(compID).stride;
	
	Pel* ref = pcCU->cs->picture->getRecoBuf(area).buf;
	
	setRefPicUsed(ref); //facilitate the access of each candidate point 
	setStride(refStride);
	
	Mv cTmpMvPred;
	cTmpMvPred.setZero();

	unsigned int uiCUPelY = area.pos().y;
	unsigned int uiCUPelX = area.pos().x;
	int blkX = 0;
	int blkY = 0;
	int iCurrY = uiCUPelY + blkY;
	int iCurrX = uiCUPelX + blkX;

	Position  ctuRsAddr = CU::getCtuXYAddr(*pcCU);
	int offsetLCUY = iCurrY - ctuRsAddr.y;
	int offsetLCUX = iCurrX - ctuRsAddr.x;

	int iYOffset, iXOffset;
	int diff;
	Pel* refCurr;

	const int regionNum = 3;
	int mvYMins[regionNum];
	int mvYMaxs[regionNum];
	int mvXMins[regionNum];
	int mvXMaxs[regionNum];
	int regionId = 0;

	//1. check the near pixels within LCU
	//above pixels in LCU
	int iTemplateSize = TMP_TEMPLATE_SIZE;
	int iBlkWidth = uiBlkWidth;
	int iBlkHeight = uiBlkHeight;
	regionId = 0;
	int iMvShift = 0;

	int iVerMin = std::max(((iTemplateSize) << iMvShift), (iCurrY - offsetLCUY - iBlkHeight + 1) << iMvShift);
	int iVerMax = (iCurrY - iBlkHeight) << iMvShift; 
	int iHorMin = std::max((iTemplateSize) << iMvShift, (iCurrX - offsetLCUX - iBlkWidth + 1) << iMvShift);
	int iHorMax = (iCurrX - iBlkWidth);

	mvXMins[regionId] = iHorMin - iCurrX;
	mvXMaxs[regionId] = iHorMax - iCurrX;
	mvYMins[regionId] = iVerMin - iCurrY;
	mvYMaxs[regionId] = iVerMax - iCurrY;

	//check within CTU pixels
	for (regionId = 0; regionId < 1; regionId++)
	{
		int mvYMin = mvYMins[regionId];
		int mvYMax = mvYMaxs[regionId];
		int mvXMin = mvXMins[regionId];
		int mvXMax = mvXMaxs[regionId];
		if (mvYMax < mvYMin || mvXMax < mvXMin)
		{
			continue;
		}

		for (iYOffset = mvYMax; iYOffset >= mvYMin; iYOffset--)
		{
			for (iXOffset = mvXMax; iXOffset >= mvXMin; iXOffset--)
			{
				refCurr = ref + iYOffset * refStride + iXOffset;
#if JVET_W0069_TMP_BOUNDARY
				diff = m_calcTemplateDiff(refCurr, refStride, tarPatch, uiPatchWidth, uiPatchHeight, pDiff, TempType);
#else
				diff = m_calcTemplateDiff(refCurr, refStride, tarPatch, uiPatchWidth, uiPatchHeight, pDiff);
#endif
				if (diff < (pDiff))
				{
					insertNode(diff, iXOffset, iYOffset, pDiff, pX, pY, pId, setId); 
				}
        if (pDiff == 0)
        {
          regionId++;
        }
			}
		}
	}

	//2. check the pixels outside CTU
	for (regionId = 0; regionId < regionNum; regionId++)
	{// this function fills in the range the template matching for pixels outside the current CTU
		clipMvIntraConstraint(pcCU, regionId, mvXMins[regionId], mvXMaxs[regionId], mvYMins[regionId], mvYMaxs[regionId], TMP_TEMPLATE_SIZE, uiBlkWidth, uiBlkHeight, iCurrY, iCurrX, offsetLCUY, offsetLCUX);
	}

	for (regionId = 0; regionId < regionNum; regionId++)
	{
		int mvYMin = mvYMins[regionId];
		int mvYMax = mvYMaxs[regionId];
		int mvXMin = mvXMins[regionId];
		int mvXMax = mvXMaxs[regionId];
		if ( mvYMax < mvYMin || mvXMax < mvXMin )
		{
			continue;
		}
		for (iYOffset = mvYMax; iYOffset >= mvYMin; iYOffset--)
		{
			for (iXOffset = mvXMax; iXOffset >= mvXMin; iXOffset--)
			{
				refCurr = ref + iYOffset * refStride + iXOffset;
#if JVET_W0069_TMP_BOUNDARY
				diff = m_calcTemplateDiff(refCurr, refStride, tarPatch, uiPatchWidth, uiPatchHeight, pDiff, TempType);
#else
				diff = m_calcTemplateDiff(refCurr, refStride, tarPatch, uiPatchWidth, uiPatchHeight, pDiff);
#endif

				if (diff < (pDiff))
				{
					insertNode(diff, iXOffset, iYOffset, pDiff, pX, pY, pId, setId);
				}

        if (pDiff == 0)
        {
          regionId = regionNum;
        }
			}
		}
	}

	m_tempLibFast.m_pX = pX;
	m_tempLibFast.m_pY = pY;
	m_tempLibFast.m_pDiff = pDiff;
	m_tempLibFast.m_pId = pId;
}
bool TrQuant::generateTMPrediction(Pel* piPred, unsigned int uiStride, unsigned int uiBlkWidth, unsigned int uiBlkHeight, int& foundCandiNum)
{
	bool bSucceedFlag = true;
	unsigned int uiPatchWidth = uiBlkWidth + TMP_TEMPLATE_SIZE;
	unsigned int uiPatchHeight = uiBlkHeight + TMP_TEMPLATE_SIZE;

	foundCandiNum = m_uiVaildCandiNum;
	if (foundCandiNum < 1)
	{
		return false;
	}

  //printf("\n Kalyan Template is matched");

	int pX = m_tempLibFast.getX();
	int pY = m_tempLibFast.getY();
	Pel* ref;
	int picStride = getStride();
	int iOffsetY, iOffsetX;
	Pel* refTarget;
	unsigned int uiHeight = uiPatchHeight - TMP_TEMPLATE_SIZE;
	unsigned int uiWidth = uiPatchWidth - TMP_TEMPLATE_SIZE;

	//the data center: we use the prediction block as the center now.
	//collect the candidates
	ref = getRefPicUsed();
	{
		iOffsetY = pY;
		iOffsetX = pX;
		refTarget = ref + iOffsetY * picStride + iOffsetX;
		for (unsigned int uiY = 0; uiY < uiHeight; uiY++)
		{
			for (unsigned int uiX = 0; uiX < uiWidth; uiX++)
			{
				piPred[uiX] = refTarget[uiX];
			}
			refTarget += picStride;
			piPred += uiStride;
		}
	}
	return bSucceedFlag;
}

#if JVET_W0069_TMP_BOUNDARY
bool TrQuant::generateTM_DC_Prediction(Pel* piPred, unsigned int uiStride, unsigned int uiBlkWidth, unsigned int uiBlkHeight, int DC_Val)
{
  bool bSucceedFlag = true;
  {
    for (unsigned int uiY = 0; uiY < uiBlkHeight; uiY++)
    {
      for (unsigned int uiX = 0; uiX < uiBlkWidth; uiX++)
      {
        piPred[uiX] = DC_Val;
      }
      piPred += uiStride;
    }
  }
	return bSucceedFlag;
}
#endif

#if JVET_W0069_TMP_BOUNDARY
int TrQuant::calcTemplateDiff(Pel* ref, unsigned int uiStride, Pel** tarPatch, unsigned int uiPatchWidth, unsigned int uiPatchHeight, int iMax, RefTemplateType TempType)
#else
int TrQuant::calcTemplateDiff(Pel* ref, unsigned int uiStride, Pel** tarPatch, unsigned int uiPatchWidth, unsigned int uiPatchHeight, int iMax)
#endif
{
  int iDiffSum = 0;
#if JVET_W0069_TMP_BOUNDARY
  Pel* refPatchRow;
  if( TempType == L_SHAPE_TEMPLATE )
  {
    refPatchRow = ref - TMP_TEMPLATE_SIZE * uiStride - TMP_TEMPLATE_SIZE;
  }
  else if( TempType == LEFT_TEMPLATE )
  {
    refPatchRow = ref - TMP_TEMPLATE_SIZE;
  }
  else if( TempType == ABOVE_TEMPLATE )
  {
    refPatchRow = ref - TMP_TEMPLATE_SIZE * uiStride;
  }
#else
  Pel* refPatchRow = ref - TMP_TEMPLATE_SIZE * uiStride - TMP_TEMPLATE_SIZE;
#endif
  Pel* tarPatchRow;

#if JVET_W0069_TMP_BOUNDARY
  if( TempType == L_SHAPE_TEMPLATE )
  {
#endif
  // horizontal difference
  for( int iY = 0; iY < TMP_TEMPLATE_SIZE; iY++ )
  {
    tarPatchRow = tarPatch[iY];
    for( int iX = 0; iX < uiPatchWidth; iX++ )
    {
      iDiffSum += abs( refPatchRow[iX] - tarPatchRow[iX] );
    }
    if( iDiffSum > iMax ) //for speeding up
    {
      return iDiffSum;
    }
    refPatchRow += uiStride;
  }

  // vertical difference
  for( int iY = TMP_TEMPLATE_SIZE; iY < uiPatchHeight; iY++ )
  {
    tarPatchRow = tarPatch[iY];
    for( int iX = 0; iX < TMP_TEMPLATE_SIZE; iX++ )
    {
      iDiffSum += abs( refPatchRow[iX] - tarPatchRow[iX] );
    }
    if( iDiffSum > iMax ) //for speeding up
    {
      return iDiffSum;
    }
    refPatchRow += uiStride;
  }
#if JVET_W0069_TMP_BOUNDARY
	}
  else if (TempType == ABOVE_TEMPLATE)
  {
    // top  template difference
    for (int iY = 0; iY < TMP_TEMPLATE_SIZE; iY++)
    {
      tarPatchRow = tarPatch[iY];
      for (int iX = 0; iX < uiPatchWidth - TMP_TEMPLATE_SIZE; iX++)
      {
        iDiffSum += abs(refPatchRow[iX] - tarPatchRow[iX]);
      }
      if (iDiffSum > iMax) //for speeding up
      {
        return iDiffSum;
      }
      refPatchRow += uiStride;
    }
  }
  else if (TempType == LEFT_TEMPLATE)
  {
	  // left template difference
	  for (int iY = TMP_TEMPLATE_SIZE; iY < uiPatchHeight; iY++)
	  {
		  tarPatchRow = tarPatch[iY];
		  for (int iX = 0; iX < TMP_TEMPLATE_SIZE; iX++)
		  {
			  iDiffSum += abs(refPatchRow[iX] - tarPatchRow[iX]);
		  }
		  if (iDiffSum > iMax) //for speeding up
		  {
			  return iDiffSum;
		  }
		  refPatchRow += uiStride;
	  }
  }
#endif

  return iDiffSum;
}
#endif

uint32_t TrQuant::getLFNSTIntraMode( int wideAngPredMode )
{
  uint32_t intraMode;

  if( wideAngPredMode < 0 )
  {
    intraMode = ( uint32_t ) ( wideAngPredMode + ( NUM_EXT_LUMA_MODE >> 1 ) + NUM_LUMA_MODE );
  }
  else if( wideAngPredMode >= NUM_LUMA_MODE )
  {
    intraMode = ( uint32_t ) ( wideAngPredMode + ( NUM_EXT_LUMA_MODE >> 1 ) );
  }
  else
  {
    intraMode = ( uint32_t ) wideAngPredMode;
  }

  return intraMode;
}

bool TrQuant::getTransposeFlag( uint32_t intraMode )
{
  return ( ( intraMode >= NUM_LUMA_MODE ) && ( intraMode >= ( NUM_LUMA_MODE + ( NUM_EXT_LUMA_MODE >> 1 ) ) ) ) ||
         ( ( intraMode <  NUM_LUMA_MODE ) && ( intraMode >  DIA_IDX ) );
}

void TrQuant::xInvLfnst( const TransformUnit &tu, const ComponentID compID )
{
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  const int maxLog2TrDynamicRange = tu.cs->sps->getMaxLog2TrDynamicRange(toChannelType(compID));
#endif
  const CompArea& area     = tu.blocks[ compID ];
  const uint32_t  width    = area.width;
  const uint32_t  height   = area.height;
  const uint32_t  lfnstIdx = tu.cu->lfnstIdx;
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  if( lfnstIdx && tu.mtsIdx[compID] != MTS_SKIP && (tu.cu->isSepTree() ? true : isLuma(compID)) )
#else
  if (lfnstIdx && tu.mtsIdx[compID] != MTS_SKIP && (CS::isDualITree(*tu.cs) ? true : isLuma(compID)))
#endif
  {
#if JVET_W0119_LFNST_EXTENSION
    const bool whge4         = PU::getUseLFNST16( width, height );
    const bool whge3         = PU::getUseLFNST8 ( width, height );
    int widthIdx             = gp_sizeIdxInfo->idxFrom(  width );
    int heightIdx            = gp_sizeIdxInfo->idxFrom( height );
    const ScanElement * scan = whge4 ? g_coefTopLeftDiagScan16x16[ widthIdx ] : ( whge3 ? g_coefTopLeftDiagScan8x8[ widthIdx ] : g_scanOrder[ SCAN_GROUPED_4x4 ][ SCAN_DIAG ][ widthIdx ][ heightIdx ] );
#else
    const bool whge3 = width >= 8 && height >= 8;
    const ScanElement * scan = whge3 ? g_coefTopLeftDiagScan8x8[ gp_sizeIdxInfo->idxFrom( width ) ] : g_scanOrder[ SCAN_GROUPED_4x4 ][ SCAN_DIAG ][ gp_sizeIdxInfo->idxFrom( width ) ][ gp_sizeIdxInfo->idxFrom( height ) ];
#endif
    uint32_t intraMode = PU::getFinalIntraMode( *tu.cs->getPU( area.pos(), toChannelType( compID ) ), toChannelType( compID ) );

#if JVET_W0123_TIMD_FUSION
    if( compID != COMPONENT_Y && PU::isLMCMode( tu.cs->getPU( area.pos(), toChannelType( compID ) )->intraDir[ toChannelType( compID ) ] ) )
#else
    if( PU::isLMCMode( tu.cs->getPU( area.pos(), toChannelType( compID ) )->intraDir[ toChannelType( compID ) ] ) )
#endif
    {
      intraMode = PU::getCoLocatedIntraLumaMode( *tu.cs->getPU( area.pos(), toChannelType( compID ) ) );
    }
    if (PU::isMIP(*tu.cs->getPU(area.pos(), toChannelType(compID)), toChannelType(compID)))
    {
      intraMode = PLANAR_IDX;
    }
#if JVET_V0130_INTRA_TMP
    if( PU::isTmp( *tu.cs->getPU( area.pos(), toChannelType( compID ) ), toChannelType( compID ) ) )
    {
      intraMode = PLANAR_IDX;
    }
#endif
#if JVET_W0123_TIMD_FUSION
    if (tu.cu->timd && compID == COMPONENT_Y)
    {
      intraMode = MAP131TO67(intraMode);
    }
#endif
    CHECK( intraMode >= NUM_INTRA_MODE - 1, "Invalid intra mode" );

#if EXTENDED_LFNST || JVET_W0119_LFNST_EXTENSION
    if (lfnstIdx < 4)
#else
    if( lfnstIdx < 3 )
#endif
    {
      intraMode = getLFNSTIntraMode( PU::getWideAngle( tu, intraMode, compID ) );
#if RExt__DECODER_DEBUG_TOOL_STATISTICS
      CodingStatistics::IncrementStatisticTool( CodingStatisticsClassType { STATS__TOOL_LFNST, width, height, compID } );
#endif
      bool          transposeFlag   = getTransposeFlag( intraMode );
#if JVET_W0119_LFNST_EXTENSION
      const int     sbSize          = whge4 ? 16 : ( whge3 ? 8 : 4 );
#else
      const int     sbSize          = whge3 ? 8 : 4;
#endif
#if !EXTENDED_LFNST && !JVET_W0119_LFNST_EXTENSION
      bool          tu4x4Flag       = ( width == 4 && height == 4 );
      bool          tu8x8Flag       = ( width == 8 && height == 8 );
#endif
      TCoeff*       lfnstTemp;
      TCoeff*       coeffTemp;
      int           y;
      lfnstTemp   = m_tempInMatrix;   // inverse spectral rearrangement
      coeffTemp   = m_tempCoeff;
      TCoeff *dst = lfnstTemp;

      const ScanElement *scanPtr = scan;
#if JVET_W0119_LFNST_EXTENSION
      int numLfnstCoeff = whge4 ? L16H : ( whge3 ? L8H : 16 );
      for( y = 0; y < numLfnstCoeff; y++ )
#else
#if EXTENDED_LFNST
      const int nSamples = sbSize * sbSize;
      for( y = 0; y < nSamples; y++ )
#else
      for (y = 0; y < 16; y++)
#endif
#endif
      {
        *dst++ = coeffTemp[ scanPtr->idx ];
        scanPtr++;
      }
#if JVET_W0119_LFNST_EXTENSION
      int zeroOutSize = PU::getLFNSTMatrixDim( width, height );
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
      invLfnstNxN( m_tempInMatrix, m_tempOutMatrix, g_lfnstLut[ intraMode ], lfnstIdx - 1, sbSize, zeroOutSize, maxLog2TrDynamicRange );
#else
      invLfnstNxN( m_tempInMatrix, m_tempOutMatrix, g_lfnstLut[ intraMode ], lfnstIdx - 1, sbSize, zeroOutSize );
#endif
#else
#if EXTENDED_LFNST
      const int trSize = whge3 ? 64 : 16;
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
      invLfnstNxN( m_tempInMatrix, m_tempOutMatrix, g_lfnstLut[ intraMode ], lfnstIdx - 1, sbSize, trSize, maxLog2TrDynamicRange );
#else
      invLfnstNxN( m_tempInMatrix, m_tempOutMatrix, g_lfnstLut[ intraMode ], lfnstIdx - 1, sbSize, trSize );
#endif          
#else
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
      invLfnstNxN( m_tempInMatrix, m_tempOutMatrix, g_lfnstLut[ intraMode ], lfnstIdx - 1, sbSize, ( tu4x4Flag || tu8x8Flag ) ? 8 : 16, maxLog2TrDynamicRange );
#else
      invLfnstNxN( m_tempInMatrix, m_tempOutMatrix, g_lfnstLut[ intraMode ], lfnstIdx - 1, sbSize, ( tu4x4Flag || tu8x8Flag ) ? 8 : 16 );
#endif          
#endif
#endif
          lfnstTemp = m_tempOutMatrix; // inverse spectral rearrangement

      if (transposeFlag)
      {
        if (sbSize == 4)
        {
          for (y = 0; y < 4; y++)
          {
            coeffTemp[0] = lfnstTemp[0];
            coeffTemp[1] = lfnstTemp[4];
            coeffTemp[2] = lfnstTemp[8];
            coeffTemp[3] = lfnstTemp[12];
            lfnstTemp++;
            coeffTemp += width;
          }
        }
#if JVET_W0119_LFNST_EXTENSION
        else if( sbSize == 8 )
#else
        else   // ( sbSize == 8 )
#endif
        {
          for (y = 0; y < 8; y++)
          {
#if EXTENDED_LFNST || JVET_W0119_LFNST_EXTENSION
            coeffTemp[0] = lfnstTemp[0];
            coeffTemp[1] = lfnstTemp[8];
            coeffTemp[2] = lfnstTemp[16];
            coeffTemp[3] = lfnstTemp[24];
            coeffTemp[4] = lfnstTemp[32];
            coeffTemp[5] = lfnstTemp[40];
            coeffTemp[6] = lfnstTemp[48];
            coeffTemp[7] = lfnstTemp[56];
#else
            coeffTemp[0] = lfnstTemp[0];
            coeffTemp[1] = lfnstTemp[8];
            coeffTemp[2] = lfnstTemp[16];
            coeffTemp[3] = lfnstTemp[24];
            if (y < 4)
            {
              coeffTemp[4] = lfnstTemp[32];
              coeffTemp[5] = lfnstTemp[36];
              coeffTemp[6] = lfnstTemp[40];
              coeffTemp[7] = lfnstTemp[44];
            }
#endif
            lfnstTemp++;
            coeffTemp += width;
          }
        }
#if JVET_W0119_LFNST_EXTENSION
        else // (sbSize == 16)
        {
          for( y = 0; y < 12; y++ )
          {
            coeffTemp[ 0 ] = lfnstTemp[  0 ];  coeffTemp[ 1 ] = lfnstTemp[ 12 ];
            coeffTemp[ 2 ] = lfnstTemp[ 24 ];  coeffTemp[ 3 ] = lfnstTemp[ 36 ];

            if( y < 8 )
            {
              coeffTemp[ 4 ] = lfnstTemp[ 48 ];  coeffTemp[ 5 ] = lfnstTemp[ 56 ];
              coeffTemp[ 6 ] = lfnstTemp[ 64 ];  coeffTemp[ 7 ] = lfnstTemp[ 72 ];
            }

            if( y < 4 )
            {
              coeffTemp[  8 ] = lfnstTemp[ 80 ];  coeffTemp[  9 ] = lfnstTemp[ 84 ];
              coeffTemp[ 10 ] = lfnstTemp[ 88 ];  coeffTemp[ 11 ] = lfnstTemp[ 92 ];
            }
            lfnstTemp++;
            coeffTemp += width;
          }
        }
#endif
      }
      else
      {
#if JVET_W0119_LFNST_EXTENSION
        if( sbSize == 16 )
        {
          for( y = 0; y < 12; y++ )
          {
            uint32_t uiStride = ( y < 4 ) ? 12 : ( ( y < 8 ) ? 8 : 4 );
            ::memcpy( coeffTemp, lfnstTemp, uiStride * sizeof( TCoeff ) );
            lfnstTemp += uiStride;
            coeffTemp += width;
          }
        }
        else
        {
#endif
        for (y = 0; y < sbSize; y++)
        {
#if EXTENDED_LFNST || JVET_W0119_LFNST_EXTENSION
          uint32_t uiStride = sbSize;
#else
          uint32_t uiStride = (y < 4) ? sbSize : 4;
#endif
          ::memcpy(coeffTemp, lfnstTemp, uiStride * sizeof(TCoeff));
          lfnstTemp += uiStride;
          coeffTemp += width;
        }
#if JVET_W0119_LFNST_EXTENSION
        }
#endif
      }
    }
  }
}

void TrQuant::xFwdLfnst( const TransformUnit &tu, const ComponentID compID, const bool loadTr )
{
  const CompArea& area     = tu.blocks[ compID ];
  const uint32_t  width    = area.width;
  const uint32_t  height   = area.height;
  const uint32_t  lfnstIdx = tu.cu->lfnstIdx;
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  if( lfnstIdx && tu.mtsIdx[compID] != MTS_SKIP && (tu.cu->isSepTree() ? true : isLuma(compID)) )
#else
  if (lfnstIdx && tu.mtsIdx[compID] != MTS_SKIP && (CS::isDualITree(*tu.cs) ? true : isLuma(compID)))
#endif
  {
#if JVET_W0119_LFNST_EXTENSION
    const bool whge4         = PU::getUseLFNST16( width, height );  // width >= 16 && height >= 16;
    const bool whge3         = PU::getUseLFNST8( width, height );
    int widthIdx             = gp_sizeIdxInfo->idxFrom( width );
    int heightIdx            = gp_sizeIdxInfo->idxFrom( height );
    const ScanElement * scan = whge4 ? g_coefTopLeftDiagScan16x16[ widthIdx ] : ( whge3 ? g_coefTopLeftDiagScan8x8[ widthIdx ] : g_scanOrder[ SCAN_GROUPED_4x4 ][ SCAN_DIAG ][ widthIdx ][ heightIdx ] );
#else
    const bool whge3 = width >= 8 && height >= 8;
    const ScanElement * scan = whge3 ? g_coefTopLeftDiagScan8x8[ gp_sizeIdxInfo->idxFrom( width ) ] : g_scanOrder[ SCAN_GROUPED_4x4 ][ SCAN_DIAG ][ gp_sizeIdxInfo->idxFrom( width ) ][ gp_sizeIdxInfo->idxFrom( height ) ];
#endif
    uint32_t intraMode = PU::getFinalIntraMode( *tu.cs->getPU( area.pos(), toChannelType( compID ) ), toChannelType( compID ) );

#if JVET_W0123_TIMD_FUSION
    if( compID != COMPONENT_Y && PU::isLMCMode( tu.cs->getPU( area.pos(), toChannelType( compID ) )->intraDir[ toChannelType( compID ) ] ) )
#else
    if( PU::isLMCMode( tu.cs->getPU( area.pos(), toChannelType( compID ) )->intraDir[ toChannelType( compID ) ] ) )
#endif
    {
      intraMode = PU::getCoLocatedIntraLumaMode( *tu.cs->getPU( area.pos(), toChannelType( compID ) ) );
    }
    if (PU::isMIP(*tu.cs->getPU(area.pos(), toChannelType(compID)), toChannelType(compID)))
    {
      intraMode = PLANAR_IDX;
    }
#if JVET_V0130_INTRA_TMP
    if( PU::isTmp( *tu.cs->getPU( area.pos(), toChannelType( compID ) ), toChannelType( compID ) ) )
    {
      intraMode = PLANAR_IDX;
    }
#endif
#if JVET_W0123_TIMD_FUSION
    if (tu.cu->timd && compID == COMPONENT_Y)
    {
      intraMode = MAP131TO67(intraMode);
    }
#endif
    CHECK( intraMode >= NUM_INTRA_MODE - 1, "Invalid intra mode" );

#if EXTENDED_LFNST || JVET_W0119_LFNST_EXTENSION
    if ( lfnstIdx < 4 )
#else
    if( lfnstIdx < 3 )
#endif
    {
      intraMode = getLFNSTIntraMode( PU::getWideAngle( tu, intraMode, compID ) );

      bool            transposeFlag   = getTransposeFlag( intraMode );
#if JVET_W0119_LFNST_EXTENSION
      const int       sbSize          = whge4 ? 16 : ( whge3 ? 8 : 4 );
#else
      const int       sbSize          = whge3 ? 8 : 4;
#endif
#if !EXTENDED_LFNST && !JVET_W0119_LFNST_EXTENSION
      bool            tu4x4Flag       = ( width == 4 && height == 4 );
      bool            tu8x8Flag       = ( width == 8 && height == 8 );
#endif
      TCoeff*         lfnstTemp;
      TCoeff*         coeffTemp;
      TCoeff *        tempCoeff = loadTr ? m_mtsCoeffs[tu.mtsIdx[compID]] : m_tempCoeff;

      int y;
      lfnstTemp = m_tempInMatrix;   // forward low frequency non-separable transform
      coeffTemp = tempCoeff;

      if (transposeFlag)
      {
        if (sbSize == 4)
        {
          for (y = 0; y < 4; y++)
          {
            lfnstTemp[0]  = coeffTemp[0];
            lfnstTemp[4]  = coeffTemp[1];
            lfnstTemp[8]  = coeffTemp[2];
            lfnstTemp[12] = coeffTemp[3];
            lfnstTemp++;
            coeffTemp += width;
          }
        }
#if JVET_W0119_LFNST_EXTENSION
        else if( sbSize == 8 )
#else
        else   // ( sbSize == 8 )
#endif
        {
          for (y = 0; y < 8; y++)
          {
#if EXTENDED_LFNST || JVET_W0119_LFNST_EXTENSION
            lfnstTemp[  0 ] = coeffTemp[ 0 ];
            lfnstTemp[  8 ] = coeffTemp[ 1 ];
            lfnstTemp[ 16 ] = coeffTemp[ 2 ];
            lfnstTemp[ 24 ] = coeffTemp[ 3 ];
            lfnstTemp[ 32 ] = coeffTemp[ 4 ];
            lfnstTemp[ 40 ] = coeffTemp[ 5 ];
            lfnstTemp[ 48 ] = coeffTemp[ 6 ];
            lfnstTemp[ 56 ] = coeffTemp[ 7 ];
#else
            lfnstTemp[0]  = coeffTemp[0];
            lfnstTemp[8]  = coeffTemp[1];
            lfnstTemp[16] = coeffTemp[2];
            lfnstTemp[24] = coeffTemp[3];
            if (y < 4)
            {
              lfnstTemp[32] = coeffTemp[4];
              lfnstTemp[36] = coeffTemp[5];
              lfnstTemp[40] = coeffTemp[6];
              lfnstTemp[44] = coeffTemp[7];
            }
#endif
            lfnstTemp++;
            coeffTemp += width;
          }
        }
#if JVET_W0119_LFNST_EXTENSION
        else // (sbSize == 16)
        {
          for( y = 0; y < 12; y++ )
          {
            lfnstTemp[  0 ] = coeffTemp[ 0 ]; lfnstTemp[ 12 ] = coeffTemp[ 1 ];
            lfnstTemp[ 24 ] = coeffTemp[ 2 ]; lfnstTemp[ 36 ] = coeffTemp[ 3 ];

            if( y < 8 )
            {
              lfnstTemp[ 48 ] = coeffTemp[ 4 ];  lfnstTemp[ 56 ] = coeffTemp[ 5 ];
              lfnstTemp[ 64 ] = coeffTemp[ 6 ];  lfnstTemp[ 72 ] = coeffTemp[ 7 ];
            }

            if( y < 4 )
            {
              lfnstTemp[ 80 ] = coeffTemp[  8 ];  lfnstTemp[ 84 ] = coeffTemp[  9 ];
              lfnstTemp[ 88 ] = coeffTemp[ 10 ];  lfnstTemp[ 92 ] = coeffTemp[ 11 ];
            }
            lfnstTemp++;
            coeffTemp += width;
          }
        }
#endif
      }
      else
      {
#if JVET_W0119_LFNST_EXTENSION
        if( sbSize == 16 )
        {
          for( y = 0; y < 16; y++ )
          {
            uint32_t uiStride = ( y < 4 ) ? 12 : ( ( y < 8 ) ? 8 : 4 );
            ::memcpy( lfnstTemp, coeffTemp, uiStride * sizeof( TCoeff ) );
            lfnstTemp += uiStride;
            coeffTemp += width;
          }
        }
        else
        {
#endif
        for( y = 0; y < sbSize; y++ )
        {
#if EXTENDED_LFNST || JVET_W0119_LFNST_EXTENSION
          uint32_t uiStride = sbSize;
#else
          uint32_t uiStride = ( y < 4 ) ? sbSize : 4;
#endif
          ::memcpy( lfnstTemp, coeffTemp, uiStride * sizeof( TCoeff ) );
          lfnstTemp += uiStride;
          coeffTemp += width;
        }
#if JVET_W0119_LFNST_EXTENSION
        }
#endif
      }

#if JVET_W0119_LFNST_EXTENSION
      int zeroOutSize = PU::getLFNSTMatrixDim( width, height );

      fwdLfnstNxN( m_tempInMatrix, m_tempOutMatrix, g_lfnstLut[ intraMode ], lfnstIdx - 1, sbSize, zeroOutSize );
#else
#if EXTENDED_LFNST
      const int trSize = whge3 ? 64 : 16;
      fwdLfnstNxN( m_tempInMatrix, m_tempOutMatrix, g_lfnstLut[ intraMode ], lfnstIdx - 1, sbSize, trSize );
#else
      fwdLfnstNxN( m_tempInMatrix, m_tempOutMatrix, g_lfnstLut[ intraMode ], lfnstIdx - 1, sbSize, ( tu4x4Flag || tu8x8Flag ) ? 8 : 16 );
#endif
#endif
      lfnstTemp = m_tempOutMatrix; // forward spectral rearrangement
      coeffTemp = tempCoeff;
          
      const ScanElement *scanPtr = scan;

#if JVET_W0119_LFNST_EXTENSION
      int lfnstCoeffNum = ( sbSize > 8 ) ? L16W_ZO : ( ( sbSize > 4 ) ? L8W_ZO : 16 );
#else
#if EXTENDED_LFNST
      int lfnstCoeffNum = sbSize * sbSize;
#else
      int lfnstCoeffNum = ( sbSize == 4 ) ? sbSize * sbSize : 48;
#endif
#endif

      for (y = 0; y < lfnstCoeffNum; y++)
      {
        coeffTemp[scanPtr->idx] = *lfnstTemp++;
        scanPtr++;
      }
    }
  }
}


void TrQuant::invTransformNxN( TransformUnit &tu, const ComponentID &compID, PelBuf &pResi, const QpParam &cQP )
{
  const CompArea &area    = tu.blocks[compID];
  const uint32_t uiWidth      = area.width;
  const uint32_t uiHeight     = area.height;

  CHECK( uiWidth > tu.cs->sps->getMaxTbSize() || uiHeight > tu.cs->sps->getMaxTbSize(), "Maximal allowed transformation size exceeded!" );
  CoeffBuf tempCoeff = CoeffBuf(m_tempCoeff, area);
  xDeQuant(tu, tempCoeff, compID, cQP);

  DTRACE_COEFF_BUF(D_TCOEFF, tempCoeff, tu, tu.cu->predMode, compID);

  if (tu.cs->sps->getUseLFNST())
  {
    xInvLfnst(tu, compID);
  }

  if (tu.mtsIdx[compID] == MTS_SKIP)
  {
    xITransformSkip(tempCoeff, pResi, tu, compID);
  }
  else
  {
    xIT(tu, compID, tempCoeff, pResi);
  }

  //DTRACE_BLOCK_COEFF(tu.getCoeffs(compID), tu, tu.cu->predMode, compID);
  DTRACE_PEL_BUF( D_RESIDUALS, pResi, tu, tu.cu->predMode, compID);
}

std::pair<int64_t,int64_t> TrQuant::fwdTransformICT( const TransformUnit &tu, const PelBuf &resCb, const PelBuf &resCr, PelBuf &resC1, PelBuf &resC2, int jointCbCr )
{
  CHECK( Size(resCb) != Size(resCr), "resCb and resCr have different sizes" );
  CHECK( Size(resCb) != Size(resC1), "resCb and resC1 have different sizes" );
  CHECK( Size(resCb) != Size(resC2), "resCb and resC2 have different sizes" );
  return (*m_fwdICT[ TU::getICTMode(tu, jointCbCr) ])( resCb, resCr, resC1, resC2 );
}

void TrQuant::invTransformICT( const TransformUnit &tu, PelBuf &resCb, PelBuf &resCr )
{
  CHECK( Size(resCb) != Size(resCr), "resCb and resCr have different sizes" );
  (*m_invICT[ TU::getICTMode(tu) ])( resCb, resCr );
}

std::vector<int> TrQuant::selectICTCandidates( const TransformUnit &tu, CompStorage* resCb, CompStorage* resCr )
{
  CHECK( !resCb[0].valid() || !resCr[0].valid(), "standard components are not valid" );

  if( !CU::isIntra( *tu.cu ) )
  {
    int cbfMask = 3;
    resCb[cbfMask].create( tu.blocks[COMPONENT_Cb] );
    resCr[cbfMask].create( tu.blocks[COMPONENT_Cr] );
    fwdTransformICT( tu, resCb[0], resCr[0], resCb[cbfMask], resCr[cbfMask], cbfMask );
    std::vector<int> cbfMasksToTest;
    cbfMasksToTest.push_back( cbfMask );
    return cbfMasksToTest;
  }

  std::pair<int64_t,int64_t> pairDist[4];
  for( int cbfMask = 0; cbfMask < 4; cbfMask++ )
  {
    if( cbfMask )
    {
      CHECK( resCb[cbfMask].valid() || resCr[cbfMask].valid(), "target components for cbfMask=" << cbfMask << " are already present" );
      resCb[cbfMask].create( tu.blocks[COMPONENT_Cb] );
      resCr[cbfMask].create( tu.blocks[COMPONENT_Cr] );
    }
    pairDist[cbfMask] = fwdTransformICT( tu, resCb[0], resCr[0], resCb[cbfMask], resCr[cbfMask], cbfMask );
  }

  std::vector<int> cbfMasksToTest;
  int64_t minDist1  = std::min<int64_t>( pairDist[0].first, pairDist[0].second );
  int64_t minDist2  = std::numeric_limits<int64_t>::max();
  int     cbfMask1  = 0;
  int     cbfMask2  = 0;
  for( int cbfMask : { 1, 2, 3 } )
  {
    if( pairDist[cbfMask].first < minDist1 )
    {
      cbfMask2  = cbfMask1; minDist2  = minDist1;
      cbfMask1  = cbfMask;  minDist1  = pairDist[cbfMask1].first;
    }
    else if( pairDist[cbfMask].first < minDist2 )
    {
      cbfMask2  = cbfMask;  minDist2  = pairDist[cbfMask2].first;
    }
  }
  if( cbfMask1 )
  {
    cbfMasksToTest.push_back( cbfMask1 );
  }
  if( cbfMask2 && ( ( minDist2 < (9*minDist1)/8 ) || ( !cbfMask1 && minDist2 < (3*minDist1)/2 ) ) )
  {
    cbfMasksToTest.push_back( cbfMask2 );
  }

  return cbfMasksToTest;
}



// ------------------------------------------------------------------------------------------------
// Logical transform
// ------------------------------------------------------------------------------------------------

void TrQuant::getTrTypes(const TransformUnit tu, const ComponentID compID, int &trTypeHor, int &trTypeVer)
{
  const bool isExplicitMTS = (CU::isIntra(*tu.cu) ? tu.cs->sps->getUseIntraMTS() : tu.cs->sps->getUseInterMTS() && CU::isInter(*tu.cu)) && isLuma(compID);
  const bool isImplicitMTS = CU::isIntra(*tu.cu) && tu.cs->sps->getUseImplicitMTS() && isLuma(compID) && tu.cu->lfnstIdx == 0 && tu.cu->mipFlag == 0;
  const bool isISP = CU::isIntra(*tu.cu) && tu.cu->ispMode && isLuma(compID);
  const bool isSBT = CU::isInter(*tu.cu) && tu.cu->sbtInfo && isLuma(compID);

  trTypeHor = DCT2;
  trTypeVer = DCT2;

  if (isISP && tu.cu->lfnstIdx)
  {
    return;
  }

  if (!tu.cs->sps->getUseMTS())
  {
    return;
  }

#if JVET_V0130_INTRA_TMP
  if (isImplicitMTS || isISP || tu.cu->tmpFlag)
#else
  if (isImplicitMTS || isISP)
#endif
  {
    int  width = tu.blocks[compID].width;
    int  height = tu.blocks[compID].height;
    bool widthDstOk = width >= 4 && width <= 16;
    bool heightDstOk = height >= 4 && height <= 16;

    if (widthDstOk)
    {
      trTypeHor = DST7;
    }
    if (heightDstOk)
    {
      trTypeVer = DST7;
    }
    return;
  }


  if (isSBT)
  {
    uint8_t sbtIdx = tu.cu->getSbtIdx();
    uint8_t sbtPos = tu.cu->getSbtPos();

    if( sbtIdx == SBT_VER_HALF || sbtIdx == SBT_VER_QUAD )
    {
      assert( tu.lwidth() <= MTS_INTER_MAX_CU_SIZE );
      if( tu.lheight() > MTS_INTER_MAX_CU_SIZE )
      {
        trTypeHor = trTypeVer = DCT2;
      }
      else
      {
        if (sbtPos == SBT_POS0)
        {
          trTypeHor = DCT8;
          trTypeVer = DST7;
        }
        else
        {
          trTypeHor = DST7;
          trTypeVer = DST7;
        }
      }
    }
    else
    {
      assert( tu.lheight() <= MTS_INTER_MAX_CU_SIZE );
      if( tu.lwidth() > MTS_INTER_MAX_CU_SIZE )
      {
        trTypeHor = trTypeVer = DCT2;
      }
      else
      {
        if (sbtPos == SBT_POS0)
        {
          trTypeHor = DST7;
          trTypeVer = DCT8;
        }
        else
        {
          trTypeHor = DST7;
          trTypeVer = DST7;
        }
      }
    }
    return;
  }

  if (isExplicitMTS)
  {
#if JVET_W0103_INTRA_MTS
    if (tu.mtsIdx[compID] > MTS_SKIP && CU::isIntra(*tu.cu))
    {
      CHECK(compID != COMPONENT_Y, " MTS activated for chroma");
      uint32_t width = tu.blocks[compID].width;
      uint32_t height = tu.blocks[compID].height;
      int TrIdx = (tu.mtsIdx[compID] - MTS_DST7_DST7);
      CHECK(width < 4 || height < 4, "width < 4 || height < 4 for MTS");
      uint8_t nSzIdxW = std::min(3, (floorLog2(width) - 2));
      uint8_t nSzIdxH = std::min(3, (floorLog2(height) - 2));
      const CompArea& area = tu.blocks[compID];
      int predMode = PU::getFinalIntraMode(*tu.cs->getPU(area.pos(), toChannelType(compID)), toChannelType(compID));
#if JVET_W0123_TIMD_FUSION
      if (tu.cu->timd && compID == COMPONENT_Y)
      {
        predMode = MAP131TO67(predMode);
      }
#endif
      int ucMode;
      int nMdIdx;
      bool isTrTransposed = false;
      if (tu.cu->mipFlag) //MIP is treated as planar.
      {
        ucMode = 0;
        nMdIdx = 35;
        isTrTransposed = (tu.cs->getPU(area.pos(), toChannelType(compID)))->mipTransposedFlag;
      }
      else
      {
        ucMode = predMode; //"ucMode" is the signaled Mode.
        predMode = PU::getWideAngle(tu, (uint32_t)predMode, compID);
        CHECK(predMode < -(NUM_EXT_LUMA_MODE >> 1) && predMode >= NUM_LUMA_MODE + (NUM_EXT_LUMA_MODE >> 1), "luma mode out of range");
        predMode = (predMode < 0) ? 2 : (predMode >= NUM_LUMA_MODE) ? 66 : predMode;
        nMdIdx = predMode > DIA_IDX ? (NUM_LUMA_MODE + 1 - predMode) : predMode;
        isTrTransposed = (predMode > DIA_IDX) ? true : false;
      }
      uint8_t nSzIdx = isTrTransposed ? (nSzIdxH * 4 + nSzIdxW) : (nSzIdxW * 4 + nSzIdxH);
      CHECK(nSzIdx >= 16, "nSzIdx >= 16");
      CHECK(nMdIdx >= 36, "nMdIdx >= 36");
      uint8_t nTrSet = g_aucIpmToTrSet[nSzIdx][nMdIdx];
      CHECK(nTrSet >= 80, "nTrSet >= 80");
      trTypeVer = g_aucTrIdxToTr[g_aucTrSet[nTrSet][TrIdx]][predMode > DIA_IDX ? 1 : 0];
      trTypeHor = g_aucTrIdxToTr[g_aucTrSet[nTrSet][TrIdx]][predMode > DIA_IDX ? 0 : 1];
      predMode = ucMode; //to Check IDTR criteria, signaled mode should be used to check the difference
      if (TrIdx == 3 && width <= 16 && height <= 16)
      {
        if (abs(predMode - HOR_IDX) <= g_aiIdLut[floorLog2(width) - 2][floorLog2(height) - 2])
        {
          trTypeVer = IDTR;
        }
        if (abs(predMode - VER_IDX) <= g_aiIdLut[floorLog2(width) - 2][floorLog2(height) - 2])
        {
          trTypeHor = IDTR;
        }
      }
    }
    else
#endif
    if (tu.mtsIdx[compID] > MTS_SKIP)
    {
      int indHor = (tu.mtsIdx[compID] - MTS_DST7_DST7) & 1;
      int indVer = (tu.mtsIdx[compID] - MTS_DST7_DST7) >> 1;
      trTypeHor = indHor ? DCT8 : DST7;
      trTypeVer = indVer ? DCT8 : DST7;
    }
  }
}

void TrQuant::xT( const TransformUnit &tu, const ComponentID &compID, const CPelBuf &resi, CoeffBuf &dstCoeff, const int width, const int height )
{
  const unsigned maxLog2TrDynamicRange  = tu.cs->sps->getMaxLog2TrDynamicRange( toChannelType( compID ) );
  const unsigned bitDepth               = tu.cs->sps->getBitDepth(              toChannelType( compID ) );
  const int      TRANSFORM_MATRIX_SHIFT = g_transformMatrixShift[TRANSFORM_FORWARD];
  const uint32_t transformWidthIndex    = floorLog2(width ) - 1;  // nLog2WidthMinus1, since transform start from 2-point
  const uint32_t transformHeightIndex   = floorLog2(height) - 1;  // nLog2HeightMinus1, since transform start from 2-point

  int trTypeHor = DCT2;
  int trTypeVer = DCT2;

  getTrTypes ( tu, compID, trTypeHor, trTypeVer );
#if TU_256
  int  skipWidth  =  width  > JVET_C0024_ZERO_OUT_TH ? width  - JVET_C0024_ZERO_OUT_TH : 0;
  int  skipHeight =  height > JVET_C0024_ZERO_OUT_TH ? height - JVET_C0024_ZERO_OUT_TH : 0;
#else
  int  skipWidth  = ( trTypeHor != DCT2 && width  == 32 ) ? 16 : width  > JVET_C0024_ZERO_OUT_TH ? width  - JVET_C0024_ZERO_OUT_TH : 0;
  int  skipHeight = ( trTypeVer != DCT2 && height == 32 ) ? 16 : height > JVET_C0024_ZERO_OUT_TH ? height - JVET_C0024_ZERO_OUT_TH : 0;
#endif

#if EXTENDED_LFNST
  if( tu.cs->sps->getUseLFNST() && tu.cu->lfnstIdx && width >= 4 && height >= 4)
  {
    const bool whge3 = width >= 8 && height >= 8;
    const int lfnst_threshold = whge3 ? 8 : 4;
    skipWidth  = width  - lfnst_threshold;
    skipHeight = height - lfnst_threshold;
  }
#else
  if( tu.cs->sps->getUseLFNST() && tu.cu->lfnstIdx )
  {
    if( (width == 4 && height > 4) || (width > 4 && height == 4) )
    {
      skipWidth  = width  - 4;
      skipHeight = height - 4;
    }
#if JVET_W0119_LFNST_EXTENSION
    else if( width >= 16 && height >= 16 )
    {
      skipWidth  = width  - 16;
      skipHeight = height - 16;
    }
#endif
    else if( (width >= 8 && height >= 8) )
    {
      skipWidth  = width  - 8;
      skipHeight = height - 8;
    }
  }
#endif

#if RExt__DECODER_DEBUG_TOOL_STATISTICS
  if ( trTypeHor != DCT2 )
  {
    CodingStatistics::IncrementStatisticTool( CodingStatisticsClassType{ STATS__TOOL_EMT, uint32_t( width ), uint32_t( height ), compID } );
  }
#endif

  ALIGN_DATA( MEMORY_ALIGN_DEF_SIZE, TCoeff block[MAX_TB_SIZEY * MAX_TB_SIZEY] );

  const Pel *resiBuf    = resi.buf;
  const int  resiStride = resi.stride;

  for( int y = 0; y < height; y++ )
  {
    for( int x = 0; x < width; x++ )
    {
      block[( y * width ) + x] = resiBuf[( y * resiStride ) + x];
    }
  }

  if( width > 1 && height > 1 ) // 2-D transform
  {
    const int      shift_1st              = ((floorLog2(width )) + bitDepth + TRANSFORM_MATRIX_SHIFT) - maxLog2TrDynamicRange + COM16_C806_TRANS_PREC;
    const int      shift_2nd              =  (floorLog2(height))            + TRANSFORM_MATRIX_SHIFT                          + COM16_C806_TRANS_PREC;
    CHECK( shift_1st < 0, "Negative shift" );
    CHECK( shift_2nd < 0, "Negative shift" );
    TCoeff *tmp = (TCoeff *) alloca(width * height * sizeof(TCoeff));

    fastFwdTrans[trTypeHor][transformWidthIndex](block, tmp, shift_1st, height, 0, skipWidth);
    fastFwdTrans[trTypeVer][transformHeightIndex](tmp, dstCoeff.buf, shift_2nd, width, skipWidth, skipHeight);
  }
  else if( height == 1 ) //1-D horizontal transform
  {
    const int      shift              = ((floorLog2(width )) + bitDepth + TRANSFORM_MATRIX_SHIFT) - maxLog2TrDynamicRange + COM16_C806_TRANS_PREC;
    CHECK( shift < 0, "Negative shift" );
    CHECKD( ( transformWidthIndex < 0 ), "There is a problem with the width." );
    fastFwdTrans[trTypeHor][transformWidthIndex]( block, dstCoeff.buf, shift, 1, 0, skipWidth );
  }
  else //if (iWidth == 1) //1-D vertical transform
  {
    int shift = ( ( floorLog2(height) ) + bitDepth + TRANSFORM_MATRIX_SHIFT ) - maxLog2TrDynamicRange + COM16_C806_TRANS_PREC;
    CHECK( shift < 0, "Negative shift" );
    CHECKD( ( transformHeightIndex < 0 ), "There is a problem with the height." );
    fastFwdTrans[trTypeVer][transformHeightIndex]( block, dstCoeff.buf, shift, 1, 0, skipHeight );
  }
}

void TrQuant::xIT( const TransformUnit &tu, const ComponentID &compID, const CCoeffBuf &pCoeff, PelBuf &pResidual )
{
  const int      width                  = pCoeff.width;
  const int      height                 = pCoeff.height;
  const unsigned maxLog2TrDynamicRange  = tu.cs->sps->getMaxLog2TrDynamicRange( toChannelType( compID ) );
  const unsigned bitDepth               = tu.cs->sps->getBitDepth(              toChannelType( compID ) );
  const int      TRANSFORM_MATRIX_SHIFT = g_transformMatrixShift[TRANSFORM_INVERSE];
  const TCoeff   clipMinimum            = -( 1 << maxLog2TrDynamicRange );
  const TCoeff   clipMaximum            =  ( 1 << maxLog2TrDynamicRange ) - 1;
  const uint32_t transformWidthIndex    = floorLog2(width ) - 1;                                // nLog2WidthMinus1, since transform start from 2-point
  const uint32_t transformHeightIndex   = floorLog2(height) - 1;                                // nLog2HeightMinus1, since transform start from 2-point

  int trTypeHor = DCT2;
  int trTypeVer = DCT2;

  getTrTypes ( tu, compID, trTypeHor, trTypeVer );
#if TU_256
  int skipWidth  =  width  > JVET_C0024_ZERO_OUT_TH ? width  - JVET_C0024_ZERO_OUT_TH : 0;
  int skipHeight =  height > JVET_C0024_ZERO_OUT_TH ? height - JVET_C0024_ZERO_OUT_TH : 0;
#else
  int skipWidth  = ( trTypeHor != DCT2 && width  == 32 ) ? 16 : width  > JVET_C0024_ZERO_OUT_TH ? width  - JVET_C0024_ZERO_OUT_TH : 0;
  int skipHeight = ( trTypeVer != DCT2 && height == 32 ) ? 16 : height > JVET_C0024_ZERO_OUT_TH ? height - JVET_C0024_ZERO_OUT_TH : 0;
#endif

#if EXTENDED_LFNST
  if (tu.cs->sps->getUseLFNST() && tu.cu->lfnstIdx && width >= 4 && height >= 4)
  {
    const bool whge3 = width >= 8 && height >= 8;
    const int lfnst_threshold = whge3 ? 8 : 4;
    skipWidth = width - lfnst_threshold;
    skipHeight = height - lfnst_threshold;
  }
#else
  if( tu.cs->sps->getUseLFNST() && tu.cu->lfnstIdx )
  {
    if( (width == 4 && height > 4) || (width > 4 && height == 4) )
    {
      skipWidth  = width  - 4;
      skipHeight = height - 4;
    }
#if JVET_W0119_LFNST_EXTENSION
    else if( ( width >= 16 && height >= 16 ) )
    {
      skipWidth  = width  - 16;
      skipHeight = height - 16;
    }
#endif
    else if( (width >= 8 && height >= 8) )
    {
      skipWidth  = width  - 8;
      skipHeight = height - 8;
    }
  }
#endif

  TCoeff *block = ( TCoeff * ) alloca( width * height * sizeof( TCoeff ) );

  if( width > 1 && height > 1 ) //2-D transform
  {
    const int      shift_1st              =   TRANSFORM_MATRIX_SHIFT + 1 + COM16_C806_TRANS_PREC; // 1 has been added to shift_1st at the expense of shift_2nd
    const int      shift_2nd              = ( TRANSFORM_MATRIX_SHIFT + maxLog2TrDynamicRange - 1 ) - bitDepth + COM16_C806_TRANS_PREC;
    CHECK( shift_1st < 0, "Negative shift" );
    CHECK( shift_2nd < 0, "Negative shift" );
    TCoeff *tmp = ( TCoeff * ) alloca( width * height * sizeof( TCoeff ) );
  fastInvTrans[trTypeVer][transformHeightIndex](pCoeff.buf, tmp, shift_1st, width, skipWidth, skipHeight, clipMinimum, clipMaximum);
  fastInvTrans[trTypeHor][transformWidthIndex] (tmp,      block, shift_2nd, height,         0, skipWidth, clipMinimum, clipMaximum);
  }
  else if( width == 1 ) //1-D vertical transform
  {
    int shift = ( TRANSFORM_MATRIX_SHIFT + maxLog2TrDynamicRange - 1 ) - bitDepth + COM16_C806_TRANS_PREC;
    CHECK( shift < 0, "Negative shift" );
    CHECK( ( transformHeightIndex < 0 ), "There is a problem with the height." );
    fastInvTrans[trTypeVer][transformHeightIndex]( pCoeff.buf, block, shift + 1, 1, 0, skipHeight, clipMinimum, clipMaximum );
  }
  else //if(iHeight == 1) //1-D horizontal transform
  {
    const int      shift              = ( TRANSFORM_MATRIX_SHIFT + maxLog2TrDynamicRange - 1 ) - bitDepth + COM16_C806_TRANS_PREC;
    CHECK( shift < 0, "Negative shift" );
    CHECK( ( transformWidthIndex < 0 ), "There is a problem with the width." );
    fastInvTrans[trTypeHor][transformWidthIndex]( pCoeff.buf, block, shift + 1, 1, 0, skipWidth, clipMinimum, clipMaximum );
  }

  Pel *resiBuf    = pResidual.buf;
  int  resiStride = pResidual.stride;

  for( int y = 0; y < height; y++ )
  {
    for( int x = 0; x < width; x++ )
    {
      resiBuf[( y * resiStride ) + x] = Pel( block[( y * width ) + x] );
    }
  }
}

/** Wrapper function between HM interface and core NxN transform skipping
 */
void TrQuant::xITransformSkip(const CCoeffBuf     &pCoeff,
                                    PelBuf        &pResidual,
                              const TransformUnit &tu,
                              const ComponentID   &compID)
{
  const CompArea &area      = tu.blocks[compID];
  const int width           = area.width;
  const int height          = area.height;

  const TCoeff *coeff = pCoeff.buf;
  for (uint32_t y = 0; y < height; y++)
  {
    for (uint32_t x = 0; x < width; x++)
    {
      pResidual.at(x, y) = coeff[x];
    }
    coeff += pCoeff.stride;
  }
}

void TrQuant::xQuant(TransformUnit &tu, const ComponentID &compID, const CCoeffBuf &pSrc, TCoeff &uiAbsSum, const QpParam &cQP, const Ctx& ctx)
{
  m_quant->quant( tu, compID, pSrc, uiAbsSum, cQP, ctx );
}

void TrQuant::transformNxN( TransformUnit& tu, const ComponentID& compID, const QpParam& cQP, std::vector<TrMode>* trModes, const int maxCand )
{
        CodingStructure &cs = *tu.cs;
  const CompArea &rect      = tu.blocks[compID];
  const uint32_t width      = rect.width;
  const uint32_t height     = rect.height;

  const CPelBuf  resiBuf    = cs.getResiBuf(rect);

  CHECK( cs.sps->getMaxTbSize() < width, "Unsupported transformation size" );

  int pos = 0;
  std::vector<TrCost> trCosts;
  std::vector<TrMode>::iterator it = trModes->begin();
#if TU_256
  const double facBB[] = { 1.2, 1.3, 1.3, 1.4, 1.5, 1.5, 1.5 };
#else
  const double facBB[] = { 1.2, 1.3, 1.3, 1.4, 1.5 };
#endif
  while( it != trModes->end() )
  {
    tu.mtsIdx[compID] = it->first;
    CoeffBuf tempCoeff( m_mtsCoeffs[tu.mtsIdx[compID]], rect);
    if( tu.noResidual )
    {
      int sumAbs = 0;
      trCosts.push_back( TrCost( sumAbs, pos++ ) );
      it++;
      continue;
    }

    if ( tu.mtsIdx[compID] == MTS_SKIP )
    {
      xTransformSkip( tu, compID, resiBuf, tempCoeff.buf );
    }
    else
    {
      xT( tu, compID, resiBuf, tempCoeff, width, height );
    }

#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT_VS
    TCoeff sumAbs = 0;
#else
    int sumAbs = 0;
#endif
    for( int pos = 0; pos < width*height; pos++ )
    {
      sumAbs += abs( tempCoeff.buf[pos] );
    }

    double scaleSAD=1.0;
    if ( tu.mtsIdx[compID] == MTS_SKIP && ((floorLog2(width) + floorLog2(height)) & 1) == 1)
    {
      scaleSAD=1.0/1.414213562; // compensate for not scaling transform skip coefficients by 1/sqrt(2)
    }
    if (tu.mtsIdx[compID] == MTS_SKIP)
    {
      int trShift = getTransformShift(tu.cu->slice->getSPS()->getBitDepth(toChannelType(compID)), rect.size(),
                                      tu.cu->slice->getSPS()->getMaxLog2TrDynamicRange(toChannelType(compID)));
      scaleSAD *= pow(2, trShift);
    }

#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT_VS
    trCosts.push_back( TrCost( int(std::min<double>(sumAbs*scaleSAD, std::numeric_limits<int>::max())), pos++ ) );
#else
    trCosts.push_back( TrCost( int(sumAbs*scaleSAD), pos++ ) );
#endif
    it++;
  }

  int numTests = 0;
  std::vector<TrCost>::iterator itC = trCosts.begin();
  const double fac   = facBB[std::max(0, floorLog2(std::max(width, height)) - 2)];
  const double thr   = fac * trCosts.begin()->first;
  const double thrTS = trCosts.begin()->first;
  while( itC != trCosts.end() )
  {
    const bool testTr = itC->first <= ( itC->second == 1 ? thrTS : thr ) && numTests <= maxCand;
    trModes->at( itC->second ).second = testTr;
    numTests += testTr;
    itC++;
  }
}
#if JVET_W0103_INTRA_MTS
// does transform for MTS candidates and return absSum of unquant Coeffs.
uint64_t TrQuant::transformNxN(TransformUnit& tu)
{
  CHECK(!tu.cu->mtsFlag, "mtsFlag should be on for selection");
  CodingStructure &cs = *tu.cs;
  const CompArea &rect = tu.blocks[COMPONENT_Y];
  const uint32_t uiWidth = rect.width;
  const uint32_t uiHeight = rect.height;

  const CPelBuf resiBuf = cs.getResiBuf(rect);
  CoeffBuf tempCoeff(m_mtsCoeffs[tu.mtsIdx[0]], rect);
  xT(tu, COMPONENT_Y, resiBuf, tempCoeff, uiWidth, uiHeight);


  const TCoeff *dstCoeffBuf = tempCoeff.buf;
  const int  dstCoeffStride = tempCoeff.stride;
  uint64_t coeffAbsSum = 0;

  for (int y = 0; y < uiHeight; y++)
  {
    for (int x = 0; x < uiWidth; x++)
    {
      coeffAbsSum += abs(dstCoeffBuf[(y * dstCoeffStride) + x]);
    }
  }
  return coeffAbsSum;
}
#endif
void TrQuant::transformNxN( TransformUnit& tu, const ComponentID& compID, const QpParam& cQP, TCoeff& uiAbsSum, const Ctx& ctx, const bool loadTr )
{
        CodingStructure &cs = *tu.cs;
  const SPS &sps            = *cs.sps;
  const CompArea &rect      = tu.blocks[compID];
  const uint32_t uiWidth        = rect.width;
  const uint32_t uiHeight       = rect.height;

  const CPelBuf resiBuf     = cs.getResiBuf(rect);

  if( tu.noResidual )
  {
    uiAbsSum = 0;
    TU::setCbfAtDepth( tu, compID, tu.depth, uiAbsSum > 0 );
    return;
  }

  if ((tu.cu->bdpcmMode && isLuma(compID)) || (!isLuma(compID) && tu.cu->bdpcmModeChroma))
  {
    tu.mtsIdx[compID] = MTS_SKIP;
  }

  uiAbsSum = 0;

  // transform and quantize
  CHECK(cs.sps->getMaxTbSize() < uiWidth, "Unsupported transformation size");

  CoeffBuf tempCoeff(loadTr ? m_mtsCoeffs[tu.mtsIdx[compID]] : m_tempCoeff, rect);

  DTRACE_PEL_BUF(D_RESIDUALS, resiBuf, tu, tu.cu->predMode, compID);

  if (!loadTr)
  {
    if (tu.mtsIdx[compID] == MTS_SKIP)
    {
      xTransformSkip(tu, compID, resiBuf, tempCoeff.buf);
    }
    else
    {
      xT(tu, compID, resiBuf, tempCoeff, uiWidth, uiHeight);
    }
  }

  if (sps.getUseLFNST())
  {
    xFwdLfnst(tu, compID, loadTr);
  }

  DTRACE_COEFF_BUF(D_TCOEFF, tempCoeff, tu, tu.cu->predMode, compID);

  xQuant(tu, compID, tempCoeff, uiAbsSum, cQP, ctx);

  DTRACE_COEFF_BUF(D_TCOEFF, tu.getCoeffs(compID), tu, tu.cu->predMode, compID);

  // set coded block flag (CBF)
  TU::setCbfAtDepth (tu, compID, tu.depth, uiAbsSum > 0);
}

void TrQuant::xTransformSkip(const TransformUnit &tu, const ComponentID &compID, const CPelBuf &resi, TCoeff* psCoeff)
{
  const CompArea &rect = tu.blocks[compID];
  const uint32_t width = rect.width;
  const uint32_t height = rect.height;
  const Pel *pelResi = resi.buf;

  for (uint32_t y = 0, coefficientIndex = 0; y < height; y++)
  {
    for (uint32_t x = 0; x < width; x++, coefficientIndex++)
    {
      psCoeff[coefficientIndex] = TCoeff(pelResi[x]);
    }
    pelResi += resi.stride;
  }
}

#if SIGN_PREDICTION
void TrQuant::predCoeffSigns(TransformUnit &tu, const ComponentID compID, const bool reshapeChroma)
{
  bool bIsJCCR = tu.jointCbCr && isChroma(compID);
  ComponentID residCompID = compID;
  bool bJccrWithCr = bIsJCCR && !(tu.jointCbCr >> 1);
  if(bJccrWithCr)
  {
    residCompID = COMPONENT_Cr;
  }
  if( !( TU::getUseSignPred( tu, residCompID ) && ( TU::getCbf( tu, compID ) || bIsJCCR ) ) )
  {
    return;
  }

  if( bIsJCCR && compID == COMPONENT_Cr )
  {
    return;
  }

  auto setCoeffSign = [](CoeffBuf &buff, uint32_t signMask, std::vector<Position> &pos) -> void
  {
    for( int i = 0, j = (int)pos.size() - 1; i < pos.size(); ++i, j-- )
    {
      bool bit_value = ( signMask >> j ) & 0x1;
      TCoeff &coeff = buff.at( pos[i] );
      coeff = std::abs( coeff ) * ( bit_value ? -1 : 1 );
    }
  };
  auto extractCoeffSign = [](CoeffBuf &buff, uint32_t &signMask, std::vector<Position> &pos) -> void
  {
    signMask = 0;
    for( int i = 0, j = (int)pos.size() - 1; i < pos.size(); ++i, j-- )
    {
      uint32_t coeffSign = buff.at( pos[i] ) < 0 ? 1 : 0;
      signMask |= coeffSign << j;
    }
  };

  auto createTemplate = [this,tu](ComponentID comp, uint32_t width, uint32_t height, uint32_t mtsIdx) -> void
  {
    // This is the function used to generate template values stored in g_initRomSignPred[]
    TCoeff *memCoeff = (TCoeff *)xMalloc(TCoeff, width*height);
    Pel      *memTmpResid = (Pel *)xMalloc(Pel,    width*height);
    CoeffBuf coeff(memCoeff, width, height);
    PelBuf   resi(memTmpResid, width, height);
    coeff.fill(0);

    int8_t  *pTemplate = (int8_t *)xMalloc(int8_t, (width + height - 1) * SIGN_PRED_FREQ_RANGE*SIGN_PRED_FREQ_RANGE);
    AreaBuf<int8_t> templateBuf(pTemplate, (width + height - 1), SIGN_PRED_FREQ_RANGE*SIGN_PRED_FREQ_RANGE);

    Position prev(0,0);
    int8_t *templ = templateBuf.buf;

    for( int j = 0; j < SIGN_PRED_FREQ_RANGE*SIGN_PRED_FREQ_RANGE; ++j)
    {
      Position curr(j%SIGN_PRED_FREQ_RANGE, j/SIGN_PRED_FREQ_RANGE);
      coeff.at(prev) = 0;
      coeff.at(curr) = 1 << SIGN_PRED_SHIFT;

      xIT( tu, comp, coeff, resi);

      Pel* pelResi = resi.bufAt(0, height - 1);

      for (uint32_t i = 0; i < height; i++)
      {
        templ[i] = (int8_t)(*pelResi);
        pelResi -= resi.stride;
      }

      pelResi = resi.buf;

      for (uint32_t i = 0; i < width; i++)
      {
        templ[i + height - 1] = (int8_t)pelResi[i];
      }

      templ += templateBuf.stride;
      prev = curr;
    }

    int log2Width = floorLog2(width);
    int log2Height = floorLog2(height);
    g_resiBorderTemplate[log2Width-2][log2Height-2][mtsIdx] = templateBuf.buf;

    xFree(memCoeff);
    xFree(memTmpResid);
  };

  

  const QpParam cQP( tu, residCompID );
  CodingStructure &cs = *tu.cs;
  PelBuf recoBuf = cs.picture->getRecoBuf(tu.blocks[residCompID]);
  PelBuf predBuf = cs.getPredBuf(tu.blocks[residCompID]);
  Pel        predResiBorder[2 * SIGN_PRED_MAX_BS];
  Pel        predResiTemplate[2 * SIGN_PRED_MAX_BS]{0};
  Pel        predResiTemplateReshape[2 * SIGN_PRED_MAX_BS]{0};
  TU::predBorderResi(tu.blocks[residCompID], recoBuf, predBuf, residCompID, tu.blocks[residCompID].width, tu.blocks[residCompID].height, predResiBorder, (1 << (tu.cs->sps->getBitDepth(toChannelType(residCompID)) - 1)));

  const uint32_t     uiWidth  = tu.blocks[residCompID].width;
  const uint32_t     uiHeight = tu.blocks[residCompID].height;
  PelBuf     bufResiTemplate(predResiTemplate, uiWidth + uiHeight - 1, 1);
  PelBuf     bufResiTemplateReshape(predResiTemplateReshape, uiWidth + uiHeight - 1, 1);
  int trHor, trVer;
  getTrTypes(tu, residCompID, trHor, trVer);
#if JVET_W0103_INTRA_MTS
  int actualTrIdx = trHor * NUM_TRANS_TYPE + trVer;
#else
  int actualTrIdx = trHor * 3 + trVer;
#endif
  int log2Width = floorLog2(uiWidth);
  int log2Height = floorLog2(uiHeight);
  if(!g_resiBorderTemplate[log2Width-2][log2Height-2][actualTrIdx])
  {
    createTemplate(residCompID, uiWidth, uiHeight, actualTrIdx);
  }
  AreaBuf<const int8_t> templateNormalizedBuf( g_resiBorderTemplate[log2Width-2][log2Height-2][actualTrIdx], uiWidth+uiHeight-1, SIGN_PRED_FREQ_RANGE*SIGN_PRED_FREQ_RANGE);
  PelBuf templateBuf(m_signPredTemplate, uiWidth+uiHeight-1, SIGN_PRED_FREQ_RANGE*SIGN_PRED_FREQ_RANGE);

  std::vector<Position> predSignsXY;
  DepQuant::getPredictedSigns( tu, residCompID, predSignsXY );
  int32_t numPredSigns = (int32_t)predSignsXY.size();

  if( !numPredSigns )
  {
    return;
  }

  uint32_t bufferSigns;
  CoeffBuf bufTmpQuant = CoeffBuf(m_tempCoeff, tu.blocks[residCompID]);
  PelBuf   piResi(m_tempSignPredResid, uiWidth, uiHeight);
  PelBuf   piResiCr(m_tempSignPredResid + uiWidth*uiHeight, uiWidth, uiHeight);

  CoeffBuf quantedCoeffBuff = tu.getCoeffs(residCompID);
  CoeffBuf bufSigns         = tu.getCoeffSigns(residCompID);

  extractCoeffSign(quantedCoeffBuff, bufferSigns, predSignsXY);
  setCoeffSign(quantedCoeffBuff, 0, predSignsXY);

  xDeQuant( tu, bufTmpQuant, residCompID, cQP );

  setCoeffSign(quantedCoeffBuff, bufferSigns, predSignsXY);

  TCoeff *tmpQuant = bufTmpQuant.bufAt( 0, 0 );
  int tmpQuantStride = bufTmpQuant.stride;

  for( auto xy : predSignsXY )
  {
    int pos = xy.y * SIGN_PRED_FREQ_RANGE + xy.x;

    Pel *templ = templateBuf.bufAt( 0, pos );
    const int8_t *templNorm = templateNormalizedBuf.bufAt( 0, pos );
    TCoeff temp = tmpQuant[xy.y * tmpQuantStride + xy.x];

    if( temp )
    {
      for( auto j = 0; j < uiWidth + uiHeight - 1; ++j )
      {
        templ[j] = ( temp * templNorm[j] + SIGN_PRED_OFFSET ) >> SIGN_PRED_SHIFT;
      }
    }
    else
    {
      std::memset( templ, 0, sizeof( *templ ) * ( uiWidth + uiHeight - 1 ) );
    }
  }

  if( bJccrWithCr )
  {
    xIT( tu, COMPONENT_Cr, bufTmpQuant, piResiCr );
  }
  else
  {
    xIT( tu, residCompID, bufTmpQuant, piResi );
  }

  Pel *pelResi = piResi.buf;
  int resiStride = piResi.stride;

  if( bIsJCCR )
  {
    invTransformICT( tu, piResi, piResiCr );
    if( bJccrWithCr )
    {
      pelResi = piResiCr.buf;
      resiStride = piResiCr.stride;
    }
  }

  Pel *resiTemplate = bufResiTemplate.buf + uiHeight - 1;

  for( int i = 1; i < uiWidth; i++ )
  {
    resiTemplate[i] += pelResi[i];
  }

  for( int i = 0; i < uiHeight; i++ )
  {
    resiTemplate[-i] += pelResi[0];
    pelResi += resiStride;
  }

  int signPrev = 0;

  for (uint32_t idx = 0; idx < (1 << numPredSigns); idx++)
  {
    const int32_t signCurr = (idx ^ (idx >> 1)); // Gray code

    if(idx)
    {
      if( reshapeChroma )
      {
        std::swap( bufResiTemplateReshape, bufResiTemplate );
      }
      int uiBit;
      int uiXor = signCurr ^ signPrev;
      for (uiBit = 0; uiBit < numPredSigns; uiBit++, uiXor >>= 1)
      {
        if (uiXor & 1)
        {
          break;
        }
      }

      bool signModifyTo = (signCurr >> uiBit) & 0x1;
      int  predSignIdx = numPredSigns - uiBit - 1;
      int pos_idx = predSignsXY[predSignIdx].y * SIGN_PRED_FREQ_RANGE + predSignsXY[predSignIdx].x;
      Pel *templ = templateBuf.bufAt(0, pos_idx);

      for (uint32_t i = 0; i < uiHeight + uiWidth - 1; i++)
      {
        bufResiTemplate.buf[i] += templ[i] * (signModifyTo ? -2 : 2);
      }
    }

    if(reshapeChroma)
    {
      bufResiTemplateReshape.copyFrom(bufResiTemplate);
      bufResiTemplate.scaleSignal(tu.getChromaAdj(), 0, tu.cu->cs->slice->clpRng(residCompID));
    }

    /*Compute cost of modificiation*/
    signPrev = signCurr;

    const bool firstBlock = !tu.blocks[residCompID].x && !tu.blocks[residCompID].y;

#if ENABLE_SIMD_SIGN_PREDICTION
    uint32_t cost = 0;

    if( tu.blocks[residCompID].x || firstBlock )
    {
      cost += m_computeSAD( predResiBorder, bufResiTemplate.buf, uiHeight );
    }

    if( tu.blocks[residCompID].y || firstBlock )
    {
      cost += m_computeSAD( predResiBorder + uiHeight, bufResiTemplate.buf + uiHeight - 1, uiWidth );
    }
#else
    uint32_t cost = 0;

    const Pel *pRef = predResiBorder;

    Pel *pCurr = bufResiTemplate.buf;

    if( tu.blocks[residCompID].x || firstBlock )
    {
      for( uint32_t i = 0; i < uiHeight; i++ )
      {
        cost += abs( pRef[i] - pCurr[i] );
      }
    }

    pRef += uiHeight;
    pCurr += uiHeight - 1;

    if( tu.blocks[residCompID].y || firstBlock )
    {
      for( uint32_t i = 0; i < uiWidth; i++ )
      {
        cost += abs( pRef[i] - pCurr[i] );
      }
    }
#endif

    m_aiSignPredCost[signCurr] = cost;
  }


  uint8_t realSign, resiSign = 1;
  uint32_t *pcCost            = m_aiSignPredCost;
  uint32_t  numSignsToProcess = numPredSigns;
  uint32_t  min_idx           = 0;

  for (uint32_t idx = 0; idx < numPredSigns; idx++)
  {
    Position xyPos = predSignsXY[idx];

    // Find predicted sign value
    if( resiSign )
    {
      uint32_t min_cost = -1;
      for (uint32_t c = 0; c < (1 << numSignsToProcess); c++)
      {
        if (pcCost[c] < min_cost)
        {
          min_cost = pcCost[c];
          min_idx  = c;
        }
      }
    }

    const uint8_t predSign = (min_idx >> (numSignsToProcess - 1)) & 1;

    numSignsToProcess--;

    if( tu.cs->pcv->isEncoder )
    {
      realSign = quantedCoeffBuff.at(xyPos) > 0 ? 0 : 1;
      resiSign = predSign ^ realSign;
    }
    else
    {
      resiSign = quantedCoeffBuff.at(xyPos) > 0 ? 0 : 1;
      realSign = predSign ^ resiSign;
      if (predSign)
      {
        quantedCoeffBuff.at(xyPos) = -quantedCoeffBuff.at(xyPos);
      }
    }

    bufSigns.at( xyPos ) = predSign ? SIGN_PRED_NEGATIVE : SIGN_PRED_POSITIVE;

    if( realSign )
    {
      pcCost += (uint32_t)(1 << numSignsToProcess);
    }
  }
}

#if ENABLE_SIMD_SIGN_PREDICTION
inline uint32_t TrQuant::xComputeSAD( const Pel* ref, const Pel* cur, const int size )
{
  uint32_t dist = 0;
  for( uint32_t i = 0; i < size; i++ )
  {
    dist += abs( ref[i] - cur[i] );
  }

  return dist;
}
#endif
#endif

//! \}
