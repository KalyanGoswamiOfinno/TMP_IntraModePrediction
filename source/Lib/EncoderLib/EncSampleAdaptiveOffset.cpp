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

/**
 \file     EncSampleAdaptiveOffset.cpp
 \brief       estimation part of sample adaptive offset class
 */
#include "EncSampleAdaptiveOffset.h"

#include "CommonLib/UnitTools.h"
#include "CommonLib/dtrace_codingstruct.h"
#include "CommonLib/dtrace_buffer.h"
#include "CommonLib/CodingStructure.h"
#if JVET_V0094_BILATERAL_FILTER
#include "CommonLib/BilateralFilter.h"
#endif

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//! \ingroup EncoderLib
//! \{


#define SAOCtx(c) SubCtx( Ctx::Sao, c )

#if JVET_W0066_CCSAO
#include <algorithm>

struct SetIdxCount
{
  uint8_t  setIdx;
  uint16_t count;
};

struct CtbCost
{
  int16_t pos;
  double  cost;
};

bool compareSetIdxCount(SetIdxCount a, SetIdxCount b) { return a.count > b.count; }

bool compareCtbCost(CtbCost a, CtbCost b) { return a.cost < b.cost; }
#endif

//! rounding with IBDI
inline double xRoundIbdi2(int bitDepth, double x)
{
#if FULL_NBIT
  return ((x) >= 0 ? ((int)((x) + 0.5)) : ((int)((x) -0.5)));
#else
  if (DISTORTION_PRECISION_ADJUSTMENT(bitDepth) == 0)
    return ((x) >= 0 ? ((int)((x) + 0.5)) : ((int)((x) -0.5)));
  else
    return ((x) > 0) ? (int)(((int)(x) + (1 << (DISTORTION_PRECISION_ADJUSTMENT(bitDepth) - 1)))
                             / (1 << DISTORTION_PRECISION_ADJUSTMENT(bitDepth)))
                     : ((int)(((int)(x) - (1 << (DISTORTION_PRECISION_ADJUSTMENT(bitDepth) - 1)))
                              / (1 << DISTORTION_PRECISION_ADJUSTMENT(bitDepth))));
#endif
}

inline double xRoundIbdi(int bitDepth, double x)
{
  return (bitDepth > 8 ? xRoundIbdi2(bitDepth, (x)) : ((x)>=0 ? ((int)((x)+0.5)) : ((int)((x)-0.5)))) ;
}


EncSampleAdaptiveOffset::EncSampleAdaptiveOffset()
{
  m_CABACEstimator = NULL;

  ::memset( m_saoDisabledRate, 0, sizeof( m_saoDisabledRate ) );
}

EncSampleAdaptiveOffset::~EncSampleAdaptiveOffset()
{
  destroyEncData();
}

void EncSampleAdaptiveOffset::createEncData(bool isPreDBFSamplesUsed, uint32_t numCTUsPic)
{
  //statistics
  const uint32_t sizeInCtus = numCTUsPic;
  m_statData.resize( sizeInCtus );
  for(uint32_t i=0; i< sizeInCtus; i++)
  {
    m_statData[i] = new SAOStatData*[MAX_NUM_COMPONENT];
    for(uint32_t compIdx=0; compIdx < MAX_NUM_COMPONENT; compIdx++)
    {
      m_statData[i][compIdx] = new SAOStatData[NUM_SAO_NEW_TYPES];
    }
  }
  if(isPreDBFSamplesUsed)
  {
    m_preDBFstatData.resize( sizeInCtus );
    for(uint32_t i=0; i< sizeInCtus; i++)
    {
      m_preDBFstatData[i] = new SAOStatData*[MAX_NUM_COMPONENT];
      for(uint32_t compIdx=0; compIdx < MAX_NUM_COMPONENT; compIdx++)
      {
        m_preDBFstatData[i][compIdx] = new SAOStatData[NUM_SAO_NEW_TYPES];
      }
    }

  }


  for(int typeIdc=0; typeIdc < NUM_SAO_NEW_TYPES; typeIdc++)
  {
    m_skipLinesR[COMPONENT_Y ][typeIdc]= 5;
    m_skipLinesR[COMPONENT_Cb][typeIdc]= m_skipLinesR[COMPONENT_Cr][typeIdc]= 3;

    m_skipLinesB[COMPONENT_Y ][typeIdc]= 4;
    m_skipLinesB[COMPONENT_Cb][typeIdc]= m_skipLinesB[COMPONENT_Cr][typeIdc]= 2;

    if(isPreDBFSamplesUsed)
    {
      switch(typeIdc)
      {
      case SAO_TYPE_EO_0:
        {
          m_skipLinesR[COMPONENT_Y ][typeIdc]= 5;
          m_skipLinesR[COMPONENT_Cb][typeIdc]= m_skipLinesR[COMPONENT_Cr][typeIdc]= 3;

          m_skipLinesB[COMPONENT_Y ][typeIdc]= 3;
          m_skipLinesB[COMPONENT_Cb][typeIdc]= m_skipLinesB[COMPONENT_Cr][typeIdc]= 1;
        }
        break;
      case SAO_TYPE_EO_90:
        {
          m_skipLinesR[COMPONENT_Y ][typeIdc]= 4;
          m_skipLinesR[COMPONENT_Cb][typeIdc]= m_skipLinesR[COMPONENT_Cr][typeIdc]= 2;

          m_skipLinesB[COMPONENT_Y ][typeIdc]= 4;
          m_skipLinesB[COMPONENT_Cb][typeIdc]= m_skipLinesB[COMPONENT_Cr][typeIdc]= 2;
        }
        break;
      case SAO_TYPE_EO_135:
      case SAO_TYPE_EO_45:
        {
          m_skipLinesR[COMPONENT_Y ][typeIdc]= 5;
          m_skipLinesR[COMPONENT_Cb][typeIdc]= m_skipLinesR[COMPONENT_Cr][typeIdc]= 3;

          m_skipLinesB[COMPONENT_Y ][typeIdc]= 4;
          m_skipLinesB[COMPONENT_Cb][typeIdc]= m_skipLinesB[COMPONENT_Cr][typeIdc]= 2;
        }
        break;
      case SAO_TYPE_BO:
        {
          m_skipLinesR[COMPONENT_Y ][typeIdc]= 4;
          m_skipLinesR[COMPONENT_Cb][typeIdc]= m_skipLinesR[COMPONENT_Cr][typeIdc]= 2;

          m_skipLinesB[COMPONENT_Y ][typeIdc]= 3;
          m_skipLinesB[COMPONENT_Cb][typeIdc]= m_skipLinesB[COMPONENT_Cr][typeIdc]= 1;
        }
        break;
      default:
        {
          THROW("Not a supported type");
        }
      }
    }
  }

#if JVET_W0066_CCSAO
  if (m_createdEnc)
  {
    return;
  }
  m_createdEnc = true;

  for (int i = 0; i < MAX_CCSAO_SET_NUM; i++)
  {
    m_ccSaoStatData[i] = new CcSaoStatData[m_numCTUsInPic];
  }

  m_bestCcSaoControl = new uint8_t[m_numCTUsInPic];
  m_tempCcSaoControl = new uint8_t[m_numCTUsInPic];
  m_initCcSaoControl = new uint8_t[m_numCTUsInPic];

  for (int i = 0; i < MAX_CCSAO_SET_NUM; i++)
  {
    m_trainingDistortion[i] = new int64_t[m_numCTUsInPic];
  }
#endif
}

void EncSampleAdaptiveOffset::destroyEncData()
{
  for(uint32_t i=0; i< m_statData.size(); i++)
  {
    for(uint32_t compIdx=0; compIdx< MAX_NUM_COMPONENT; compIdx++)
    {
      delete[] m_statData[i][compIdx];
    }
    delete[] m_statData[i];
  }
  m_statData.clear();


  for(int i=0; i< m_preDBFstatData.size(); i++)
  {
    for(int compIdx=0; compIdx< MAX_NUM_COMPONENT; compIdx++)
    {
      delete[] m_preDBFstatData[i][compIdx];
    }
    delete[] m_preDBFstatData[i];
  }
  m_preDBFstatData.clear();

#if JVET_W0066_CCSAO
  if (!m_createdEnc)
  {
    return;
  }
  m_createdEnc = false;

  for (int i = 0; i < MAX_CCSAO_SET_NUM; i++)
  {
    if (m_ccSaoStatData[i]) { delete[] m_ccSaoStatData[i]; m_ccSaoStatData[i] = nullptr; }
  }

  if (m_bestCcSaoControl) { delete[] m_bestCcSaoControl; m_bestCcSaoControl = nullptr; }
  if (m_tempCcSaoControl) { delete[] m_tempCcSaoControl; m_tempCcSaoControl = nullptr; }
  if (m_initCcSaoControl) { delete[] m_initCcSaoControl; m_initCcSaoControl = nullptr; }

  for (int i = 0; i < MAX_CCSAO_SET_NUM; i++)
  {
    if (m_trainingDistortion[i]) { delete[] m_trainingDistortion[i]; m_trainingDistortion[i] = nullptr; }
  }
#endif
}

void EncSampleAdaptiveOffset::initCABACEstimator( CABACEncoder* cabacEncoder, CtxCache* ctxCache, Slice* pcSlice )
{
  m_CABACEstimator = cabacEncoder->getCABACEstimator( pcSlice->getSPS() );
  m_CtxCache       = ctxCache;
  m_CABACEstimator->initCtxModels( *pcSlice );
  m_CABACEstimator->resetBits();
}


void EncSampleAdaptiveOffset::SAOProcess( CodingStructure& cs, bool* sliceEnabled, const double* lambdas,
#if ENABLE_QPA
                                          const double lambdaChromaWeight,
#endif
                                          const bool bTestSAODisableAtPictureLevel, const double saoEncodingRate, const double saoEncodingRateChroma, const bool isPreDBFSamplesUsed, bool isGreedyMergeEncoding
#if JVET_V0094_BILATERAL_FILTER
                                          ,BIFCabacEst* BifCABACEstimator
#endif
                                         )
{
#if ALF_SAO_TRUE_ORG && !JVET_V0094_BILATERAL_FILTER
  PelUnitBuf org = cs.getTrueOrgBuf();
#else
  PelUnitBuf org = cs.getOrgBuf();
#endif
  PelUnitBuf res = cs.getRecoBuf();
  PelUnitBuf src = m_tempBuf;
#if !JVET_V0094_BILATERAL_FILTER
  // Moved until after the bilateral filter has been initialized
  memcpy(m_lambda, lambdas, sizeof(m_lambda));
#endif

  src.copyFrom(res);

#if JVET_V0094_BILATERAL_FILTER
  const PreCalcValues& pcv = *cs.pcv;
  BifParams& bifParams = cs.picture->getBifParam();
  int width = cs.picture->lwidth();
  int height = cs.picture->lheight();
  int block_width = pcv.maxCUWidth;
  int block_height = pcv.maxCUHeight;

  int width_in_blocks = width / block_width + (width % block_width != 0);
  int height_in_blocks = height / block_height + (height % block_height != 0);

  bifParams.numBlocks = width_in_blocks * height_in_blocks;
  bifParams.ctuOn.resize(bifParams.numBlocks);

  std::fill(bifParams.ctuOn.begin(), bifParams.ctuOn.end(), 0);

  // Currently no RDO to figure out if we should turn CTUs on or off
  bifParams.frmOn = 1;
  bifParams.allCtuOn = 1;

  if (bifParams.frmOn == 0)
    std::fill(bifParams.ctuOn.begin(), bifParams.ctuOn.end(), 0);
  else if (bifParams.allCtuOn)
      std::fill(bifParams.ctuOn.begin(), bifParams.ctuOn.end(), 1);

  //double MseNoFltFrame = 0;
  //double MseFltDefFrame = 0;
  //double MseFltDefSwitchFrame = 0;
  //int CtuIdx = 0;
#endif
  
#if JVET_V0094_BILATERAL_FILTER
  BilateralFilter bilateralFilter;
  
  // Special case when SAO = 0 and BIF = 1.
  // Just filter reconstruction and return.
  // No need to estimate SAO parameters.
  if(!cs.sps->getSAOEnabledFlag() && cs.pps->getUseBIF())
  {
    bilateralFilter.create();
    bilateralFilter.bilateralFilterPicRDOperCTU(cs, src, BifCABACEstimator); // Filters from src to res
    bilateralFilter.destroy();
    return;
  }
  memcpy(m_lambda, lambdas, sizeof(m_lambda));
#endif
  
  //collect statistics
#if JVET_V0094_BILATERAL_FILTER
  //apply BILAT to res
  if(cs.pps->getUseBIF())
  {
    bilateralFilter.create();
    bilateralFilter.bilateralFilterPicRDOperCTU(cs, src, BifCABACEstimator); // Filters from src to res
    getStatistics(m_statData, org, src, res, cs);
    bilateralFilter.destroy();
  }
  else
  {
    getStatistics(m_statData, org, src, src, cs);
  }
#else
  getStatistics(m_statData, org, src, cs);
#endif

  if(isPreDBFSamplesUsed)
  {
    addPreDBFStatistics(m_statData);
  }
  
#if JVET_V0094_BILATERAL_FILTER
  //undo BILAT on res
  if(cs.pps->getUseBIF())
    res.copyFrom(src);
#endif
  
  //slice on/off
  decidePicParams(*cs.slice, sliceEnabled, saoEncodingRate, saoEncodingRateChroma);

  //block on/off
  std::vector<SAOBlkParam> reconParams(cs.pcv->sizeInCtus);
  decideBlkParams( cs, sliceEnabled, m_statData, src, res, &reconParams[0], cs.picture->getSAO(), bTestSAODisableAtPictureLevel,
#if ENABLE_QPA
                   lambdaChromaWeight,
#endif
                   saoEncodingRate, saoEncodingRateChroma, isGreedyMergeEncoding );

  DTRACE_UPDATE(g_trace_ctx, (std::make_pair("poc", cs.slice->getPOC())));
  DTRACE_PIC_COMP(D_REC_CB_LUMA_SAO, cs, cs.getRecoBuf(), COMPONENT_Y);
  DTRACE_PIC_COMP(D_REC_CB_CHROMA_SAO, cs, cs.getRecoBuf(), COMPONENT_Cb);
  DTRACE_PIC_COMP(D_REC_CB_CHROMA_SAO, cs, cs.getRecoBuf(), COMPONENT_Cr);

  DTRACE    ( g_trace_ctx, D_CRC, "SAO" );
  DTRACE_CRC( g_trace_ctx, D_CRC, cs, cs.getRecoBuf() );

}


void EncSampleAdaptiveOffset::getPreDBFStatistics(CodingStructure& cs)
{
#if ALF_SAO_TRUE_ORG
  PelUnitBuf org = cs.getTrueOrgBuf();
#else
  PelUnitBuf org = cs.getOrgBuf();
#endif
  PelUnitBuf rec = cs.getRecoBuf();
  getStatistics(m_preDBFstatData, org,
#if JVET_V0094_BILATERAL_FILTER
                rec, rec,
#else
                rec,
#endif
                
                cs, true);
}

void EncSampleAdaptiveOffset::addPreDBFStatistics(std::vector<SAOStatData**>& blkStats)
{
  const uint32_t numCTUsPic = (uint32_t)blkStats.size();
  for(uint32_t n=0; n< numCTUsPic; n++)
  {
    for(uint32_t compIdx=0; compIdx < MAX_NUM_COMPONENT; compIdx++)
    {
      for(uint32_t typeIdc=0; typeIdc < NUM_SAO_NEW_TYPES; typeIdc++)
      {
        blkStats[n][compIdx][typeIdc] += m_preDBFstatData[n][compIdx][typeIdc];
      }
    }
  }
}

void EncSampleAdaptiveOffset::getStatistics(std::vector<SAOStatData**>& blkStats, PelUnitBuf& orgYuv, PelUnitBuf& srcYuv,
#if JVET_V0094_BILATERAL_FILTER
                                            PelUnitBuf& bifYuv,
#endif
                                            CodingStructure& cs, bool isCalculatePreDeblockSamples)
{
  bool isLeftAvail, isRightAvail, isAboveAvail, isBelowAvail, isAboveLeftAvail, isAboveRightAvail;

  const PreCalcValues& pcv = *cs.pcv;
  const int numberOfComponents = getNumberValidComponents(pcv.chrFormat);

  size_t lineBufferSize = pcv.maxCUWidth + 1;
  if (m_signLineBuf1.size() != lineBufferSize)
  {
    m_signLineBuf1.resize(lineBufferSize);
    m_signLineBuf2.resize(lineBufferSize);
  }

  int ctuRsAddr = 0;
  for( uint32_t yPos = 0; yPos < pcv.lumaHeight; yPos += pcv.maxCUHeight )
  {
    for( uint32_t xPos = 0; xPos < pcv.lumaWidth; xPos += pcv.maxCUWidth )
    {
      const uint32_t width  = (xPos + pcv.maxCUWidth  > pcv.lumaWidth)  ? (pcv.lumaWidth - xPos)  : pcv.maxCUWidth;
      const uint32_t height = (yPos + pcv.maxCUHeight > pcv.lumaHeight) ? (pcv.lumaHeight - yPos) : pcv.maxCUHeight;
      const UnitArea area( cs.area.chromaFormat, Area(xPos , yPos, width, height) );

      deriveLoopFilterBoundaryAvailibility(cs, area.Y(), isLeftAvail, isAboveAvail, isAboveLeftAvail );

      //NOTE: The number of skipped lines during gathering CTU statistics depends on the slice boundary availabilities.
      //For simplicity, here only picture boundaries are considered.

      isRightAvail      = (xPos + pcv.maxCUWidth  < pcv.lumaWidth );
      isBelowAvail      = (yPos + pcv.maxCUHeight < pcv.lumaHeight);
      isAboveRightAvail = ((yPos > 0) && (isRightAvail));

      int numHorVirBndry = 0, numVerVirBndry = 0;
      int horVirBndryPos[] = { -1,-1,-1 };
      int verVirBndryPos[] = { -1,-1,-1 };
      int horVirBndryPosComp[] = { -1,-1,-1 };
      int verVirBndryPosComp[] = { -1,-1,-1 };
      bool isCtuCrossedByVirtualBoundaries = isCrossedByVirtualBoundaries(xPos, yPos, width, height, numHorVirBndry, numVerVirBndry, horVirBndryPos, verVirBndryPos, cs.picHeader );

      for(int compIdx = 0; compIdx < numberOfComponents; compIdx++)
      {
        const ComponentID compID = ComponentID(compIdx);
        const CompArea& compArea = area.block( compID );

        int  srcStride  = srcYuv.get(compID).stride;
        Pel* srcBlk     = srcYuv.get(compID).bufAt( compArea );

        int  orgStride  = orgYuv.get(compID).stride;
        Pel* orgBlk     = orgYuv.get(compID).bufAt( compArea );

#if JVET_V0094_BILATERAL_FILTER
        int  bifStride  = bifYuv.get(compID).stride;
        Pel* bifBlk     = bifYuv.get(compID).bufAt( compArea );
#endif
        
        for (int i = 0; i < numHorVirBndry; i++)
        {
          horVirBndryPosComp[i] = (horVirBndryPos[i] >> ::getComponentScaleY(compID, area.chromaFormat)) - compArea.y;
        }
        for (int i = 0; i < numVerVirBndry; i++)
        {
          verVirBndryPosComp[i] = (verVirBndryPos[i] >> ::getComponentScaleX(compID, area.chromaFormat)) - compArea.x;
        }

        getBlkStats(compID, cs.sps->getBitDepth(toChannelType(compID)), blkStats[ctuRsAddr][compID]
                  , srcBlk, orgBlk,
#if JVET_V0094_BILATERAL_FILTER
                    bifBlk, bifStride,
#endif
                    srcStride, orgStride, compArea.width, compArea.height
                  , isLeftAvail,  isRightAvail, isAboveAvail, isBelowAvail, isAboveLeftAvail, isAboveRightAvail
                  , isCalculatePreDeblockSamples
                  , isCtuCrossedByVirtualBoundaries, horVirBndryPosComp, verVirBndryPosComp, numHorVirBndry, numVerVirBndry
                  );
      }
      ctuRsAddr++;
    }
  }
}

void EncSampleAdaptiveOffset::decidePicParams(const Slice& slice, bool* sliceEnabled, const double saoEncodingRate, const double saoEncodingRateChroma)
{
  if ( slice.getPendingRasInit() )
  { // reset
    for (int compIdx = 0; compIdx < MAX_NUM_COMPONENT; compIdx++)
    {
      for (int tempLayer = 1; tempLayer < MAX_TLAYER; tempLayer++)
      {
        m_saoDisabledRate[compIdx][tempLayer] = 0.0;
      }
    }
  }

  const int picTempLayer = slice.getDepth();

  //decide sliceEnabled[compIdx]
  const int numberOfComponents = m_numberOfComponents;
  for (int compIdx = 0; compIdx < MAX_NUM_COMPONENT; compIdx++)
  {
    sliceEnabled[compIdx] = false;
  }

  for (int compIdx = 0; compIdx < numberOfComponents; compIdx++)
  {
    // reset flags & counters
    sliceEnabled[compIdx] = true;

    if (saoEncodingRate>0.0)
    {
      if (saoEncodingRateChroma>0.0)
      {
        // decide slice-level on/off based on previous results
        if( (picTempLayer > 0)
          && (m_saoDisabledRate[compIdx][picTempLayer-1] > ((compIdx==COMPONENT_Y) ? saoEncodingRate : saoEncodingRateChroma)) )
        {
          sliceEnabled[compIdx] = false;
        }
      }
      else
      {
        // decide slice-level on/off based on previous results
        if( (picTempLayer > 0)
          && (m_saoDisabledRate[COMPONENT_Y][0] > saoEncodingRate) )
        {
          sliceEnabled[compIdx] = false;
        }
      }
    }
  }
}

int64_t EncSampleAdaptiveOffset::getDistortion(const int channelBitDepth, int typeIdc, int typeAuxInfo, int* invQuantOffset, SAOStatData& statData)
{
  int64_t dist        = 0;
  int shift = 2 * DISTORTION_PRECISION_ADJUSTMENT(channelBitDepth);

  switch(typeIdc)
  {
    case SAO_TYPE_EO_0:
    case SAO_TYPE_EO_90:
    case SAO_TYPE_EO_135:
    case SAO_TYPE_EO_45:
      {
        for (int offsetIdx=0; offsetIdx<NUM_SAO_EO_CLASSES; offsetIdx++)
        {
          dist += estSaoDist( statData.count[offsetIdx], invQuantOffset[offsetIdx], statData.diff[offsetIdx], shift);
        }
      }
      break;
    case SAO_TYPE_BO:
      {
        for (int offsetIdx=typeAuxInfo; offsetIdx<typeAuxInfo+4; offsetIdx++)
        {
          int bandIdx = offsetIdx % NUM_SAO_BO_CLASSES ;
          dist += estSaoDist( statData.count[bandIdx], invQuantOffset[bandIdx], statData.diff[bandIdx], shift);
        }
      }
      break;
    default:
      {
        THROW("Not a supported type");
      }
  }

  return dist;
}

inline int64_t EncSampleAdaptiveOffset::estSaoDist(int64_t count, int64_t offset, int64_t diffSum, int shift)
{
  return (( count*offset*offset-diffSum*offset*2 ) >> shift);
}


inline int EncSampleAdaptiveOffset::estIterOffset(int typeIdx, double lambda, int offsetInput, int64_t count, int64_t diffSum, int shift, int bitIncrease, int64_t& bestDist, double& bestCost, int offsetTh )
{
  int iterOffset, tempOffset;
  int64_t tempDist, tempRate;
  double tempCost, tempMinCost;
  int offsetOutput = 0;
  iterOffset = offsetInput;
  // Assuming sending quantized value 0 results in zero offset and sending the value zero needs 1 bit. entropy coder can be used to measure the exact rate here.
  tempMinCost = lambda;
  while (iterOffset != 0)
  {
    // Calculate the bits required for signaling the offset
    tempRate = (typeIdx == SAO_TYPE_BO) ? (abs((int)iterOffset)+2) : (abs((int)iterOffset)+1);
    if (abs((int)iterOffset)==offsetTh) //inclusive
    {
      tempRate --;
    }
    // Do the dequantization before distortion calculation
    tempOffset  = iterOffset << bitIncrease;
    tempDist    = estSaoDist( count, tempOffset, diffSum, shift);
    tempCost    = ((double)tempDist + lambda * (double) tempRate);
    if(tempCost < tempMinCost)
    {
      tempMinCost = tempCost;
      offsetOutput = iterOffset;
      bestDist = tempDist;
      bestCost = tempCost;
    }
    iterOffset = (iterOffset > 0) ? (iterOffset-1):(iterOffset+1);
  }
  return offsetOutput;
}

void EncSampleAdaptiveOffset::deriveOffsets(ComponentID compIdx, const int channelBitDepth, int typeIdc, SAOStatData& statData, int* quantOffsets, int& typeAuxInfo)
{
  int bitDepth = channelBitDepth;
  int shift = 2 * DISTORTION_PRECISION_ADJUSTMENT(bitDepth);
  int offsetTh = SampleAdaptiveOffset::getMaxOffsetQVal(channelBitDepth);  //inclusive

  ::memset(quantOffsets, 0, sizeof(int)*MAX_NUM_SAO_CLASSES);

  //derive initial offsets
  int numClasses = (typeIdc == SAO_TYPE_BO)?((int)NUM_SAO_BO_CLASSES):((int)NUM_SAO_EO_CLASSES);
  for(int classIdx=0; classIdx< numClasses; classIdx++)
  {
    if( (typeIdc != SAO_TYPE_BO) && (classIdx==SAO_CLASS_EO_PLAIN)  )
    {
      continue; //offset will be zero
    }

    if(statData.count[classIdx] == 0)
    {
      continue; //offset will be zero
    }

    quantOffsets[classIdx] =
      (int) xRoundIbdi(bitDepth, (double)(statData.diff[classIdx] << DISTORTION_PRECISION_ADJUSTMENT(bitDepth))
                                   / (double)(statData.count[classIdx] << m_offsetStepLog2[compIdx]));
    quantOffsets[classIdx] = Clip3(-offsetTh, offsetTh, quantOffsets[classIdx]);
  }

  // adjust offsets
  switch(typeIdc)
  {
    case SAO_TYPE_EO_0:
    case SAO_TYPE_EO_90:
    case SAO_TYPE_EO_135:
    case SAO_TYPE_EO_45:
      {
        int64_t classDist;
        double classCost;
        for(int classIdx=0; classIdx<NUM_SAO_EO_CLASSES; classIdx++)
        {
          if(classIdx==SAO_CLASS_EO_FULL_VALLEY && quantOffsets[classIdx] < 0)
          {
            quantOffsets[classIdx] =0;
          }
          if(classIdx==SAO_CLASS_EO_HALF_VALLEY && quantOffsets[classIdx] < 0)
          {
            quantOffsets[classIdx] =0;
          }
          if(classIdx==SAO_CLASS_EO_HALF_PEAK   && quantOffsets[classIdx] > 0)
          {
            quantOffsets[classIdx] =0;
          }
          if(classIdx==SAO_CLASS_EO_FULL_PEAK   && quantOffsets[classIdx] > 0)
          {
            quantOffsets[classIdx] =0;
          }

          if( quantOffsets[classIdx] != 0 ) //iterative adjustment only when derived offset is not zero
          {
            quantOffsets[classIdx] = estIterOffset( typeIdc, m_lambda[compIdx], quantOffsets[classIdx], statData.count[classIdx], statData.diff[classIdx], shift, m_offsetStepLog2[compIdx], classDist , classCost , offsetTh );
          }
        }

        typeAuxInfo =0;
      }
      break;
    case SAO_TYPE_BO:
      {
        int64_t  distBOClasses[NUM_SAO_BO_CLASSES];
        double costBOClasses[NUM_SAO_BO_CLASSES];
        ::memset(distBOClasses, 0, sizeof(int64_t)*NUM_SAO_BO_CLASSES);
        for(int classIdx=0; classIdx< NUM_SAO_BO_CLASSES; classIdx++)
        {
          costBOClasses[classIdx]= m_lambda[compIdx];
          if( quantOffsets[classIdx] != 0 ) //iterative adjustment only when derived offset is not zero
          {
            quantOffsets[classIdx] = estIterOffset( typeIdc, m_lambda[compIdx], quantOffsets[classIdx], statData.count[classIdx], statData.diff[classIdx], shift, m_offsetStepLog2[compIdx], distBOClasses[classIdx], costBOClasses[classIdx], offsetTh );
          }
        }

        //decide the starting band index
        double minCost = MAX_DOUBLE, cost;
        for(int band=0; band< NUM_SAO_BO_CLASSES- 4+ 1; band++)
        {
          cost  = costBOClasses[band  ];
          cost += costBOClasses[band+1];
          cost += costBOClasses[band+2];
          cost += costBOClasses[band+3];

          if(cost < minCost)
          {
            minCost = cost;
            typeAuxInfo = band;
          }
        }
        //clear those unused classes
        int clearQuantOffset[NUM_SAO_BO_CLASSES];
        ::memset(clearQuantOffset, 0, sizeof(int)*NUM_SAO_BO_CLASSES);
        for(int i=0; i< 4; i++)
        {
          int band = (typeAuxInfo+i)%NUM_SAO_BO_CLASSES;
          clearQuantOffset[band] = quantOffsets[band];
        }
        ::memcpy(quantOffsets, clearQuantOffset, sizeof(int)*NUM_SAO_BO_CLASSES);
      }
      break;
    default:
      {
        THROW("Not a supported type");
      }

  }


}

void EncSampleAdaptiveOffset::deriveModeNewRDO(const BitDepths &bitDepths, int ctuRsAddr, SAOBlkParam* mergeList[NUM_SAO_MERGE_TYPES], bool* sliceEnabled, std::vector<SAOStatData**>& blkStats, SAOBlkParam& modeParam, double& modeNormCost )
{
  double minCost, cost;
  uint64_t previousFracBits;
  const int numberOfComponents = m_numberOfComponents;

  int64_t dist[MAX_NUM_COMPONENT], modeDist[MAX_NUM_COMPONENT];
  SAOOffset testOffset[MAX_NUM_COMPONENT];
  int invQuantOffset[MAX_NUM_SAO_CLASSES];
  for(int comp=0; comp < MAX_NUM_COMPONENT; comp++)
  {
    modeDist[comp] = 0;
  }

  //pre-encode merge flags
  modeParam[COMPONENT_Y].modeIdc = SAO_MODE_OFF;
  const TempCtx ctxStartBlk   ( m_CtxCache, SAOCtx( m_CABACEstimator->getCtx() ) );
  m_CABACEstimator->sao_block_pars( modeParam, bitDepths, sliceEnabled, (mergeList[SAO_MERGE_LEFT]!= NULL), (mergeList[SAO_MERGE_ABOVE]!= NULL), true );
  const TempCtx ctxStartLuma  ( m_CtxCache, SAOCtx( m_CABACEstimator->getCtx() ) );
  TempCtx       ctxBestLuma   ( m_CtxCache );

    //------ luma --------//
  {
    const ComponentID compIdx = COMPONENT_Y;
    //"off" case as initial cost
    modeParam[compIdx].modeIdc = SAO_MODE_OFF;
    m_CABACEstimator->resetBits();
    m_CABACEstimator->sao_offset_pars( modeParam[compIdx], compIdx, sliceEnabled[compIdx], bitDepths.recon[CHANNEL_TYPE_LUMA] );
    modeDist[compIdx] = 0;
    minCost           = m_lambda[compIdx] * (FRAC_BITS_SCALE * m_CABACEstimator->getEstFracBits());
    ctxBestLuma = SAOCtx( m_CABACEstimator->getCtx() );
    if(sliceEnabled[compIdx])
    {
      for(int typeIdc=0; typeIdc< NUM_SAO_NEW_TYPES; typeIdc++)
      {
        testOffset[compIdx].modeIdc = SAO_MODE_NEW;
        testOffset[compIdx].typeIdc = typeIdc;

        //derive coded offset
        deriveOffsets(compIdx, bitDepths.recon[CHANNEL_TYPE_LUMA], typeIdc, blkStats[ctuRsAddr][compIdx][typeIdc], testOffset[compIdx].offset, testOffset[compIdx].typeAuxInfo);

        //inversed quantized offsets
        invertQuantOffsets(compIdx, typeIdc, testOffset[compIdx].typeAuxInfo, invQuantOffset, testOffset[compIdx].offset);

        //get distortion
        dist[compIdx] = getDistortion(bitDepths.recon[CHANNEL_TYPE_LUMA], testOffset[compIdx].typeIdc, testOffset[compIdx].typeAuxInfo, invQuantOffset, blkStats[ctuRsAddr][compIdx][typeIdc]);

        //get rate
        m_CABACEstimator->getCtx() = SAOCtx( ctxStartLuma );
        m_CABACEstimator->resetBits();
        m_CABACEstimator->sao_offset_pars( testOffset[compIdx], compIdx, sliceEnabled[compIdx], bitDepths.recon[CHANNEL_TYPE_LUMA] );
        double rate = FRAC_BITS_SCALE * m_CABACEstimator->getEstFracBits();
        cost = (double)dist[compIdx] + m_lambda[compIdx]*rate;
        if(cost < minCost)
        {
          minCost = cost;
          modeDist[compIdx] = dist[compIdx];
          modeParam[compIdx]= testOffset[compIdx];
          ctxBestLuma = SAOCtx( m_CABACEstimator->getCtx() );
        }
      }
    }
    m_CABACEstimator->getCtx() = SAOCtx( ctxBestLuma );
  }

  //------ chroma --------//
//"off" case as initial cost
  cost = 0;
  previousFracBits = 0;
  m_CABACEstimator->resetBits();
  for(uint32_t componentIndex = COMPONENT_Cb; componentIndex < numberOfComponents; componentIndex++)
  {
    const ComponentID component = ComponentID(componentIndex);

    modeParam[component].modeIdc = SAO_MODE_OFF;
    modeDist [component]         = 0;
    m_CABACEstimator->sao_offset_pars( modeParam[component], component, sliceEnabled[component], bitDepths.recon[CHANNEL_TYPE_CHROMA] );
    const uint64_t currentFracBits = m_CABACEstimator->getEstFracBits();
    cost += m_lambda[component] * FRAC_BITS_SCALE * (currentFracBits - previousFracBits);
    previousFracBits = currentFracBits;
  }

  minCost = cost;

  //doesn't need to store cabac status here since the whole CTU parameters will be re-encoded at the end of this function

  for(int typeIdc=0; typeIdc< NUM_SAO_NEW_TYPES; typeIdc++)
  {
    m_CABACEstimator->getCtx() = SAOCtx( ctxBestLuma );
    m_CABACEstimator->resetBits();
    previousFracBits = 0;
    cost = 0;

    for(uint32_t componentIndex = COMPONENT_Cb; componentIndex < numberOfComponents; componentIndex++)
    {
      const ComponentID component = ComponentID(componentIndex);
      if(!sliceEnabled[component])
      {
        testOffset[component].modeIdc = SAO_MODE_OFF;
        dist[component]= 0;
        continue;
      }
      testOffset[component].modeIdc = SAO_MODE_NEW;
      testOffset[component].typeIdc = typeIdc;

      //derive offset & get distortion
      deriveOffsets(component, bitDepths.recon[CHANNEL_TYPE_CHROMA], typeIdc, blkStats[ctuRsAddr][component][typeIdc], testOffset[component].offset, testOffset[component].typeAuxInfo);
      invertQuantOffsets(component, typeIdc, testOffset[component].typeAuxInfo, invQuantOffset, testOffset[component].offset);
      dist[component] = getDistortion(bitDepths.recon[CHANNEL_TYPE_CHROMA], typeIdc, testOffset[component].typeAuxInfo, invQuantOffset, blkStats[ctuRsAddr][component][typeIdc]);
      m_CABACEstimator->sao_offset_pars( testOffset[component], component, sliceEnabled[component], bitDepths.recon[CHANNEL_TYPE_CHROMA] );
      const uint64_t currentFracBits = m_CABACEstimator->getEstFracBits();
      cost += dist[component] + (m_lambda[component] * FRAC_BITS_SCALE * (currentFracBits - previousFracBits));
      previousFracBits = currentFracBits;
    }

    if(cost < minCost)
    {
      minCost = cost;
      for(uint32_t componentIndex = COMPONENT_Cb; componentIndex < numberOfComponents; componentIndex++)
      {
        modeDist[componentIndex]  = dist[componentIndex];
        modeParam[componentIndex] = testOffset[componentIndex];
      }
    }

  } // SAO_TYPE loop

  //----- re-gen rate & normalized cost----//
  modeNormCost = 0;
  for(uint32_t componentIndex = COMPONENT_Y; componentIndex < numberOfComponents; componentIndex++)
  {
    modeNormCost += (double)modeDist[componentIndex] / m_lambda[componentIndex];
  }

  m_CABACEstimator->getCtx() = SAOCtx( ctxStartBlk );
  m_CABACEstimator->resetBits();
  m_CABACEstimator->sao_block_pars( modeParam, bitDepths, sliceEnabled, (mergeList[SAO_MERGE_LEFT]!= NULL), (mergeList[SAO_MERGE_ABOVE]!= NULL), false );
  modeNormCost += FRAC_BITS_SCALE * m_CABACEstimator->getEstFracBits();
}

void EncSampleAdaptiveOffset::deriveModeMergeRDO(const BitDepths &bitDepths, int ctuRsAddr, SAOBlkParam* mergeList[NUM_SAO_MERGE_TYPES], bool* sliceEnabled, std::vector<SAOStatData**>& blkStats, SAOBlkParam& modeParam, double& modeNormCost )
{
  modeNormCost = MAX_DOUBLE;

  double cost;
  SAOBlkParam testBlkParam;
  const int numberOfComponents = m_numberOfComponents;

  const TempCtx ctxStart  ( m_CtxCache, SAOCtx( m_CABACEstimator->getCtx() ) );
  TempCtx       ctxBest   ( m_CtxCache );

  for(int mergeType=0; mergeType< NUM_SAO_MERGE_TYPES; mergeType++)
  {
    if(mergeList[mergeType] == NULL)
    {
      continue;
    }

    testBlkParam = *(mergeList[mergeType]);
    //normalized distortion
    double normDist=0;
    for(int compIdx = 0; compIdx < numberOfComponents; compIdx++)
    {
      testBlkParam[compIdx].modeIdc = SAO_MODE_MERGE;
      testBlkParam[compIdx].typeIdc = mergeType;

      SAOOffset& mergedOffsetParam = (*(mergeList[mergeType]))[compIdx];

      if( mergedOffsetParam.modeIdc != SAO_MODE_OFF)
      {
        //offsets have been reconstructed. Don't call inversed quantization function.
        normDist += (((double)getDistortion(bitDepths.recon[toChannelType(ComponentID(compIdx))], mergedOffsetParam.typeIdc, mergedOffsetParam.typeAuxInfo, mergedOffsetParam.offset, blkStats[ctuRsAddr][compIdx][mergedOffsetParam.typeIdc]))
                       /m_lambda[compIdx] );
      }
    }

    //rate
    m_CABACEstimator->getCtx() = SAOCtx( ctxStart );
    m_CABACEstimator->resetBits();
    m_CABACEstimator->sao_block_pars( testBlkParam, bitDepths, sliceEnabled, (mergeList[SAO_MERGE_LEFT]!= NULL), (mergeList[SAO_MERGE_ABOVE]!= NULL), false );
    double rate = FRAC_BITS_SCALE * m_CABACEstimator->getEstFracBits();
    cost = normDist+rate;

    if(cost < modeNormCost)
    {
      modeNormCost = cost;
      modeParam    = testBlkParam;
      ctxBest      = SAOCtx( m_CABACEstimator->getCtx() );
    }
  }
  if( modeNormCost < MAX_DOUBLE )
  {
    m_CABACEstimator->getCtx() = SAOCtx( ctxBest );
  }
}

void EncSampleAdaptiveOffset::decideBlkParams(CodingStructure& cs, bool* sliceEnabled, std::vector<SAOStatData**>& blkStats, PelUnitBuf& srcYuv, PelUnitBuf& resYuv,
                                               SAOBlkParam* reconParams, SAOBlkParam* codedParams, const bool bTestSAODisableAtPictureLevel,
#if ENABLE_QPA
                                               const double chromaWeight,
#endif
                                               const double saoEncodingRate, const double saoEncodingRateChroma, const bool isGreedymergeEncoding)

{
  const PreCalcValues& pcv = *cs.pcv;
  bool allBlksDisabled = true;
#if JVET_V0094_BILATERAL_FILTER
  if(cs.sps->getSAOEnabledFlag())
  {
    // If SAO is enabled, we should investigate the components.
    // If SAO is disabled, we should stick with allBlksDisabled=true;
#endif
  const uint32_t numberOfComponents = m_numberOfComponents;
  for(uint32_t compId = COMPONENT_Y; compId < numberOfComponents; compId++)
  {
    if (sliceEnabled[compId])
    {
      allBlksDisabled = false;
    }
  }
#if JVET_V0094_BILATERAL_FILTER
  }
#endif

#if JVET_V0094_BILATERAL_FILTER
  BilateralFilter bilateralFilter;
  bilateralFilter.create();
#endif
  
  const TempCtx ctxPicStart ( m_CtxCache, SAOCtx( m_CABACEstimator->getCtx() ) );

  SAOBlkParam modeParam;
  double minCost, modeCost;

  double minCost2 = 0;
  std::vector<SAOStatData**> groupBlkStat;
  if (isGreedymergeEncoding)
  {
    groupBlkStat.resize(cs.pcv->sizeInCtus);
    for (uint32_t k = 0; k < cs.pcv->sizeInCtus; k++)
    {
      groupBlkStat[k] = new SAOStatData*[MAX_NUM_COMPONENT];
      for (uint32_t compIdx = 0; compIdx < MAX_NUM_COMPONENT; compIdx++)
      {
        groupBlkStat[k][compIdx] = new SAOStatData[NUM_SAO_NEW_TYPES];
      }
    }
  }
  SAOBlkParam  testBlkParam;
  SAOBlkParam  groupParam;
  SAOBlkParam* tempMergeList[NUM_SAO_MERGE_TYPES] = { NULL };
  SAOBlkParam* startingMergeList[NUM_SAO_MERGE_TYPES] = { NULL };

  int     mergeCtuAddr = 1; //Ctu to be merged
  int     groupSize = 1;
  double  Cost[2] = { 0, 0 };
  TempCtx ctxBeforeMerge(m_CtxCache);
  TempCtx ctxAfterMerge(m_CtxCache);

  double totalCost = 0; // Used if bTestSAODisableAtPictureLevel==true

  int ctuRsAddr = 0;
#if ENABLE_QPA
  CHECK ((chromaWeight > 0.0) && (cs.slice->getFirstCtuRsAddrInSlice() != 0), "incompatible start CTU address, must be 0");
#endif

  for( uint32_t yPos = 0; yPos < pcv.lumaHeight; yPos += pcv.maxCUHeight )
  {
    for( uint32_t xPos = 0; xPos < pcv.lumaWidth; xPos += pcv.maxCUWidth )
    {
      const uint32_t width  = (xPos + pcv.maxCUWidth  > pcv.lumaWidth)  ? (pcv.lumaWidth - xPos)  : pcv.maxCUWidth;
      const uint32_t height = (yPos + pcv.maxCUHeight > pcv.lumaHeight) ? (pcv.lumaHeight - yPos) : pcv.maxCUHeight;
      const UnitArea area( pcv.chrFormat, Area( xPos , yPos, width, height) );

      if(allBlksDisabled)
      {
        codedParams[ctuRsAddr].reset();
#if JVET_V0094_BILATERAL_FILTER
        // In the combined filter, we cannot continue here, even if SAO is
        // turned off for all blocks, since we need to perform the bilateral
        // filter later on (see next JVET_V0094_BILATERAL_FILTER).
        if(!cs.pps->getUseBIF())
        {
          // Unless we are not using BIF. Then it is safe to continue.
          continue;
        }
#else
        continue;
#endif
      }

      const TempCtx  ctxStart ( m_CtxCache, SAOCtx( m_CABACEstimator->getCtx() ) );
      TempCtx        ctxBest  ( m_CtxCache );

      if (ctuRsAddr == (mergeCtuAddr - 1))
      {
        ctxBeforeMerge = SAOCtx(m_CABACEstimator->getCtx());
      }

      //get merge list
      SAOBlkParam* mergeList[NUM_SAO_MERGE_TYPES] = { NULL };
      getMergeList(cs, ctuRsAddr, reconParams, mergeList);

      minCost = MAX_DOUBLE;
#if ENABLE_QPA
      if (chromaWeight > 0.0) // temporarily adopt local (CTU-wise) lambdas from QPA
      {
        for (int compIdx = 0; compIdx < MAX_NUM_COMPONENT; compIdx++)
        {
          m_lambda[compIdx] = isLuma ((ComponentID)compIdx) ? cs.picture->m_uEnerHpCtu[ctuRsAddr] : cs.picture->m_uEnerHpCtu[ctuRsAddr] / chromaWeight;
        }
      }
#endif
      for(int mode=1; mode < NUM_SAO_MODES; mode++)
      {
        if( mode > 1 )
        {
          m_CABACEstimator->getCtx() = SAOCtx( ctxStart );
        }
        switch(mode)
        {
        case SAO_MODE_NEW:
          {
            deriveModeNewRDO(cs.sps->getBitDepths(), ctuRsAddr, mergeList, sliceEnabled, blkStats, modeParam, modeCost );
          }
          break;
        case SAO_MODE_MERGE:
          {
            deriveModeMergeRDO(cs.sps->getBitDepths(), ctuRsAddr, mergeList, sliceEnabled, blkStats , modeParam, modeCost );
          }
          break;
        default:
          {
            THROW( "Not a supported SAO mode." );
          }
        }

        if(modeCost < minCost)
        {
          minCost                = modeCost;
          codedParams[ctuRsAddr] = modeParam;
          ctxBest                = SAOCtx( m_CABACEstimator->getCtx() );
        }
      } //mode

      if (!isGreedymergeEncoding)
      {
      totalCost += minCost;
      }


      m_CABACEstimator->getCtx() = SAOCtx( ctxBest );

      //apply reconstructed offsets
      reconParams[ctuRsAddr] = codedParams[ctuRsAddr];
      reconstructBlkSAOParam(reconParams[ctuRsAddr], mergeList);

      if (isGreedymergeEncoding)
      {
        if (ctuRsAddr == (mergeCtuAddr - 1))
        {
          Cost[0] = minCost;  //previous
          groupSize = 1;
          getMergeList(cs, ctuRsAddr, reconParams, startingMergeList);
        }
        else if (ctuRsAddr == mergeCtuAddr)
        {
          Cost[1] = minCost;
          minCost2 = MAX_DOUBLE;
          for (int tmp = groupSize; tmp >= 0; tmp--)
          {
            for (int compIdx = 0; compIdx < MAX_NUM_COMPONENT; compIdx++)
            {
              for (int i = 0; i < NUM_SAO_NEW_TYPES; i++)
              {
                for (int j = 0; j < MAX_NUM_SAO_CLASSES; j++)
                {
                  if (tmp == groupSize)
                  {
                    groupBlkStat[ctuRsAddr][compIdx][i].count[j] = blkStats[ctuRsAddr - tmp][compIdx][i].count[j];
                    groupBlkStat[ctuRsAddr][compIdx][i].diff[j] = blkStats[ctuRsAddr - tmp][compIdx][i].diff[j];
                  }
                  else
                  {
                    groupBlkStat[ctuRsAddr][compIdx][i].count[j] += blkStats[ctuRsAddr - tmp][compIdx][i].count[j];
                    groupBlkStat[ctuRsAddr][compIdx][i].diff[j] += blkStats[ctuRsAddr - tmp][compIdx][i].diff[j];
                  }
                }
              }
            }
          }

          // Derive new offset for grouped CTUs
          m_CABACEstimator->getCtx() = SAOCtx(ctxBeforeMerge);
          deriveModeNewRDO(cs.sps->getBitDepths(), ctuRsAddr, startingMergeList, sliceEnabled, groupBlkStat, modeParam, modeCost);

          //rate for mergeLeft CTB
          testBlkParam[COMPONENT_Y].modeIdc = SAO_MODE_MERGE;
          testBlkParam[COMPONENT_Y].typeIdc = SAO_MERGE_LEFT;
          m_CABACEstimator->resetBits();
          m_CABACEstimator->sao_block_pars(testBlkParam, cs.sps->getBitDepths(), sliceEnabled, true, false, true);
          double rate = FRAC_BITS_SCALE * m_CABACEstimator->getEstFracBits();
          modeCost += rate * groupSize;
          if (modeCost < minCost2)
          {
            groupParam = modeParam;
            minCost2 = modeCost;
            ctxAfterMerge = SAOCtx(m_CABACEstimator->getCtx());
          }

          // Test merge mode for grouped CTUs
          m_CABACEstimator->getCtx() = SAOCtx(ctxStart);
          deriveModeMergeRDO(cs.sps->getBitDepths(), ctuRsAddr, startingMergeList, sliceEnabled, groupBlkStat, modeParam, modeCost);
          modeCost += rate * groupSize;
          if (modeCost < minCost2)
          {
            minCost2 = modeCost;
            groupParam = modeParam;
            ctxAfterMerge = SAOCtx(m_CABACEstimator->getCtx());
          }

          totalCost += Cost[0];
          totalCost += Cost[1];

          if ((Cost[0] + Cost[1]) > minCost2) //merge current CTU
          {
            //original merge all
            totalCost = totalCost - Cost[0] - Cost[1] + minCost2;
            codedParams[ctuRsAddr - groupSize] = groupParam;
            for (int compIdx = 0; compIdx < MAX_NUM_COMPONENT; compIdx++)
            {
              codedParams[ctuRsAddr][compIdx].modeIdc = SAO_MODE_MERGE;
              codedParams[ctuRsAddr][compIdx].typeIdc = SAO_MERGE_LEFT;
            }
            for (int i = groupSize; i >= 0; i--) //change previous results
            {
              reconParams[ctuRsAddr - i] = codedParams[ctuRsAddr - i];
              getMergeList(cs, ctuRsAddr - i, reconParams, tempMergeList);
              reconstructBlkSAOParam(reconParams[ctuRsAddr - i], tempMergeList);
            }

            mergeCtuAddr += 1;
            if (mergeCtuAddr % pcv.widthInCtus == 0) //reaching the end of a row
            {
              mergeCtuAddr += 1;
            }
            else //next CTU can be merged with current group
            {
              Cost[0] = minCost2;
              groupSize += 1;
            }
            m_CABACEstimator->getCtx() = SAOCtx(ctxAfterMerge);
          }
          else // don't merge current CTU
          {
            mergeCtuAddr += 1;
            // Current block will be the starting block for successive operations
            Cost[0] = Cost[1];
            getMergeList(cs, ctuRsAddr, reconParams, startingMergeList);
            groupSize = 1;
            m_CABACEstimator->getCtx() = SAOCtx(ctxStart);
            ctxBeforeMerge = SAOCtx(m_CABACEstimator->getCtx());
            m_CABACEstimator->getCtx() = SAOCtx(ctxBest);
            if (mergeCtuAddr% pcv.widthInCtus == 0) //reaching the end of a row
            {
              mergeCtuAddr += 1;
            }
          } //else, if(Cost[0] + Cost[1] > minCost2)
        }//else if (ctuRsAddr == mergeCtuAddr)
      }
      else
      {
#if JVET_W0066_CCSAO
      offsetCTUnoClip(area, srcYuv, resYuv, reconParams[ctuRsAddr], cs);
#if JVET_V0094_BILATERAL_FILTER
      if (cs.pps->getUseBIF())
      {
        BifParams& bifParams = cs.picture->getBifParam();
        for (auto& currCU : cs.traverseCUs(CS::getArea(cs, area, CH_L), CH_L))
        {
          for (auto& currTU : CU::traverseTUs(currCU))
          {

            bool isInter = (currCU.predMode == MODE_INTER) ? true : false;
            if (bifParams.ctuOn[ctuRsAddr] && ((TU::getCbf(currTU, COMPONENT_Y) || isInter == false) && (currTU.cu->qp > 17)) && (128 > std::max(currTU.lumaSize().width, currTU.lumaSize().height)) && ((isInter == false) || (32 > std::min(currTU.lumaSize().width, currTU.lumaSize().height))))
            {
              bilateralFilter.bilateralFilterDiamond5x5NoClip(srcYuv, resYuv, currTU.cu->qp, cs.slice->clpRng(COMPONENT_Y), currTU);
            }
          }
        }
      }
#endif
#else
#if JVET_V0094_BILATERAL_FILTER
        if(cs.pps->getUseBIF())
        {
          offsetCTUnoClip(area, srcYuv, resYuv, reconParams[ctuRsAddr], cs);
          // Avoid old slow code
          // offsetCTUonlyBIF(area, srcYuv, resYuv, reconParams[ctuRsAddr], cs);
          // and instead do the code included below.
          
          // We don't need to clip if SAO was not performed on luma.
          SAOBlkParam mySAOblkParam = cs.picture->getSAO()[ctuRsAddr];
          SAOOffset& myCtbOffset     = mySAOblkParam[0];
          BifParams& bifParams = cs.picture->getBifParam();
          
          bool clipLumaIfNoBilat = false;
          if(myCtbOffset.modeIdc != SAO_MODE_OFF)
            clipLumaIfNoBilat = true;
          
          for (auto &currCU : cs.traverseCUs(CS::getArea(cs, area, CH_L), CH_L))
          {
            for (auto &currTU : CU::traverseTUs(currCU))
            {
              
              bool isInter = (currCU.predMode == MODE_INTER) ? true : false;
              if ( bifParams.ctuOn[ctuRsAddr] && ((TU::getCbf(currTU, COMPONENT_Y) || isInter == false) && (currTU.cu->qp > 17)) && (128 > std::max(currTU.lumaSize().width, currTU.lumaSize().height)) && ((isInter == false) || (32 > std::min(currTU.lumaSize().width, currTU.lumaSize().height))))
              {
                bilateralFilter.bilateralFilterDiamond5x5(srcYuv, resYuv, currTU.cu->qp, cs.slice->clpRng(COMPONENT_Y), currTU);
              }
              else
              {
                // We don't need to clip if SAO was not performed on luma.
                if(clipLumaIfNoBilat)
                  bilateralFilter.clipNotBilaterallyFilteredBlocks(srcYuv, resYuv, cs.slice->clpRng(COMPONENT_Y), currTU);
              }
            }
          }
        }
        else
        {
          // We do not do BIF for this sequence, so we can use the regular SAO function
          offsetCTU(area, srcYuv, resYuv, reconParams[ctuRsAddr], cs);
        }
#else
      offsetCTU(area, srcYuv, resYuv, reconParams[ctuRsAddr], cs);
#endif
#endif

      }

      ctuRsAddr++;
    } //ctuRsAddr
  }

#if ENABLE_QPA
  // restore global lambdas (might be unnecessary)
  if (chromaWeight > 0.0) memcpy (m_lambda, cs.slice->getLambdas(), sizeof (m_lambda));

#endif
  //reconstruct
#if JVET_V0094_BILATERAL_FILTER
  if (isGreedymergeEncoding || (cs.pps->getUseBIF() &&allBlksDisabled) )
#else
  if (isGreedymergeEncoding)
#endif
  {
    ctuRsAddr = 0;
    for (uint32_t yPos = 0; yPos < pcv.lumaHeight; yPos += pcv.maxCUHeight)
    {
      for (uint32_t xPos = 0; xPos < pcv.lumaWidth; xPos += pcv.maxCUWidth)
      {
        const uint32_t width = (xPos + pcv.maxCUWidth > pcv.lumaWidth) ? (pcv.lumaWidth - xPos) : pcv.maxCUWidth;
        const uint32_t height = (yPos + pcv.maxCUHeight > pcv.lumaHeight) ? (pcv.lumaHeight - yPos) : pcv.maxCUHeight;

        const UnitArea area(pcv.chrFormat, Area(xPos, yPos, width, height));
        
#if JVET_V0094_BILATERAL_FILTER
        if(cs.pps->getUseBIF())
        {
          // Sorry for not using nice copy method and using ugly code instead:
          int  myResStride = resYuv.get(COMPONENT_Y).stride;
          const CompArea& myCompArea = area.block(COMPONENT_Y);
          Pel* myResBlk = resYuv.get(COMPONENT_Y).bufAt(myCompArea);
          int mySrcStride = srcYuv.get(COMPONENT_Y).stride;
          Pel* mySrcBlk = srcYuv.get(COMPONENT_Y).bufAt(myCompArea);
          
          for(int yy = 0; yy<area.lheight(); yy++)
            for(int xx = 0; xx<area.lwidth(); xx++)
              myResBlk[yy*myResStride+xx] = mySrcBlk[yy*mySrcStride+xx];
        }
#endif

#if JVET_W0066_CCSAO
        offsetCTUnoClip(area, srcYuv, resYuv, reconParams[ctuRsAddr], cs);
#if JVET_V0094_BILATERAL_FILTER
        if (cs.pps->getUseBIF())
        {
          BifParams& bifParams = cs.picture->getBifParam();
          for (auto& currCU : cs.traverseCUs(CS::getArea(cs, area, CH_L), CH_L))
          {
            for (auto& currTU : CU::traverseTUs(currCU))
            {

              bool isInter = (currCU.predMode == MODE_INTER) ? true : false;
              if (bifParams.ctuOn[ctuRsAddr] && ((TU::getCbf(currTU, COMPONENT_Y) || isInter == false) && (currTU.cu->qp > 17)) && (128 > std::max(currTU.lumaSize().width, currTU.lumaSize().height)) && ((isInter == false) || (32 > std::min(currTU.lumaSize().width, currTU.lumaSize().height))))
              {
                bilateralFilter.bilateralFilterDiamond5x5NoClip(srcYuv, resYuv, currTU.cu->qp, cs.slice->clpRng(COMPONENT_Y), currTU);
              }
            }
          }
        }
#endif
#else
#if JVET_V0094_BILATERAL_FILTER
        if(cs.pps->getUseBIF())
        {
          offsetCTUnoClip(area, srcYuv, resYuv, reconParams[ctuRsAddr], cs);
          
          // Avoid old slow code
          // offsetCTUonlyBIF(area, srcYuv, resYuv, reconParams[ctuRsAddr], cs);
          // and instead do the code included below.
          
          // We don't need to clip if SAO was not performed on luma.
          SAOBlkParam mySAOblkParam = cs.picture->getSAO()[ctuRsAddr];
          SAOOffset& myCtbOffset     = mySAOblkParam[0];
          BifParams& bifParams = cs.picture->getBifParam();
          
          bool clipLumaIfNoBilat = false;
          if(myCtbOffset.modeIdc != SAO_MODE_OFF)
            clipLumaIfNoBilat = true;
          
          for (auto &currCU : cs.traverseCUs(CS::getArea(cs, area, CH_L), CH_L))
          {
            for (auto &currTU : CU::traverseTUs(currCU))
            {
              
              bool isInter = (currCU.predMode == MODE_INTER) ? true : false;
              if ( bifParams.ctuOn[ctuRsAddr] && ((TU::getCbf(currTU, COMPONENT_Y) || isInter == false) && (currTU.cu->qp > 17)) && (128 > std::max(currTU.lumaSize().width, currTU.lumaSize().height)) && ((isInter == false) || (32 > std::min(currTU.lumaSize().width, currTU.lumaSize().height))))
              {
                bilateralFilter.bilateralFilterDiamond5x5(srcYuv, resYuv, currTU.cu->qp, cs.slice->clpRng(COMPONENT_Y), currTU);
              }
              else
              {
                // We don't need to clip if SAO was not performed on luma.
                if(clipLumaIfNoBilat)
                  bilateralFilter.clipNotBilaterallyFilteredBlocks(srcYuv, resYuv, cs.slice->clpRng(COMPONENT_Y), currTU);
              }
            }
          }
        }
        else
        {
          // We do not use BIF so we can use the regular SAO function call
          offsetCTU(area, srcYuv, resYuv, reconParams[ctuRsAddr], cs);
        }
#else
        offsetCTU(area, srcYuv, resYuv, reconParams[ctuRsAddr], cs);
#endif
#endif
        ctuRsAddr++;
      }
    }
    //delete memory
#if JVET_V0094_BILATERAL_FILTER
    if (!(cs.pps->getUseBIF()) ||  !allBlksDisabled)
    {
#endif
    for (uint32_t i = 0; i< groupBlkStat.size(); i++)
    {
      for (uint32_t compIdx = 0; compIdx< MAX_NUM_COMPONENT; compIdx++)
      {
        delete[] groupBlkStat[i][compIdx];
      }
      delete[] groupBlkStat[i];
    }
    groupBlkStat.clear();
#if JVET_V0094_BILATERAL_FILTER
    }
#endif
  }
  if (!allBlksDisabled && (totalCost >= 0) && bTestSAODisableAtPictureLevel) //SAO has not beneficial in this case - disable it
  {
    for( ctuRsAddr = 0; ctuRsAddr < pcv.sizeInCtus; ctuRsAddr++)
    {
      codedParams[ctuRsAddr].reset();
    }

    for (uint32_t componentIndex = 0; componentIndex < MAX_NUM_COMPONENT; componentIndex++)
    {
      sliceEnabled[componentIndex] = false;
    }
    m_CABACEstimator->getCtx() = SAOCtx(ctxPicStart);
  }

  EncSampleAdaptiveOffset::disabledRate( cs, reconParams, saoEncodingRate, saoEncodingRateChroma );
  
#if JVET_V0094_BILATERAL_FILTER
  bilateralFilter.destroy();
#endif
}

void EncSampleAdaptiveOffset::disabledRate( CodingStructure& cs, SAOBlkParam* reconParams, const double saoEncodingRate, const double saoEncodingRateChroma )
{
  if (saoEncodingRate > 0.0)
  {
    const PreCalcValues& pcv = *cs.pcv;
    const uint32_t numberOfComponents = getNumberValidComponents( cs.picture->chromaFormat );
    int picTempLayer = cs.slice->getDepth();
    int numCtusForSAOOff[MAX_NUM_COMPONENT];

    for (int compIdx = 0; compIdx < numberOfComponents; compIdx++)
    {
      numCtusForSAOOff[compIdx] = 0;
      for( int ctuRsAddr=0; ctuRsAddr< pcv.sizeInCtus; ctuRsAddr++)
      {
        if( reconParams[ctuRsAddr][compIdx].modeIdc == SAO_MODE_OFF)
        {
          numCtusForSAOOff[compIdx]++;
        }
      }
    }
    if (saoEncodingRateChroma > 0.0)
    {
      for (int compIdx = 0; compIdx < numberOfComponents; compIdx++)
      {
        m_saoDisabledRate[compIdx][picTempLayer] = (double)numCtusForSAOOff[compIdx]/(double)pcv.sizeInCtus;
      }
    }
    else if (picTempLayer == 0)
    {
      m_saoDisabledRate[COMPONENT_Y][0] = (double)(numCtusForSAOOff[COMPONENT_Y]+numCtusForSAOOff[COMPONENT_Cb]+numCtusForSAOOff[COMPONENT_Cr])/(double)(pcv.sizeInCtus *3);
    }
  }
}

void EncSampleAdaptiveOffset::getBlkStats(const ComponentID compIdx, const int channelBitDepth, SAOStatData* statsDataTypes
                        , Pel* srcBlk, Pel* orgBlk,
#if JVET_V0094_BILATERAL_FILTER
                          Pel* bifBlk, int bifStride,
#endif
                                          int srcStride, int orgStride, int width, int height
                        , bool isLeftAvail,  bool isRightAvail, bool isAboveAvail, bool isBelowAvail, bool isAboveLeftAvail, bool isAboveRightAvail
                        , bool isCalculatePreDeblockSamples
                        , bool isCtuCrossedByVirtualBoundaries, int horVirBndryPos[], int verVirBndryPos[], int numHorVirBndry, int numVerVirBndry
                        )
{
  int x,y, startX, startY, endX, endY, edgeType, firstLineStartX, firstLineEndX;
  int8_t signLeft, signRight, signDown;
  int64_t *diff, *count;
  Pel *srcLine, *orgLine;
#if JVET_V0094_BILATERAL_FILTER
  Pel *bifLine;
#endif
  int* skipLinesR = m_skipLinesR[compIdx];
  int* skipLinesB = m_skipLinesB[compIdx];

  for(int typeIdx=0; typeIdx< NUM_SAO_NEW_TYPES; typeIdx++)
  {
    SAOStatData& statsData= statsDataTypes[typeIdx];
    statsData.reset();

    srcLine = srcBlk;
    orgLine = orgBlk;
#if JVET_V0094_BILATERAL_FILTER
    bifLine = bifBlk;
#endif
    diff    = statsData.diff;
    count   = statsData.count;
    switch(typeIdx)
    {
    case SAO_TYPE_EO_0:
      {
        diff +=2;
        count+=2;
        endY   = (isBelowAvail) ? (height - skipLinesB[typeIdx]) : height;
        startX = (!isCalculatePreDeblockSamples) ? (isLeftAvail  ? 0 : 1)
                                                 : (isRightAvail ? (width - skipLinesR[typeIdx]) : (width - 1))
                                                 ;
        endX   = (!isCalculatePreDeblockSamples) ? (isRightAvail ? (width - skipLinesR[typeIdx]) : (width - 1))
                                                 : (isRightAvail ? width : (width - 1))
                                                 ;
        for (y=0; y<endY; y++)
        {
          signLeft = (int8_t)sgn(srcLine[startX] - srcLine[startX-1]);
          for (x=startX; x<endX; x++)
          {
            signRight =  (int8_t)sgn(srcLine[x] - srcLine[x+1]);
            if (isCtuCrossedByVirtualBoundaries && isProcessDisabled(x, y, numVerVirBndry, 0, verVirBndryPos, horVirBndryPos))
            {
              signLeft = -signRight;
              continue;
            }
            edgeType  =  signRight + signLeft;
            signLeft  = -signRight;

#if JVET_V0094_BILATERAL_FILTER
            diff [edgeType] += (orgLine[x] - bifLine[x]);
#else
            diff [edgeType] += (orgLine[x] - srcLine[x]);
#endif
            count[edgeType] ++;
          }
          srcLine  += srcStride;
          orgLine  += orgStride;
#if JVET_V0094_BILATERAL_FILTER
          bifLine  += bifStride;
#endif
        }
        if(isCalculatePreDeblockSamples)
        {
          if(isBelowAvail)
          {
            startX = isLeftAvail  ? 0 : 1;
            endX   = isRightAvail ? width : (width -1);

            for(y=0; y<skipLinesB[typeIdx]; y++)
            {
              signLeft = (int8_t)sgn(srcLine[startX] - srcLine[startX-1]);
              for (x=startX; x<endX; x++)
              {
                signRight =  (int8_t)sgn(srcLine[x] - srcLine[x+1]);
                if (isCtuCrossedByVirtualBoundaries && isProcessDisabled(x, endY + y, numVerVirBndry, 0, verVirBndryPos, horVirBndryPos))
                {
                  signLeft = -signRight;
                  continue;
                }
                edgeType  =  signRight + signLeft;
                signLeft  = -signRight;

#if JVET_V0094_BILATERAL_FILTER
                diff [edgeType] += (orgLine[x] - bifLine[x]);
#else
                diff [edgeType] += (orgLine[x] - srcLine[x]);
#endif
                count[edgeType] ++;
              }
              srcLine  += srcStride;
              orgLine  += orgStride;
#if JVET_V0094_BILATERAL_FILTER
              bifLine  += bifStride;
#endif
            }
          }
        }
      }
      break;
    case SAO_TYPE_EO_90:
      {
        diff +=2;
        count+=2;
        int8_t *signUpLine = &m_signLineBuf1[0];

        startX = (!isCalculatePreDeblockSamples) ? 0
                                                 : (isRightAvail ? (width - skipLinesR[typeIdx]) : width)
                                                 ;
        startY = isAboveAvail ? 0 : 1;
        endX   = (!isCalculatePreDeblockSamples) ? (isRightAvail ? (width - skipLinesR[typeIdx]) : width)
                                                 : width
                                                 ;
        endY   = isBelowAvail ? (height - skipLinesB[typeIdx]) : (height - 1);
        if (!isAboveAvail)
        {
          srcLine += srcStride;
          orgLine += orgStride;
#if JVET_V0094_BILATERAL_FILTER
          bifLine += bifStride;
#endif
        }

        Pel* srcLineAbove = srcLine - srcStride;
        for (x=startX; x<endX; x++)
        {
          signUpLine[x] = (int8_t)sgn(srcLine[x] - srcLineAbove[x]);
        }

        Pel* srcLineBelow;
        for (y=startY; y<endY; y++)
        {
          srcLineBelow = srcLine + srcStride;

          for (x=startX; x<endX; x++)
          {
            signDown  = (int8_t)sgn(srcLine[x] - srcLineBelow[x]);
            if (isCtuCrossedByVirtualBoundaries && isProcessDisabled(x, y, 0, numHorVirBndry, verVirBndryPos, horVirBndryPos))
            {
              signUpLine[x] = -signDown;
              continue;
            }
            edgeType  = signDown + signUpLine[x];
            signUpLine[x]= -signDown;

#if JVET_V0094_BILATERAL_FILTER
            diff [edgeType] += (orgLine[x] - bifLine[x]);
#else
            diff [edgeType] += (orgLine[x] - srcLine[x]);

#endif
            count[edgeType] ++;
          }
          srcLine += srcStride;
          orgLine += orgStride;
#if JVET_V0094_BILATERAL_FILTER
          bifLine += bifStride;
#endif
        }
        if(isCalculatePreDeblockSamples)
        {
          if(isBelowAvail)
          {
            startX = 0;
            endX   = width;

            for(y=0; y<skipLinesB[typeIdx]; y++)
            {
              srcLineBelow = srcLine + srcStride;
              srcLineAbove = srcLine - srcStride;

              for (x=startX; x<endX; x++)
              {
                if (isCtuCrossedByVirtualBoundaries && isProcessDisabled(x, y + endY, 0, numHorVirBndry, verVirBndryPos, horVirBndryPos))
                {
                  continue;
                }
                edgeType = sgn(srcLine[x] - srcLineBelow[x]) + sgn(srcLine[x] - srcLineAbove[x]);
#if JVET_V0094_BILATERAL_FILTER
                diff [edgeType] += (orgLine[x] - bifLine[x]);
#else
                diff [edgeType] += (orgLine[x] - srcLine[x]);
#endif
                count[edgeType] ++;
              }
              srcLine  += srcStride;
              orgLine  += orgStride;
#if JVET_V0094_BILATERAL_FILTER
              bifLine  += bifStride;
#endif
            }
          }
        }

      }
      break;
    case SAO_TYPE_EO_135:
      {
        diff +=2;
        count+=2;
        int8_t *signUpLine, *signDownLine, *signTmpLine;

        signUpLine  = &m_signLineBuf1[0];
        signDownLine= &m_signLineBuf2[0];

        startX = (!isCalculatePreDeblockSamples) ? (isLeftAvail  ? 0 : 1)
                                                 : (isRightAvail ? (width - skipLinesR[typeIdx]) : (width - 1))
                                                 ;

        endX   = (!isCalculatePreDeblockSamples) ? (isRightAvail ? (width - skipLinesR[typeIdx]): (width - 1))
                                                 : (isRightAvail ? width : (width - 1))
                                                 ;
        endY   = isBelowAvail ? (height - skipLinesB[typeIdx]) : (height - 1);

        //prepare 2nd line's upper sign
        Pel* srcLineBelow = srcLine + srcStride;
        for (x=startX; x<endX+1; x++)
        {
          signUpLine[x] = (int8_t)sgn(srcLineBelow[x] - srcLine[x-1]);
        }

        //1st line
        Pel* srcLineAbove = srcLine - srcStride;
        firstLineStartX = (!isCalculatePreDeblockSamples) ? (isAboveLeftAvail ? 0    : 1) : startX;
        firstLineEndX   = (!isCalculatePreDeblockSamples) ? (isAboveAvail     ? endX : 1) : endX;
        for(x=firstLineStartX; x<firstLineEndX; x++)
        {
          if (isCtuCrossedByVirtualBoundaries && isProcessDisabled(x, 0, numVerVirBndry, numHorVirBndry, verVirBndryPos, horVirBndryPos))
          {
            continue;
          }
          edgeType = sgn(srcLine[x] - srcLineAbove[x-1]) - signUpLine[x+1];
#if JVET_V0094_BILATERAL_FILTER
          diff [edgeType] += (orgLine[x] - bifLine[x]);
#else
          diff [edgeType] += (orgLine[x] - srcLine[x]);
#endif
          count[edgeType] ++;
        }
        srcLine  += srcStride;
        orgLine  += orgStride;
#if JVET_V0094_BILATERAL_FILTER
        bifLine  += bifStride;
#endif
        //middle lines
        for (y=1; y<endY; y++)
        {
          srcLineBelow = srcLine + srcStride;

          for (x=startX; x<endX; x++)
          {
            signDown = (int8_t)sgn(srcLine[x] - srcLineBelow[x+1]);
            if (isCtuCrossedByVirtualBoundaries && isProcessDisabled(x, y, numVerVirBndry, numHorVirBndry, verVirBndryPos, horVirBndryPos))
            {
              signDownLine[x + 1] = -signDown;
              continue;
            }
            edgeType = signDown + signUpLine[x];
#if JVET_V0094_BILATERAL_FILTER
            diff [edgeType] += (orgLine[x] - bifLine[x]);
#else
            diff [edgeType] += (orgLine[x] - srcLine[x]);
#endif
            count[edgeType] ++;

            signDownLine[x+1] = -signDown;
          }
          signDownLine[startX] = (int8_t)sgn(srcLineBelow[startX] - srcLine[startX-1]);

          signTmpLine  = signUpLine;
          signUpLine   = signDownLine;
          signDownLine = signTmpLine;

          srcLine += srcStride;
          orgLine += orgStride;
#if JVET_V0094_BILATERAL_FILTER
          bifLine += bifStride;
#endif
        }
        if(isCalculatePreDeblockSamples)
        {
          if(isBelowAvail)
          {
            startX = isLeftAvail  ? 0     : 1 ;
            endX   = isRightAvail ? width : (width -1);

            for(y=0; y<skipLinesB[typeIdx]; y++)
            {
              srcLineBelow = srcLine + srcStride;
              srcLineAbove = srcLine - srcStride;

              for (x=startX; x< endX; x++)
              {
                if (isCtuCrossedByVirtualBoundaries && isProcessDisabled(x, y + endY, numVerVirBndry, numHorVirBndry, verVirBndryPos, horVirBndryPos))
                {
                  continue;
                }
                edgeType = sgn(srcLine[x] - srcLineBelow[x+1]) + sgn(srcLine[x] - srcLineAbove[x-1]);
#if JVET_V0094_BILATERAL_FILTER
                diff [edgeType] += (orgLine[x] - bifLine[x]);
#else
                diff [edgeType] += (orgLine[x] - srcLine[x]);
#endif
                count[edgeType] ++;
              }
              srcLine  += srcStride;
              orgLine  += orgStride;
#if JVET_V0094_BILATERAL_FILTER
              bifLine  += bifStride;
#endif
            }
          }
        }
      }
      break;
    case SAO_TYPE_EO_45:
      {
        diff +=2;
        count+=2;
        int8_t *signUpLine = &m_signLineBuf1[1];

        startX = (!isCalculatePreDeblockSamples) ? (isLeftAvail  ? 0 : 1)
                                                 : (isRightAvail ? (width - skipLinesR[typeIdx]) : (width - 1))
                                                 ;
        endX   = (!isCalculatePreDeblockSamples) ? (isRightAvail ? (width - skipLinesR[typeIdx]) : (width - 1))
                                                 : (isRightAvail ? width : (width - 1))
                                                 ;
        endY   = isBelowAvail ? (height - skipLinesB[typeIdx]) : (height - 1);

        //prepare 2nd line upper sign
        Pel* srcLineBelow = srcLine + srcStride;
        for (x=startX-1; x<endX; x++)
        {
          signUpLine[x] = (int8_t)sgn(srcLineBelow[x] - srcLine[x+1]);
        }


        //first line
        Pel* srcLineAbove = srcLine - srcStride;
        firstLineStartX = (!isCalculatePreDeblockSamples) ? (isAboveAvail ? startX : endX)
                                                          : startX
                                                          ;
        firstLineEndX   = (!isCalculatePreDeblockSamples) ? ((!isRightAvail && isAboveRightAvail) ? width : endX)
                                                          : endX
                                                          ;
        for(x=firstLineStartX; x<firstLineEndX; x++)
        {
          if (isCtuCrossedByVirtualBoundaries && isProcessDisabled(x, 0, numVerVirBndry, numHorVirBndry, verVirBndryPos, horVirBndryPos))
          {
            continue;
          }
          edgeType = sgn(srcLine[x] - srcLineAbove[x+1]) - signUpLine[x-1];
#if JVET_V0094_BILATERAL_FILTER
          diff [edgeType] += (orgLine[x] - bifLine[x]);
#else
          diff [edgeType] += (orgLine[x] - srcLine[x]);
#endif
          count[edgeType] ++;
        }

        srcLine += srcStride;
        orgLine += orgStride;
#if JVET_V0094_BILATERAL_FILTER
        bifLine += bifStride;
#endif
        //middle lines
        for (y=1; y<endY; y++)
        {
          srcLineBelow = srcLine + srcStride;

          for(x=startX; x<endX; x++)
          {
            signDown = (int8_t)sgn(srcLine[x] - srcLineBelow[x-1]);
            if (isCtuCrossedByVirtualBoundaries && isProcessDisabled(x, y, numVerVirBndry, numHorVirBndry, verVirBndryPos, horVirBndryPos))
            {
              signUpLine[x - 1] = -signDown;
              continue;
            }
            edgeType = signDown + signUpLine[x];
#if JVET_V0094_BILATERAL_FILTER
            diff [edgeType] += (orgLine[x] - bifLine[x]);
#else
            diff [edgeType] += (orgLine[x] - srcLine[x]);
#endif
            count[edgeType] ++;

            signUpLine[x-1] = -signDown;
          }
          signUpLine[endX-1] = (int8_t)sgn(srcLineBelow[endX-1] - srcLine[endX]);
          srcLine  += srcStride;
          orgLine  += orgStride;
#if JVET_V0094_BILATERAL_FILTER
          bifLine  += bifStride;
#endif
        }
        if(isCalculatePreDeblockSamples)
        {
          if(isBelowAvail)
          {
            startX = isLeftAvail  ? 0     : 1 ;
            endX   = isRightAvail ? width : (width -1);

            for(y=0; y<skipLinesB[typeIdx]; y++)
            {
              srcLineBelow = srcLine + srcStride;
              srcLineAbove = srcLine - srcStride;

              for (x=startX; x<endX; x++)
              {
                if (isCtuCrossedByVirtualBoundaries && isProcessDisabled(x, y + endY, numVerVirBndry, numHorVirBndry, verVirBndryPos, horVirBndryPos))
                {
                  continue;
                }
                edgeType = sgn(srcLine[x] - srcLineBelow[x-1]) + sgn(srcLine[x] - srcLineAbove[x+1]);
#if JVET_V0094_BILATERAL_FILTER
                diff [edgeType] += (orgLine[x] - bifLine[x]);
#else
                diff [edgeType] += (orgLine[x] - srcLine[x]);
#endif
                count[edgeType] ++;
              }
              srcLine  += srcStride;
              orgLine  += orgStride;
#if JVET_V0094_BILATERAL_FILTER
              bifLine  += bifStride;
#endif
            }
          }
        }
      }
      break;
    case SAO_TYPE_BO:
      {
        startX = (!isCalculatePreDeblockSamples)?0
                                                :( isRightAvail?(width- skipLinesR[typeIdx]):width)
                                                ;
        endX   = (!isCalculatePreDeblockSamples)?(isRightAvail ? (width - skipLinesR[typeIdx]) : width )
                                                :width
                                                ;
        endY = isBelowAvail ? (height- skipLinesB[typeIdx]) : height;
        int shiftBits = channelBitDepth - NUM_SAO_BO_CLASSES_LOG2;
        for (y=0; y< endY; y++)
        {
          for (x=startX; x< endX; x++)
          {

            int bandIdx= srcLine[x] >> shiftBits;
#if JVET_V0094_BILATERAL_FILTER
            diff [bandIdx] += (orgLine[x] - bifLine[x]);
#else
            diff [bandIdx] += (orgLine[x] - srcLine[x]);
#endif
            count[bandIdx] ++;
          }
          srcLine += srcStride;
          orgLine += orgStride;
#if JVET_V0094_BILATERAL_FILTER
          bifLine += bifStride;
#endif
        }
        if(isCalculatePreDeblockSamples)
        {
          if(isBelowAvail)
          {
            startX = 0;
            endX   = width;

            for(y= 0; y< skipLinesB[typeIdx]; y++)
            {
              for (x=startX; x< endX; x++)
              {
                int bandIdx= srcLine[x] >> shiftBits;
#if JVET_V0094_BILATERAL_FILTER
                diff [bandIdx] += (orgLine[x] - bifLine[x]);
#else
                diff [bandIdx] += (orgLine[x] - srcLine[x]);
#endif
                count[bandIdx] ++;
              }
              srcLine  += srcStride;
              orgLine  += orgStride;
#if JVET_V0094_BILATERAL_FILTER
              bifLine  += bifStride;
#endif
            }

          }
        }
      }
      break;
    default:
      {
        THROW("Not a supported SAO type");
      }
    }
  }
}

#if JVET_W0066_CCSAO
void EncSampleAdaptiveOffset::CCSAOProcess(CodingStructure& cs, const double* lambdas, const int intraPeriod)
{
  PelUnitBuf orgYuv = cs.getOrgBuf(); 
  PelUnitBuf dstYuv = cs.getRecoBuf();
  PelUnitBuf srcYuv = m_ccSaoBuf.getBuf( cs.area );
  srcYuv.extendBorderPel( MAX_CCSAO_FILTER_LENGTH >> 1 );
  m_intraPeriod = intraPeriod;

  setupCcSaoLambdas(cs, lambdas);

  if (cs.slice->getSPS()->getCCSAOEnabledFlag())
  {
    const TempCtx ctxStartCcSao(m_CtxCache, SubCtx(Ctx::CcSaoControlIdc, m_CABACEstimator->getCtx()));
    m_CABACEstimator->getCtx() = SubCtx(Ctx::CcSaoControlIdc, ctxStartCcSao); deriveCcSao(cs, COMPONENT_Y,  orgYuv, srcYuv, dstYuv);
    m_CABACEstimator->getCtx() = SubCtx(Ctx::CcSaoControlIdc, ctxStartCcSao); deriveCcSao(cs, COMPONENT_Cb, orgYuv, srcYuv, dstYuv);
    m_CABACEstimator->getCtx() = SubCtx(Ctx::CcSaoControlIdc, ctxStartCcSao); deriveCcSao(cs, COMPONENT_Cr, orgYuv, srcYuv, dstYuv);
    applyCcSao(cs, *cs.pcv, srcYuv, dstYuv);
  }
}

void EncSampleAdaptiveOffset::setupCcSaoLambdas(CodingStructure& cs, const double* lambdas)
{
  m_lambda[COMPONENT_Y ] = m_picWidth * m_picHeight <= 416 * 240 
                         ? lambdas[COMPONENT_Y ] * 4.0 
                         : lambdas[COMPONENT_Y ];
  m_lambda[COMPONENT_Cb] = lambdas[COMPONENT_Cb];
  m_lambda[COMPONENT_Cr] = lambdas[COMPONENT_Cr];
}

void EncSampleAdaptiveOffset::deriveCcSao(CodingStructure& cs, const ComponentID compID, const CPelUnitBuf& orgYuv, const CPelUnitBuf& srcYuv, const CPelUnitBuf& dstYuv)
{
  double bestCost = 0;
  double tempCost = 0;
  double bestCostS[MAX_CCSAO_SET_NUM + 1] = { 0 };

  double bestCostG[17] = { 0 };
  int    classNumG[17] = { 0 };
  int    stageNum = m_intraPeriod == 1 ? MAX_CCSAO_CLASS_NUM / 4 : MAX_CCSAO_CLASS_NUM / 16;
  for (int stage = 1; stage <= stageNum; stage++)
    classNumG[stage] = stage * (MAX_CCSAO_CLASS_NUM / stageNum);

  m_bestCcSaoParam.reset();
  memset(m_bestCcSaoControl, 0, sizeof(uint8_t) * m_numCTUsInPic);

  for (int setNum = 1; setNum <= MAX_CCSAO_SET_NUM; setNum++)
  {
    if (setNum > 1)
    {
      getCcSaoStatistics(cs, compID, orgYuv, srcYuv, dstYuv, m_ccSaoStatData, m_bestCcSaoParam);
    }
    setupInitCcSaoParam(cs, compID, setNum, m_trainingDistortion, m_ccSaoStatData, m_ccSaoStatFrame,
                        m_initCcSaoParam, m_bestCcSaoParam, m_initCcSaoControl, m_bestCcSaoControl);

    for (int stage = 1; stage <= stageNum; stage++)
    {
      for (int bandNumY = 1; bandNumY <= MAX_CCSAO_BAND_NUM_Y; bandNumY++)
      for (int bandNumU = 1; bandNumU <= MAX_CCSAO_BAND_NUM_U; bandNumU++)
      for (int bandNumV = 1; bandNumV <= MAX_CCSAO_BAND_NUM_V; bandNumV++)
      for (int candPosY = 0; candPosY <  MAX_CCSAO_CAND_POS_Y && bandNumY > 1; candPosY++)
      {
        if (bandNumY < bandNumU || bandNumY < bandNumV)
          continue;

        int classNum = bandNumY * bandNumU * bandNumV;
        if (classNum > MAX_CCSAO_CLASS_NUM)
          continue;

        if (classNum <= classNumG[stage - 1] || classNum > classNumG[stage])
          continue;

        setupTempCcSaoParam(cs, compID, setNum, candPosY, bandNumY, bandNumU, bandNumV, m_tempCcSaoParam, m_initCcSaoParam, m_tempCcSaoControl, m_initCcSaoControl);
        getCcSaoStatistics(cs, compID, orgYuv, srcYuv, dstYuv, m_ccSaoStatData, m_tempCcSaoParam);
        deriveCcSaoRDO(cs, compID, m_trainingDistortion, m_ccSaoStatData, m_ccSaoStatFrame,
                       m_bestCcSaoParam, m_tempCcSaoParam, m_bestCcSaoControl, m_tempCcSaoControl, bestCost, tempCost);
      }

      bestCostG[stage] = bestCost;
      if (bestCostG[stage] >= bestCostG[stage - 1])
        break;
    }

    bestCostS[setNum] = bestCost;
    if (bestCostS[setNum] >= bestCostS[setNum - 1])
      break;
  }

  bool oneBlockFiltered = false;
  for (int ctbIdx = 0; m_bestCcSaoParam.setNum > 0 && ctbIdx < m_numCTUsInPic; ctbIdx++)
  {
    if (m_bestCcSaoControl[ctbIdx])
    {
      oneBlockFiltered = true;
      break;
    }
  }
  
  m_ccSaoComParam.reset(compID);
  memset(m_ccSaoControl[compID], 0, sizeof(uint8_t) * m_numCTUsInPic);

  m_ccSaoComParam.enabled[compID] = oneBlockFiltered;
  if (oneBlockFiltered)
  {
    CcSaoEncParam storedBestCcSaoParam = m_bestCcSaoParam;
    memcpy(m_tempCcSaoControl, m_bestCcSaoControl, sizeof(uint8_t) * m_numCTUsInPic);

    int setNum = 0;
    for (int setIdx = 0; setIdx < MAX_CCSAO_SET_NUM; setIdx++)
    {
      uint8_t setIdc = m_bestCcSaoParam.mapIdxToIdc[setIdx];
      if (m_bestCcSaoParam.setEnabled[setIdx])
      {
        for (int ctbIdx = 0; ctbIdx < m_numCTUsInPic; ctbIdx++)
        {
          if (m_tempCcSaoControl[ctbIdx] == (setIdx + 1) )
          {
            m_bestCcSaoControl[ctbIdx] = setIdc;
          }
        }
        m_bestCcSaoParam.candPos[setIdc - 1][COMPONENT_Y ] = storedBestCcSaoParam.candPos[setIdx][COMPONENT_Y ];
        m_bestCcSaoParam.bandNum[setIdc - 1][COMPONENT_Y ] = storedBestCcSaoParam.bandNum[setIdx][COMPONENT_Y ];
        m_bestCcSaoParam.bandNum[setIdc - 1][COMPONENT_Cb] = storedBestCcSaoParam.bandNum[setIdx][COMPONENT_Cb];
        m_bestCcSaoParam.bandNum[setIdc - 1][COMPONENT_Cr] = storedBestCcSaoParam.bandNum[setIdx][COMPONENT_Cr];
        memcpy(m_bestCcSaoParam.offset[setIdc - 1], storedBestCcSaoParam.offset[setIdx], sizeof(storedBestCcSaoParam.offset[setIdx]));
        setNum++;
      }
      m_bestCcSaoParam.setEnabled[setIdx] = setIdx < m_bestCcSaoParam.setNum ? true : false;
    }
    CHECK(setNum != m_bestCcSaoParam.setNum, "Number of sets enabled != setNum");

    m_ccSaoComParam.setNum [compID] = m_bestCcSaoParam.setNum;

    for ( int setIdx = 0; setIdx < m_bestCcSaoParam.setNum; setIdx++ )
    {
      m_ccSaoComParam.setEnabled[compID][setIdx]               = m_bestCcSaoParam.setEnabled[setIdx];
      m_ccSaoComParam.candPos   [compID][setIdx][COMPONENT_Y ] = m_bestCcSaoParam.candPos   [setIdx][COMPONENT_Y ];
      m_ccSaoComParam.bandNum   [compID][setIdx][COMPONENT_Y ] = m_bestCcSaoParam.bandNum   [setIdx][COMPONENT_Y ];
      m_ccSaoComParam.bandNum   [compID][setIdx][COMPONENT_Cb] = m_bestCcSaoParam.bandNum   [setIdx][COMPONENT_Cb];
      m_ccSaoComParam.bandNum   [compID][setIdx][COMPONENT_Cr] = m_bestCcSaoParam.bandNum   [setIdx][COMPONENT_Cr];
      memcpy(m_ccSaoComParam.offset[compID][setIdx], m_bestCcSaoParam.offset[setIdx], sizeof(m_bestCcSaoParam.offset[setIdx]));
    }
    memcpy(m_ccSaoControl[compID], m_bestCcSaoControl, sizeof(uint8_t) * m_numCTUsInPic);
  }
}

void EncSampleAdaptiveOffset::setupInitCcSaoParam(CodingStructure& cs, const ComponentID compID, const int setNum, int64_t* trainingDistortion[MAX_CCSAO_SET_NUM]
                                                , CcSaoStatData* blkStats[MAX_CCSAO_SET_NUM], CcSaoStatData frameStats[MAX_CCSAO_SET_NUM]
                                                , CcSaoEncParam& initCcSaoParam, CcSaoEncParam& bestCcSaoParam
                                                , uint8_t* initCcSaoControl, uint8_t* bestCcSaoControl)
{
  initCcSaoParam.reset();
  memset(initCcSaoControl, 0, sizeof(uint8_t) * m_numCTUsInPic);

  if (setNum == 1)
  {
    std::fill_n(initCcSaoControl, m_numCTUsInPic, 1);
    return;
  }

  for (int setIdx = 0; setIdx < MAX_CCSAO_SET_NUM; setIdx++)
  {
    if (bestCcSaoParam.setEnabled[setIdx])
    {
      getCcSaoFrameStats(compID, setIdx, bestCcSaoControl, blkStats, frameStats);
      getCcSaoDistortion(compID, setIdx, blkStats, bestCcSaoParam.offset, trainingDistortion);
    }
  }

  initCcSaoParam = bestCcSaoParam;

  int ctbCntOn = 0;
  CtbCost *ctbCost = new CtbCost[m_numCTUsInPic];

  for (int ctbIdx = 0; ctbIdx < m_numCTUsInPic; ctbIdx++)
  {
    int64_t dist = 0;

    if (bestCcSaoControl[ctbIdx])
    {
      int setIdx = bestCcSaoControl[ctbIdx] - 1;
      dist = trainingDistortion[setIdx][ctbIdx];
      ctbCntOn++;
    }

    ctbCost[ctbIdx].pos = ctbIdx;
    ctbCost[ctbIdx].cost = (double)dist;
  }

  std::stable_sort(ctbCost, ctbCost + m_numCTUsInPic, compareCtbCost);

  for (int ctbIdx = 0; ctbIdx < m_numCTUsInPic; ctbIdx++)
  {
    int ctbPos = ctbCost[ctbIdx].pos;
    if (ctbIdx < ctbCntOn)
    {
      if (ctbIdx * 2 > ctbCntOn)
      {
        initCcSaoControl[ctbPos] = setNum;
      }
      else
      {
        initCcSaoControl[ctbPos] = bestCcSaoControl[ctbPos];
      }
    }
    else
    {
      initCcSaoControl[ctbPos] = 0;
    }
  }
  
  delete[] ctbCost;
  ctbCost = nullptr;
}

void EncSampleAdaptiveOffset::setupTempCcSaoParam(CodingStructure& cs, const ComponentID compID, const int setNum
                                                , const int candPosY, const int bandNumY, const int bandNumU, const int bandNumV
                                                , CcSaoEncParam& tempCcSaoParam, CcSaoEncParam& initCcSaoParam
                                                , uint8_t* tempCcSaoControl, uint8_t* initCcSaoControl)
{
  tempCcSaoParam.reset();
  memset(tempCcSaoControl, 0, sizeof(uint8_t) * m_numCTUsInPic);

  tempCcSaoParam = initCcSaoParam;;
  memcpy(tempCcSaoControl, initCcSaoControl, sizeof(uint8_t) * m_numCTUsInPic);

  tempCcSaoParam.setNum = setNum;
  tempCcSaoParam.setEnabled[setNum - 1] = true;
  tempCcSaoParam.candPos   [setNum - 1][COMPONENT_Y ] = candPosY;
  tempCcSaoParam.bandNum   [setNum - 1][COMPONENT_Y ] = bandNumY;
  tempCcSaoParam.bandNum   [setNum - 1][COMPONENT_Cb] = bandNumU;
  tempCcSaoParam.bandNum   [setNum - 1][COMPONENT_Cr] = bandNumV;

  for (int setIdx = 0; setIdx <= setNum; setIdx++)
  {
    tempCcSaoParam.mapIdxToIdc[setIdx] = setIdx < setNum ? setIdx + 1 : 0;
  }
}

void EncSampleAdaptiveOffset::getCcSaoStatistics(CodingStructure& cs, const ComponentID compID
                                               , const CPelUnitBuf& orgYuv, const CPelUnitBuf& srcYuv, const CPelUnitBuf& dstYuv
                                               , CcSaoStatData* blkStats[MAX_CCSAO_SET_NUM], const CcSaoEncParam& ccSaoParam)
{
  bool isLeftAvail, isRightAvail, isAboveAvail, isBelowAvail, isAboveLeftAvail, isAboveRightAvail;

  const PreCalcValues& pcv = *cs.pcv;

  int ctuRsAddr = 0;
  for( uint32_t yPos = 0; yPos < pcv.lumaHeight; yPos += pcv.maxCUHeight )
  {
    for( uint32_t xPos = 0; xPos < pcv.lumaWidth; xPos += pcv.maxCUWidth )
    {
      const uint32_t width  = (xPos + pcv.maxCUWidth  > pcv.lumaWidth)  ? (pcv.lumaWidth - xPos)  : pcv.maxCUWidth;
      const uint32_t height = (yPos + pcv.maxCUHeight > pcv.lumaHeight) ? (pcv.lumaHeight - yPos) : pcv.maxCUHeight;
      const UnitArea area( cs.area.chromaFormat, Area(xPos , yPos, width, height) );

      deriveLoopFilterBoundaryAvailibility(cs, area.Y(), isLeftAvail, isAboveAvail, isAboveLeftAvail );

      //NOTE: The number of skipped lines during gathering CTU statistics depends on the slice boundary availabilities.
      //For simplicity, here only picture boundaries are considered.

      isRightAvail      = (xPos + pcv.maxCUWidth  < pcv.lumaWidth );
      isBelowAvail      = (yPos + pcv.maxCUHeight < pcv.lumaHeight);
      isAboveRightAvail = ((yPos > 0) && (isRightAvail));

      for (int setIdx = 0; setIdx < MAX_CCSAO_SET_NUM; setIdx++)
      {
        blkStats[setIdx][ctuRsAddr].reset();
        if (!ccSaoParam.setEnabled[setIdx])
          continue;

        const CompArea   &compArea   = area.block(compID);
        const int         srcStrideY = srcYuv.get(COMPONENT_Y ).stride;
        const int         srcStrideU = srcYuv.get(COMPONENT_Cb).stride;
        const int         srcStrideV = srcYuv.get(COMPONENT_Cr).stride;
        const int         dstStride  = dstYuv.get(compID      ).stride;
        const int         orgStride  = orgYuv.get(compID      ).stride;
        const Pel        *srcBlkY    = srcYuv.get(COMPONENT_Y ).bufAt(area.block(COMPONENT_Y ));
        const Pel        *srcBlkU    = srcYuv.get(COMPONENT_Cb).bufAt(area.block(COMPONENT_Cb));
        const Pel        *srcBlkV    = srcYuv.get(COMPONENT_Cr).bufAt(area.block(COMPONENT_Cr));
        const Pel        *dstBlk     = dstYuv.get(compID      ).bufAt(compArea);
        const Pel        *orgBlk     = orgYuv.get(compID      ).bufAt(compArea);

        const uint16_t    candPosY   = ccSaoParam.candPos[setIdx][COMPONENT_Y ];
        const uint16_t    bandNumY   = ccSaoParam.bandNum[setIdx][COMPONENT_Y ];
        const uint16_t    bandNumU   = ccSaoParam.bandNum[setIdx][COMPONENT_Cb];
        const uint16_t    bandNumV   = ccSaoParam.bandNum[setIdx][COMPONENT_Cr];

        getCcSaoBlkStats(compID, cs.area.chromaFormat, cs.sps->getBitDepth(toChannelType(compID))
                       , setIdx, blkStats, ctuRsAddr
                       , candPosY, bandNumY, bandNumU, bandNumV
                       , srcBlkY, srcBlkU, srcBlkV, orgBlk, dstBlk
                       , srcStrideY, srcStrideU, srcStrideV, orgStride, dstStride, compArea.width, compArea.height
                       , isLeftAvail, isRightAvail, isAboveAvail, isBelowAvail, isAboveLeftAvail, isAboveRightAvail
                       );
      }
      ctuRsAddr++;
    }
  }
}

void EncSampleAdaptiveOffset::getCcSaoBlkStats(const ComponentID compID, const ChromaFormat chromaFormat, const int bitDepth
                                             , const int setIdx, CcSaoStatData* blkStats[MAX_CCSAO_SET_NUM], const int ctuRsAddr
                                             , const uint16_t candPosY
                                             , const uint16_t bandNumY, const uint16_t bandNumU, const uint16_t bandNumV
                                             , const Pel* srcY, const Pel* srcU, const Pel* srcV, const Pel* org, const Pel* dst
                                             , const int srcStrideY, const int srcStrideU, const int srcStrideV, const int orgStride, const int dstStride
                                             , const int width, const int height
                                             , bool isLeftAvail, bool isRightAvail, bool isAboveAvail, bool isBelowAvail, bool isAboveLeftAvail, bool isAboveRightAvail
                                             )
{
  const int candPosYX = g_ccSaoCandPosX[COMPONENT_Y][candPosY];
  const int candPosYY = g_ccSaoCandPosY[COMPONENT_Y][candPosY];

  switch(compID)
  {
  case COMPONENT_Y:
    {
      for (int y = 0; y < height; y++)
      {
        for (int x = 0; x < width; x++)
        {
          const Pel *colY = srcY +  x  + srcStrideY * candPosYY + candPosYX;
          const Pel *colU = srcU + (x >> 1);
          const Pel *colV = srcV + (x >> 1);

          const int bandY    = (*colY * bandNumY) >> bitDepth;
          const int bandU    = (*colU * bandNumU) >> bitDepth;
          const int bandV    = (*colV * bandNumV) >> bitDepth;
          const int bandIdx  = bandY * bandNumU * bandNumV
                             + bandU * bandNumV
                             + bandV;
          const int classIdx = bandIdx;

          blkStats[setIdx][ctuRsAddr].diff [classIdx] += org[x] - dst[x];
          blkStats[setIdx][ctuRsAddr].count[classIdx]++;
        }

        srcY += srcStrideY;
        srcU += srcStrideU * (y & 0x1);
        srcV += srcStrideV * (y & 0x1);
        org  += orgStride;
        dst  += dstStride;
      }
    }
    break;
  case COMPONENT_Cb:
  case COMPONENT_Cr:
    {
      for (int y = 0; y < height; y++)
      {
        for (int x = 0; x < width; x++)
        {
          const Pel *colY = srcY + (x << 1) + srcStrideY * candPosYY + candPosYX;
          const Pel *colU = srcU + x;
          const Pel *colV = srcV + x;

          const int bandY    = (*colY * bandNumY) >> bitDepth;
          const int bandU    = (*colU * bandNumU) >> bitDepth;
          const int bandV    = (*colV * bandNumV) >> bitDepth;
          const int bandIdx  = bandY * bandNumU * bandNumV
                             + bandU * bandNumV
                             + bandV;
          const int classIdx = bandIdx;

          blkStats[setIdx][ctuRsAddr].diff [classIdx] += org[x] - dst[x];
          blkStats[setIdx][ctuRsAddr].count[classIdx]++;
        }

        srcY += srcStrideY << 1;
        srcU += srcStrideU;
        srcV += srcStrideV;
        org  += orgStride;
        dst  += dstStride;
      }
    }
    break;
  default:
    {
      THROW("Not a supported CCSAO compID\n");
    }
  }
}

void EncSampleAdaptiveOffset::getCcSaoFrameStats(const ComponentID compID, const int setIdx, const uint8_t* ccSaoControl
                                               , CcSaoStatData* blkStats[MAX_CCSAO_SET_NUM], CcSaoStatData frameStats[MAX_CCSAO_SET_NUM])
{
  frameStats[setIdx].reset();
  int setIdc = setIdx + 1;

  for (int ctbIdx = 0; ctbIdx < m_numCTUsInPic; ctbIdx++)
  {
    if (ccSaoControl[ctbIdx] == setIdc)
    {
      frameStats[setIdx] += blkStats[setIdx][ctbIdx];
    }
  }
}

inline int EncSampleAdaptiveOffset::estCcSaoIterOffset(const double lambda, const int offsetInput, const int64_t count, const int64_t diffSum, const int shift, const int bitIncrease, int64_t& bestDist, double& bestCost, const int offsetTh)
{
  int iterOffset, tempOffset;
  int64_t tempDist, tempRate;
  double tempCost, tempMinCost;
  int offsetOutput = 0;
  iterOffset = offsetInput;
  // Assuming sending quantized value 0 results in zero offset and sending the value zero needs 1 bit. entropy coder can be used to measure the exact rate here.
  tempMinCost = lambda;
  while (iterOffset != 0)
  {
    // Calculate the bits required for signaling the offset
    tempRate = lengthUvlc(abs(iterOffset)) + (iterOffset == 0 ? 0 : 1);

    // Do the dequantization before distortion calculation
    tempOffset = iterOffset << bitIncrease;
    tempDist = estSaoDist(count, tempOffset, diffSum, shift);
    tempCost = ((double)tempDist + lambda * (double)tempRate);
    if (tempCost < tempMinCost)
    {
      tempMinCost = tempCost;
      offsetOutput = iterOffset;
      bestDist = tempDist;
      bestCost = tempCost;
    }
    iterOffset = (iterOffset > 0) ? (iterOffset - 1) : (iterOffset + 1);
  }
  return offsetOutput;
}

void EncSampleAdaptiveOffset::deriveCcSaoOffsets(const ComponentID compID, const int bitDepth, const int setIdx
                                               , CcSaoStatData frameStats[MAX_CCSAO_SET_NUM]
                                               , short offset[MAX_CCSAO_SET_NUM][MAX_CCSAO_CLASS_NUM])
{
  int quantOffsets[MAX_CCSAO_CLASS_NUM] = { 0 };

  for(int k = 0; k < MAX_CCSAO_CLASS_NUM; k++)
  {
    if(frameStats[setIdx].count[k] == 0)
      continue;

    quantOffsets[k] =
      (int) xRoundIbdi(bitDepth, (double)(frameStats[setIdx].diff [k] << DISTORTION_PRECISION_ADJUSTMENT(bitDepth))
                               / (double)(frameStats[setIdx].count[k]));
    quantOffsets[k] = Clip3(-MAX_CCSAO_OFFSET_THR, MAX_CCSAO_OFFSET_THR, quantOffsets[k]);
  }

  int64_t dist[MAX_CCSAO_CLASS_NUM] = { 0 };
  double  cost[MAX_CCSAO_CLASS_NUM] = { 0 };
  for (int k = 0; k < MAX_CCSAO_CLASS_NUM; k++)
  {
    cost[k] = m_lambda[compID];
    if (quantOffsets[k] != 0)
    {
      quantOffsets[k] = estCcSaoIterOffset(m_lambda[compID], quantOffsets[k], frameStats[setIdx].count[k], frameStats[setIdx].diff[k], 0, 0, dist[k], cost[k], MAX_CCSAO_OFFSET_THR);
    }
  }

  for (int k = 0; k < MAX_CCSAO_CLASS_NUM; k++)
  {
    CHECK(quantOffsets[k] < -MAX_CCSAO_OFFSET_THR || quantOffsets[k] > MAX_CCSAO_OFFSET_THR, "Exceeded valid range for CCSAO offset");
    offset[setIdx][k] = quantOffsets[k];
  }
}

void EncSampleAdaptiveOffset::getCcSaoDistortion(const ComponentID compID, const int setIdx, CcSaoStatData* blkStats[MAX_CCSAO_SET_NUM]
                                               , short offset[MAX_CCSAO_SET_NUM][MAX_CCSAO_CLASS_NUM]
                                               , int64_t* trainingDistortion[MAX_CCSAO_SET_NUM])
{
  ::memset(trainingDistortion[setIdx], 0, sizeof(int64_t) * m_numCTUsInPic);

  for (int ctbIdx = 0; ctbIdx < m_numCTUsInPic; ctbIdx++)
  {
    for (int k = 0; k < MAX_CCSAO_CLASS_NUM; k++)
    {
      trainingDistortion[setIdx][ctbIdx]
        += estSaoDist(blkStats[setIdx][ctbIdx].count[k], offset[setIdx][k], blkStats[setIdx][ctbIdx].diff[k], 0);
    }
  }
}

void EncSampleAdaptiveOffset::determineCcSaoControlIdc(CodingStructure& cs, const ComponentID compID 
                                                     , const int ctuWidthC, const int ctuHeightC, const int picWidthC, const int picHeightC
                                                     , CcSaoEncParam& ccSaoParam, uint8_t* ccSaoControl
                                                     , int64_t* trainingDistorsion[MAX_CCSAO_SET_NUM]
                                                     , int64_t& curTotalDist, double& curTotalRate)
{
  bool setEnabled[MAX_CCSAO_SET_NUM];
  std::fill_n(setEnabled, MAX_CCSAO_SET_NUM, false);

  SetIdxCount setIdxCount[MAX_CCSAO_SET_NUM];
  for (int i = 0; i < MAX_CCSAO_SET_NUM; i++)
  {
    setIdxCount[i].setIdx = i;
    setIdxCount[i].count  = 0;
  }

  double prevRate = curTotalRate;

  TempCtx ctxInitial(m_CtxCache);
  TempCtx ctxBest(m_CtxCache);
  TempCtx ctxStart(m_CtxCache);
  ctxInitial = SubCtx(Ctx::CcSaoControlIdc, m_CABACEstimator->getCtx());
  ctxBest    = SubCtx(Ctx::CcSaoControlIdc, m_CABACEstimator->getCtx());

  int ctbIdx = 0;
  for (int yCtb = 0; yCtb < picHeightC; yCtb += ctuHeightC)
  {
    for (int xCtb = 0; xCtb < picWidthC; xCtb += ctuWidthC)
    {
      int64_t  bestDist   = MAX_INT;
      double   bestRate   = MAX_DOUBLE;
      double   bestCost   = MAX_DOUBLE;
      uint8_t  bestSetIdc = 0;
      uint8_t  bestSetIdx = 0;

      m_CABACEstimator->getCtx() = ctxBest;
      ctxStart                   = SubCtx(Ctx::CcSaoControlIdc, m_CABACEstimator->getCtx());

      for (int setIdx = 0; setIdx <= MAX_CCSAO_SET_NUM; setIdx++)
      {
        if (setIdx < MAX_CCSAO_SET_NUM && !ccSaoParam.setEnabled[setIdx])
          continue;

        uint8_t setIdc = ccSaoParam.mapIdxToIdc[setIdx];
        m_CABACEstimator->getCtx() = ctxStart;
        m_CABACEstimator->resetBits();
        const Position lumaPos = Position({ xCtb << getComponentScaleX(compID, cs.pcv->chrFormat),
                                            yCtb << getComponentScaleY(compID, cs.pcv->chrFormat) });
        m_CABACEstimator->codeCcSaoControlIdc(setIdc, cs, compID, ctbIdx, ccSaoControl, lumaPos, ccSaoParam.setNum);
        
        int64_t dist = setIdx == MAX_CCSAO_SET_NUM ? 0 : trainingDistorsion[setIdx][ctbIdx];
        double  rate = FRAC_BITS_SCALE * m_CABACEstimator->getEstFracBits();
        double  cost = rate * m_lambda[compID] + dist;

        if (cost < bestCost)
        {
          bestCost   = cost;
          bestRate   = rate;
          bestDist   = dist;
          bestSetIdc = setIdc;
          bestSetIdx = setIdx;
          ctxBest = SubCtx(Ctx::CcSaoControlIdc, m_CABACEstimator->getCtx());
          ccSaoControl[ctbIdx] = setIdx == MAX_CCSAO_SET_NUM ? 0 : setIdx + 1;
        }
      }
      if (bestSetIdc != 0)
      {
        setEnabled [bestSetIdx] = true;
        setIdxCount[bestSetIdx].count++;
      }
      curTotalRate += bestRate;
      curTotalDist += bestDist;
      ctbIdx++;
    }
  }

  std::copy_n(setEnabled, MAX_CCSAO_SET_NUM, ccSaoParam.setEnabled);

  std::stable_sort(setIdxCount, setIdxCount + MAX_CCSAO_SET_NUM, compareSetIdxCount);

  int setIdc = 1;
  ccSaoParam.setNum = 0;
  for (SetIdxCount &s : setIdxCount)
  {
    int setIdx = s.setIdx;
    if (ccSaoParam.setEnabled[setIdx])
    {
      ccSaoParam.mapIdxToIdc[setIdx] = setIdc;
      ccSaoParam.setNum++;
      setIdc++;
    }
  }

  curTotalRate = prevRate;
  m_CABACEstimator->getCtx() = ctxInitial;
  m_CABACEstimator->resetBits();
  ctbIdx = 0;
  for (int yCtb = 0; yCtb < picHeightC; yCtb += ctuHeightC)
  {
    for (int xCtb = 0; xCtb < picWidthC; xCtb += ctuWidthC)
    {
      const int setIdxPlus1 = ccSaoControl[ctbIdx];
      const Position lumaPos = Position({ xCtb << getComponentScaleX(compID, cs.pcv->chrFormat), 
                                          yCtb << getComponentScaleY(compID, cs.pcv->chrFormat) });

      m_CABACEstimator->codeCcSaoControlIdc(setIdxPlus1 == 0 ? 0 : ccSaoParam.mapIdxToIdc[setIdxPlus1 - 1],
                                            cs, compID, ctbIdx, ccSaoControl, lumaPos, ccSaoParam.setNum);
      ctbIdx++;
    }
  }
  curTotalRate += FRAC_BITS_SCALE*m_CABACEstimator->getEstFracBits();

  // restore for next iteration
  m_CABACEstimator->getCtx() = ctxInitial;
}

int EncSampleAdaptiveOffset::lengthUvlc(int uiCode)
{
  int uiLength = 1;
  int uiTemp = ++uiCode;

  CHECK(!uiTemp, "Integer overflow");

  while (1 != uiTemp)
  {
    uiTemp >>= 1;
    uiLength += 2;
  }
  // Take care of cases where uiLength > 32
  return (uiLength >> 1) + ((uiLength + 1) >> 1);
}

int EncSampleAdaptiveOffset::getCcSaoParamRate(const ComponentID compID, const CcSaoEncParam& ccSaoParam)
{
  int bits = 0;

  if (ccSaoParam.setNum > 0 )
  {
    bits += lengthUvlc(ccSaoParam.setNum - 1);

    int signaledSetNum = 0;
    for (int setIdx = 0; setIdx < MAX_CCSAO_SET_NUM; setIdx++)
    {
      if (ccSaoParam.setEnabled[setIdx])
      {
        bits += MAX_CCSAO_CAND_POS_Y_BITS;
        bits += MAX_CCSAO_BAND_NUM_Y_BITS;
        bits += MAX_CCSAO_BAND_NUM_U_BITS;
        bits += MAX_CCSAO_BAND_NUM_V_BITS;

        int classNum = ccSaoParam.bandNum[setIdx][COMPONENT_Y ]
                     * ccSaoParam.bandNum[setIdx][COMPONENT_Cb]
                     * ccSaoParam.bandNum[setIdx][COMPONENT_Cr];
        for (int i = 0; i < classNum; i++)
        {
          bits += lengthUvlc(abs(ccSaoParam.offset[setIdx][i])) + (ccSaoParam.offset[setIdx][i] == 0 ? 0 : 1);
        }
        signaledSetNum++;
      }
    }
    CHECK(signaledSetNum != ccSaoParam.setNum, "Number of sets signaled not the same as indicated");
  }
  return bits;
}

void EncSampleAdaptiveOffset::deriveCcSaoRDO(CodingStructure& cs, const ComponentID compID, int64_t* trainingDistortion[MAX_CCSAO_SET_NUM]
                                           , CcSaoStatData* blkStats[MAX_CCSAO_SET_NUM], CcSaoStatData frameStats[MAX_CCSAO_SET_NUM]
                                           , CcSaoEncParam& bestCcSaoParam, CcSaoEncParam& tempCcSaoParam
                                           , uint8_t* bestCcSaoControl, uint8_t* tempCcSaoControl
                                           , double& bestCost, double& tempCost)
{
  const int scaleX          = getComponentScaleX(compID, cs.pcv->chrFormat);
  const int scaleY          = getComponentScaleY(compID, cs.pcv->chrFormat);
  const int ctuWidthC       = cs.pcv->maxCUWidth  >> scaleX;
  const int ctuHeightC      = cs.pcv->maxCUHeight >> scaleY;
  const int picWidthC       = cs.pcv->lumaWidth   >> scaleX;
  const int picHeightC      = cs.pcv->lumaHeight  >> scaleY;
  const int maxTrainingIter = 15;

  const TempCtx ctxStartCcSaoControlFlag  ( m_CtxCache, SubCtx( Ctx::CcSaoControlIdc, m_CABACEstimator->getCtx() ) );

  int    trainingIter = 0;
  bool   keepTraining = true;
  bool   improved = false;
  double prevCost = MAX_DOUBLE;
  while (keepTraining)
  {
    improved = false;

    for (int setIdx = 0; setIdx < MAX_CCSAO_SET_NUM; setIdx++)
    {
      if (tempCcSaoParam.setEnabled[setIdx])
      {
        getCcSaoFrameStats(compID, setIdx, tempCcSaoControl, blkStats, frameStats);
        deriveCcSaoOffsets(compID, cs.sps->getBitDepth(toChannelType(compID)), setIdx, frameStats, tempCcSaoParam.offset);
        getCcSaoDistortion(compID, setIdx, blkStats, tempCcSaoParam.offset, trainingDistortion);
      }
    }

    m_CABACEstimator->getCtx() = ctxStartCcSaoControlFlag;

    int64_t curTotalDist = 0;
    double  curTotalRate = 0;
    determineCcSaoControlIdc(cs, compID, ctuWidthC, ctuHeightC, picWidthC, picHeightC,
                             tempCcSaoParam, tempCcSaoControl, trainingDistortion,
                             curTotalDist, curTotalRate);

    if (tempCcSaoParam.setNum > 0)
    {
      curTotalRate += getCcSaoParamRate(compID, tempCcSaoParam);
      tempCost = curTotalRate * m_lambda[compID] + curTotalDist;

      if (tempCost < prevCost)
      {
        prevCost = tempCost;
        improved = true;
      }

      if (tempCost < bestCost)
      {
        bestCost = tempCost;
        bestCcSaoParam = tempCcSaoParam;
        memcpy(bestCcSaoControl, tempCcSaoControl, sizeof(uint8_t) * m_numCTUsInPic);
      }
    }

    trainingIter++;
    if (!improved || trainingIter > maxTrainingIter)
    {
      keepTraining = false;
    }
  }
}
#endif

void EncSampleAdaptiveOffset::deriveLoopFilterBoundaryAvailibility(CodingStructure& cs, const Position &pos, bool& isLeftAvail, bool& isAboveAvail, bool& isAboveLeftAvail) const
{
  bool isLoopFiltAcrossSlicePPS = cs.pps->getLoopFilterAcrossSlicesEnabledFlag();
  bool isLoopFiltAcrossTilePPS = cs.pps->getLoopFilterAcrossTilesEnabledFlag();

  const int width = cs.pcv->maxCUWidth;
  const int height = cs.pcv->maxCUHeight;
  const CodingUnit* cuCurr = cs.getCU(pos, CH_L);
  const CodingUnit* cuLeft = cs.getCU(pos.offset(-width, 0), CH_L);
  const CodingUnit* cuAbove = cs.getCU(pos.offset(0, -height), CH_L);
  const CodingUnit* cuAboveLeft = cs.getCU(pos.offset(-width, -height), CH_L);

  if (!isLoopFiltAcrossSlicePPS)
  {
    isLeftAvail      = (cuLeft == NULL)      ? false : CU::isSameTile(*cuCurr, *cuLeft);
    isAboveAvail     = (cuAbove == NULL)     ? false : CU::isSameTile(*cuCurr, *cuAbove);
    isAboveLeftAvail = (cuAboveLeft == NULL) ? false : CU::isSameTile(*cuCurr, *cuAboveLeft);
  }
  else
  {
    isLeftAvail      = (cuLeft != NULL);
    isAboveAvail     = (cuAbove != NULL);
    isAboveLeftAvail = (cuAboveLeft != NULL);
  }

  if (!isLoopFiltAcrossTilePPS)
  {
    isLeftAvail      = (!isLeftAvail)      ? false : CU::isSameTile(*cuCurr, *cuLeft);
    isAboveAvail     = (!isAboveAvail)     ? false : CU::isSameTile(*cuCurr, *cuAbove);
    isAboveLeftAvail = (!isAboveLeftAvail) ? false : CU::isSameTile(*cuCurr, *cuAboveLeft);
  }

  const SubPic& curSubPic = cs.pps->getSubPicFromCU(*cuCurr);
  if (!curSubPic.getloopFilterAcrossEnabledFlag())
  {
    isLeftAvail      = (!isLeftAvail)      ? false : CU::isSameSubPic(*cuCurr, *cuLeft);
    isAboveAvail     = (!isAboveAvail)     ? false : CU::isSameSubPic(*cuCurr, *cuAbove);
    isAboveLeftAvail = (!isAboveLeftAvail) ? false : CU::isSameSubPic(*cuCurr, *cuAboveLeft);
  }
}

//! \}
