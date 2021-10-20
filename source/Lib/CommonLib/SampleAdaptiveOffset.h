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

/** \file     SampleAdaptiveOffset.h
    \brief    sample adaptive offset class (header)
*/

#ifndef __SAMPLEADAPTIVEOFFSET__
#define __SAMPLEADAPTIVEOFFSET__

#include "CommonDef.h"
#include "Unit.h"
#include "Reshape.h"
#if JVET_V0094_BILATERAL_FILTER
#include "BilateralFilter.h"
#endif
//! \ingroup CommonLib
//! \{


// ====================================================================================================================
// Constants
// ====================================================================================================================

#define MAX_SAO_TRUNCATED_BITDEPTH     10

// ====================================================================================================================
// Class definition
// ====================================================================================================================

template <typename T> int sgn(T val)
{
  return (T(0) < val) - (val < T(0));
}

class SampleAdaptiveOffset
{
public:
  SampleAdaptiveOffset();
  virtual ~SampleAdaptiveOffset();
  void SAOProcess( CodingStructure& cs, SAOBlkParam* saoBlkParams
                   );
  void create( int picWidth, int picHeight, ChromaFormat format, uint32_t maxCUWidth, uint32_t maxCUHeight, uint32_t maxCUDepth, uint32_t lumaBitShift, uint32_t chromaBitShift );
  void destroy();
  static int getMaxOffsetQVal(const int channelBitDepth) { return (1<<(std::min<int>(channelBitDepth,MAX_SAO_TRUNCATED_BITDEPTH)-5))-1; } //Table 9-32, inclusive
  void setReshaper(Reshape * p) { m_pcReshape = p; }
#if JVET_V0094_BILATERAL_FILTER
  BilateralFilter m_bilateralFilter;
#endif
#if JVET_W0066_CCSAO
  void CCSAOProcess(CodingStructure& cs);
  CcSaoComParam& getCcSaoComParam() { return m_ccSaoComParam; }
  uint8_t* getCcSaoControlIdc(const ComponentID compID) { return m_ccSaoControl[compID]; }
  PelStorage& getCcSaoBuf() { return m_ccSaoBuf; }
  void jointClipSaoBifCcSao(CodingStructure& cs);
#endif
protected:
  void deriveLoopFilterBoundaryAvailibility(CodingStructure& cs, const Position &pos,
    bool& isLeftAvail,
    bool& isRightAvail,
    bool& isAboveAvail,
    bool& isBelowAvail,
    bool& isAboveLeftAvail,
    bool& isAboveRightAvail,
    bool& isBelowLeftAvail,
    bool& isBelowRightAvail
    ) const;

  void offsetBlock(const int channelBitDepth, const ClpRng& clpRng, int typeIdx, int* offset, const Pel* srcBlk, Pel* resBlk, int srcStride, int resStride,  int width, int height
                  , bool isLeftAvail, bool isRightAvail, bool isAboveAvail, bool isBelowAvail, bool isAboveLeftAvail, bool isAboveRightAvail, bool isBelowLeftAvail, bool isBelowRightAvail
                  , bool isCtuCrossedByVirtualBoundaries, int horVirBndryPos[], int verVirBndryPos[], int numHorVirBndry, int numVerVirBndry
    );
#if JVET_V0094_BILATERAL_FILTER
  void offsetBlockBIFonly(const int channelBitDepth, const ClpRng& clpRng, int typeIdx, int* offset, const Pel* srcBlk, Pel* resBlk, int srcStride, int resStride,  int width, int height
                          , bool isLeftAvail, bool isRightAvail, bool isAboveAvail, bool isBelowAvail, bool isAboveLeftAvail, bool isAboveRightAvail, bool isBelowLeftAvail, bool isBelowRightAvail, CodingStructure &cs, const UnitArea &area);
  void offsetBlockNoClip(const int channelBitDepth, const ClpRng& clpRng, int typeIdx, int* offset
                         , const Pel* srcBlk, Pel* resBlk, int srcStride, int resStride,  int width, int height
                         , bool isLeftAvail,  bool isRightAvail, bool isAboveAvail, bool isBelowAvail, bool isAboveLeftAvail, bool isAboveRightAvail, bool isBelowLeftAvail, bool isBelowRightAvail);
#endif
  void invertQuantOffsets(ComponentID compIdx, int typeIdc, int typeAuxInfo, int* dstOffsets, int* srcOffsets);
  void reconstructBlkSAOParam(SAOBlkParam& recParam, SAOBlkParam* mergeList[NUM_SAO_MERGE_TYPES]);
  int  getMergeList(CodingStructure& cs, int ctuRsAddr, SAOBlkParam* blkParams, SAOBlkParam* mergeList[NUM_SAO_MERGE_TYPES]);
  void offsetCTU(const UnitArea& area, const CPelUnitBuf& src, PelUnitBuf& res, SAOBlkParam& saoblkParam, CodingStructure& cs);
#if JVET_V0094_BILATERAL_FILTER
  void offsetCTUnoClip( const UnitArea& area, const CPelUnitBuf& src, PelUnitBuf& res, SAOBlkParam& saoblkParam, CodingStructure& cs);
#if JVET_W0066_CCSAO
  void clipCTU(CodingStructure& cs, PelUnitBuf& dstYuv, const UnitArea& area, const ComponentID compID);
#endif
  void offsetCTUonlySAO(const UnitArea& area, const CPelUnitBuf& src, PelUnitBuf& res, SAOBlkParam& saoblkParam, CodingStructure& cs);
#endif
  void xReconstructBlkSAOParams(CodingStructure& cs, SAOBlkParam* saoBlkParams);
  bool isCrossedByVirtualBoundaries(const int xPos, const int yPos, const int width, const int height, int& numHorVirBndry, int& numVerVirBndry, int horVirBndryPos[], int verVirBndryPos[], const PicHeader* picHeader);
  inline bool isProcessDisabled(int xPos, int yPos, int numVerVirBndry, int numHorVirBndry, int verVirBndryPos[], int horVirBndryPos[])
  {
    bool bDisabledFlag = false;
    for (int i = 0; i < numVerVirBndry; i++)
    {
      if ((xPos == verVirBndryPos[i]) || (xPos == verVirBndryPos[i] - 1))
      {
        bDisabledFlag = true;
        break;
      }
    }
    for (int i = 0; i < numHorVirBndry; i++)
    {
      if ((yPos == horVirBndryPos[i]) || (yPos == horVirBndryPos[i] - 1))
      {
        bDisabledFlag = true;
        break;
      }
    }
    return bDisabledFlag;
  }
  Reshape* m_pcReshape;
#if JVET_W0066_CCSAO
  void applyCcSao(CodingStructure &cs, const PreCalcValues& pcv, const CPelUnitBuf& srcYuv, PelUnitBuf& dstYuv);
  void offsetCTUCcSao(CodingStructure& cs, const UnitArea& area, const CPelUnitBuf& srcYuv, PelUnitBuf& dstYuv, const int ctuRsAddr);
  void offsetCTUCcSaoNoClip(CodingStructure& cs, const UnitArea& area, const CPelUnitBuf& srcYuv, PelUnitBuf& dstYuv, const int ctuRsAddr);
  void offsetBlockCcSao(const ComponentID compID, const int bitDepth, const ClpRng& clpRng
                      , const uint16_t candPosY
                      , const uint16_t bandNumY, const uint16_t bandNumU, const uint16_t bandNumV
                      , const short* offset
                      , const Pel* srcY, const Pel* srcU, const Pel* srcV, Pel* dst
                      , const int srcStrideY, const int srcStrideU, const int srcStrideV, const int dstStride
                      , const int width, const int height
                      , bool isLeftAvail, bool isRightAvail, bool isAboveAvail, bool isBelowAvail, bool isAboveLeftAvail, bool isAboveRightAvail, bool isBelowLeftAvail, bool isBelowRightAvail
                       );
  void offsetBlockCcSaoNoClip(const ComponentID compID, const int bitDepth, const ClpRng& clpRng
                            , const uint16_t candPosY
                            , const uint16_t bandNumY, const uint16_t bandNumU, const uint16_t bandNumV
                            , const short* offset
                            , const Pel* srcY, const Pel* srcU, const Pel* srcV, Pel* dst
                            , const int srcStrideY, const int srcStrideU, const int srcStrideV, const int dstStride
                            , const int width, const int height
                            , bool isLeftAvail, bool isRightAvail, bool isAboveAvail, bool isBelowAvail, bool isAboveLeftAvail, bool isAboveRightAvail, bool isBelowLeftAvail, bool isBelowRightAvail
                             );
#endif
protected:
  uint32_t m_offsetStepLog2[MAX_NUM_COMPONENT]; //offset step
  PelStorage m_tempBuf;
  uint32_t m_numberOfComponents;

  std::vector<int8_t> m_signLineBuf1;
  std::vector<int8_t> m_signLineBuf2;
#if JVET_W0066_CCSAO
  bool                m_created = false;
  PelStorage          m_ccSaoBuf;
  int                 m_picWidth;
  int                 m_picHeight;
  int                 m_maxCUWidth;
  int                 m_maxCUHeight;
  int                 m_numCTUsInWidth;
  int                 m_numCTUsInHeight;
  int                 m_numCTUsInPic;

  CcSaoComParam       m_ccSaoComParam;
  uint8_t*            m_ccSaoControl[MAX_NUM_COMPONENT];
#endif
private:
  bool m_picSAOEnabled[MAX_NUM_COMPONENT];
};

//! \}
#endif

