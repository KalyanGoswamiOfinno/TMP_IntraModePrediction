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

/** \file     IntraPrediction.h
    \brief    prediction class (header)
*/

#ifndef __INTRAPREDICTION__
#define __INTRAPREDICTION__


// Include files
#include "Unit.h"
#include "Buffer.h"
#include "Picture.h"
#if JVET_W0123_TIMD_FUSION
#include "RdCost.h"
#endif

#include "MatrixIntraPrediction.h"

//! \ingroup CommonLib
//! \{

// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// prediction class
enum PredBuf
{
  PRED_BUF_UNFILTERED = 0,
  PRED_BUF_FILTERED   = 1,
  NUM_PRED_BUF        = 2
};

static const uint32_t MAX_INTRA_FILTER_DEPTHS=8;

class IntraPrediction
{
#if MMLM
public:
  bool m_encPreRDRun;
#endif
protected:
  Pel      m_refBuffer[MAX_NUM_COMPONENT][NUM_PRED_BUF][(MAX_CU_SIZE * 2 + 1 + MAX_REF_LINE_IDX) * 2];
  uint32_t m_refBufferStride[MAX_NUM_COMPONENT];

private:

#if !MERGE_ENC_OPT
  Pel* m_yuvExt2[MAX_NUM_COMPONENT][4];
  int  m_yuvExtSize2;
#endif

  static const uint8_t m_aucIntraFilter[MAX_INTRA_FILTER_DEPTHS];
#if JVET_W0123_TIMD_FUSION
  static const uint8_t m_aucIntraFilterExt[MAX_INTRA_FILTER_DEPTHS];
  RdCost* m_timdSatdCost;
#endif
#if LMS_LINEAR_MODEL
  unsigned m_auShiftLM[32]; // Table for substituting division operation by multiplication
#endif
  struct IntraPredParam //parameters of Intra Prediction
  {
    bool refFilterFlag;
    bool applyPDPC;
#if GRAD_PDPC
    bool useGradPDPC;
#endif
    bool isModeVer;
    int  multiRefIndex;
    int  intraPredAngle;
    int  absInvAngle;
    bool interpolationFlag;
    int  angularScale;

    // clang-format off
    IntraPredParam()
      : refFilterFlag(false)
      , applyPDPC(false)
#if GRAD_PDPC
      , useGradPDPC(false)
#endif
      , isModeVer(false)
      , multiRefIndex(-1)
      , intraPredAngle(std::numeric_limits<int>::max())
      , absInvAngle(std::numeric_limits<int>::max())
      , interpolationFlag(false)
      , angularScale(-1)
    // clang-format on
    {
    }
  };

  IntraPredParam m_ipaParam;

  Pel* m_piTemp;
  Pel* m_pMdlmTemp; // for MDLM mode
  MatrixIntraPrediction m_matrixIntraPred;



protected:
  ChromaFormat  m_currChromaFormat;

  int m_topRefLength;
  int m_leftRefLength;
  ScanElement* m_scanOrder;
  bool         m_bestScanRotationMode;
  std::vector<PelStorage>   m_tempBuffer;

  // prediction
  void xPredIntraPlanar           ( const CPelBuf &pSrc, PelBuf &pDst );
  void xPredIntraDc               ( const CPelBuf &pSrc, PelBuf &pDst, const ChannelType channelType, const bool enableBoundaryFilter = true );
#if JVET_W0123_TIMD_FUSION
  void xPredIntraAng              ( const CPelBuf &pSrc, PelBuf &pDst, const ChannelType channelType, const ClpRng& clpRng, const bool bExtIntraDir);
#else
  void xPredIntraAng              ( const CPelBuf &pSrc, PelBuf &pDst, const ChannelType channelType, const ClpRng& clpRng);
#endif

  void initPredIntraParams        ( const PredictionUnit & pu,  const CompArea compArea, const SPS& sps );

  static bool isIntegerSlope(const int absAng) { return (0 == (absAng & 0x1F)); }
#if JVET_W0123_TIMD_FUSION
  static bool isIntegerSlopeExt(const int absAng) { return (0 == (absAng & 0x3F)); }
#endif

  void xPredIntraBDPCM            ( const CPelBuf &pSrc, PelBuf &pDst, const uint32_t dirMode, const ClpRng& clpRng );
  Pel  xGetPredValDc              ( const CPelBuf &pSrc, const Size &dstSize );
#if JVET_V0130_INTRA_TMP && !JVET_W0069_TMP_BOUNDARY
  bool isRefTemplateAvailable(CodingUnit& cu, CompArea& area);
#endif

  void xFillReferenceSamples      ( const CPelBuf &recoBuf,      Pel* refBufUnfiltered, const CompArea &area, const CodingUnit &cu );
  void xFilterReferenceSamples(const Pel *refBufUnfiltered, Pel *refBufFiltered, const CompArea &area, const SPS &sps,
                               int multiRefIdx
  );

  static int getModifiedWideAngle         ( int width, int height, int predMode );
#if JVET_W0123_TIMD_FUSION
  static int getWideAngleExt      ( int width, int height, int predMode );
#endif
  void setReferenceArrayLengths   ( const CompArea &area );

  void destroy                    ();
#if LMS_LINEAR_MODEL
  void xPadMdlmTemplateSample(Pel*pSrc, Pel*pCur, int cWidth, int cHeight, int existSampNum, int targetSampNum);
#if MMLM
  void xGetLMParameters_LMS(const PredictionUnit &pu, const ComponentID compID, const CompArea& chromaArea, int& a, int& b, int& iShift, int &a2, int &b2, int &iShift2, int &yThres);
#else
  void xGetLMParameters_LMS(const PredictionUnit &pu, const ComponentID compID, const CompArea& chromaArea, int& a, int& b, int& iShift);
#endif
#else
#if MMLM
  void xGetLMParameters(const PredictionUnit &pu, const ComponentID compID, const CompArea& chromaArea, int& a, int& b, int& iShift, int &a2, int &b2, int &iShift2, int &yThres);
#else
  void xGetLMParameters(const PredictionUnit &pu, const ComponentID compID, const CompArea& chromaArea, int& a, int& b, int& iShift);
#endif
#endif
#if LMS_LINEAR_MODEL && MMLM
  struct MMLM_parameter
  {
    int a;
    int b;
    int shift;
  };
  int xCalcLMParametersGeneralized(int x, int y, int xx, int xy, int count, int bitDepth, int &a, int &b, int &iShift);
  int xLMSampleClassifiedTraining(int count, int mean, int meanC, int LumaSamples[], int ChrmSamples[], int bitDepth, MMLM_parameter parameters[]);
#endif

public:
  IntraPrediction();
  virtual ~IntraPrediction();

#if JVET_W0069_TMP_BOUNDARY
  RefTemplateType GetRefTemplateType(CodingUnit& cu, CompArea& area);
#endif

  void init                       (ChromaFormat chromaFormatIDC, const unsigned bitDepthY);
#if ENABLE_DIMD
  static void deriveDimdMode      (const CPelBuf &recoBuf, const CompArea &area, CodingUnit &cu);
#endif
#if JVET_W0123_TIMD_FUSION
  void xIntraPredTimdHorVerPdpc   (Pel* pDsty,const int dstStride, Pel* refSide, const int width, const int height, int xOffset, int yOffset, int scale, const Pel* refMain, const ClpRng& clpRng);
  void xPredTimdIntraPlanar       (const CPelBuf &pSrc, Pel* pDst, int iDstStride, int width, int height, TEMPLATE_TYPE eTempType, int iTemplateWidth , int iTemplateHeight);
  void xPredTimdIntraDc           ( const PredictionUnit &pu, const CPelBuf &pSrc, Pel* pDst, int iDstStride, int iWidth, int iHeight, TEMPLATE_TYPE eTempType, int iTemplateWidth , int iTemplateHeight);
  void xPredTimdIntraAng          ( const CPelBuf &pSrc, const ClpRng& clpRng, Pel* pTrueDst, int iDstStride, int iWidth, int iHeight, TEMPLATE_TYPE eTempType, int iTemplateWidth , int iTemplateHeight, uint32_t dirMode);
  void xIntraPredTimdAngLuma(Pel* pDstBuf, const ptrdiff_t dstStride, Pel* refMain, int width, int height, int deltaPos, int intraPredAngle, const ClpRng& clpRng, int xOffset, int yOffset);
  void xIntraPredTimdPlanarDcPdpc (const CPelBuf &pSrc, Pel* pDst, int iDstStride, int width, int height, TEMPLATE_TYPE eTempType, int iTemplateWidth , int iTemplateHeight);
  void xIntraPredTimdAngPdpc(Pel* pDsty,const int dstStride,Pel* refSide,const int width,const int height, int xOffset, int yOffset, int scale, int invAngle);
  void xFillTimdReferenceSamples  ( const CPelBuf &recoBuf, Pel* refBufUnfiltered, const CompArea &area, const CodingUnit &cu, int iTemplateWidth, int iTemplateHeight );
  Pel  xGetPredTimdValDc          ( const CPelBuf &pSrc, const Size &dstSize, TEMPLATE_TYPE eTempType, int iTempHeight, int iTempWidth );
  void initPredTimdIntraParams    (const PredictionUnit & pu, const CompArea area, int dirMode);
  void predTimdIntraAng           ( const ComponentID compId, const PredictionUnit &pu, uint32_t uiDirMode, Pel* pPred, uint32_t uiStride, uint32_t iWidth, uint32_t iHeight, TEMPLATE_TYPE eTempType, int32_t iTemplateWidth, int32_t iTemplateHeight);
  int deriveTimdMode       ( const CPelBuf &recoBuf, const CompArea &area, CodingUnit &cu );
  void initTimdIntraPatternLuma   (const CodingUnit &cu, const CompArea &area, int iTemplateWidth, int iTemplateHeight, uint32_t uiRefWidth, uint32_t uiRefHeight);
#if GRAD_PDPC
  void xIntraPredTimdAngGradPdpc  (Pel* pDsty, const int dstStride, Pel* refMain, Pel* refSide, const int width, const int height, int xOffset, int yOffset, int scale, int deltaPos, int intraPredAngle, const ClpRng& clpRng);
#endif
#endif
  // Angular Intra
  void predIntraAng               ( const ComponentID compId, PelBuf &piPred, const PredictionUnit &pu);
  Pel *getPredictorPtr(const ComponentID compId)
  {
    return m_refBuffer[compId][m_ipaParam.refFilterFlag ? PRED_BUF_FILTERED : PRED_BUF_UNFILTERED];
  }

  // Cross-component Chroma
  void predIntraChromaLM(const ComponentID compID, PelBuf &piPred, const PredictionUnit &pu, const CompArea& chromaArea, int intraDir);
  void xGetLumaRecPixels(const PredictionUnit &pu, CompArea chromaArea);
  /// set parameters from CU data for accessing intra data
  void initIntraPatternChType     (const CodingUnit &cu, const CompArea &area, const bool forceRefFilterFlag = false); // use forceRefFilterFlag to get both filtered and unfiltered buffers
  void initIntraPatternChTypeISP  (const CodingUnit& cu, const CompArea& area, PelBuf& piReco, const bool forceRefFilterFlag = false); // use forceRefFilterFlag to get both filtered and unfiltered buffers

  // Matrix-based intra prediction
  void initIntraMip               (const PredictionUnit &pu, const CompArea &area);
  void predIntraMip               (const ComponentID compId, PelBuf &piPred, const PredictionUnit &pu);

  template<bool lmcs>
  void geneWeightedPred           ( const ComponentID compId, PelBuf& pred, const PredictionUnit &pu, const PelBuf& interPred, const PelBuf& intraPred, const Pel* pLUT = nullptr );
  void geneIntrainterPred         (const CodingUnit &cu, PelStorage& pred);
  void reorderPLT                 (CodingStructure& cs, Partitioner& partitioner, ComponentID compBegin, uint32_t numComp);
#if !MERGE_ENC_OPT
  void geneWeightedPred           (const ComponentID compId, PelBuf &pred, const PredictionUnit &pu, Pel *srcBuf);
  Pel* getPredictorPtr2           (const ComponentID compID, uint32_t idx) { return m_yuvExt2[compID][idx]; }
  void switchBuffer               (const PredictionUnit &pu, ComponentID compID, PelBuf srcBuff, Pel *dst);
#endif
};
#if ENABLE_DIMD
int  buildHistogram(const Pel *pReco, int iStride, uint32_t uiHeight, uint32_t uiWidth, int* piHistogram, int direction, int bw, int bh);
#endif
//! \}

#endif // __INTRAPREDICTION__
