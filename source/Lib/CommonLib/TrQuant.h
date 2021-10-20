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

/** \file     TrQuant.h
    \brief    transform and quantization class (header)
*/

#ifndef __TRQUANT__
#define __TRQUANT__

#include <array>
#include "CommonDef.h"
#include "Unit.h"
#include "ChromaFormat.h"
#include "Contexts.h"
#include "ContextModelling.h"

#include "UnitPartitioner.h"
#include "Quant.h"

#include "DepQuant.h"
//! \ingroup CommonLib
//! \{

typedef void FwdTrans(const TCoeff*, TCoeff*, int, int, int, int);
typedef void InvTrans(const TCoeff*, TCoeff*, int, int, int, int, const TCoeff, const TCoeff);



#if JVET_V0130_INTRA_TMP
extern unsigned int g_uiDepth2Width[5];
extern unsigned int g_uiDepth2MaxCandiNum[5];

class TempLibFast
{
public:
	int m_pX;    //offset X
	int m_pY;    //offset Y
	int m_pXInteger;    //offset X for integer pixel search
	int m_pYInteger;    //offset Y for integer pixel search
	int m_pDiffInteger;
	int getXInteger() { return m_pXInteger; }
	int getYInteger() { return m_pYInteger; }
	int getDiffInteger() { return m_pDiffInteger; }
	short m_pIdInteger; //frame id
	short getIdInteger() { return m_pIdInteger; }
	int m_pDiff; //mse
	short m_pId; //frame id
	

	TempLibFast();
	~TempLibFast();
	//void init();
	int getX() { return m_pX; }
	int getY() { return m_pY; }
	int getDiff() { return m_pDiff; }
	short getId() { return m_pId; }
	/*void initDiff(unsigned int uiPatchSize, int bitDepth);
	void initDiff(unsigned int uiPatchSize, int bitDepth, int iCandiNumber);*/
	void initTemplateDiff(unsigned int uiPatchWidth, unsigned int uiPatchHeight, unsigned int uiBlkWidth, unsigned int uiBlkHeight, int bitDepth);
	int m_diffMax;
	int getDiffMax() { return m_diffMax; }
};


typedef short TrainDataType;
#endif


// ====================================================================================================================
// Class definition
// ====================================================================================================================


/// transform and quantization class
class TrQuant
{
public:
  TrQuant();
  ~TrQuant();

  // initialize class
  void init      (
                    const Quant* otherQuant,
                    const uint32_t uiMaxTrSize,
                    const bool bUseRDOQ,
                    const bool bUseRDOQTS,
#if T0196_SELECTIVE_RDOQ
                    const bool useSelectiveRDOQ,
#endif
                    const bool bEnc
  );
  void getTrTypes(const TransformUnit tu, const ComponentID compID, int &trTypeHor, int &trTypeVer);

#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  void fwdLfnstNxN( TCoeff* src, TCoeff* dst, const uint32_t mode, const uint32_t index, const uint32_t size, int zeroOutSize );
  void invLfnstNxN( TCoeff* src, TCoeff* dst, const uint32_t mode, const uint32_t index, const uint32_t size, int zeroOutSize, const int maxLog2TrDynamicRange );
#else
  void fwdLfnstNxN( int* src, int* dst, const uint32_t mode, const uint32_t index, const uint32_t size, int zeroOutSize );
  void invLfnstNxN( int* src, int* dst, const uint32_t mode, const uint32_t index, const uint32_t size, int zeroOutSize );
#endif
#if JVET_V0130_INTRA_TMP
#if JVET_W0069_TMP_BOUNDARY
  int (*m_calcTemplateDiff)(Pel* ref, unsigned int uiStride, Pel** tarPatch, unsigned int uiPatchWidth, unsigned int uiPatchHeight, int iMax, RefTemplateType TempType);
  static int calcTemplateDiff(Pel* ref, unsigned int uiStride, Pel** tarPatch, unsigned int uiPatchWidth, unsigned int uiPatchHeight, int iMax, RefTemplateType TempType);
#else
  int ( *m_calcTemplateDiff )(Pel* ref, unsigned int uiStride, Pel** tarPatch, unsigned int uiPatchWidth, unsigned int uiPatchHeight, int iMax);
  static int calcTemplateDiff(Pel* ref, unsigned int uiStride, Pel** tarPatch, unsigned int uiPatchWidth, unsigned int uiPatchHeight, int iMax);
#endif
  Pel** getTargetPatch(unsigned int uiDepth)       { return m_pppTarPatch[uiDepth]; }
  Pel* getRefPicUsed()                             { return m_refPicUsed; }
  void setRefPicUsed(Pel* ref)                     { m_refPicUsed = ref; }
  unsigned int getStride()                         { return m_uiPicStride; }
  void         setStride(unsigned int uiPicStride) { m_uiPicStride = uiPicStride; }

#if JVET_W0069_TMP_BOUNDARY
  void searchCandidateFromOnePicIntra(CodingUnit* pcCU, Pel** tarPatch, unsigned int uiPatchWidth, unsigned int uiPatchHeight, unsigned int setId, RefTemplateType TempType);
  void candidateSearchIntra(CodingUnit* pcCU, unsigned int uiBlkWidth, unsigned int uiBlkHeight, RefTemplateType TempType);
#else
  void searchCandidateFromOnePicIntra(CodingUnit* pcCU, Pel** tarPatch, unsigned int uiPatchWidth, unsigned int uiPatchHeight, unsigned int setId);
  void candidateSearchIntra(CodingUnit* pcCU, unsigned int uiBlkWidth, unsigned int uiBlkHeight);
#endif
  bool generateTMPrediction(Pel* piPred, unsigned int uiStride, unsigned int uiBlkWidth, unsigned int uiBlkHeight, int& foundCandiNum);
#if JVET_W0069_TMP_BOUNDARY
  bool generateTM_DC_Prediction(Pel* piPred, unsigned int uiStride, unsigned int uiBlkWidth, unsigned int uiBlkHeight, int DC_Val);
#endif
#if JVET_W0069_TMP_BOUNDARY
  void getTargetTemplate(CodingUnit* pcCU, unsigned int uiBlkWidth, unsigned int uiBlkHeight, RefTemplateType TempType);
  int  getVaildCandiNumTNP() { return m_uiVaildCandiNum; } // Modified by Kalyan
  TempLibFast getTempLibFast() { return m_tempLibFast; }            // Modified by Kalyan
#else
  void getTargetTemplate(CodingUnit* pcCU, unsigned int uiBlkWidth, unsigned int uiBlkHeight);
#endif
#endif

  uint32_t getLFNSTIntraMode( int wideAngPredMode );
  bool     getTransposeFlag ( uint32_t intraMode  );

protected:

  void xFwdLfnst( const TransformUnit &tu, const ComponentID compID, const bool loadTr = false );
  void xInvLfnst( const TransformUnit &tu, const ComponentID compID );

public:

  void invTransformNxN  (TransformUnit &tu, const ComponentID &compID, PelBuf &pResi, const QpParam &cQPs);
  void transformNxN     ( TransformUnit& tu, const ComponentID& compID, const QpParam& cQP, std::vector<TrMode>* trModes, const int maxCand );
  void transformNxN     ( TransformUnit& tu, const ComponentID& compID, const QpParam& cQP, TCoeff& uiAbsSum, const Ctx& ctx, const bool loadTr = false );
#if JVET_W0103_INTRA_MTS
  uint64_t transformNxN(TransformUnit& tu);
#endif

  void transformSkipQuantOneSample(TransformUnit &tu, const ComponentID &compID, const TCoeff &resiDiff, TCoeff &coeff,    const uint32_t &uiPos, const QpParam &cQP, const bool bUseHalfRoundingPoint);
  void invTrSkipDeQuantOneSample  (TransformUnit &tu, const ComponentID &compID, const TCoeff &pcCoeff,  Pel &reconSample, const uint32_t &uiPos, const QpParam &cQP);

  void                        invTransformICT     ( const TransformUnit &tu, PelBuf &resCb, PelBuf &resCr );
  std::pair<int64_t,int64_t>  fwdTransformICT     ( const TransformUnit &tu, const PelBuf &resCb, const PelBuf &resCr, PelBuf& resC1, PelBuf& resC2, int jointCbCr = -1 );
  std::vector<int>            selectICTCandidates ( const TransformUnit &tu, CompStorage* resCb, CompStorage* resCr );

#if RDOQ_CHROMA_LAMBDA
  void   setLambdas  ( const double lambdas[MAX_NUM_COMPONENT] )   { m_quant->setLambdas( lambdas ); }
  void   selectLambda( const ComponentID compIdx )                 { m_quant->selectLambda( compIdx ); }
  void   getLambdas  ( double (&lambdas)[MAX_NUM_COMPONENT]) const { m_quant->getLambdas( lambdas ); }
#endif
  void   setLambda   ( const double dLambda )                      { m_quant->setLambda( dLambda ); }
  double getLambda   () const                                      { return m_quant->getLambda(); }

  DepQuant* getQuant() { return m_quant; }
  void   lambdaAdjustColorTrans(bool forward) { m_quant->lambdaAdjustColorTrans(forward); }
  void   resetStore() { m_quant->resetStore(); }
#if SIGN_PREDICTION
  void predCoeffSigns( TransformUnit &tu, const ComponentID compID, const bool reshapeChroma );
#endif

#if ENABLE_SPLIT_PARALLELISM
  void    copyState( const TrQuant& other );
#endif

#if SIGN_PREDICTION
  enum SIGN_PRED_TYPE
  {
      SIGN_PRED_BYPASS   = 0,
      SIGN_PRED_POSITIVE = 1,
      SIGN_PRED_NEGATIVE = 2,
      SIGN_PRED_HIDDEN   = 3,
  };
  uint32_t m_aiSignPredCost[1 << SIGN_PRED_MAX_NUM];
#endif

protected:
  TCoeff   m_tempCoeff[MAX_TB_SIZEY * MAX_TB_SIZEY];
#if JVET_V0130_INTRA_TMP
  int          m_uiPartLibSize;
  TempLibFast  m_tempLibFast;
  Pel*         m_refPicUsed;
  Picture*     m_refPicBuf;
  unsigned int m_uiPicStride;
  unsigned int m_uiVaildCandiNum;
  Pel***       m_pppTarPatch;
#endif
#if SIGN_PREDICTION
  Pel      m_tempSignPredResid[SIGN_PRED_MAX_BS * SIGN_PRED_MAX_BS * 2]{0};
  Pel      m_signPredTemplate[SIGN_PRED_FREQ_RANGE*SIGN_PRED_FREQ_RANGE*SIGN_PRED_MAX_BS*2];
#endif

private:
  DepQuant *m_quant;          //!< Quantizer
  TCoeff    m_mtsCoeffs[NUM_TRAFO_MODES_MTS][MAX_TB_SIZEY * MAX_TB_SIZEY];
#if JVET_W0119_LFNST_EXTENSION
  TCoeff   m_tempInMatrix [ L16W_ZO ];
  TCoeff   m_tempOutMatrix[ L16W_ZO ];
#else
#if EXTENDED_LFNST
  TCoeff   m_tempInMatrix [ 64 ];
  TCoeff   m_tempOutMatrix[ 64 ];
#else
  TCoeff   m_tempInMatrix [ 48 ];
  TCoeff   m_tempOutMatrix[ 48 ];
#endif
#endif
  static const int maxAbsIctMode = 3;
  void                      (*m_invICTMem[1+2*maxAbsIctMode])(PelBuf&,PelBuf&);
  std::pair<int64_t,int64_t>(*m_fwdICTMem[1+2*maxAbsIctMode])(const PelBuf&,const PelBuf&,PelBuf&,PelBuf&);
  void                      (**m_invICT)(PelBuf&,PelBuf&);
  std::pair<int64_t,int64_t>(**m_fwdICT)(const PelBuf&,const PelBuf&,PelBuf&,PelBuf&);

  std::array<std::array<FwdTrans*, g_numTransformMatrixSizes>, NUM_TRANS_TYPE> fastFwdTrans;
  std::array<std::array<InvTrans*, g_numTransformMatrixSizes>, NUM_TRANS_TYPE> fastInvTrans;

#if TRANSFORM_SIMD_OPT
  static std::array<std::array<const TMatrixCoeff*, g_numTransformMatrixSizes>, NUM_TRANS_TYPE> m_forwardTransformKernels;
  static std::array<std::array<const TMatrixCoeff*, g_numTransformMatrixSizes>, NUM_TRANS_TYPE> m_inverseTransformKernels;
#endif

  // forward Transform
  void xT               (const TransformUnit &tu, const ComponentID &compID, const CPelBuf &resi, CoeffBuf &dstCoeff, const int width, const int height);

  // skipping Transform
  void xTransformSkip   (const TransformUnit &tu, const ComponentID &compID, const CPelBuf &resi, TCoeff* psCoeff);

  // quantization
  void xQuant           (TransformUnit &tu, const ComponentID &compID, const CCoeffBuf &pSrc, TCoeff &uiAbsSum, const QpParam &cQP, const Ctx& ctx);

  // dequantization
  void xDeQuant( const TransformUnit &tu,
                       CoeffBuf      &dstCoeff,
                 const ComponentID   &compID,
                 const QpParam       &cQP      );

  // inverse transform
  void xIT     ( const TransformUnit &tu, const ComponentID &compID, const CCoeffBuf &pCoeff, PelBuf &pResidual );

  // inverse skipping transform
  void xITransformSkip(
                 const CCoeffBuf     &plCoef,
                       PelBuf        &pResidual,
                 const TransformUnit &tu,
                 const ComponentID   &component);

  void xGetCoeffEnergy(
                       TransformUnit  &tu,
                 const ComponentID    &compID,
                 const CoeffBuf       &coeffs,
                       double*        diagRatio,
                       double*        horVerRatio );

#if ENABLE_SIMD_SIGN_PREDICTION
  uint32_t( *m_computeSAD ) ( const Pel* ref, const Pel* cur, const int size );
  static inline uint32_t xComputeSAD( const Pel* ref, const Pel* cur, const int size );
#endif

#if TRANSFORM_SIMD_OPT
  template <TransType transType, int transSize>
  static void fastForwardTransform_SIMD( const TCoeff *block, TCoeff *coeff, int shift, int line, int iSkipLine, int iSkipLine2 );

  template <TransType transType, int transSize>
  static void fastInverseTransform_SIMD( const TCoeff *coeff, TCoeff *block, int shift, int line, int iSkipLine, int iSkipLine2, const TCoeff outputMinimum, const TCoeff outputMaximum );
#endif

#if ENABLE_SIMD_SIGN_PREDICTION || TRANSFORM_SIMD_OPT || ENABLE_SIMD_TMP
#ifdef TARGET_SIMD_X86
  void    initTrQuantX86();
  template <X86_VEXT vext>
  void    _initTrQuantX86();
#endif
#endif

};// END CLASS DEFINITION TrQuant

//! \}

#endif // __TRQUANT__
