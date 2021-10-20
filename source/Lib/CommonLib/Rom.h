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

/** \file     Rom.h
    \brief    global variables & functions (header)
*/

#ifndef __ROM__
#define __ROM__

#include "CommonDef.h"
#include "Common.h"

#include <stdio.h>
#include <iostream>


//! \ingroup CommonLib
//! \{

// ====================================================================================================================
// Initialize / destroy functions
// ====================================================================================================================

void         initROM();
void         destroyROM();

// ====================================================================================================================
// Data structure related table & variable
// ====================================================================================================================


#if SIGN_PREDICTION
extern const int8_t * g_resiBorderTemplate[6][6][NUM_TRANS_TYPE*NUM_TRANS_TYPE];
#endif

// flexible conversion from relative to absolute index
struct ScanElement
{
  uint32_t idx;
  uint16_t x;
  uint16_t y;
};

extern       uint32_t   g_log2SbbSize[MAX_CU_DEPTH + 1][MAX_CU_DEPTH + 1][2];
extern ScanElement
  *g_scanOrder[SCAN_NUMBER_OF_GROUP_TYPES][SCAN_NUMBER_OF_TYPES][MAX_CU_SIZE / 2 + 1][MAX_CU_SIZE / 2 + 1];
#if JVET_W0119_LFNST_EXTENSION
extern       ScanElement   g_coefTopLeftDiagScan8x8  [ MAX_CU_DEPTH + 1 ][  64 ];
extern       ScanElement   g_coefTopLeftDiagScan16x16[ MAX_CU_DEPTH + 1 ][ 256 ];
#else
extern       ScanElement   g_coefTopLeftDiagScan8x8[ MAX_CU_SIZE / 2 + 1 ][ 64 ];
#endif

extern const int g_quantScales   [2/*0=4^n blocks, 1=2*4^n blocks*/][SCALING_LIST_REM_NUM];          // Q(QP%6)
extern const int g_invQuantScales[2/*0=4^n blocks, 1=2*4^n blocks*/][SCALING_LIST_REM_NUM];          // IQ(QP%6)
#if TU_256
static const int g_numTransformMatrixSizes = 8;
#else
static const int g_numTransformMatrixSizes = 6;
#endif
#if RExt__HIGH_PRECISION_FORWARD_TRANSFORM
static const int g_transformMatrixShift[TRANSFORM_NUMBER_OF_DIRECTIONS] = { 14, 6 };
#else
static const int g_transformMatrixShift[TRANSFORM_NUMBER_OF_DIRECTIONS] = {  6, 6 };
#endif

#if LMS_LINEAR_MODEL || TRANSFORM_SIMD_OPT
extern int8_t          g_aucLog2[MAX_CU_SIZE + 1];
#endif
#if LMS_LINEAR_MODEL
extern int8_t          g_aucNextLog2[MAX_CU_SIZE + 1];
extern int8_t          g_aucPrevLog2[MAX_CU_SIZE + 1];
#endif
#if MMLM && !LMS_LINEAR_MODEL
extern int8_t          g_aucPrevLog2[MAX_CU_SIZE + 1];
#endif
// ====================================================================================================================
// Scanning order & context mapping table
// ====================================================================================================================

#if TU_256
extern const uint32_t   g_uiGroupIdx[];
#else
extern const uint32_t   g_uiGroupIdx[ MAX_TB_SIZEY ];
#endif
extern const uint32_t   g_uiMinInGroup[ LAST_SIGNIFICANT_GROUPS ];
extern const uint32_t   g_auiGoRiceParsCoeff     [ 32 ];
inline uint32_t g_auiGoRicePosCoeff0(int st, uint32_t ricePar)
{
#if TCQ_8STATES
	return (1 + (st & 1)) << ricePar;
#else
  return (st < 2 ? 1 : 2) << ricePar;
#endif
}

#if TCQ_8STATES
extern const uint64_t g_stateTransTab[3];
#endif
// ====================================================================================================================
// Intra prediction table
// ====================================================================================================================

extern const uint8_t  g_aucIntraModeNumFast_UseMPM_2D[7 - MIN_CU_LOG2 + 1][7 - MIN_CU_LOG2 + 1];
extern const uint8_t  g_aucIntraModeNumFast_UseMPM   [MAX_CU_DEPTH];
extern const uint8_t  g_aucIntraModeNumFast_NotUseMPM[MAX_CU_DEPTH];

extern const uint8_t  g_chroma422IntraAngleMappingTable[NUM_INTRA_MODE];


// ====================================================================================================================
// Mode-Dependent DST Matrices
// ====================================================================================================================


extern const TMatrixCoeff g_trCoreDCT2P2  [TRANSFORM_NUMBER_OF_DIRECTIONS][  2][  2];
extern const TMatrixCoeff g_trCoreDCT2P4  [TRANSFORM_NUMBER_OF_DIRECTIONS][  4][  4];
extern const TMatrixCoeff g_trCoreDCT2P8  [TRANSFORM_NUMBER_OF_DIRECTIONS][  8][  8];
extern const TMatrixCoeff g_trCoreDCT2P16 [TRANSFORM_NUMBER_OF_DIRECTIONS][ 16][ 16];
extern const TMatrixCoeff g_trCoreDCT2P32 [TRANSFORM_NUMBER_OF_DIRECTIONS][ 32][ 32];
extern const TMatrixCoeff g_trCoreDCT2P64 [TRANSFORM_NUMBER_OF_DIRECTIONS][ 64][ 64];
#if TU_256
extern const TMatrixCoeff g_trCoreDCT2P128[TRANSFORM_NUMBER_OF_DIRECTIONS][128][128];
#endif

extern const TMatrixCoeff g_trCoreDCT8P4  [TRANSFORM_NUMBER_OF_DIRECTIONS][  4][  4];
extern const TMatrixCoeff g_trCoreDCT8P8  [TRANSFORM_NUMBER_OF_DIRECTIONS][  8][  8];
extern const TMatrixCoeff g_trCoreDCT8P16 [TRANSFORM_NUMBER_OF_DIRECTIONS][ 16][ 16];
extern const TMatrixCoeff g_trCoreDCT8P32 [TRANSFORM_NUMBER_OF_DIRECTIONS][ 32][ 32];
#if TU_256
extern const TMatrixCoeff g_trCoreDCT8P64 [TRANSFORM_NUMBER_OF_DIRECTIONS][ 64][ 64];
extern const TMatrixCoeff g_trCoreDCT8P128[TRANSFORM_NUMBER_OF_DIRECTIONS][128][128];
#endif

extern const TMatrixCoeff g_trCoreDST7P4  [TRANSFORM_NUMBER_OF_DIRECTIONS][  4][  4];
extern const TMatrixCoeff g_trCoreDST7P8  [TRANSFORM_NUMBER_OF_DIRECTIONS][  8][  8];
extern const TMatrixCoeff g_trCoreDST7P16 [TRANSFORM_NUMBER_OF_DIRECTIONS][ 16][ 16];
extern const TMatrixCoeff g_trCoreDST7P32 [TRANSFORM_NUMBER_OF_DIRECTIONS][ 32][ 32];
#if TU_256
extern const TMatrixCoeff g_trCoreDST7P64 [TRANSFORM_NUMBER_OF_DIRECTIONS][ 64][ 64];
extern const TMatrixCoeff g_trCoreDST7P128[TRANSFORM_NUMBER_OF_DIRECTIONS][128][128];

extern TMatrixCoeff g_trCoreDCT2P256[256][256];
extern TMatrixCoeff g_trCoreDCT8P256[256][256];
extern TMatrixCoeff g_trCoreDST7P256[256][256];
#endif
#if JVET_W0103_INTRA_MTS
extern TMatrixCoeff g_aiTr2[NUM_TRANS_TYPE][2][2];
extern TMatrixCoeff g_aiTr4[NUM_TRANS_TYPE][4][4];
extern TMatrixCoeff g_aiTr8[NUM_TRANS_TYPE][8][8];
extern TMatrixCoeff g_aiTr16[NUM_TRANS_TYPE][16][16];
extern TMatrixCoeff g_aiTr32[NUM_TRANS_TYPE][32][32];
extern TMatrixCoeff g_aiTr64[NUM_TRANS_TYPE][64][64];
extern TMatrixCoeff g_aiTr128[NUM_TRANS_TYPE][128][128];
extern TMatrixCoeff g_aiTr256[NUM_TRANS_TYPE][256][256];

extern const uint8_t g_aucIpmToTrSet[16][36];
extern const uint8_t g_aucTrSet[80][4];
extern const int8_t  g_aiIdLut[3][3];
extern const uint8_t g_aucTrIdxToTr[25][2];
#endif
#if JVET_W0119_LFNST_EXTENSION
extern const     int8_t   g_lfnst16x16[ 35 ][ 3 ][ L16H ][ L16W ];
extern const     int8_t   g_lfnst8x8  [ 35 ][ 3 ][  L8H ][  L8W ];
extern const     int8_t   g_lfnst4x4  [ 35 ][ 3 ][   16 ][   16 ];
#else
#if EXTENDED_LFNST
extern const     int8_t   g_lfnst8x8[ 35 ][ 3 ][ 64 ][ 64 ];
extern const     int8_t   g_lfnst4x4[ 35 ][ 3 ][ 16 ][ 16 ];
#else
extern const     int8_t   g_lfnst8x8[ 4 ][ 2 ][ 16 ][ 48 ];
extern const     int8_t   g_lfnst4x4[ 4 ][ 2 ][ 16 ][ 16 ];
#endif
#endif

#if EXTENDED_LFNST || JVET_W0119_LFNST_EXTENSION
extern const     uint8_t  g_lfnstLut[ 95 ];
#else
extern const     uint8_t  g_lfnstLut[ NUM_INTRA_MODE + NUM_EXT_LUMA_MODE - 1 ];
#endif

// ====================================================================================================================
// Misc.
// ====================================================================================================================
extern SizeIndexInfo* gp_sizeIdxInfo;

extern const int       g_ictModes[2][4];

inline bool is34( const SizeType& size )
{
  return ( size & ( ( int64_t ) 1 << ( floorLog2(size) - 1 ) ) );
}

inline bool is58( const SizeType& size )
{
  return ( size & ( ( int64_t ) 1 << ( floorLog2(size) - 2 ) ) );
}

inline bool isNonLog2BlockSize( const Size& size )
{
  return ( ( 1 << floorLog2(size.width) ) != size.width ) || ( ( 1 << floorLog2(size.height) ) != size.height );
}

inline bool isNonLog2Size( const SizeType& size )
{
  return ( ( 1 << floorLog2(size) ) != size );
}

extern UnitScale     g_miScaling; // scaling object for motion scaling

/*! Sophisticated Trace-logging */
#if ENABLE_TRACING
#include "dtrace.h"
extern CDTrace* g_trace_ctx;
#endif

const char* nalUnitTypeToString(NalUnitType type);

extern const char *MatrixType   [SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM];
extern const char *MatrixType_DC[SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM];

extern const int g_quantTSDefault4x4   [4*4];
extern const int g_quantIntraDefault8x8[8*8];
extern const int g_quantInterDefault8x8[8*8];

extern const uint32_t g_scalingListSize [SCALING_LIST_SIZE_NUM];
extern const uint32_t g_scalingListSizeX[SCALING_LIST_SIZE_NUM];
extern const uint32_t g_scalingListId[SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM];

extern MsgLevel g_verbosity;
#if MMLM
extern int g_aiLMDivTableLow[];
extern int g_aiLMDivTableHigh[];
#endif

extern const int8_t g_BcwLog2WeightBase;
extern const int8_t g_BcwWeightBase;
extern const int8_t g_BcwWeights[BCW_NUM];
extern const int8_t g_BcwSearchOrder[BCW_NUM];
extern       int8_t g_BcwCodingOrder[BCW_NUM];
extern       int8_t g_BcwParsingOrder[BCW_NUM];

class CodingStructure;
int8_t getBcwWeight(uint8_t bcwIdx, uint8_t uhRefFrmList);
void resetBcwCodingOrder(bool bRunDecoding, const CodingStructure &cs);
uint32_t deriveWeightIdxBits(uint8_t bcwIdx);

constexpr uint8_t g_tbMax[257] = { 0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7,
7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8 };

#if SIGN_PREDICTION
extern const int8_t    g_initRomSignPred[];
extern const int32_t   g_aucIdxTrCombo[5][5][5][NUM_INTRA_MODE];
extern const int32_t   g_aucNumTrCombo[5][5];
#endif

//! \}


extern bool g_mctsDecCheckEnabled;

class  Mv;
#if CTU_256
extern Mv   g_reusedUniMVs        [MAX_CU_SIZE >> MIN_CU_LOG2][MAX_CU_SIZE >> MIN_CU_LOG2][MAX_CU_DEPTH + 1][MAX_CU_DEPTH + 1][2][33];
extern bool g_isReusedUniMVsFilled[MAX_CU_SIZE >> MIN_CU_LOG2][MAX_CU_SIZE >> MIN_CU_LOG2][MAX_CU_DEPTH + 1][MAX_CU_DEPTH + 1];
#else
extern Mv   g_reusedUniMVs[32][32][8][8][2][33];
extern bool g_isReusedUniMVsFilled[32][32][8][8];
#endif
#if INTER_LIC
extern Mv   g_reusedUniMVsLIC        [MAX_CU_SIZE >> MIN_CU_LOG2][MAX_CU_SIZE >> MIN_CU_LOG2][MAX_CU_DEPTH + 1][MAX_CU_DEPTH + 1][2][33];
extern bool g_isReusedUniMVsFilledLIC[MAX_CU_SIZE >> MIN_CU_LOG2][MAX_CU_SIZE >> MIN_CU_LOG2][MAX_CU_DEPTH + 1][MAX_CU_DEPTH + 1];
#endif

extern uint16_t g_paletteQuant[57];
extern uint8_t g_paletteRunTopLut[5];
extern uint8_t g_paletteRunLeftLut[5];

#if CTU_256
const int g_IBCBufferSize = 256 * 128 * 4;
#else
const int g_IBCBufferSize = 256 * 128;
#endif

void initGeoTemplate();
extern int16_t** g_GeoParams;
extern int16_t*  g_globalGeoWeights   [GEO_NUM_PRESTORED_MASK];
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
extern Pel*      g_globalGeoEncSADmask[GEO_NUM_PRESTORED_MASK];
#else
extern int16_t*  g_globalGeoEncSADmask[GEO_NUM_PRESTORED_MASK];
#endif
extern int16_t   g_weightOffset       [GEO_NUM_PARTITION_MODE][GEO_NUM_CU_SIZE][GEO_NUM_CU_SIZE][2];
extern int8_t    g_angle2mask         [GEO_NUM_ANGLES];
extern int8_t    g_Dis[GEO_NUM_ANGLES];
extern int8_t    g_angle2mirror[GEO_NUM_ANGLES];
#if JVET_W0097_GPM_MMVD_TM && TM_MRG
extern uint8_t g_geoTmShape[2][GEO_NUM_ANGLES];
#endif
#if MULTI_HYP_PRED
extern const int g_addHypWeight[MULTI_HYP_PRED_NUM_WEIGHTS];
#endif
#if JVET_W0066_CCSAO
extern const int8_t g_ccSaoCandPosX[MAX_NUM_LUMA_COMP][MAX_CCSAO_CAND_POS_Y];
extern const int8_t g_ccSaoCandPosY[MAX_NUM_LUMA_COMP][MAX_CCSAO_CAND_POS_Y];
#endif
#endif  //__TCOMROM__

