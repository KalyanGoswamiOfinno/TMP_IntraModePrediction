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

/** \file     AdaptiveLoopFilter.cpp
    \brief    adaptive loop filter class
*/

#include "AdaptiveLoopFilter.h"
#if ALF_IMPROVEMENT
#include "AlfFixedFilters.h"
#endif
#include "CodingStructure.h"
#include "Picture.h"
#include <array>
#include <cmath>

constexpr int AdaptiveLoopFilter::AlfNumClippingValues[];

AdaptiveLoopFilter::AdaptiveLoopFilter()
  : m_classifier( nullptr )
{
  for (size_t i = 0; i < NUM_DIRECTIONS; i++)
  {
    m_laplacian[i] = m_laplacianPtr[i];
    for (size_t j = 0; j < sizeof(m_laplacianPtr[i]) / sizeof(m_laplacianPtr[i][0]); j++)
    {
      m_laplacianPtr[i][j] = m_laplacianData[i][j];
    }
  }

#if ALF_IMPROVEMENT
  int ind = 0;
  memset(m_mappingDir, 0, sizeof(m_mappingDir));
  for (int i = 0; i < NUM_DIR_FIX; i++)
  {
    for (int j = 0; j <= i; j++)
    {
      m_mappingDir[i][j] = ind;
      ind++;
    }
  }
#endif

  for( int compIdx = 0; compIdx < MAX_NUM_COMPONENT; compIdx++ )
  {
    m_ctuEnableFlag[compIdx] = nullptr;
    m_ctuAlternative[compIdx] = nullptr;
  }
  m_filterCcAlf = filterBlkCcAlf<CC_ALF>;
  m_filter5x5Blk = filterBlk<ALF_FILTER_5>;
  m_filter7x7Blk = filterBlk<ALF_FILTER_7>;
#if ALF_IMPROVEMENT
  m_filter9x9Blk = filterBlk<ALF_FILTER_9>;
  m_filter9x9BlkExt = filterBlk<ALF_FILTER_9_EXT>;
  m_filterBlkExt = filterBlk<ALF_FILTER_EXT>;
  m_filter13x13Blk = fixedFiltering;
  m_deriveClassificationLaplacian = deriveClassificationLaplacian;
  m_deriveClassificationLaplacianBig = deriveClassificationLaplacianBig;
  m_calcClass0 = calcClass;
  m_calcClass1 = calcClass;
  m_fixFilterResult = nullptr;
#else
  m_deriveClassificationBlk = deriveClassificationBlk;
#endif

#if ENABLE_SIMD_OPT_ALF
#ifdef TARGET_SIMD_X86
  initAdaptiveLoopFilterX86();
#endif
#endif
}

bool AdaptiveLoopFilter::isCrossedByVirtualBoundaries( const CodingStructure& cs, const int xPos, const int yPos, const int width, const int height, bool& clipTop, bool& clipBottom, bool& clipLeft, bool& clipRight, int& numHorVirBndry, int& numVerVirBndry, int horVirBndryPos[], int verVirBndryPos[], int& rasterSliceAlfPad )
{
  clipTop = false; clipBottom = false; clipLeft = false; clipRight = false;
  numHorVirBndry = 0; numVerVirBndry = 0;
  const PPS*   pps = cs.pps;
  const PicHeader* picHeader = cs.picHeader;

  if( picHeader->getVirtualBoundariesPresentFlag() )
  {
    for( int i = 0; i < picHeader->getNumHorVirtualBoundaries(); i++ )
    {
      if( picHeader->getVirtualBoundariesPosY(i) == yPos )
      {
        clipTop = true;
      }
      else if( picHeader->getVirtualBoundariesPosY(i) == yPos + height )
      {
        clipBottom = true;
      }
      else if( yPos < picHeader->getVirtualBoundariesPosY(i) && picHeader->getVirtualBoundariesPosY(i) < yPos + height )
      {
        horVirBndryPos[numHorVirBndry++] = picHeader->getVirtualBoundariesPosY(i);
      }
    }
    for( int i = 0; i < picHeader->getNumVerVirtualBoundaries(); i++ )
    {
      if( picHeader->getVirtualBoundariesPosX(i) == xPos )
      {
        clipLeft = true;
      }
      else if( picHeader->getVirtualBoundariesPosX(i) == xPos + width )
      {
        clipRight = true;
      }
      else if( xPos < picHeader->getVirtualBoundariesPosX(i) && picHeader->getVirtualBoundariesPosX(i) < xPos + width )
      {
        verVirBndryPos[numVerVirBndry++] = picHeader->getVirtualBoundariesPosX(i);
      }
    }
  }

  const Slice& slice = *(cs.slice);
  int   ctuSize = slice.getSPS()->getCTUSize();
  const Position currCtuPos(xPos, yPos);
  const CodingUnit *currCtu = cs.getCU(currCtuPos, CHANNEL_TYPE_LUMA);
  const SubPic& curSubPic = slice.getPPS()->getSubPicFromPos(currCtuPos);
  bool loopFilterAcrossSubPicEnabledFlag = curSubPic.getloopFilterAcrossEnabledFlag();
  //top
  if (yPos >= ctuSize && clipTop == false)
  {
    const Position prevCtuPos(xPos, yPos - ctuSize);
    const CodingUnit *prevCtu = cs.getCU(prevCtuPos, CHANNEL_TYPE_LUMA);
    if ((!pps->getLoopFilterAcrossSlicesEnabledFlag() && !CU::isSameSlice(*currCtu, *prevCtu)) ||
        (!pps->getLoopFilterAcrossTilesEnabledFlag()  && !CU::isSameTile(*currCtu,  *prevCtu))
      || (!loopFilterAcrossSubPicEnabledFlag && !CU::isSameSubPic(*currCtu, *prevCtu))
      )
    {
      clipTop = true;
    }
  }

  //bottom
  if (yPos + ctuSize < cs.pcv->lumaHeight && clipBottom == false)
  {
    const Position nextCtuPos(xPos, yPos + ctuSize);
    const CodingUnit *nextCtu = cs.getCU(nextCtuPos, CHANNEL_TYPE_LUMA);
    if ((!pps->getLoopFilterAcrossSlicesEnabledFlag() && !CU::isSameSlice(*currCtu, *nextCtu)) ||
        (!pps->getLoopFilterAcrossTilesEnabledFlag()  && !CU::isSameTile(*currCtu,  *nextCtu))
      || (!loopFilterAcrossSubPicEnabledFlag && !CU::isSameSubPic(*currCtu, *nextCtu))
      )
    {
      clipBottom = true;
    }
  }

  //left
  if (xPos >= ctuSize && clipLeft == false)
  {
    const Position prevCtuPos(xPos - ctuSize, yPos);
    const CodingUnit *prevCtu = cs.getCU(prevCtuPos, CHANNEL_TYPE_LUMA);
    if ((!pps->getLoopFilterAcrossSlicesEnabledFlag() && !CU::isSameSlice(*currCtu, *prevCtu)) ||
        (!pps->getLoopFilterAcrossTilesEnabledFlag()  && !CU::isSameTile(*currCtu,  *prevCtu))
      || (!loopFilterAcrossSubPicEnabledFlag && !CU::isSameSubPic(*currCtu, *prevCtu))
      )
    {
      clipLeft = true;
    }
  }

  //right
  if (xPos + ctuSize < cs.pcv->lumaWidth && clipRight == false)
  {
    const Position nextCtuPos(xPos + ctuSize, yPos);
    const CodingUnit *nextCtu = cs.getCU(nextCtuPos, CHANNEL_TYPE_LUMA);
    if ((!pps->getLoopFilterAcrossSlicesEnabledFlag() && !CU::isSameSlice(*currCtu, *nextCtu)) ||
        (!pps->getLoopFilterAcrossTilesEnabledFlag()  && !CU::isSameTile(*currCtu,  *nextCtu))
      || (!loopFilterAcrossSubPicEnabledFlag && !CU::isSameSubPic(*currCtu, *nextCtu))
      )
    {
      clipRight = true;
    }
  }

  rasterSliceAlfPad = 0;
  if ( !clipTop && !clipLeft )
  {
    //top-left CTU
    if ( xPos >= ctuSize && yPos >= ctuSize )
    {
      const Position prevCtuPos( xPos - ctuSize, yPos - ctuSize );
      const CodingUnit *prevCtu = cs.getCU( prevCtuPos, CHANNEL_TYPE_LUMA );
      if ( !pps->getLoopFilterAcrossSlicesEnabledFlag() && !CU::isSameSlice( *currCtu, *prevCtu ) )
      {
        rasterSliceAlfPad = 1;
      }
    }
  }

  if ( !clipBottom && !clipRight )
  {
    //bottom-right CTU
    if ( xPos + ctuSize < cs.pcv->lumaWidth && yPos + ctuSize < cs.pcv->lumaHeight )
    {
      const Position nextCtuPos( xPos + ctuSize, yPos + ctuSize );
      const CodingUnit *nextCtu = cs.getCU( nextCtuPos, CHANNEL_TYPE_LUMA );
      if ( !pps->getLoopFilterAcrossSlicesEnabledFlag() && !CU::isSameSlice( *currCtu, *nextCtu ) )
      {
        rasterSliceAlfPad += 2;
      }
    }
  }

  return numHorVirBndry > 0 || numVerVirBndry > 0 || clipTop || clipBottom || clipLeft || clipRight || rasterSliceAlfPad;
}
#if !ALF_IMPROVEMENT
const int AdaptiveLoopFilter::m_fixedFilterSetCoeff[ALF_FIXED_FILTER_NUM][MAX_NUM_ALF_LUMA_COEFF] =
{
  { 0,   0,   2,  -3,   1,  -4,   1,   7,  -1,   1,  -1,   5, 0 },
  { 0,   0,   0,   0,   0,  -1,   0,   1,   0,   0,  -1,   2, 0 },
  { 0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0, 0 },
  { 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -1,   1, 0 },
  { 2,   2,  -7,  -3,   0,  -5,  13,  22,  12,  -3,  -3,  17,  0 },
  { -1,   0,   6,  -8,   1,  -5,   1,  23,   0,   2,  -5,  10,  0 },
  { 0,   0,  -1,  -1,   0,  -1,   2,   1,   0,   0,  -1,   4, 0 },
  { 0,   0,   3, -11,   1,   0,  -1,  35,   5,   2,  -9,   9,  0 },
  { 0,   0,   8,  -8,  -2,  -7,   4,   4,   2,   1,  -1,  25,  0 },
  { 0,   0,   1,  -1,   0,  -3,   1,   3,  -1,   1,  -1,   3, 0 },
  { 0,   0,   3,  -3,   0,  -6,   5,  -1,   2,   1,  -4,  21,  0 },
  { -7,   1,   5,   4,  -3,   5,  11,  13,  12,  -8,  11,  12,  0 },
  { -5,  -3,   6,  -2,  -3,   8,  14,  15,   2,  -7,  11,  16,  0 },
  { 2,  -1,  -6,  -5,  -2,  -2,  20,  14,  -4,   0,  -3,  25,  0 },
  { 3,   1,  -8,  -4,   0,  -8,  22,   5,  -3,   2, -10,  29,  0 },
  { 2,   1,  -7,  -1,   2, -11,  23,  -5,   0,   2, -10,  29,  0 },
  { -6,  -3,   8,   9,  -4,   8,   9,   7,  14,  -2,   8,   9,  0 },
  { 2,   1,  -4,  -7,   0,  -8,  17,  22,   1,  -1,  -4,  23,  0 },
  { 3,   0,  -5,  -7,   0,  -7,  15,  18,  -5,   0,  -5,  27,  0 },
  { 2,   0,   0,  -7,   1, -10,  13,  13,  -4,   2,  -7,  24,  0 },
  { 3,   3, -13,   4,  -2,  -5,   9,  21,  25,  -2,  -3,  12,  0 },
  { -5,  -2,   7,  -3,  -7,   9,   8,   9,  16,  -2,  15,  12,  0 },
  { 0,  -1,   0,  -7,  -5,   4,  11,  11,   8,  -6,  12,  21,  0 },
  { 3,  -2,  -3,  -8,  -4,  -1,  16,  15,  -2,  -3,   3,  26,  0 },
  { 2,   1,  -5,  -4,  -1,  -8,  16,   4,  -2,   1,  -7,  33,  0 },
  { 2,   1,  -4,  -2,   1, -10,  17,  -2,   0,   2, -11,  33,  0 },
  { 1,  -2,   7, -15, -16,  10,   8,   8,  20,  11,  14,  11,  0 },
  { 2,   2,   3, -13, -13,   4,   8,  12,   2,  -3,  16,  24,  0 },
  { 1,   4,   0,  -7,  -8,  -4,   9,   9,  -2,  -2,   8,  29,  0 },
  { 1,   1,   2,  -4,  -1,  -6,   6,   3,  -1,  -1,  -3,  30,  0 },
  { -7,   3,   2,  10,  -2,   3,   7,  11,  19,  -7,   8,  10, 0 },
  { 0,  -2,  -5,  -3,  -2,   4,  20,  15,  -1,  -3,  -1,  22,  0 },
  { 3,  -1,  -8,  -4,  -1,  -4,  22,   8,  -4,   2,  -8,  28,  0 },
  { 0,   3, -14,   3,   0,   1,  19,  17,   8,  -3,  -7,  20,  0 },
  { 0,   2,  -1,  -8,   3,  -6,   5,  21,   1,   1,  -9,  13,  0 },
  { -4,  -2,   8,  20,  -2,   2,   3,   5,  21,   4,   6,   1, 0 },
  { 2,  -2,  -3,  -9,  -4,   2,  14,  16,   3,  -6,   8,  24,  0 },
  { 2,   1,   5, -16,  -7,   2,   3,  11,  15,  -3,  11,  22,  0 },
  { 1,   2,   3, -11,  -2,  -5,   4,   8,   9,  -3,  -2,  26,  0 },
  { 0,  -1,  10,  -9,  -1,  -8,   2,   3,   4,   0,   0,  29,  0 },
  { 1,   2,   0,  -5,   1,  -9,   9,   3,   0,   1,  -7,  20,  0 },
  { -2,   8,  -6,  -4,   3,  -9,  -8,  45,  14,   2, -13,   7, 0 },
  { 1,  -1,  16, -19,  -8,  -4,  -3,   2,  19,   0,   4,  30,  0 },
  { 1,   1,  -3,   0,   2, -11,  15,  -5,   1,   2,  -9,  24,  0 },
  { 0,   1,  -2,   0,   1,  -4,   4,   0,   0,   1,  -4,   7,  0 },
  { 0,   1,   2,  -5,   1,  -6,   4,  10,  -2,   1,  -4,  10,  0 },
  { 3,   0,  -3,  -6,  -2,  -6,  14,   8,  -1,  -1,  -3,  31,  0 },
  { 0,   1,   0,  -2,   1,  -6,   5,   1,   0,   1,  -5,  13,  0 },
  { 3,   1,   9, -19, -21,   9,   7,   6,  13,   5,  15,  21,  0 },
  { 2,   4,   3, -12, -13,   1,   7,   8,   3,   0,  12,  26,  0 },
  { 3,   1,  -8,  -2,   0,  -6,  18,   2,  -2,   3, -10,  23,  0 },
  { 1,   1,  -4,  -1,   1,  -5,   8,   1,  -1,   2,  -5,  10,  0 },
  { 0,   1,  -1,   0,   0,  -2,   2,   0,   0,   1,  -2,   3,  0 },
  { 1,   1,  -2,  -7,   1,  -7,  14,  18,   0,   0,  -7,  21,  0 },
  { 0,   1,   0,  -2,   0,  -7,   8,   1,  -2,   0,  -3,  24,  0 },
  { 0,   1,   1,  -2,   2, -10,  10,   0,  -2,   1,  -7,  23,  0 },
  { 0,   2,   2, -11,   2,  -4,  -3,  39,   7,   1, -10,   9,  0 },
  { 1,   0,  13, -16,  -5,  -6,  -1,   8,   6,   0,   6,  29,  0 },
  { 1,   3,   1,  -6,  -4,  -7,   9,   6,  -3,  -2,   3,  33,  0 },
  { 4,   0, -17,  -1,  -1,   5,  26,   8,  -2,   3, -15,  30,  0 },
  { 0,   1,  -2,   0,   2,  -8,  12,  -6,   1,   1,  -6,  16,  0 },
  { 0,   0,   0,  -1,   1,  -4,   4,   0,   0,   0,  -3,  11,  0 },
  { 0,   1,   2,  -8,   2,  -6,   5,  15,   0,   2,  -7,   9,  0 },
  { 1,  -1,  12, -15,  -7,  -2,   3,   6,   6,  -1,   7,  30,  0 },
};
const int AdaptiveLoopFilter::m_classToFilterMapping[NUM_FIXED_FILTER_SETS][MAX_NUM_ALF_CLASSES] =
{
  { 8,   2,   2,   2,   3,   4,  53,   9,   9,  52,   4,   4,   5,   9,   2,   8,  10,   9,   1,   3,  39,  39,  10,   9,  52 },
  { 11,  12,  13,  14,  15,  30,  11,  17,  18,  19,  16,  20,  20,   4,  53,  21,  22,  23,  14,  25,  26,  26,  27,  28,  10 },
  { 16,  12,  31,  32,  14,  16,  30,  33,  53,  34,  35,  16,  20,   4,   7,  16,  21,  36,  18,  19,  21,  26,  37,  38,  39 },
  { 35,  11,  13,  14,  43,  35,  16,   4,  34,  62,  35,  35,  30,  56,   7,  35,  21,  38,  24,  40,  16,  21,  48,  57,  39 },
  { 11,  31,  32,  43,  44,  16,   4,  17,  34,  45,  30,  20,  20,   7,   5,  21,  22,  46,  40,  47,  26,  48,  63,  58,  10 },
  { 12,  13,  50,  51,  52,  11,  17,  53,  45,   9,  30,   4,  53,  19,   0,  22,  23,  25,  43,  44,  37,  27,  28,  10,  55 },
  { 30,  33,  62,  51,  44,  20,  41,  56,  34,  45,  20,  41,  41,  56,   5,  30,  56,  38,  40,  47,  11,  37,  42,  57,   8 },
  { 35,  11,  23,  32,  14,  35,  20,   4,  17,  18,  21,  20,  20,  20,   4,  16,  21,  36,  46,  25,  41,  26,  48,  49,  58 },
  { 12,  31,  59,  59,   3,  33,  33,  59,  59,  52,   4,  33,  17,  59,  55,  22,  36,  59,  59,  60,  22,  36,  59,  25,  55 },
  { 31,  25,  15,  60,  60,  22,  17,  19,  55,  55,  20,  20,  53,  19,  55,  22,  46,  25,  43,  60,  37,  28,  10,  55,  52 },
  { 12,  31,  32,  50,  51,  11,  33,  53,  19,  45,  16,   4,   4,  53,   5,  22,  36,  18,  25,  43,  26,  27,  27,  28,  10 },
  { 5,   2,  44,  52,   3,   4,  53,  45,   9,   3,   4,  56,   5,   0,   2,   5,  10,  47,  52,   3,  63,  39,  10,   9,  52 },
  { 12,  34,  44,  44,   3,  56,  56,  62,  45,   9,  56,  56,   7,   5,   0,  22,  38,  40,  47,  52,  48,  57,  39,  10,   9 },
  { 35,  11,  23,  14,  51,  35,  20,  41,  56,  62,  16,  20,  41,  56,   7,  16,  21,  38,  24,  40,  26,  26,  42,  57,  39 },
  { 33,  34,  51,  51,  52,  41,  41,  34,  62,   0,  41,  41,  56,   7,   5,  56,  38,  38,  40,  44,  37,  42,  57,  39,  10 },
  { 16,  31,  32,  15,  60,  30,   4,  17,  19,  25,  22,  20,   4,  53,  19,  21,  22,  46,  25,  55,  26,  48,  63,  58,  55 },
};
#endif
void AdaptiveLoopFilter::applyCcAlfFilter(CodingStructure &cs, ComponentID compID, const PelBuf &dstBuf,
                                          const PelUnitBuf &recYuvExt, uint8_t *filterControl,
                                          const short filterSet[MAX_NUM_CC_ALF_FILTERS][MAX_NUM_CC_ALF_CHROMA_COEFF],
                                          const int   selectedFilterIdx)
{
  bool clipTop = false, clipBottom = false, clipLeft = false, clipRight = false;
  int  numHorVirBndry = 0, numVerVirBndry = 0;
  int  horVirBndryPos[] = { 0, 0, 0 };
  int  verVirBndryPos[] = { 0, 0, 0 };

  int ctuIdx = 0;
  for( int yPos = 0; yPos < m_picHeight; yPos += m_maxCUHeight )
  {
    for( int xPos = 0; xPos < m_picWidth; xPos += m_maxCUWidth )
    {
      int filterIdx =
        (filterControl == nullptr)
          ? selectedFilterIdx
          : filterControl[(yPos >> cs.pcv->maxCUHeightLog2) * cs.pcv->widthInCtus + (xPos >> cs.pcv->maxCUWidthLog2)];
      bool skipFiltering = (filterControl != nullptr && filterIdx == 0) ? true : false;
      if (!skipFiltering)
      {
        if (filterControl != nullptr)
        {
          filterIdx--;
        }

        const int16_t *filterCoeff = filterSet[filterIdx];

        const int width        = (xPos + m_maxCUWidth > m_picWidth) ? (m_picWidth - xPos) : m_maxCUWidth;
        const int height       = (yPos + m_maxCUHeight > m_picHeight) ? (m_picHeight - yPos) : m_maxCUHeight;
        const int chromaScaleX = getComponentScaleX(compID, m_chromaFormat);
        const int chromaScaleY = getComponentScaleY(compID, m_chromaFormat);

        int rasterSliceAlfPad = 0;
        if (isCrossedByVirtualBoundaries(cs, xPos, yPos, width, height, clipTop, clipBottom, clipLeft, clipRight,
                                         numHorVirBndry, numVerVirBndry, horVirBndryPos, verVirBndryPos,
                                         rasterSliceAlfPad))
        {
          int yStart = yPos;
          for (int i = 0; i <= numHorVirBndry; i++)
          {
            const int  yEnd   = i == numHorVirBndry ? yPos + height : horVirBndryPos[i];
            const int  h      = yEnd - yStart;
            const bool clipT  = (i == 0 && clipTop) || (i > 0) || (yStart == 0);
            const bool clipB  = (i == numHorVirBndry && clipBottom) || (i < numHorVirBndry) || (yEnd == m_picHeight);
            int        xStart = xPos;
            for (int j = 0; j <= numVerVirBndry; j++)
            {
              const int  xEnd  = j == numVerVirBndry ? xPos + width : verVirBndryPos[j];
              const int  w     = xEnd - xStart;
              const bool clipL = (j == 0 && clipLeft) || (j > 0) || (xStart == 0);
              const bool clipR = (j == numVerVirBndry && clipRight) || (j < numVerVirBndry) || (xEnd == m_picWidth);
              const int  wBuf  = w + (clipL ? 0 : MAX_ALF_PADDING_SIZE) + (clipR ? 0 : MAX_ALF_PADDING_SIZE);
              const int  hBuf  = h + (clipT ? 0 : MAX_ALF_PADDING_SIZE) + (clipB ? 0 : MAX_ALF_PADDING_SIZE);
              PelUnitBuf buf   = m_tempBuf2.subBuf(UnitArea(cs.area.chromaFormat, Area(0, 0, wBuf, hBuf)));
              buf.copyFrom(recYuvExt.subBuf(
                UnitArea(cs.area.chromaFormat, Area(xStart - (clipL ? 0 : MAX_ALF_PADDING_SIZE),
                                                    yStart - (clipT ? 0 : MAX_ALF_PADDING_SIZE), wBuf, hBuf))));
              // pad top-left unavailable samples for raster slice
              if (xStart == xPos && yStart == yPos && (rasterSliceAlfPad & 1))
              {
                buf.padBorderPel(MAX_ALF_PADDING_SIZE, 1);
              }

              // pad bottom-right unavailable samples for raster slice
              if (xEnd == xPos + width && yEnd == yPos + height && (rasterSliceAlfPad & 2))
              {
                buf.padBorderPel(MAX_ALF_PADDING_SIZE, 2);
              }
              buf.extendBorderPel(MAX_ALF_PADDING_SIZE);
              buf = buf.subBuf(UnitArea(
                cs.area.chromaFormat, Area(clipL ? 0 : MAX_ALF_PADDING_SIZE, clipT ? 0 : MAX_ALF_PADDING_SIZE, w, h)));

              const Area blkSrc(0, 0, w, h);

              const Area blkDst(xStart >> chromaScaleX, yStart >> chromaScaleY, w >> chromaScaleX, h >> chromaScaleY);
#if ALF_IMPROVEMENT
              m_filterCcAlf( dstBuf, buf, blkDst, blkSrc, compID, filterCoeff, m_clpRngs, cs );
#else
              m_filterCcAlf(dstBuf, buf, blkDst, blkSrc, compID, filterCoeff, m_clpRngs, cs, m_alfVBLumaCTUHeight, m_alfVBLumaPos);
#endif

              xStart = xEnd;
            }

            yStart = yEnd;
          }
        }
        else
        {
          const UnitArea area(m_chromaFormat, Area(xPos, yPos, width, height));

          Area blkDst(xPos >> chromaScaleX, yPos >> chromaScaleY, width >> chromaScaleX, height >> chromaScaleY);
          Area blkSrc(xPos, yPos, width, height);
#if ALF_IMPROVEMENT
          m_filterCcAlf( dstBuf, recYuvExt, blkDst, blkSrc, compID, filterCoeff, m_clpRngs, cs );
#else
          m_filterCcAlf(dstBuf, recYuvExt, blkDst, blkSrc, compID, filterCoeff, m_clpRngs, cs, m_alfVBLumaCTUHeight, m_alfVBLumaPos);
#endif
        }
      }
      ctuIdx++;
    }
  }
}

void AdaptiveLoopFilter::ALFProcess(CodingStructure& cs)
{

  // set clipping range
  m_clpRngs = cs.slice->getClpRngs();

  // set CTU enable flags
  for( int compIdx = 0; compIdx < MAX_NUM_COMPONENT; compIdx++ )
  {
    m_ctuEnableFlag[compIdx] = cs.picture->getAlfCtuEnableFlag( compIdx );
    m_ctuAlternative[compIdx] = cs.picture->getAlfCtuAlternativeData( compIdx );
  }
  short* alfCtuFilterIndex = nullptr;
  uint32_t lastSliceIdx = 0xFFFFFFFF;

  PelUnitBuf recYuv = cs.getRecoBuf();
  m_tempBuf.copyFrom( recYuv );
  PelUnitBuf tmpYuv = m_tempBuf.getBuf( cs.area );
#if ALF_IMPROVEMENT
  tmpYuv.extendBorderPel( MAX_FILTER_LENGTH_FIXED >> 1 );
#else
  tmpYuv.extendBorderPel( MAX_ALF_FILTER_LENGTH >> 1 );
#endif

  const PreCalcValues& pcv = *cs.pcv;

  int ctuIdx = 0;
  bool clipTop = false, clipBottom = false, clipLeft = false, clipRight = false;
  int numHorVirBndry = 0, numVerVirBndry = 0;
  int horVirBndryPos[] = { 0, 0, 0 };
  int verVirBndryPos[] = { 0, 0, 0 };

#if ALF_IMPROVEMENT
  int fixedFilterSetIdx = cs.slice->getTileGroupAlfFixedFilterSetIdx();
#endif

  for( int yPos = 0; yPos < pcv.lumaHeight; yPos += pcv.maxCUHeight )
  {
    for( int xPos = 0; xPos < pcv.lumaWidth; xPos += pcv.maxCUWidth )
    {
      // get first CU in CTU
      const CodingUnit *cu = cs.getCU( Position(xPos, yPos), CHANNEL_TYPE_LUMA );

      // skip this CTU if ALF is disabled
      if (!cu->slice->getTileGroupAlfEnabledFlag(COMPONENT_Y) && !cu->slice->getTileGroupAlfEnabledFlag(COMPONENT_Cb) && !cu->slice->getTileGroupAlfEnabledFlag(COMPONENT_Cr))
      {
        ctuIdx++;
        continue;
      }

      // reload ALF APS each time the slice changes during raster scan filtering
      if(ctuIdx == 0 || lastSliceIdx != cu->slice->getSliceID() || alfCtuFilterIndex==nullptr)
      {
        cs.slice = cu->slice;
        reconstructCoeffAPSs(cs, true, cu->slice->getTileGroupAlfEnabledFlag(COMPONENT_Cb) || cu->slice->getTileGroupAlfEnabledFlag(COMPONENT_Cr), false);
        alfCtuFilterIndex = cu->slice->getPic()->getAlfCtbFilterIndex();
        m_ccAlfFilterParam = cu->slice->m_ccAlfFilterParam;
      }
      lastSliceIdx = cu->slice->getSliceID();

      const int width = ( xPos + pcv.maxCUWidth > pcv.lumaWidth ) ? ( pcv.lumaWidth - xPos ) : pcv.maxCUWidth;
      const int height = ( yPos + pcv.maxCUHeight > pcv.lumaHeight ) ? ( pcv.lumaHeight - yPos ) : pcv.maxCUHeight;
      bool ctuEnableFlag = m_ctuEnableFlag[COMPONENT_Y][ctuIdx];
      for( int compIdx = 1; compIdx < MAX_NUM_COMPONENT; compIdx++ )
      {
        ctuEnableFlag |= m_ctuEnableFlag[compIdx][ctuIdx] > 0;
        if (cu->slice->m_ccAlfFilterParam.ccAlfFilterEnabled[compIdx - 1])
        {
          ctuEnableFlag |= m_ccAlfFilterControl[compIdx - 1][ctuIdx] > 0;
        }
      }
      int rasterSliceAlfPad = 0;
      if( ctuEnableFlag && isCrossedByVirtualBoundaries( cs, xPos, yPos, width, height, clipTop, clipBottom, clipLeft, clipRight, numHorVirBndry, numVerVirBndry, horVirBndryPos, verVirBndryPos, rasterSliceAlfPad ) )
      {
        int yStart = yPos;
        for( int i = 0; i <= numHorVirBndry; i++ )
        {
          const int yEnd = i == numHorVirBndry ? yPos + height : horVirBndryPos[i];
          const int h = yEnd - yStart;
          const bool clipT = ( i == 0 && clipTop ) || ( i > 0 ) || ( yStart == 0 );
          const bool clipB = ( i == numHorVirBndry && clipBottom ) || ( i < numHorVirBndry ) || ( yEnd == pcv.lumaHeight );
          int xStart = xPos;
          for( int j = 0; j <= numVerVirBndry; j++ )
          {
            const int xEnd = j == numVerVirBndry ? xPos + width : verVirBndryPos[j];
            const int w = xEnd - xStart;
            const bool clipL = ( j == 0 && clipLeft ) || ( j > 0 ) || ( xStart == 0 );
            const bool clipR = ( j == numVerVirBndry && clipRight ) || ( j < numVerVirBndry ) || ( xEnd == pcv.lumaWidth );
            const int wBuf = w + ( clipL ? 0 : MAX_ALF_PADDING_SIZE ) + ( clipR ? 0 : MAX_ALF_PADDING_SIZE );
            const int hBuf = h + ( clipT ? 0 : MAX_ALF_PADDING_SIZE ) + ( clipB ? 0 : MAX_ALF_PADDING_SIZE );
            PelUnitBuf buf = m_tempBuf2.subBuf( UnitArea( cs.area.chromaFormat, Area( 0, 0, wBuf, hBuf ) ) );
            buf.copyFrom( tmpYuv.subBuf( UnitArea( cs.area.chromaFormat, Area( xStart - ( clipL ? 0 : MAX_ALF_PADDING_SIZE ), yStart - ( clipT ? 0 : MAX_ALF_PADDING_SIZE ), wBuf, hBuf ) ) ) );
            // pad top-left unavailable samples for raster slice
            if( xStart == xPos && yStart == yPos && ( rasterSliceAlfPad & 1 ) )
            {
              buf.padBorderPel( MAX_ALF_PADDING_SIZE, 1 );
            }

            // pad bottom-right unavailable samples for raster slice
            if( xEnd == xPos + width && yEnd == yPos + height && ( rasterSliceAlfPad & 2 ) )
            {
              buf.padBorderPel( MAX_ALF_PADDING_SIZE, 2 );
            }
            buf.extendBorderPel( MAX_ALF_PADDING_SIZE );
            buf = buf.subBuf( UnitArea( cs.area.chromaFormat, Area( clipL ? 0 : MAX_ALF_PADDING_SIZE, clipT ? 0 : MAX_ALF_PADDING_SIZE, w, h ) ) );

            if( m_ctuEnableFlag[COMPONENT_Y][ctuIdx] )
            {
              const Area blkSrc( 0, 0, w, h );
              const Area blkDst( xStart, yStart, w, h );
              short filterSetIndex = alfCtuFilterIndex[ctuIdx];
              deriveClassification( m_classifier, buf.get( COMPONENT_Y ), blkDst, blkSrc 
#if ALF_IMPROVEMENT
              , cs, filterSetIndex < NUM_FIXED_FILTER_SETS ? filterSetIndex : -1
#endif
              );

              short *coeff;
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
              Pel *clip;
#else
              short *clip;
#endif
#if ALF_IMPROVEMENT
              if( filterSetIndex < NUM_FIXED_FILTER_SETS )
              {
                copyFixedFilterResults(recYuv, blkDst, COMPONENT_Y, m_fixFilterResult, fixedFilterSetIdx, filterSetIndex);
              }
              else
              {
                uint8_t alt_num = m_ctuAlternative[COMPONENT_Y][ctuIdx];
                coeff = m_coeffApsLuma[filterSetIndex - NUM_FIXED_FILTER_SETS][alt_num];
                clip = m_clippApsLuma[filterSetIndex - NUM_FIXED_FILTER_SETS][alt_num];
                AlfFilterType filterTypeCtb = m_filterTypeApsLuma[filterSetIndex - NUM_FIXED_FILTER_SETS];
                alfFiltering( m_classifier, recYuv, buf, blkDst, blkSrc, COMPONENT_Y, coeff, clip, m_clpRngs.comp[COMPONENT_Y], cs, filterTypeCtb, m_fixFilterResult, fixedFilterSetIdx );
              }
#else
              if( filterSetIndex >= NUM_FIXED_FILTER_SETS )
              {
                coeff = m_coeffApsLuma[filterSetIndex - NUM_FIXED_FILTER_SETS];
                clip = m_clippApsLuma[filterSetIndex - NUM_FIXED_FILTER_SETS];
              }
              else
              {
                coeff = m_fixedFilterSetCoeffDec[filterSetIndex];
                clip = m_clipDefault;
              }
              m_filter7x7Blk( m_classifier, recYuv, buf, blkDst, blkSrc, COMPONENT_Y, coeff, clip, m_clpRngs.comp[COMPONENT_Y], cs, m_alfVBLumaCTUHeight, m_alfVBLumaPos );
#endif
          }

            for( int compIdx = 1; compIdx < MAX_NUM_COMPONENT; compIdx++ )
            {
              ComponentID compID = ComponentID( compIdx );
              const int chromaScaleX = getComponentScaleX( compID, tmpYuv.chromaFormat );
              const int chromaScaleY = getComponentScaleY( compID, tmpYuv.chromaFormat );

              if( m_ctuEnableFlag[compIdx][ctuIdx] )
              {
                const Area blkSrc( 0, 0, w >> chromaScaleX, h >> chromaScaleY );
                const Area blkDst( xStart >> chromaScaleX, yStart >> chromaScaleY, w >> chromaScaleX, h >> chromaScaleY );
                uint8_t alt_num = m_ctuAlternative[compIdx][ctuIdx];
#if ALF_IMPROVEMENT
                alfFiltering( m_classifier, recYuv, buf, blkDst, blkSrc, compID, m_chromaCoeffFinal[alt_num], m_chromaClippFinal[alt_num], m_clpRngs.comp[compIdx], cs, m_filterTypeApsChroma, nullptr, -1 );
#else
                m_filter5x5Blk(m_classifier, recYuv, buf, blkDst, blkSrc, compID, m_chromaCoeffFinal[alt_num], m_chromaClippFinal[alt_num], m_clpRngs.comp[compIdx], cs, m_alfVBChmaCTUHeight, m_alfVBChmaPos );
#endif
              }
              if (cu->slice->m_ccAlfFilterParam.ccAlfFilterEnabled[compIdx - 1])
              {
                const int filterIdx = m_ccAlfFilterControl[compIdx - 1][ctuIdx];

                if (filterIdx != 0)
                {
                  const Area blkSrc(0, 0, w, h);
                  Area blkDst(xStart >> chromaScaleX, yStart >> chromaScaleY, w >> chromaScaleX, h >> chromaScaleY);

                  const int16_t *filterCoeff = m_ccAlfFilterParam.ccAlfCoeff[compIdx - 1][filterIdx - 1];
#if ALF_IMPROVEMENT
                  m_filterCcAlf( recYuv.get( compID ), buf, blkDst, blkSrc, compID, filterCoeff, m_clpRngs, cs );
#else
                  m_filterCcAlf(recYuv.get(compID), buf, blkDst, blkSrc, compID, filterCoeff, m_clpRngs, cs, m_alfVBLumaCTUHeight, m_alfVBLumaPos);
#endif
                }
              }
            }

            xStart = xEnd;
          }

          yStart = yEnd;
        }
      }
      else
      {
        const UnitArea area( cs.area.chromaFormat, Area( xPos, yPos, width, height ) );
        if( m_ctuEnableFlag[COMPONENT_Y][ctuIdx] )
        {
          Area blk( xPos, yPos, width, height );
          short filterSetIndex = alfCtuFilterIndex[ctuIdx];
          deriveClassification( m_classifier, tmpYuv.get( COMPONENT_Y ), blk, blk 
#if ALF_IMPROVEMENT
            , cs, filterSetIndex < NUM_FIXED_FILTER_SETS ? filterSetIndex : -1
#endif
          );          
          short *coeff;
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
          Pel *clip;
#else
          short *clip;
#endif
#if ALF_IMPROVEMENT
          if( filterSetIndex < NUM_FIXED_FILTER_SETS )
          {
            copyFixedFilterResults( recYuv, blk, COMPONENT_Y, m_fixFilterResult, fixedFilterSetIdx, filterSetIndex );
          }
          else
          {
            uint8_t alt_num = m_ctuAlternative[COMPONENT_Y][ctuIdx];
            coeff = m_coeffApsLuma[filterSetIndex - NUM_FIXED_FILTER_SETS][alt_num];
            clip = m_clippApsLuma[filterSetIndex - NUM_FIXED_FILTER_SETS][alt_num];
            AlfFilterType filterTypeCtb = m_filterTypeApsLuma[filterSetIndex - NUM_FIXED_FILTER_SETS];
            alfFiltering( m_classifier, recYuv, tmpYuv, blk, blk, COMPONENT_Y, coeff, clip, m_clpRngs.comp[COMPONENT_Y], cs, filterTypeCtb, m_fixFilterResult, fixedFilterSetIdx );
          }
#else
          if( filterSetIndex >= NUM_FIXED_FILTER_SETS )
          {
            coeff = m_coeffApsLuma[filterSetIndex - NUM_FIXED_FILTER_SETS];
            clip = m_clippApsLuma[filterSetIndex - NUM_FIXED_FILTER_SETS];
          }
          else
          {
            coeff = m_fixedFilterSetCoeffDec[filterSetIndex];
            clip = m_clipDefault;
          }
          m_filter7x7Blk( m_classifier, recYuv, tmpYuv, blk, blk, COMPONENT_Y, coeff, clip, m_clpRngs.comp[COMPONENT_Y], cs, m_alfVBLumaCTUHeight, m_alfVBLumaPos );
#endif
        }

        for( int compIdx = 1; compIdx < MAX_NUM_COMPONENT; compIdx++ )
        {
          ComponentID compID = ComponentID( compIdx );
          const int chromaScaleX = getComponentScaleX( compID, tmpYuv.chromaFormat );
          const int chromaScaleY = getComponentScaleY( compID, tmpYuv.chromaFormat );

          if (m_ctuEnableFlag[compIdx][ctuIdx])
          {
            Area    blk(xPos >> chromaScaleX, yPos >> chromaScaleY, width >> chromaScaleX, height >> chromaScaleY);
            uint8_t alt_num = m_ctuAlternative[compIdx][ctuIdx];

#if ALF_IMPROVEMENT
            alfFiltering( m_classifier, recYuv, tmpYuv, blk, blk, compID, m_chromaCoeffFinal[alt_num], m_chromaClippFinal[alt_num], m_clpRngs.comp[compIdx], cs, m_filterTypeApsChroma  , nullptr, -1 );
#else
            m_filter5x5Blk( m_classifier, recYuv, tmpYuv, blk, blk, compID, m_chromaCoeffFinal[alt_num], m_chromaClippFinal[alt_num], m_clpRngs.comp[compIdx], cs, m_alfVBChmaCTUHeight, m_alfVBChmaPos );
#endif
          }
          if (cu->slice->m_ccAlfFilterParam.ccAlfFilterEnabled[compIdx - 1])
          {
            const int filterIdx = m_ccAlfFilterControl[compIdx - 1][ctuIdx];

            if (filterIdx != 0)
            {
              Area blkDst(xPos >> chromaScaleX, yPos >> chromaScaleY, width >> chromaScaleX, height >> chromaScaleY);
              Area blkSrc(xPos, yPos, width, height);

              const int16_t *filterCoeff = m_ccAlfFilterParam.ccAlfCoeff[compIdx - 1][filterIdx - 1];

#if ALF_IMPROVEMENT
              m_filterCcAlf( recYuv.get( compID ), tmpYuv, blkDst, blkSrc, compID, filterCoeff, m_clpRngs, cs );
#else
              m_filterCcAlf(recYuv.get(compID), tmpYuv, blkDst, blkSrc, compID, filterCoeff, m_clpRngs, cs, m_alfVBLumaCTUHeight, m_alfVBLumaPos);
#endif
            }
          }
        }
      }
      ctuIdx++;
    }
  }
}

void AdaptiveLoopFilter::reconstructCoeffAPSs(CodingStructure& cs, bool luma, bool chroma, bool isRdo)
{
  //luma
  APS** aps = cs.slice->getAlfAPSs();
  AlfParam alfParamTmp;
  APS* curAPS;
  if (luma)
  {
    for (int i = 0; i < cs.slice->getTileGroupNumAps(); i++)
    {
      int apsIdx = cs.slice->getTileGroupApsIdLuma()[i];
      curAPS = aps[apsIdx];
      CHECK(curAPS == NULL, "invalid APS");
      alfParamTmp = curAPS->getAlfAPSParam();
      reconstructCoeff(alfParamTmp, CHANNEL_TYPE_LUMA, isRdo, true);
      memcpy(m_coeffApsLuma[i], m_coeffFinal, sizeof(m_coeffFinal));
      memcpy(m_clippApsLuma[i], m_clippFinal, sizeof(m_clippFinal));
#if ALF_IMPROVEMENT
      m_filterTypeApsLuma[i] = alfParamTmp.filterType[CHANNEL_TYPE_LUMA];
      m_numLumaAltAps[i] = alfParamTmp.numAlternativesLuma;
#endif
    }
  }

  //chroma
  if (chroma)
  {
    int apsIdxChroma = cs.slice->getTileGroupApsIdChroma();
    curAPS = aps[apsIdxChroma];
    m_alfParamChroma = &curAPS->getAlfAPSParam();
    alfParamTmp = *m_alfParamChroma;
    reconstructCoeff(alfParamTmp, CHANNEL_TYPE_CHROMA, isRdo, true);
#if ALF_IMPROVEMENT
    m_filterTypeApsChroma = alfParamTmp.filterType[CHANNEL_TYPE_CHROMA];
#endif
  }
}

void AdaptiveLoopFilter::reconstructCoeff( AlfParam& alfParam, ChannelType channel, const bool isRdo, const bool isRedo )
{
  int factor = isRdo ? 0 : (1 << (m_NUM_BITS - 1));
#if ALF_IMPROVEMENT
  AlfFilterType filterType = alfParam.filterType[channel];
  int numClasses = isLuma(channel) ? MAX_NUM_ALF_CLASSES : 1;
  int numCoeff = m_filterShapes[channel][m_filterTypeToStatIndex[channel][filterType]].numCoeff;
#else
  AlfFilterType filterType = isLuma( channel ) ? ALF_FILTER_7 : ALF_FILTER_5;
  int numClasses = isLuma( channel ) ? MAX_NUM_ALF_CLASSES : 1;
  int numCoeff = filterType == ALF_FILTER_5 ? 7 : 13;
#endif
  int numCoeffMinus1 = numCoeff - 1;
#if ALF_IMPROVEMENT
  const int numAlts = isLuma( channel ) ? alfParam.numAlternativesLuma : alfParam.numAlternativesChroma;
#else
  const int numAlts = isLuma( channel ) ? 1 : alfParam.numAlternativesChroma;
#endif

  for( int altIdx = 0; altIdx < numAlts; ++ altIdx )
  {
#if ALF_IMPROVEMENT
    int numFilters = isLuma( channel ) ? alfParam.numLumaFilters[altIdx] : 1;
    short* coeff = isLuma( channel ) ? alfParam.lumaCoeff[altIdx] : alfParam.chromaCoeff[altIdx];
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
    Pel* clipp = isLuma( channel ) ? alfParam.lumaClipp[altIdx] : alfParam.chromaClipp[altIdx];
#else
    short* clipp = isLuma( channel ) ? alfParam.lumaClipp[altIdx] : alfParam.chromaClipp[altIdx];
#endif
#else
    int numFilters = isLuma( channel ) ? alfParam.numLumaFilters : 1;
    short* coeff = isLuma( channel ) ? alfParam.lumaCoeff : alfParam.chromaCoeff[altIdx];
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
    Pel* clipp = isLuma( channel ) ? alfParam.lumaClipp : alfParam.chromaClipp[altIdx];
#else
    short* clipp = isLuma( channel ) ? alfParam.lumaClipp : alfParam.chromaClipp[altIdx];
#endif
#endif

    for( int filterIdx = 0; filterIdx < numFilters; filterIdx++ )
    {
      coeff[filterIdx* MAX_NUM_ALF_LUMA_COEFF + numCoeffMinus1] = factor;
    }

    if( isChroma( channel ) )
    {
      for( int coeffIdx = 0; coeffIdx < numCoeffMinus1; ++coeffIdx )
      {
        m_chromaCoeffFinal[altIdx][coeffIdx] = coeff[coeffIdx];
#if ALF_IMPROVEMENT
        int clipIdx = alfParam.nonLinearFlag[channel][altIdx] ? clipp[coeffIdx] : 0;
#else
        int clipIdx = alfParam.nonLinearFlag[channel] ? clipp[coeffIdx] : 0;
#endif
        m_chromaClippFinal[altIdx][coeffIdx] = isRdo ? clipIdx : m_alfClippingValues[channel][clipIdx];
      }
      m_chromaCoeffFinal[altIdx][numCoeffMinus1] = factor;
      m_chromaClippFinal[altIdx][numCoeffMinus1] = isRdo ? 0 : m_alfClippingValues[channel][0];
      continue;
    }
    for( int classIdx = 0; classIdx < numClasses; classIdx++ )
    {
#if ALF_IMPROVEMENT
      int filterIdx = alfParam.filterCoeffDeltaIdx[altIdx][classIdx];
      CHECK(!(filterIdx >= 0 && filterIdx < alfParam.numLumaFilters[altIdx]), "Bad coeff delta idx in ALF");
      for (int coeffIdx = 0; coeffIdx < numCoeffMinus1; ++coeffIdx)
      {
        m_coeffFinal[altIdx][classIdx * MAX_NUM_ALF_LUMA_COEFF + coeffIdx] = coeff[filterIdx * MAX_NUM_ALF_LUMA_COEFF + coeffIdx];
      }
      m_coeffFinal[altIdx][classIdx* MAX_NUM_ALF_LUMA_COEFF + numCoeffMinus1] = factor;
      m_clippFinal[altIdx][classIdx* MAX_NUM_ALF_LUMA_COEFF + numCoeffMinus1] = isRdo ? 0 : m_alfClippingValues[channel][0];
      for (int coeffIdx = 0; coeffIdx < numCoeffMinus1; ++coeffIdx)
      {
        int clipIdx = alfParam.nonLinearFlag[channel][altIdx] ? clipp[filterIdx * MAX_NUM_ALF_LUMA_COEFF + coeffIdx] : 0;
        CHECK(!(clipIdx >= 0 && clipIdx < MaxAlfNumClippingValues), "Bad clip idx in ALF");
        m_clippFinal[altIdx][classIdx * MAX_NUM_ALF_LUMA_COEFF + coeffIdx] = isRdo ? clipIdx : m_alfClippingValues[channel][clipIdx];
      }
      m_clippFinal[altIdx][classIdx* MAX_NUM_ALF_LUMA_COEFF + numCoeffMinus1] = isRdo ? 0 : m_alfClippingValues[channel][0];
#else
      int filterIdx = alfParam.filterCoeffDeltaIdx[classIdx];

      CHECK(!(filterIdx >= 0 && filterIdx < alfParam.numLumaFilters), "Bad coeff delta idx in ALF");
      for (int coeffIdx = 0; coeffIdx < numCoeffMinus1; ++coeffIdx)
      {
        m_coeffFinal[classIdx * MAX_NUM_ALF_LUMA_COEFF + coeffIdx] = coeff[filterIdx * MAX_NUM_ALF_LUMA_COEFF + coeffIdx];
      }
      m_coeffFinal[classIdx* MAX_NUM_ALF_LUMA_COEFF + numCoeffMinus1] = factor;
      m_clippFinal[classIdx* MAX_NUM_ALF_LUMA_COEFF + numCoeffMinus1] = isRdo ? 0 : m_alfClippingValues[channel][0];
      for( int coeffIdx = 0; coeffIdx < numCoeffMinus1; ++coeffIdx )
      {
        int clipIdx = alfParam.nonLinearFlag[channel] ? clipp[filterIdx * MAX_NUM_ALF_LUMA_COEFF + coeffIdx] : 0;
        CHECK(!(clipIdx >= 0 && clipIdx < MaxAlfNumClippingValues), "Bad clip idx in ALF");
        m_clippFinal[classIdx * MAX_NUM_ALF_LUMA_COEFF + coeffIdx] = isRdo ? clipIdx : m_alfClippingValues[channel][clipIdx];
      }
      m_clippFinal[classIdx* MAX_NUM_ALF_LUMA_COEFF + numCoeffMinus1] = isRdo ? 0 : m_alfClippingValues[channel][0];
#endif
    }
  }
}

void AdaptiveLoopFilter::create(const int picWidth, const int picHeight, const ChromaFormat format, const int maxCUWidth, const int maxCUHeight, const int maxCUDepth, const int inputBitDepth[MAX_NUM_CHANNEL_TYPE])
{
  destroy();
  std::memcpy(m_inputBitDepth, inputBitDepth, sizeof(m_inputBitDepth));
  m_picWidth = picWidth;
  m_picHeight = picHeight;
  m_maxCUWidth = maxCUWidth;
  m_maxCUHeight = maxCUHeight;
  m_maxCUDepth = maxCUDepth;
  m_chromaFormat = format;

  m_numCTUsInWidth = (m_picWidth / m_maxCUWidth) + ((m_picWidth % m_maxCUWidth) ? 1 : 0);
  m_numCTUsInHeight = (m_picHeight / m_maxCUHeight) + ((m_picHeight % m_maxCUHeight) ? 1 : 0);
  m_numCTUsInPic = m_numCTUsInHeight * m_numCTUsInWidth;
  m_filterShapesCcAlf[0].push_back(AlfFilterShape(size_CC_ALF));
  m_filterShapesCcAlf[1].push_back(AlfFilterShape(size_CC_ALF));

#if ALF_IMPROVEMENT
  m_filterShapes[CHANNEL_TYPE_LUMA].push_back(AlfFilterShape(5));
  m_filterShapes[CHANNEL_TYPE_LUMA].push_back(AlfFilterShape(7));
  m_filterShapes[CHANNEL_TYPE_LUMA].push_back(AlfFilterShape(9));
  m_filterShapes[CHANNEL_TYPE_LUMA].push_back(AlfFilterShape(size_ALF_FILTER_9_EXT));
  m_filterShapes[CHANNEL_TYPE_LUMA].push_back(AlfFilterShape(size_ALF_FILTER_EXT));
  m_filterShapes[CHANNEL_TYPE_CHROMA].push_back(AlfFilterShape(5));
  m_filterShapes[CHANNEL_TYPE_CHROMA].push_back(AlfFilterShape(7));
  m_filterShapes[CHANNEL_TYPE_CHROMA].push_back(AlfFilterShape(9));
  memset(m_filterTypeTest, 0, sizeof(m_filterTypeTest));
  memset(m_filterTypeToStatIndex, 0, sizeof(m_filterTypeToStatIndex));
  m_filterTypeToStatIndex[CHANNEL_TYPE_LUMA][ALF_FILTER_5] = 0;
  m_filterTypeToStatIndex[CHANNEL_TYPE_LUMA][ALF_FILTER_7] = 1;
  m_filterTypeToStatIndex[CHANNEL_TYPE_LUMA][ALF_FILTER_9] = 2;
  m_filterTypeToStatIndex[CHANNEL_TYPE_LUMA][ALF_FILTER_9_EXT] = 3;
  m_filterTypeToStatIndex[CHANNEL_TYPE_LUMA][ALF_FILTER_EXT] = 4;
  m_filterTypeToStatIndex[CHANNEL_TYPE_CHROMA][ALF_FILTER_5] = 0;
  m_filterTypeToStatIndex[CHANNEL_TYPE_CHROMA][ALF_FILTER_7] = 1;
  m_filterTypeToStatIndex[CHANNEL_TYPE_CHROMA][ALF_FILTER_9] = 2;
  m_filterTypeTest[CHANNEL_TYPE_LUMA][ALF_FILTER_9_EXT] = true;
  m_filterTypeTest[CHANNEL_TYPE_LUMA][ALF_FILTER_EXT] = false;
  m_filterTypeTest[CHANNEL_TYPE_CHROMA][ALF_FILTER_9] = true;
#else
  m_filterShapes[CHANNEL_TYPE_LUMA].push_back(AlfFilterShape(7));
  m_filterShapes[CHANNEL_TYPE_CHROMA].push_back(AlfFilterShape(5));
  m_alfVBLumaPos = m_maxCUHeight - ALF_VB_POS_ABOVE_CTUROW_LUMA;
  m_alfVBChmaPos = (m_maxCUHeight >> ((m_chromaFormat == CHROMA_420) ? 1 : 0)) - ALF_VB_POS_ABOVE_CTUROW_CHMA;

  m_alfVBLumaCTUHeight = m_maxCUHeight;
  m_alfVBChmaCTUHeight = (m_maxCUHeight >> ((m_chromaFormat == CHROMA_420) ? 1 : 0));
#endif

  static_assert(AlfNumClippingValues[CHANNEL_TYPE_LUMA] > 0, "AlfNumClippingValues[CHANNEL_TYPE_LUMA] must be at least one");
  m_alfClippingValues[CHANNEL_TYPE_LUMA][0] = 1 << m_inputBitDepth[CHANNEL_TYPE_LUMA];
  int shiftLuma = m_inputBitDepth[CHANNEL_TYPE_LUMA] - 8;
  for (int i = 1; i < AlfNumClippingValues[CHANNEL_TYPE_LUMA]; ++i)
  {
    m_alfClippingValues[CHANNEL_TYPE_LUMA][i] = 1 << (7 - 2 * i + shiftLuma);
  }
  static_assert(AlfNumClippingValues[CHANNEL_TYPE_CHROMA] > 0, "AlfNumClippingValues[CHANNEL_TYPE_CHROMA] must be at least one");
  m_alfClippingValues[CHANNEL_TYPE_CHROMA][0] = 1 << m_inputBitDepth[CHANNEL_TYPE_CHROMA];
  int shiftChroma = m_inputBitDepth[CHANNEL_TYPE_CHROMA] - 8;
  for (int i = 1; i < AlfNumClippingValues[CHANNEL_TYPE_CHROMA]; ++i)
  {
    m_alfClippingValues[CHANNEL_TYPE_CHROMA][i] = 1 << (7 - 2 * i + shiftChroma);
  }
  if (m_created)
  {
    return;
  }

  m_tempBuf.destroy();
  // NOTE: make border 1 sample wider to avoid out-of-bounds memory access in SIMD code (simdDeriveClassificationBlk
  // function)
#if ALF_IMPROVEMENT
  m_tempBuf.create(format, Area(0, 0, picWidth, picHeight), maxCUWidth, MAX_FILTER_LENGTH_FIXED, 0, false);
#else
  m_tempBuf.create(format, Area(0, 0, picWidth, picHeight), maxCUWidth, (MAX_ALF_FILTER_LENGTH + 1) >> 1, 0, false);
#endif
  m_tempBuf2.destroy();
  m_tempBuf2.create( format, Area( 0, 0, maxCUWidth + (MAX_ALF_PADDING_SIZE << 1), maxCUHeight + (MAX_ALF_PADDING_SIZE << 1) ), maxCUWidth, MAX_ALF_PADDING_SIZE, 0, false );

#if ALF_IMPROVEMENT
  int numFixedFilters = EXT_LENGTH << 1;
  //destroy
  if (m_fixFilterResult)
  {
    for (int i = 0; i < numFixedFilters; i++)
    {
      if (m_fixFilterResult[i])
      {
        //for( int j = 0; j < m_picWidth; j++ )
        for (int j = 0; j < m_picHeight; j++)
        {
          if (m_fixFilterResult[i][j])
          {
            delete[] m_fixFilterResult[i][j];
            m_fixFilterResult[i][j] = nullptr;
          }
        }
        delete[] m_fixFilterResult[i];
        m_fixFilterResult[i] = nullptr;
      }
    }
    delete[] m_fixFilterResult;
    m_fixFilterResult = nullptr;
  }
  //create
  m_fixFilterResult = new Pel**[numFixedFilters];
  for (int i = 0; i < numFixedFilters; i++)
  {
    m_fixFilterResult[i] = new Pel*[picHeight];
    for (int j = 0; j < picHeight; j++)
    {
      m_fixFilterResult[i][j] = new Pel[picWidth];
    }
  }
#endif

  // Classification
  if ( m_classifier == nullptr )
  {
    m_classifier = new AlfClassifier*[picHeight];
    m_classifier[0] = new AlfClassifier[picWidth * picHeight];

    for (int i = 1; i < picHeight; i++)
    {
      m_classifier[i] = m_classifier[0] + i * picWidth;
    }
  }
#if !ALF_IMPROVEMENT
  for (int filterSetIndex = 0; filterSetIndex < NUM_FIXED_FILTER_SETS; filterSetIndex++)
  {
    for (int classIdx = 0; classIdx < MAX_NUM_ALF_CLASSES; classIdx++)
    {
      int fixedFilterIdx = m_classToFilterMapping[filterSetIndex][classIdx];
      for (int i = 0; i < MAX_NUM_ALF_LUMA_COEFF - 1; i++)
      {
       m_fixedFilterSetCoeffDec[filterSetIndex][classIdx * MAX_NUM_ALF_LUMA_COEFF + i] = m_fixedFilterSetCoeff[fixedFilterIdx][i];
      }
      m_fixedFilterSetCoeffDec[filterSetIndex][classIdx * MAX_NUM_ALF_LUMA_COEFF + MAX_NUM_ALF_LUMA_COEFF - 1] = (1 << (m_NUM_BITS - 1));
    }
  }
  for (int i = 0; i < MAX_NUM_ALF_LUMA_COEFF * MAX_NUM_ALF_CLASSES; i++)
  {
    m_clipDefault[i] = m_alfClippingValues[CHANNEL_TYPE_LUMA][0];
  }
#endif
  m_created = true;

  m_ccAlfFilterControl[0] = new uint8_t[m_numCTUsInPic];
  m_ccAlfFilterControl[1] = new uint8_t[m_numCTUsInPic];
}

void AdaptiveLoopFilter::destroy()
{
  if (!m_created)
  {
    return;
  }

  if( m_classifier )
  {
    delete[] m_classifier[0];
    delete[] m_classifier;
    m_classifier = nullptr;
  }

#if ALF_IMPROVEMENT
  if( m_fixFilterResult )
  {
    int numFixedFilters = EXT_LENGTH << 1;
    //for( int i = 0; i < m_picHeight; i++ )
    for( int i = 0; i < numFixedFilters; i++ )
    {
      if( m_fixFilterResult[i] )
      {
        //for( int j = 0; j < m_picWidth; j++ )
        for( int j = 0; j < m_picHeight; j++ )
        {
          if( m_fixFilterResult[i][j] )
          {
            delete[] m_fixFilterResult[i][j];
            m_fixFilterResult[i][j] = nullptr;
          }
        }
        delete[] m_fixFilterResult[i];
        m_fixFilterResult[i] = nullptr;
      }
    }
    delete[] m_fixFilterResult;
    m_fixFilterResult = nullptr;
  }
#endif

  m_tempBuf.destroy();
  m_tempBuf2.destroy();
  m_filterShapes[CHANNEL_TYPE_LUMA].clear();
  m_filterShapes[CHANNEL_TYPE_CHROMA].clear();
  m_created = false;

  m_filterShapesCcAlf[0].clear();
  m_filterShapesCcAlf[1].clear();
  if ( m_ccAlfFilterControl[0] )
  {
    delete [] m_ccAlfFilterControl[0];
    m_ccAlfFilterControl[0] = nullptr;
  }

  if ( m_ccAlfFilterControl[1] )
  {
    delete [] m_ccAlfFilterControl[1];
    m_ccAlfFilterControl[1] = nullptr;
  }
}

#if ALF_IMPROVEMENT
void  AdaptiveLoopFilter::copyFixedFilterResults( const PelUnitBuf &recDst, const Area &blkDst, ComponentID compId, Pel*** fixedFilterResults, const int fixedFilterSetIdx, const int classifierIdx )
{
  PelBuf dstLuma = recDst.get( compId );
  const int dstStride = dstLuma.stride;
  Pel* dst = dstLuma.buf + blkDst.y * dstStride + blkDst.x;
  int filterIndex = classifierIdx * NUM_FIXED_FILTER_SETS + fixedFilterSetIdx;
  for( int y = 0; y < blkDst.height; y++ )
  {
    for( int x = 0; x < blkDst.width; x++ )
    {
      dst[x] = fixedFilterResults[filterIndex][y + blkDst.y][x + blkDst.x];
    }
    dst += dstStride;
  }
}

void  AdaptiveLoopFilter::alfFiltering( AlfClassifier **classifier, const PelUnitBuf &recDst, const CPelUnitBuf &recSrc, const Area &blkDst, const Area &blk, const ComponentID compId, const short *filterSet, 
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  const Pel *fClipSet,
#else
  const short *fClipSet, 
#endif
  const ClpRng &clpRng, CodingStructure &cs, AlfFilterType filterType, Pel*** fixedFilterResults, int fixedFilterSetIdx )
{
  if (filterType == ALF_FILTER_9)
  {
    m_filter9x9Blk( classifier, recDst, recSrc, blkDst, blk, compId, filterSet, fClipSet, clpRng, cs, nullptr, fixedFilterSetIdx );
  }
  else if (filterType == ALF_FILTER_7)
  {
    m_filter7x7Blk( classifier, recDst, recSrc, blkDst, blk, compId, filterSet, fClipSet, clpRng, cs, nullptr, fixedFilterSetIdx );
  }
  else if (filterType == ALF_FILTER_5)
  {
    m_filter5x5Blk( classifier, recDst, recSrc, blkDst, blk, compId, filterSet, fClipSet, clpRng, cs, nullptr, fixedFilterSetIdx );
  }
  else if (filterType == ALF_FILTER_9_EXT)
  {
    m_filter9x9BlkExt( classifier, recDst, recSrc, blkDst, blk, compId, filterSet, fClipSet, clpRng, cs, fixedFilterResults, fixedFilterSetIdx );
  }
  else if (filterType == ALF_FILTER_EXT)
  {
    m_filterBlkExt( classifier, recDst, recSrc, blkDst, blk, compId, filterSet, fClipSet, clpRng, cs, fixedFilterResults, fixedFilterSetIdx );
  }
  else
  {
    CHECK(1, "unsupported filter type");
  }
}
#endif

void AdaptiveLoopFilter::deriveClassification( AlfClassifier** classifier, const CPelBuf& srcLuma, const Area& blkDst, const Area& blk 
#if ALF_IMPROVEMENT
  , CodingStructure &cs, const int classifierIdx
#endif
)
{
#if ALF_IMPROVEMENT
  if( cs.slice->getCuQpDeltaSubdiv() )
  {
    UnitArea curArea( cs.area.chromaFormat, blkDst );
    for( auto &currCU : cs.traverseCUs( curArea, CH_L ) )
    {
      deriveClassificationAndFixFilterResultsBlk( classifier, m_fixFilterResult, srcLuma, Area(currCU.lumaPos().x, currCU.lumaPos().y, currCU.lwidth(), currCU.lheight()), Area(currCU.lumaPos().x, currCU.lumaPos().y, currCU.lwidth(), currCU.lheight()), m_inputBitDepth[CHANNEL_TYPE_LUMA], cs, m_clpRngs.comp[COMPONENT_Y], m_alfClippingValues[CHANNEL_TYPE_LUMA], currCU.qp, cs.slice->getTileGroupAlfFixedFilterSetIdx(), m_mappingDir, m_laplacian, classifierIdx );
    }
  }
  else
  {
    int height = blk.pos().y + blk.height;
    int width = blk.pos().x + blk.width;
    for( int i = blk.pos().y; i < height; i += m_CLASSIFICATION_BLK_SIZE )
    {
      int nHeight = std::min( i + m_CLASSIFICATION_BLK_SIZE, height ) - i;

      for( int j = blk.pos().x; j < width; j += m_CLASSIFICATION_BLK_SIZE )
      {
        int nWidth = std::min( j + m_CLASSIFICATION_BLK_SIZE, width ) - j;
        deriveClassificationAndFixFilterResultsBlk( classifier, m_fixFilterResult, srcLuma, Area(j - blk.pos().x + blkDst.pos().x, i - blk.pos().y + blkDst.pos().y, nWidth, nHeight), Area(j, i, nWidth, nHeight), m_inputBitDepth[CHANNEL_TYPE_LUMA], cs, m_clpRngs.comp[COMPONENT_Y], m_alfClippingValues[CHANNEL_TYPE_LUMA], cs.slice->getSliceQp(), cs.slice->getTileGroupAlfFixedFilterSetIdx(), m_mappingDir, m_laplacian, classifierIdx );
      }
    }
  }
#else
  int height = blk.pos().y + blk.height;
  int width = blk.pos().x + blk.width;

  for( int i = blk.pos().y; i < height; i += m_CLASSIFICATION_BLK_SIZE )
  {
    int nHeight = std::min( i + m_CLASSIFICATION_BLK_SIZE, height ) - i;

    for( int j = blk.pos().x; j < width; j += m_CLASSIFICATION_BLK_SIZE )
    {
      int nWidth = std::min( j + m_CLASSIFICATION_BLK_SIZE, width ) - j;
      m_deriveClassificationBlk(classifier, m_laplacian, srcLuma, Area( j - blk.pos().x + blkDst.pos().x, i - blk.pos().y + blkDst.pos().y, nWidth, nHeight ), Area(j, i, nWidth, nHeight), m_inputBitDepth[CHANNEL_TYPE_LUMA] + 4, m_alfVBLumaCTUHeight, m_alfVBLumaPos);
    }
  }
#endif
}

#if ALF_IMPROVEMENT
int AdaptiveLoopFilter::assignAct( int avg_varPrec, int shift, int noAct )
{
  if (noAct == 5)
  {
    int th[] = { 0, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4 };
    return th[std::min(avg_varPrec >> shift, 15)];
  }
  else
  {
    int avg_var = avg_varPrec >> (shift - 4);
    static const int th[193] = { 0,  1,  2,  3,  4,  4,  5,  5,  6,  6,  6,  6,  7,  7,  7,  7,  8,  8,  8,  8,  8,  8,  8,  8,  9,  9,  9,  9,  9,  9,  9,  9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 15 };
    return th[std::min(avg_var, 192)];
  }
}

void AdaptiveLoopFilter::calcClass(AlfClassifier **classifier, const Area &blkDst, const Area &curBlk, int dirWindSize, int classDir, int noDir, int noAct, int bitDepth, int subBlkSize, int mappingDir[NUM_DIR_FIX][NUM_DIR_FIX], uint32_t **laplacian[NUM_DIRECTIONS])
{
  int shift = 9 + bitDepth;
  const static int multTab[] = { 5628, 1407, 624, 351, 225, 156 };
  int mult = multTab[dirWindSize];

  for( int i = 0; i < curBlk.height; i += subBlkSize)
  {
    for( int j = 0; j < curBlk.width; j += subBlkSize)
    {
      int sum = 0, sumV = 0, sumH = 0, sumD0 = 0, sumD1 = 0;
      if( dirWindSize == 1 )
      {
        int yOffset = (i >> 1) + 2;
        int xOffset = (j >> 1) + 2;
        sumV = laplacian[VER][yOffset][xOffset];
        sumH = laplacian[HOR][yOffset][xOffset];
        sumD0 = laplacian[DIAG0][yOffset][xOffset];
        sumD1 = laplacian[DIAG1][yOffset][xOffset];
      }
      else if (dirWindSize == 5)
      {
        int yOffset = i >> 1;
        int xOffset = j >> 1;
        sumV = laplacian[VER][yOffset][xOffset];
        sumH = laplacian[HOR][yOffset][xOffset];
        sumD0 = laplacian[DIAG0][yOffset][xOffset];
        sumD1 = laplacian[DIAG1][yOffset][xOffset];
      }
      else
      {
        CHECK(1, "other window sizes are not supported");
      }
      sum += sumV + sumH;
      int direction = 0, mainDirection, secondaryDirection, dirTempHV = 0, dirTempD = 0, hv1, hv0, d1, d0, hvd1, hvd0;
      if (sumV > sumH)
      {
        hv1 = sumV;
        hv0 = sumH;
        dirTempHV = 1;
      }
      else
      {
        hv1 = sumH;
        hv0 = sumV;
        dirTempHV = 3;
      }
      if (sumD0 > sumD1)
      {
        d1 = sumD0;
        d0 = sumD1;
        dirTempD = 0;
      }
      else
      {
        d1 = sumD1;
        d0 = sumD0;
        dirTempD = 2;
      }

      if ( ( uint64_t )d1 * hv0 > ( uint64_t )d0 * hv1 )
      {
        mainDirection = dirTempD;
        secondaryDirection = dirTempHV;
        hvd1 = d1;
        hvd0 = d0;
      }
      else
      {
        mainDirection = dirTempHV;
        secondaryDirection = dirTempD;
        hvd1 = hv1;
        hvd0 = hv0;
      }
      //activity
      int activity = 0;
      if (noAct == 5)
      {
        const int th[] = { 0, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4 };
        activity = th[std::min((mult * sum) >> shift, 15)];
      }
      else
      {
        const int th[193] = { 0,  1,  2,  3,  4,  4,  5,  5,  6,  6,  6,  6,  7,  7,  7,  7,  8,  8,  8,  8,  8,  8,  8,  8,  9,  9,  9,  9,  9,  9,  9,  9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 15 };
        activity = th[std::min((mult * sum) >> (shift - 4), 192)];
      }

      //direction
      if (classDir == 0)
      {
        int directionStrength = 0;
        if (hvd1 * 2 > 9 * hvd0)
        {
          directionStrength = 2;
        }
        else if (hvd1 > 2 * hvd0)
        {
          directionStrength = 1;
        }

        if (directionStrength)
        {
          direction = ((mainDirection & 0x1) << 1) + directionStrength;
        }
      }
      else
      {
        int edgeStrengthHV = 0, edgeStrengthD = 0;
        if (hv1 > 8 * hv0)
        {
          edgeStrengthHV = 6;
        }
        else if (hv1 * 2 > 9 * hv0)
        {
          edgeStrengthHV = 5;
        }
        else if (hv1 > 3 * hv0)
        {
          edgeStrengthHV = 4;
        }
        else if (hv1 > 2 * hv0)
        {
          edgeStrengthHV = 3;
        }
        else if (hv1 * 2 > 3 * hv0)
        {
          edgeStrengthHV = 2;
        }
        else if (hv1 * 4 > 5 * hv0)
        {
          edgeStrengthHV = 1;
        }

        if (d1 > 8 * d0)
        {
          edgeStrengthD = 6;
        }
        else if (d1 * 2 > 9 * d0)
        {
          edgeStrengthD = 5;
        }
        else if (d1 > 3 * d0)
        {
          edgeStrengthD = 4;
        }
        else if (d1 > 2 * d0)
        {
          edgeStrengthD = 3;
        }
        else if (d1 * 2 > 3 * d0)
        {
          edgeStrengthD = 2;
        }
        else if (d1 * 4 > 5 * d0)
        {
          edgeStrengthD = 1;
        }

        if ( ( uint64_t )hv1 * d0 >  ( uint64_t )hv0 * d1)
        {
          direction = mappingDir[edgeStrengthHV][edgeStrengthD];
        }
        else
        {
          direction = 28 + mappingDir[edgeStrengthD][edgeStrengthHV];
        }
      }

      int actDirInd = 0;
      if ( classDir == 0 )
      {
        actDirInd = noDir * activity + direction;
      }
      else
      {
        actDirInd = noDir * ( noDir + 1 ) * activity + direction;
      }

      static const int transposeTable[8] = { 0, 1, 0, 2, 2, 3, 1, 3 };
      int transposeIdx = transposeTable[mainDirection * 2 + (secondaryDirection >> 1)];

      for( int ii = curBlk.y + i; ii < curBlk.y + i + subBlkSize; ii++ )
      {
        for( int jj = curBlk.x + j; jj < curBlk.x + j + subBlkSize; jj++ )
        {
          classifier[ii][jj] = (actDirInd << 2) + transposeIdx;
        }
      }
    }
  }
}

void AdaptiveLoopFilter::fixedFiltering( AlfClassifier **classifier, const CPelBuf &srcLuma, const Area& curBlk, Pel ***fixedFilterResults, int picWidth, const int fixedFiltInd, const short classIndFixed[NUM_CLASSES_FIX], int fixedFiltQpInd, int dirWindSize, const ClpRng &clpRng, const Pel clippingValues[4])
{
  const int shift = m_NUM_BITS_FIXED_FILTER - 1;
  const int offset = 1 << (shift - 1);
  const int srcStride = srcLuma.stride;
  int posY = curBlk.y;
  int posX = curBlk.x;
  const Pel *pImgYPad0 = srcLuma.buf + posY * srcStride + posX;
  const Pel *pImg0, *pImg1, *pImg2, *pImg3, *pImg4, *pImg5, *pImg6, *pImg7, *pImg8, *pImg9, *pImg10, *pImg11, *pImg12;
  const int clsSizeY = 2;
  const int clsSizeX = 2;
  const int srcStride2 = srcStride * clsSizeY;

  for( int i = 0; i < curBlk.height; i+= clsSizeY)
  {
    for (int j = 0; j < curBlk.width; j += clsSizeX)
    {
      int classIdx = classifier[posY + i][posX + j] >> 2;
      int transposeIdx = classifier[posY + i][posX + j] & 0x3;
      int filterIdx = classIndFixed[classIdx];

      std::array<short, FIX_FILTER_NUM_COEFF> filterCoeff;
      std::array<short, FIX_FILTER_NUM_COEFF> filterClipp;

      const short* coeff = m_filterCoeffFixed[fixedFiltQpInd][dirWindSize][filterIdx];
      const short* clipp = m_clippingFixed[fixedFiltQpInd][dirWindSize][filterIdx];

      if (transposeIdx == 0)
      {
        filterCoeff = { coeff[0], coeff[1], coeff[2], coeff[3], coeff[4], coeff[5], coeff[6], coeff[7], coeff[8], coeff[9], coeff[10], coeff[11], coeff[12], coeff[13], coeff[14], coeff[15], coeff[16], coeff[17], coeff[18], coeff[19], coeff[20], coeff[21], coeff[22], coeff[23], coeff[24], coeff[25], coeff[26], coeff[27], coeff[28], coeff[29], coeff[30], coeff[31], coeff[32], coeff[33], coeff[34], coeff[35], coeff[36], coeff[37], coeff[38], coeff[39], coeff[40], coeff[41] };
        filterClipp = { clipp[0], clipp[1], clipp[2], clipp[3], clipp[4], clipp[5], clipp[6], clipp[7], clipp[8], clipp[9], clipp[10], clipp[11], clipp[12], clipp[13], clipp[14], clipp[15], clipp[16], clipp[17], clipp[18], clipp[19], clipp[20], clipp[21], clipp[22], clipp[23], clipp[24], clipp[25], clipp[26], clipp[27], clipp[28], clipp[29], clipp[30], clipp[31], clipp[32], clipp[33], clipp[34], clipp[35], clipp[36], clipp[37], clipp[38], clipp[39], clipp[40], clipp[41] };
      }
      else if (transposeIdx == 1)
      {
        filterCoeff = { coeff[36], coeff[25], coeff[37], coeff[35], coeff[16], coeff[26], coeff[38], coeff[34], coeff[24], coeff[9], coeff[17], coeff[27], coeff[39], coeff[33], coeff[23], coeff[15], coeff[4], coeff[10], coeff[18], coeff[28], coeff[40], coeff[32], coeff[22], coeff[14], coeff[8], coeff[1], coeff[5], coeff[11], coeff[19], coeff[29], coeff[41], coeff[31], coeff[21], coeff[13], coeff[7], coeff[3], coeff[0], coeff[2], coeff[6], coeff[12], coeff[20], coeff[30] };
        filterClipp = { clipp[36], clipp[25], clipp[37], clipp[35], clipp[16], clipp[26], clipp[38], clipp[34], clipp[24], clipp[9], clipp[17], clipp[27], clipp[39], clipp[33], clipp[23], clipp[15], clipp[4], clipp[10], clipp[18], clipp[28], clipp[40], clipp[32], clipp[22], clipp[14], clipp[8], clipp[1], clipp[5], clipp[11], clipp[19], clipp[29], clipp[41], clipp[31], clipp[21], clipp[13], clipp[7], clipp[3], clipp[0], clipp[2], clipp[6], clipp[12], clipp[20], clipp[30] };
      }
      else if (transposeIdx == 2)
      {
        filterCoeff = { coeff[0], coeff[3], coeff[2], coeff[1], coeff[8], coeff[7], coeff[6], coeff[5], coeff[4], coeff[15], coeff[14], coeff[13], coeff[12], coeff[11], coeff[10], coeff[9], coeff[24], coeff[23], coeff[22], coeff[21], coeff[20], coeff[19], coeff[18], coeff[17], coeff[16], coeff[35], coeff[34], coeff[33], coeff[32], coeff[31], coeff[30], coeff[29], coeff[28], coeff[27], coeff[26], coeff[25], coeff[36], coeff[37], coeff[38], coeff[39], coeff[40], coeff[41] };
        filterClipp = { clipp[0], clipp[3], clipp[2], clipp[1], clipp[8], clipp[7], clipp[6], clipp[5], clipp[4], clipp[15], clipp[14], clipp[13], clipp[12], clipp[11], clipp[10], clipp[9], clipp[24], clipp[23], clipp[22], clipp[21], clipp[20], clipp[19], clipp[18], clipp[17], clipp[16], clipp[35], clipp[34], clipp[33], clipp[32], clipp[31], clipp[30], clipp[29], clipp[28], clipp[27], clipp[26], clipp[25], clipp[36], clipp[37], clipp[38], clipp[39], clipp[40], clipp[41] };
      }
      else
      {
        filterCoeff = { coeff[36], coeff[35], coeff[37], coeff[25], coeff[24], coeff[34], coeff[38], coeff[26], coeff[16], coeff[15], coeff[23], coeff[33], coeff[39], coeff[27], coeff[17], coeff[9], coeff[8], coeff[14], coeff[22], coeff[32], coeff[40], coeff[28], coeff[18], coeff[10], coeff[4], coeff[3], coeff[7], coeff[13], coeff[21], coeff[31], coeff[41], coeff[29], coeff[19], coeff[11], coeff[5], coeff[1], coeff[0], coeff[2], coeff[6], coeff[12], coeff[20], coeff[30] };
        filterClipp = { clipp[36], clipp[35], clipp[37], clipp[25], clipp[24], clipp[34], clipp[38], clipp[26], clipp[16], clipp[15], clipp[23], clipp[33], clipp[39], clipp[27], clipp[17], clipp[9], clipp[8], clipp[14], clipp[22], clipp[32], clipp[40], clipp[28], clipp[18], clipp[10], clipp[4], clipp[3], clipp[7], clipp[13], clipp[21], clipp[31], clipp[41], clipp[29], clipp[19], clipp[11], clipp[5], clipp[1], clipp[0], clipp[2], clipp[6], clipp[12], clipp[20], clipp[30] };
      }

      for (int k = 0; k < 42; k++)
      {
        filterClipp[k] = clippingValues[filterClipp[k]];
      }

      for (int ii = 0; ii < clsSizeY; ii++)
      {
        pImg0 = pImgYPad0 + ii * srcStride + j;
        for (int jj = 0; jj < clsSizeX; jj++)
        {
          pImg1 = pImg0 - srcStride;
          pImg3 = pImg1 - srcStride;
          pImg5 = pImg3 - srcStride;
          pImg7 = pImg5 - srcStride;
          pImg9 = pImg7 - srcStride;
          pImg11 = pImg9 - srcStride;
          pImg2 = pImg0 + srcStride;
          pImg4 = pImg2 + srcStride;
          pImg6 = pImg4 + srcStride;
          pImg8 = pImg6 + srcStride;
          pImg10 = pImg8 + srcStride;
          pImg12 = pImg10 + srcStride;
          Pel curr = pImg0[0];

          int sum = 0;
          sum += filterCoeff[0] * (clipALF(filterClipp[0], curr, pImg11[+0], pImg12[+0]));

          sum += filterCoeff[1] * (clipALF(filterClipp[1], curr, pImg9[-1], pImg10[+1]));
          sum += filterCoeff[2] * (clipALF(filterClipp[2], curr, pImg9[+0], pImg10[+0]));
          sum += filterCoeff[3] * (clipALF(filterClipp[3], curr, pImg9[+1], pImg10[-1]));

          sum += filterCoeff[4] * (clipALF(filterClipp[4], curr, pImg7[-2], pImg8[+2]));
          sum += filterCoeff[5] * (clipALF(filterClipp[5], curr, pImg7[-1], pImg8[+1]));
          sum += filterCoeff[6] * (clipALF(filterClipp[6], curr, pImg7[+0], pImg8[+0]));
          sum += filterCoeff[7] * (clipALF(filterClipp[7], curr, pImg7[+1], pImg8[-1]));
          sum += filterCoeff[8] * (clipALF(filterClipp[8], curr, pImg7[+2], pImg8[-2]));

          sum += filterCoeff[9] * (clipALF(filterClipp[9], curr, pImg5[-3], pImg6[+3]));
          sum += filterCoeff[10] * (clipALF(filterClipp[10], curr, pImg5[-2], pImg6[+2]));
          sum += filterCoeff[11] * (clipALF(filterClipp[11], curr, pImg5[-1], pImg6[+1]));
          sum += filterCoeff[12] * (clipALF(filterClipp[12], curr, pImg5[+0], pImg6[+0]));
          sum += filterCoeff[13] * (clipALF(filterClipp[13], curr, pImg5[+1], pImg6[-1]));
          sum += filterCoeff[14] * (clipALF(filterClipp[14], curr, pImg5[+2], pImg6[-2]));
          sum += filterCoeff[15] * (clipALF(filterClipp[15], curr, pImg5[+3], pImg6[-3]));

          sum += filterCoeff[16] * (clipALF(filterClipp[16], curr, pImg3[-4], pImg4[+4]));
          sum += filterCoeff[17] * (clipALF(filterClipp[17], curr, pImg3[-3], pImg4[+3]));
          sum += filterCoeff[18] * (clipALF(filterClipp[18], curr, pImg3[-2], pImg4[+2]));
          sum += filterCoeff[19] * (clipALF(filterClipp[19], curr, pImg3[-1], pImg4[+1]));
          sum += filterCoeff[20] * (clipALF(filterClipp[20], curr, pImg3[+0], pImg4[+0]));
          sum += filterCoeff[21] * (clipALF(filterClipp[21], curr, pImg3[+1], pImg4[-1]));
          sum += filterCoeff[22] * (clipALF(filterClipp[22], curr, pImg3[+2], pImg4[-2]));
          sum += filterCoeff[23] * (clipALF(filterClipp[23], curr, pImg3[+3], pImg4[-3]));
          sum += filterCoeff[24] * (clipALF(filterClipp[24], curr, pImg3[+4], pImg4[-4]));

          sum += filterCoeff[25] * (clipALF(filterClipp[25], curr, pImg1[-5], pImg2[+5]));
          sum += filterCoeff[26] * (clipALF(filterClipp[26], curr, pImg1[-4], pImg2[+4]));
          sum += filterCoeff[27] * (clipALF(filterClipp[27], curr, pImg1[-3], pImg2[+3]));
          sum += filterCoeff[28] * (clipALF(filterClipp[28], curr, pImg1[-2], pImg2[+2]));
          sum += filterCoeff[29] * (clipALF(filterClipp[29], curr, pImg1[-1], pImg2[+1]));
          sum += filterCoeff[30] * (clipALF(filterClipp[30], curr, pImg1[+0], pImg2[+0]));
          sum += filterCoeff[31] * (clipALF(filterClipp[31], curr, pImg1[+1], pImg2[-1]));
          sum += filterCoeff[32] * (clipALF(filterClipp[32], curr, pImg1[+2], pImg2[-2]));
          sum += filterCoeff[33] * (clipALF(filterClipp[33], curr, pImg1[+3], pImg2[-3]));
          sum += filterCoeff[34] * (clipALF(filterClipp[34], curr, pImg1[+4], pImg2[-4]));
          sum += filterCoeff[35] * (clipALF(filterClipp[35], curr, pImg1[+5], pImg2[-5]));

          sum += filterCoeff[36] * (clipALF(filterClipp[36], curr, pImg0[-6], pImg0[+6]));
          sum += filterCoeff[37] * (clipALF(filterClipp[37], curr, pImg0[-5], pImg0[+5]));
          sum += filterCoeff[38] * (clipALF(filterClipp[38], curr, pImg0[-4], pImg0[+4]));
          sum += filterCoeff[39] * (clipALF(filterClipp[39], curr, pImg0[-3], pImg0[+3]));
          sum += filterCoeff[40] * (clipALF(filterClipp[40], curr, pImg0[-2], pImg0[+2]));
          sum += filterCoeff[41] * (clipALF(filterClipp[41], curr, pImg0[-1], pImg0[+1]));

          sum = (sum + offset) >> shift;
          sum += curr;

          fixedFilterResults[fixedFiltInd][posY + ii+ i][posX + jj + j] = ClipPel(sum, clpRng);

          pImg0++;
        }//jj
      }//ii
    }//j
    pImgYPad0 += srcStride2;
  }
}

void AdaptiveLoopFilter::deriveClassificationLaplacianBig(const Area &curBlk, uint32_t **laplacian[NUM_DIRECTIONS])
{
  int fl2 = ALF_CLASSIFIER_FL << 1;
  for (int i = 0; i < (curBlk.height + fl2) >> 1; i++)
  {
    for (int j = 4; j <= (curBlk.width + fl2 - 2) >> 1; j++)
    {
      int jM2 = j - 2;
      int jM4 = j - 4;
      laplacian[VER][i][jM4] += laplacian[VER][i][jM2] + laplacian[VER][i][j];
      laplacian[HOR][i][jM4] += laplacian[HOR][i][jM2] + laplacian[HOR][i][j];
      laplacian[DIAG0][i][jM4] += laplacian[DIAG0][i][jM2] + laplacian[DIAG0][i][j];
      laplacian[DIAG1][i][jM4] += laplacian[DIAG1][i][jM2] + laplacian[DIAG1][i][j];
    }
    if (i >= 4)
    {
      int iM4 = i - 4;
      int iM2 = i - 2;
      for (int j = 0; j <= (curBlk.width + fl2 - 11) >> 1; j++)
      {
        laplacian[VER][iM4][j] += laplacian[VER][iM2][j] + laplacian[VER][i][j];
        laplacian[HOR][iM4][j] += laplacian[HOR][iM2][j] + laplacian[HOR][i][j];
        laplacian[DIAG0][iM4][j] += laplacian[DIAG0][iM2][j] + laplacian[DIAG0][i][j];
        laplacian[DIAG1][iM4][j] += laplacian[DIAG1][iM2][j] + laplacian[DIAG1][i][j];
      }
    }
  }
}

void AdaptiveLoopFilter::deriveClassificationLaplacian(const CPelBuf &srcLuma, const Area &blkDst, const Area &blk, uint32_t **laplacian[NUM_DIRECTIONS])
{
  int fl = ALF_CLASSIFIER_FL;
  int fl2 = 2 * fl;
  const Pel* src = srcLuma.buf;
  int stride = srcLuma.stride;
  int yOffset = (blk.pos().y - fl) * stride;
  const Pel *src0 = &src[yOffset - stride];
  const Pel *src1 = &src[yOffset];
  const Pel *src2 = &src[yOffset + stride];
  int stride2 = 2 * stride;
  const Pel *src3 = &src[yOffset + stride2];
  for (int i = 0; i < blk.height + fl2; i += 2)
  {
    for (int j = 0; j < blk.width + fl2; j += 2)
    {
      int xOffset = blk.pos().x - fl + j;
      const Pel *pY = src1 + xOffset;
      const Pel *pYup = src0 + xOffset;
      const Pel *pYdn = src2 + xOffset;
      const Pel *pYdn2 = src3 + xOffset;
      Pel y0 = pY[0] << 1;
      Pel y01 = pY[1] << 1;
      Pel y1 = pYdn[0] << 1;
      Pel y11 = pYdn[1] << 1;
      laplacian[VER][i >> 1][j >> 1] = abs(y0 - pYup[0] - pYdn[0]) + abs(y01 - pYup[1] - pYdn[1]) + abs(y1 - pY[0] - pYdn2[0]) + abs(y11 - pY[1] - pYdn2[1]);
      laplacian[HOR][i >> 1][j >> 1] = abs(y0 - pY[-1] - pY[1]) + abs(y01 - pY[0] - pY[2]) + abs(y1 - pYdn[-1] - pYdn[1]) + abs(y11 - pYdn[0] - pYdn[2]);
      laplacian[DIAG0][i >> 1][j >> 1] = abs(y0 - pYup[-1] - pYdn[1]) + abs(y01 - pYup[0] - pYdn[2]) + abs(y1 - pY[-1] - pYdn2[1]) + abs(y11 - pY[0] - pYdn2[2]);
      laplacian[DIAG1][i >> 1][j >> 1] = abs(y0 - pYup[1] - pYdn[-1]) + abs(y01 - pYup[2] - pYdn[0]) + abs(y1 - pY[1] - pYdn2[-1]) + abs(y11 - pY[2] - pYdn2[0]);
    }
    src0 = src0 + stride2;
    src1 = src1 + stride2;
    src2 = src2 + stride2;
    src3 = src3 + stride2;
  }

  for (int i = 0; i < (blk.height + fl2) >> 1; i++)
  {
    for (int j = 1; j < (blk.width + fl2) >> 1; j++)
    {
      int jM1 = j - 1;
      laplacian[VER][i][jM1] = laplacian[VER][i][jM1] + laplacian[VER][i][j];
      laplacian[HOR][i][jM1] = laplacian[HOR][i][jM1] + laplacian[HOR][i][j];
      laplacian[DIAG0][i][jM1] = laplacian[DIAG0][i][jM1] + laplacian[DIAG0][i][j];
      laplacian[DIAG1][i][jM1] = laplacian[DIAG1][i][jM1] + laplacian[DIAG1][i][j];
    }
    if (i > 0)
    {
      int iM1 = i - 1;
      for (int j = 0; j < (blk.width + fl2) >> 1; j++)
      {
        laplacian[VER][iM1][j] = laplacian[VER][iM1][j] + laplacian[VER][i][j];
        laplacian[HOR][iM1][j] = laplacian[HOR][iM1][j] + laplacian[HOR][i][j];
        laplacian[DIAG0][iM1][j] = laplacian[DIAG0][iM1][j] + laplacian[DIAG0][i][j];
        laplacian[DIAG1][iM1][j] = laplacian[DIAG1][iM1][j] + laplacian[DIAG1][i][j];
      }
    }
  }
}

void AdaptiveLoopFilter::deriveClassificationAndFixFilterResultsBlk( AlfClassifier **classifier, Pel ***fixedFilterResults, const CPelBuf &srcLuma, const Area &blkDst, const Area &blk, const int bits, CodingStructure& cs, const ClpRng &clpRng, const Pel clippingValues[4], int qp, int fixedFilterSetIdx, int mappingDir[NUM_DIR_FIX][NUM_DIR_FIX], uint32_t **laplacian[NUM_DIRECTIONS], const int classifierIdx )
{
  m_deriveClassificationLaplacian(srcLuma, blkDst, blk, laplacian);
  
  int fixedFiltSetIndMin = (qp - 22) / 4;
  int fixedFiltSetIndMax = (qp - 18) / 4;
  if (fixedFiltSetIndMin <= 0)
  {
    fixedFiltSetIndMin = 0;
    fixedFiltSetIndMax = 1;
  }

  if (fixedFiltSetIndMax > (NUM_SETS_FIXED_FILTERS - 1))
  {
    fixedFiltSetIndMin = NUM_SETS_FIXED_FILTERS - 2;
    fixedFiltSetIndMax = NUM_SETS_FIXED_FILTERS - 1;
  }

  int targetFixedFilterSetInd = fixedFilterSetIdx == 0 ? fixedFiltSetIndMin : (fixedFilterSetIdx == 1 ? fixedFiltSetIndMax : -1);
      
  int fixedFiltInd = 0;
  for (int dirWindSize = 0; dirWindSize < NUM_CLASSIFIER; dirWindSize++)
  {
    bool reuse = false;
    for (int fixedFiltSetInd = fixedFiltSetIndMin; fixedFiltSetInd <= fixedFiltSetIndMax; fixedFiltSetInd++)
    {
      if (dirWindSize == classifierIdx || classifierIdx == -1)
      {
        if (fixedFiltSetInd == targetFixedFilterSetInd || targetFixedFilterSetInd == -1)
        {
          if (!reuse)
          {
            if (dirWindSize > 0)
            {
              m_deriveClassificationLaplacianBig(blkDst, laplacian);
            }
            m_calcClass0(classifier, blkDst, blkDst, usedWindowIdx[dirWindSize], 1, NUM_DIR_FIX, NUM_ACT_FIX, bits, 2, mappingDir, laplacian);

            reuse = true;
          }
          //fixed filtering
          m_filter13x13Blk(classifier, srcLuma, blkDst, fixedFilterResults, m_picWidth, fixedFiltInd, m_classIdnFixedFilter[fixedFiltSetInd][dirWindSize], fixedFiltSetInd, dirWindSize, clpRng, clippingValues);
        }
      }
      fixedFiltInd++;
    }
  }
  if( classifierIdx == -1 )
  {
    m_calcClass1( classifier, blkDst, Area(blkDst.pos().x, blkDst.pos().y, blkDst.width, blkDst.height), 5, 0, 5, 5, bits, 2, mappingDir, laplacian );
  }
}
#else
void AdaptiveLoopFilter::deriveClassificationBlk(AlfClassifier **classifier, int **laplacian[NUM_DIRECTIONS], const CPelBuf &srcLuma, const Area &blkDst, const Area &blk, const int shift, const int vbCTUHeight, int vbPos)
{
  CHECK((vbCTUHeight & (vbCTUHeight - 1)) != 0, "vbCTUHeight must be a power of 2");
  static const int th[16] = { 0, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4 };
  const int stride = srcLuma.stride;
  const Pel* src = srcLuma.buf;
  const int maxActivity = 15;

  int fl = 2;
  int flP1 = fl + 1;
  int fl2 = 2 * fl;

  int mainDirection, secondaryDirection, dirTempHV, dirTempD;

  int pixY;
  int height = blk.height + fl2;
  int width = blk.width + fl2;
  int posX = blk.pos().x;
  int posY = blk.pos().y;
  int startHeight = posY - flP1;

  for( int i = 0; i < height; i += 2 )
  {
    int yoffset = ( i + 1 + startHeight ) * stride - flP1;
    const Pel *src0 = &src[yoffset - stride];
    const Pel *src1 = &src[yoffset];
    const Pel *src2 = &src[yoffset + stride];
    const Pel *src3 = &src[yoffset + stride * 2];

    const int y = blkDst.pos().y - 2 + i;
    if (y > 0 && (y & (vbCTUHeight - 1)) == vbPos - 2)
    {
      src3 = &src[yoffset + stride];
    }
    else if (y > 0 && (y & (vbCTUHeight - 1)) == vbPos)
    {
      src0 = &src[yoffset];
    }
    int* pYver = laplacian[VER][i];
    int* pYhor = laplacian[HOR][i];
    int* pYdig0 = laplacian[DIAG0][i];
    int* pYdig1 = laplacian[DIAG1][i];

    for( int j = 0; j < width; j += 2 )
    {
      pixY = j + 1 + posX;
      const Pel *pY = src1 + pixY;
      const Pel* pYdown = src0 + pixY;
      const Pel* pYup = src2 + pixY;
      const Pel* pYup2 = src3 + pixY;

      const Pel y0 = pY[0] << 1;
      const Pel yup1 = pYup[1] << 1;

      pYver[j] = abs( y0 - pYdown[0] - pYup[0] ) + abs( yup1 - pY[1] - pYup2[1] );
      pYhor[j] = abs( y0 - pY[1] - pY[-1] ) + abs( yup1 - pYup[2] - pYup[0] );
      pYdig0[j] = abs( y0 - pYdown[-1] - pYup[1] ) + abs( yup1 - pY[0] - pYup2[2] );
      pYdig1[j] = abs( y0 - pYup[-1] - pYdown[1] ) + abs( yup1 - pYup2[0] - pY[2] );

      if( j > 4 && ( j - 6 ) % 4 == 0 )
      {
        int jM6 = j - 6;
        int jM4 = j - 4;
        int jM2 = j - 2;
        pYver[jM6] += pYver[jM4] + pYver[jM2] + pYver[j];
        pYhor[jM6] += pYhor[jM4] + pYhor[jM2] + pYhor[j];
        pYdig0[jM6] += pYdig0[jM4] + pYdig0[jM2] + pYdig0[j];
        pYdig1[jM6] += pYdig1[jM4] + pYdig1[jM2] + pYdig1[j];
      }
    }
  }
  // classification block size
  const int clsSizeY = 4;
  const int clsSizeX = 4;

  for( int i = 0; i < blk.height; i += clsSizeY )
  {
    int* pYver = laplacian[VER][i];
    int* pYver2 = laplacian[VER][i + 2];
    int* pYver4 = laplacian[VER][i + 4];
    int* pYver6 = laplacian[VER][i + 6];

    int* pYhor = laplacian[HOR][i];
    int* pYhor2 = laplacian[HOR][i + 2];
    int* pYhor4 = laplacian[HOR][i + 4];
    int* pYhor6 = laplacian[HOR][i + 6];

    int* pYdig0 = laplacian[DIAG0][i];
    int* pYdig02 = laplacian[DIAG0][i + 2];
    int* pYdig04 = laplacian[DIAG0][i + 4];
    int* pYdig06 = laplacian[DIAG0][i + 6];

    int* pYdig1 = laplacian[DIAG1][i];
    int* pYdig12 = laplacian[DIAG1][i + 2];
    int* pYdig14 = laplacian[DIAG1][i + 4];
    int* pYdig16 = laplacian[DIAG1][i + 6];

    for( int j = 0; j < blk.width; j += clsSizeX )
    {
      int sumV = 0; int sumH = 0; int sumD0 = 0; int sumD1 = 0;

      if (((i + blkDst.pos().y) % vbCTUHeight) == (vbPos - 4))
      {
        sumV = pYver[j] + pYver2[j] + pYver4[j];
        sumH = pYhor[j] + pYhor2[j] + pYhor4[j];
        sumD0 = pYdig0[j] + pYdig02[j] + pYdig04[j];
        sumD1 = pYdig1[j] + pYdig12[j] + pYdig14[j];
      }
      else if (((i + blkDst.pos().y) % vbCTUHeight) == vbPos)
      {
        sumV = pYver2[j] + pYver4[j] + pYver6[j];
        sumH = pYhor2[j] + pYhor4[j] + pYhor6[j];
        sumD0 = pYdig02[j] + pYdig04[j] + pYdig06[j];
        sumD1 = pYdig12[j] + pYdig14[j] + pYdig16[j];
      }
      else
      {
        sumV = pYver[j] + pYver2[j] + pYver4[j] + pYver6[j];
        sumH = pYhor[j] + pYhor2[j] + pYhor4[j] + pYhor6[j];
        sumD0 = pYdig0[j] + pYdig02[j] + pYdig04[j] + pYdig06[j];
        sumD1 = pYdig1[j] + pYdig12[j] + pYdig14[j] + pYdig16[j];
      }

      int tempAct = sumV + sumH;
      int activity = 0;

      const int y = (i + blkDst.pos().y) & (vbCTUHeight - 1);
      if (y == vbPos - 4 || y == vbPos)
      {
        activity = (Pel)Clip3<int>(0, maxActivity, (tempAct * 96) >> shift);
      }
      else
      {
        activity = (Pel)Clip3<int>(0, maxActivity, (tempAct * 64) >> shift);
      }

      int classIdx = th[activity];

      int hv1, hv0, d1, d0, hvd1, hvd0;

      if( sumV > sumH )
      {
        hv1 = sumV;
        hv0 = sumH;
        dirTempHV = 1;
      }
      else
      {
        hv1 = sumH;
        hv0 = sumV;
        dirTempHV = 3;
      }
      if( sumD0 > sumD1 )
      {
        d1 = sumD0;
        d0 = sumD1;
        dirTempD = 0;
      }
      else
      {
        d1 = sumD1;
        d0 = sumD0;
        dirTempD = 2;
      }
      if( (uint32_t)d1 * (uint32_t)hv0 > (uint32_t)hv1 * (uint32_t)d0 )
      {
        hvd1 = d1;
        hvd0 = d0;
        mainDirection = dirTempD;
        secondaryDirection = dirTempHV;
      }
      else
      {
        hvd1 = hv1;
        hvd0 = hv0;
        mainDirection = dirTempHV;
        secondaryDirection = dirTempD;
      }

      int directionStrength = 0;
      if( hvd1 > 2 * hvd0 )
      {
        directionStrength = 1;
      }
      if( hvd1 * 2 > 9 * hvd0 )
      {
        directionStrength = 2;
      }

      if( directionStrength )
      {
        classIdx += ( ( ( mainDirection & 0x1 ) << 1 ) + directionStrength ) * 5;
      }

      static const int transposeTable[8] = { 0, 1, 0, 2, 2, 3, 1, 3 };
      int transposeIdx = transposeTable[mainDirection * 2 + ( secondaryDirection >> 1 )];
      int yOffset = i + blkDst.pos().y;
      int xOffset = j + blkDst.pos().x;

      AlfClassifier *cl0 = classifier[yOffset] + xOffset;
      AlfClassifier *cl1 = classifier[yOffset + 1] + xOffset;
      AlfClassifier *cl2 = classifier[yOffset + 2] + xOffset;
      AlfClassifier *cl3 = classifier[yOffset + 3] + xOffset;
      cl0[0] = cl0[1] = cl0[2] = cl0[3] = cl1[0] = cl1[1] = cl1[2] = cl1[3] = cl2[0] = cl2[1] = cl2[2] = cl2[3] = cl3[0] = cl3[1] = cl3[2] = cl3[3] = AlfClassifier( classIdx, transposeIdx );
    }
  }
}
#endif

template<AlfFilterType filtType>
void AdaptiveLoopFilter::filterBlk(AlfClassifier **classifier, const PelUnitBuf &recDst, const CPelUnitBuf &recSrc,
                                   const Area &blkDst, const Area &blk, const ComponentID compId,
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
                                   const short *filterSet, const Pel *fClipSet, const ClpRng &clpRng,
#else
                                   const short *filterSet, const short *fClipSet, const ClpRng &clpRng,
#endif
                                   CodingStructure &cs, 
#if ALF_IMPROVEMENT
                                   Pel ***fixedFilterResults, int fixedFilterSetIdx)
#else
                                   const int vbCTUHeight, int vbPos)
#endif
{
#if !ALF_IMPROVEMENT
  CHECK((vbCTUHeight & (vbCTUHeight - 1)) != 0, "vbCTUHeight must be a power of 2");
#endif

  const bool bChroma = isChroma( compId );
#if !ALF_IMPROVEMENT
  if( bChroma )
  {
    CHECK( filtType != 0, "Chroma needs to have filtType == 0" );
  }
#endif

  const CPelBuf srcLuma = recSrc.get( compId );
  PelBuf dstLuma = recDst.get( compId );

  const int srcStride = srcLuma.stride;
  const int dstStride = dstLuma.stride;

  const int startHeight = blk.y;
  const int endHeight = blk.y + blk.height;
  const int startWidth = blk.x;
  const int endWidth = blk.x + blk.width;

  const Pel* src = srcLuma.buf;
  Pel* dst = dstLuma.buf + blkDst.y * dstStride;

#if ALF_IMPROVEMENT
  const Pel *pImgYPad0, *pImgYPad1, *pImgYPad2, *pImgYPad3, *pImgYPad4, *pImgYPad5, *pImgYPad6, *pImgYPad7, *pImgYPad8;
  const Pel *pImg0, *pImg1, *pImg2, *pImg3, *pImg4, *pImg5, *pImg6, *pImg7, *pImg8;
#else
  const Pel *pImgYPad0, *pImgYPad1, *pImgYPad2, *pImgYPad3, *pImgYPad4, *pImgYPad5, *pImgYPad6;
  const Pel *pImg0, *pImg1, *pImg2, *pImg3, *pImg4, *pImg5, *pImg6;
#endif

  const short *coef = filterSet;
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  const Pel *clip = fClipSet;
#else
  const short *clip = fClipSet;
#endif
  const int shift = m_NUM_BITS - 1;

  const int offset = 1 << ( shift - 1 );

  int transposeIdx = 0;
#if ALF_IMPROVEMENT
  const int clsSizeY = 2;
  const int clsSizeX = 2;
#else
  const int clsSizeY = 4;
  const int clsSizeX = 4;
#endif


  CHECK( startHeight % clsSizeY, "Wrong startHeight in filtering" );
  CHECK( startWidth % clsSizeX, "Wrong startWidth in filtering" );
  CHECK( ( endHeight - startHeight ) % clsSizeY, "Wrong endHeight in filtering" );
  CHECK( ( endWidth - startWidth ) % clsSizeX, "Wrong endWidth in filtering" );

  AlfClassifier *pClass = nullptr;

  int dstStride2 = dstStride * clsSizeY;
  int srcStride2 = srcStride * clsSizeY;

  std::array<int, MAX_NUM_ALF_LUMA_COEFF> filterCoeff;
  std::array<int, MAX_NUM_ALF_LUMA_COEFF> filterClipp;

  pImgYPad0 = src + startHeight * srcStride + startWidth;
  pImgYPad1 = pImgYPad0 + srcStride;
  pImgYPad2 = pImgYPad0 - srcStride;
  pImgYPad3 = pImgYPad1 + srcStride;
  pImgYPad4 = pImgYPad2 - srcStride;
  pImgYPad5 = pImgYPad3 + srcStride;
  pImgYPad6 = pImgYPad4 - srcStride;
#if ALF_IMPROVEMENT
  pImgYPad7 = pImgYPad5 + srcStride;
  pImgYPad8 = pImgYPad6 - srcStride;
#endif

  Pel* pRec0 = dst + blkDst.x;
  Pel* pRec1 = pRec0 + dstStride;

  for( int i = 0; i < endHeight - startHeight; i += clsSizeY )
  {
    if( !bChroma )
    {
      pClass = classifier[blkDst.y + i] + blkDst.x;
    }

    for( int j = 0; j < endWidth - startWidth; j += clsSizeX )
    {
      if( !bChroma )
      {
        AlfClassifier& cl = pClass[j];
#if ALF_IMPROVEMENT
        transposeIdx = cl & 0x3;
        int classIdx = cl >> 2;
        coef = filterSet + classIdx * MAX_NUM_ALF_LUMA_COEFF;
        clip = fClipSet + classIdx * MAX_NUM_ALF_LUMA_COEFF;
#else
        transposeIdx = cl.transposeIdx;
        coef = filterSet + cl.classIdx * MAX_NUM_ALF_LUMA_COEFF;
        clip = fClipSet + cl.classIdx * MAX_NUM_ALF_LUMA_COEFF;
#endif
      }

#if ALF_IMPROVEMENT
      if( filtType == ALF_FILTER_9 || filtType == ALF_FILTER_9_EXT )
      {
        if( transposeIdx == 1 )
        {
          filterCoeff = { coef[16], coef[9], coef[17], coef[15], coef[4], coef[10], coef[18], coef[14], coef[8], coef[1], coef[5], coef[11], coef[19], coef[13], coef[7], coef[3], coef[0], coef[2], coef[6], coef[12], coef[20] };
          filterClipp = { clip[16], clip[9], clip[17], clip[15], clip[4], clip[10], clip[18], clip[14], clip[8], clip[1], clip[5], clip[11], clip[19], clip[13], clip[7], clip[3], clip[0], clip[2], clip[6], clip[12], clip[20] };
        }
        else if( transposeIdx == 2 )
        {
          filterCoeff = { coef[0], coef[3], coef[2], coef[1], coef[8], coef[7], coef[6], coef[5], coef[4], coef[15], coef[14] , coef[13], coef[12], coef[11], coef[10] , coef[9], coef[16], coef[17], coef[18], coef[19], coef[20] };
          filterClipp = { clip[0], clip[3], clip[2], clip[1], clip[8], clip[7], clip[6], clip[5], clip[4], clip[15], clip[14] , clip[13], clip[12], clip[11], clip[10] , clip[9], clip[16], clip[17], clip[18], clip[19], clip[20] };
        }
        else if( transposeIdx == 3 )
        {
          filterCoeff = { coef[16], coef[15], coef[17], coef[9], coef[8], coef[14], coef[18], coef[10], coef[4], coef[3], coef[7], coef[13], coef[19], coef[11], coef[5], coef[1], coef[0], coef[2], coef[6], coef[12], coef[20] };
          filterClipp = { clip[16], clip[15], clip[17], clip[9], clip[8], clip[14], clip[18], clip[10], clip[4], clip[3], clip[7], clip[13], clip[19], clip[11], clip[5], clip[1], clip[0], clip[2], clip[6], clip[12], clip[20] };
        }
        else
        {
          filterCoeff = { coef[0], coef[1], coef[2], coef[3], coef[4], coef[5], coef[6], coef[7], coef[8], coef[9], coef[10], coef[11], coef[12], coef[13], coef[14], coef[15], coef[16], coef[17], coef[18], coef[19], coef[20] };
          filterClipp = { clip[0], clip[1], clip[2], clip[3], clip[4], clip[5], clip[6], clip[7], clip[8], clip[9], clip[10], clip[11], clip[12], clip[13], clip[14], clip[15], clip[16], clip[17], clip[18], clip[19], clip[20] };
        }
      }
      else
#endif
      if( filtType == ALF_FILTER_7 )
      {
        if( transposeIdx == 1 )
        {
          filterCoeff = { coef[9], coef[4], coef[10], coef[8], coef[1], coef[5], coef[11], coef[7], coef[3], coef[0], coef[2], coef[6], coef[12] };
          filterClipp = { clip[9], clip[4], clip[10], clip[8], clip[1], clip[5], clip[11], clip[7], clip[3], clip[0], clip[2], clip[6], clip[12] };
        }
        else if( transposeIdx == 2 )
        {
          filterCoeff = { coef[0], coef[3], coef[2], coef[1], coef[8], coef[7], coef[6], coef[5], coef[4], coef[9], coef[10], coef[11], coef[12] };
          filterClipp = { clip[0], clip[3], clip[2], clip[1], clip[8], clip[7], clip[6], clip[5], clip[4], clip[9], clip[10], clip[11], clip[12] };
        }
        else if( transposeIdx == 3 )
        {
          filterCoeff = { coef[9], coef[8], coef[10], coef[4], coef[3], coef[7], coef[11], coef[5], coef[1], coef[0], coef[2], coef[6], coef[12] };
          filterClipp = { clip[9], clip[8], clip[10], clip[4], clip[3], clip[7], clip[11], clip[5], clip[1], clip[0], clip[2], clip[6], clip[12] };
        }
        else
        {
          filterCoeff = { coef[0], coef[1], coef[2], coef[3], coef[4], coef[5], coef[6], coef[7], coef[8], coef[9], coef[10], coef[11], coef[12] };
          filterClipp = { clip[0], clip[1], clip[2], clip[3], clip[4], clip[5], clip[6], clip[7], clip[8], clip[9], clip[10], clip[11], clip[12] };
        }
      }
      else
      {
        if( transposeIdx == 1 )
        {
          filterCoeff = { coef[4], coef[1], coef[5], coef[3], coef[0], coef[2], coef[6] };
          filterClipp = { clip[4], clip[1], clip[5], clip[3], clip[0], clip[2], clip[6] };
        }
        else if( transposeIdx == 2 )
        {
          filterCoeff = { coef[0], coef[3], coef[2], coef[1], coef[4], coef[5], coef[6] };
          filterClipp = { clip[0], clip[3], clip[2], clip[1], clip[4], clip[5], clip[6] };
        }
        else if( transposeIdx == 3 )
        {
          filterCoeff = { coef[4], coef[3], coef[5], coef[1], coef[0], coef[2], coef[6] };
          filterClipp = { clip[4], clip[3], clip[5], clip[1], clip[0], clip[2], clip[6] };
        }
        else
        {
          filterCoeff = { coef[0], coef[1], coef[2], coef[3], coef[4], coef[5], coef[6] };
          filterClipp = { clip[0], clip[1], clip[2], clip[3], clip[4], clip[5], clip[6] };
        }
      }

      for( int ii = 0; ii < clsSizeY; ii++ )
      {
        pImg0 = pImgYPad0 + j + ii * srcStride;
        pImg1 = pImgYPad1 + j + ii * srcStride;
        pImg2 = pImgYPad2 + j + ii * srcStride;
        pImg3 = pImgYPad3 + j + ii * srcStride;
        pImg4 = pImgYPad4 + j + ii * srcStride;
        pImg5 = pImgYPad5 + j + ii * srcStride;
        pImg6 = pImgYPad6 + j + ii * srcStride;
#if ALF_IMPROVEMENT
        pImg7 = pImgYPad7 + j + ii * srcStride;
        pImg8 = pImgYPad8 + j + ii * srcStride;
#endif

        pRec1 = pRec0 + j + ii * dstStride;
#if !ALF_IMPROVEMENT
        const int yVb = (blkDst.y + i + ii) & (vbCTUHeight - 1);
        if (yVb < vbPos && (yVb >= vbPos - (bChroma ? 2 : 4)))   // above
        {
          pImg1 = (yVb == vbPos - 1) ? pImg0 : pImg1;
          pImg3 = (yVb >= vbPos - 2) ? pImg1 : pImg3;
          pImg5 = (yVb >= vbPos - 3) ? pImg3 : pImg5;

          pImg2 = (yVb == vbPos - 1) ? pImg0 : pImg2;
          pImg4 = (yVb >= vbPos - 2) ? pImg2 : pImg4;
          pImg6 = (yVb >= vbPos - 3) ? pImg4 : pImg6;
        }
        else if (yVb >= vbPos && (yVb <= vbPos + (bChroma ? 1 : 3)))   // bottom
        {
          pImg2 = (yVb == vbPos) ? pImg0 : pImg2;
          pImg4 = (yVb <= vbPos + 1) ? pImg2 : pImg4;
          pImg6 = (yVb <= vbPos + 2) ? pImg4 : pImg6;

          pImg1 = (yVb == vbPos) ? pImg0 : pImg1;
          pImg3 = (yVb <= vbPos + 1) ? pImg1 : pImg3;
          pImg5 = (yVb <= vbPos + 2) ? pImg3 : pImg5;
        }
        bool isNearVBabove = yVb < vbPos && (yVb >= vbPos - 1);
        bool isNearVBbelow = yVb >= vbPos && (yVb <= vbPos);
#endif

        for( int jj = 0; jj < clsSizeX; jj++ )
        {
          int32_t   sum  = 0;
          const Pel curr = pImg0[+0];
#if ALF_IMPROVEMENT
          if( filtType == ALF_FILTER_9 || filtType == ALF_FILTER_9_EXT )
          {
            sum += filterCoeff[0] * ( clipALF(filterClipp[0], curr, pImg7[+0], pImg8[+0]) );

            sum += filterCoeff[1] * ( clipALF(filterClipp[1], curr, pImg5[+1], pImg6[-1]) );
            sum += filterCoeff[2] * ( clipALF(filterClipp[2], curr, pImg5[+0], pImg6[+0]) );
            sum += filterCoeff[3] * ( clipALF(filterClipp[3], curr, pImg5[-1], pImg6[+1]) );

            sum += filterCoeff[4] * ( clipALF(filterClipp[4], curr, pImg3[+2], pImg4[-2]) );
            sum += filterCoeff[5] * ( clipALF(filterClipp[5], curr, pImg3[+1], pImg4[-1]) );
            sum += filterCoeff[6] * ( clipALF(filterClipp[6], curr, pImg3[+0], pImg4[+0]) );
            sum += filterCoeff[7] * ( clipALF(filterClipp[7], curr, pImg3[-1], pImg4[+1]) );
            sum += filterCoeff[8] * ( clipALF(filterClipp[8], curr, pImg3[-2], pImg4[+2]) );

            sum += filterCoeff[9] * ( clipALF(filterClipp[9], curr, pImg1[+3], pImg2[-3]) );
            sum += filterCoeff[10] * ( clipALF(filterClipp[10], curr, pImg1[+2], pImg2[-2]) );
            sum += filterCoeff[11] * ( clipALF(filterClipp[11], curr, pImg1[+1], pImg2[-1]) );
            sum += filterCoeff[12] * ( clipALF(filterClipp[12], curr, pImg1[+0], pImg2[+0]) );
            sum += filterCoeff[13] * ( clipALF(filterClipp[13], curr, pImg1[-1], pImg2[+1]) );
            sum += filterCoeff[14] * ( clipALF(filterClipp[14], curr, pImg1[-2], pImg2[+2]) );
            sum += filterCoeff[15] * ( clipALF(filterClipp[15], curr, pImg1[-3], pImg2[+3]) );

            sum += filterCoeff[16] * ( clipALF(filterClipp[16], curr, pImg0[+4], pImg0[-4]) );
            sum += filterCoeff[17] * ( clipALF(filterClipp[17], curr, pImg0[+3], pImg0[-3]) );
            sum += filterCoeff[18] * ( clipALF(filterClipp[18], curr, pImg0[+2], pImg0[-2]) );
            sum += filterCoeff[19] * ( clipALF(filterClipp[19], curr, pImg0[+1], pImg0[-1]) );
          }
          else
#endif
          if( filtType == ALF_FILTER_7 )
          {
            sum += filterCoeff[0] * ( clipALF(filterClipp[0], curr, pImg5[+0], pImg6[+0]) );

            sum += filterCoeff[1] * ( clipALF(filterClipp[1], curr, pImg3[+1], pImg4[-1]) );
            sum += filterCoeff[2] * ( clipALF(filterClipp[2], curr, pImg3[+0], pImg4[+0]) );
            sum += filterCoeff[3] * ( clipALF(filterClipp[3], curr, pImg3[-1], pImg4[+1]) );

            sum += filterCoeff[4] * ( clipALF(filterClipp[4], curr, pImg1[+2], pImg2[-2]) );
            sum += filterCoeff[5] * ( clipALF(filterClipp[5], curr, pImg1[+1], pImg2[-1]) );
            sum += filterCoeff[6] * ( clipALF(filterClipp[6], curr, pImg1[+0], pImg2[+0]) );
            sum += filterCoeff[7] * ( clipALF(filterClipp[7], curr, pImg1[-1], pImg2[+1]) );
            sum += filterCoeff[8] * ( clipALF(filterClipp[8], curr, pImg1[-2], pImg2[+2]) );

            sum += filterCoeff[9] * ( clipALF(filterClipp[9], curr, pImg0[+3], pImg0[-3]) );
            sum += filterCoeff[10] * ( clipALF(filterClipp[10], curr, pImg0[+2], pImg0[-2]) );
            sum += filterCoeff[11] * ( clipALF(filterClipp[11], curr, pImg0[+1], pImg0[-1]) );
          }
          else
#if ALF_IMPROVEMENT
          if( filtType == ALF_FILTER_5 )
#endif
          {
            sum += filterCoeff[0] * ( clipALF(filterClipp[0], curr, pImg3[+0], pImg4[+0]) );

            sum += filterCoeff[1] * ( clipALF(filterClipp[1], curr, pImg1[+1], pImg2[-1]) );
            sum += filterCoeff[2] * ( clipALF(filterClipp[2], curr, pImg1[+0], pImg2[+0]) );
            sum += filterCoeff[3] * ( clipALF(filterClipp[3], curr, pImg1[-1], pImg2[+1]) );

            sum += filterCoeff[4] * ( clipALF(filterClipp[4], curr, pImg0[+2], pImg0[-2]) );
            sum += filterCoeff[5] * ( clipALF(filterClipp[5], curr, pImg0[+1], pImg0[-1]) );
          }
#if ALF_IMPROVEMENT
          if( filtType == ALF_FILTER_9_EXT )
          {
            for( int k = 0; k < EXT_LENGTH; k++ )
            {
              sum += coef[k + 20] * ( clipALF(clip[k + 20], curr, fixedFilterResults[k * EXT_LENGTH + fixedFilterSetIdx][blkDst.y + i + ii][blkDst.x + j + jj]) );
            }
          }
          else if( filtType == ALF_FILTER_EXT )
          {
            for( int k = 0; k < EXT_LENGTH; k++ )
            {
              sum += coef[k] * ( clipALF(clip[k], curr, fixedFilterResults[k * EXT_LENGTH + fixedFilterSetIdx][blkDst.y + i + ii][blkDst.x + j + jj] ) );
            }
          }
          sum = ( sum + offset ) >> shift;
#else
          if (!(isNearVBabove || isNearVBbelow))
          {
            sum = (sum + offset) >> shift;
          }
          else
          {
            sum = (sum + (1 << ((shift + 3) - 1))) >> (shift + 3);
          }
#endif
          sum += curr;
          pRec1[jj] = ClipPel( sum, clpRng );

          pImg0++;
          pImg1++;
          pImg2++;
          pImg3++;
          pImg4++;
          pImg5++;
          pImg6++;
#if ALF_IMPROVEMENT
          pImg7++;
          pImg8++;
#endif
        }
      }
    }

    pRec0 += dstStride2;
    pRec1 += dstStride2;

    pImgYPad0 += srcStride2;
    pImgYPad1 += srcStride2;
    pImgYPad2 += srcStride2;
    pImgYPad3 += srcStride2;
    pImgYPad4 += srcStride2;
    pImgYPad5 += srcStride2;
    pImgYPad6 += srcStride2;
#if ALF_IMPROVEMENT
    pImgYPad7 += srcStride2;
    pImgYPad8 += srcStride2;
#endif
  }
}

template<AlfFilterType filtTypeCcAlf>
#if ALF_IMPROVEMENT
void AdaptiveLoopFilter::filterBlkCcAlf(const PelBuf &dstBuf, const CPelUnitBuf &recSrc, const Area &blkDst, const Area &blkSrc, const ComponentID compId, const int16_t *filterCoeff, const ClpRngs &clpRngs, CodingStructure &cs)
#else
void AdaptiveLoopFilter::filterBlkCcAlf(const PelBuf &dstBuf, const CPelUnitBuf &recSrc, const Area &blkDst, const Area &blkSrc, const ComponentID compId, const int16_t *filterCoeff, const ClpRngs &clpRngs, CodingStructure &cs, int vbCTUHeight, int vbPos)
#endif
{
#if !ALF_IMPROVEMENT
  CHECK(1 << floorLog2(vbCTUHeight) != vbCTUHeight, "Not a power of 2");
#endif

  CHECK(!isChroma(compId), "Must be chroma");

  const SPS*     sps           = cs.slice->getSPS();
  ChromaFormat nChromaFormat   = sps->getChromaFormatIdc();
  const int clsSizeY           = 4;
  const int clsSizeX           = 4;
  const int      startHeight        = blkDst.y;
  const int      endHeight          = blkDst.y + blkDst.height;
  const int      startWidth         = blkDst.x;
  const int      endWidth           = blkDst.x + blkDst.width;
  const int scaleX             = getComponentScaleX(compId, nChromaFormat);
  const int scaleY             = getComponentScaleY(compId, nChromaFormat);

  CHECK( startHeight % clsSizeY, "Wrong startHeight in filtering" );
  CHECK( startWidth % clsSizeX, "Wrong startWidth in filtering" );
  CHECK( ( endHeight - startHeight ) % clsSizeY, "Wrong endHeight in filtering" );
  CHECK( ( endWidth - startWidth ) % clsSizeX, "Wrong endWidth in filtering" );

  CPelBuf     srcBuf     = recSrc.get(COMPONENT_Y);
  const int   lumaStride = srcBuf.stride;
  const Pel * lumaPtr    = srcBuf.buf + blkSrc.y * lumaStride + blkSrc.x;

  const int   chromaStride = dstBuf.stride;
  Pel *       chromaPtr    = dstBuf.buf + blkDst.y * chromaStride + blkDst.x;

  for( int i = 0; i < endHeight - startHeight; i += clsSizeY )
  {
    for( int j = 0; j < endWidth - startWidth; j += clsSizeX )
    {
      for( int ii = 0; ii < clsSizeY; ii++ )
      {
        int row       = ii;
        int col       = j;
        Pel *srcSelf  = chromaPtr + col + row * chromaStride;

        int offset1 = lumaStride;
        int offset2 = -lumaStride;
        int offset3 = 2 * lumaStride;
        row <<= scaleY;
        col <<= scaleX;
        const Pel *srcCross = lumaPtr + col + row * lumaStride;

#if !ALF_IMPROVEMENT
        int pos = ((startHeight + i + ii) << scaleY) & (vbCTUHeight - 1);
        if (scaleY == 0 && (pos == vbPos || pos == vbPos + 1))
        {
          continue;
        }
        if (pos == (vbPos - 2) || pos == (vbPos + 1))
        {
          offset3 = offset1;
        }
        else if (pos == (vbPos - 1) || pos == vbPos)
        {
          offset1 = 0;
          offset2 = 0;
          offset3 = 0;
        }
#endif

        for (int jj = 0; jj < clsSizeX; jj++)
        {
          const int jj2     = (jj << scaleX);
          const int offset0 = 0;

          int sum = 0;
          const Pel currSrcCross = srcCross[offset0 + jj2];
          sum += filterCoeff[0] * (srcCross[offset2 + jj2    ] - currSrcCross);
          sum += filterCoeff[1] * (srcCross[offset0 + jj2 - 1] - currSrcCross);
          sum += filterCoeff[2] * (srcCross[offset0 + jj2 + 1] - currSrcCross);
          sum += filterCoeff[3] * (srcCross[offset1 + jj2 - 1] - currSrcCross);
          sum += filterCoeff[4] * (srcCross[offset1 + jj2    ] - currSrcCross);
          sum += filterCoeff[5] * (srcCross[offset1 + jj2 + 1] - currSrcCross);
          sum += filterCoeff[6] * (srcCross[offset3 + jj2    ] - currSrcCross);

          sum = (sum + ((1 << m_scaleBits ) >> 1)) >> m_scaleBits;
          const int offset = 1 << clpRngs.comp[compId].bd >> 1;
          sum = ClipPel(sum + offset, clpRngs.comp[compId]) - offset;
          sum += srcSelf[jj];
          srcSelf[jj] = ClipPel(sum, clpRngs.comp[compId]);
        }
      }
    }

    chromaPtr += chromaStride * clsSizeY;

    lumaPtr += lumaStride * clsSizeY << getComponentScaleY(compId, nChromaFormat);
  }
}

