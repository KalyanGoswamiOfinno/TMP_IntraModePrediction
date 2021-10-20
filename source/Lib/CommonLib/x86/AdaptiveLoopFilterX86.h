/* The copyright in this software is being made available under the BSD
 * License, included below. This software may be subject to other third party
 * and contributor rights, including patent rights, and no such rights are
 * granted under this license.
 *
 * Copyright (c) 2010-2017, ITU/ISO/IEC
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

#pragma once

#include "CommonDefX86.h"
#include "../AdaptiveLoopFilter.h"
#if ALF_IMPROVEMENT
#include "AlfFixedFilters.h"
#endif

#ifdef TARGET_SIMD_X86
#if defined _MSC_VER
#include <tmmintrin.h>
#else
#include <x86intrin.h>
#endif

constexpr uint16_t sh(int x)
{
  return 0x0202 * (x & 7) + 0x0100 + 0x1010 * (x & 8);
}

#if ALF_IMPROVEMENT
static const uint16_t shuffleTab5[4][8] = {
  {
    sh(0), sh(1), sh(2), sh(3), sh(4), sh(5), sh(6), sh(7),
  },
  {
    sh(4), sh(1), sh(5), sh(3), sh(0), sh(2), sh(6), sh(7),
  },
  {
    sh(0), sh(3), sh(2), sh(1), sh(4), sh(5), sh(6), sh(7),
  },
  {
    sh(4), sh(3), sh(5), sh(1), sh(0), sh(2), sh(6), sh(7) ,
  },
};

constexpr int p0[42] =
{
  0, 2, 6, 12, 20, 30,
  36, 37, 38, 39, 40, 41,
  1, 4, 9, 16, 25, 11, 18, 27,
  3, 8, 15, 24, 35, 13, 22, 33,
  5, 10, 17, 26, 19, 28, 29, 31,
  7, 14, 23, 34, 21, 32
};

static std::array<std::array<std::array<std::array<short, 42>, NUM_FIXED_FILTERS>, NUM_CLASSIFIER>, NUM_SETS_FIXED_FILTERS> packedDataFixedFilters;
#else
template<X86_VEXT vext>
static void simdDeriveClassificationBlk(AlfClassifier **classifier, int **laplacian[NUM_DIRECTIONS], const CPelBuf &srcLuma, const Area &blkDst, const Area &blk, const int shift, const int vbCTUHeight, int vbPos)
{
  CHECK((blk.height & 7) != 0, "Block height must be a multiple of 8");
  CHECK((blk.width & 7) != 0, "Block width must be a multiple of 8");
  CHECK((vbCTUHeight & (vbCTUHeight - 1)) != 0, "vbCTUHeight must be a power of 2");

  const size_t imgStride = srcLuma.stride;
  const Pel *  srcExt    = srcLuma.buf;

  const int imgHExtended = blk.height + 4;
  const int imgWExtended = blk.width + 4;

  const int posX = blk.pos().x;
  const int posY = blk.pos().y;

  // 18x40 array
  uint16_t colSums[(AdaptiveLoopFilter::m_CLASSIFICATION_BLK_SIZE + 4) >> 1]
                  [AdaptiveLoopFilter::m_CLASSIFICATION_BLK_SIZE + 8];

  for (int i = 0; i < imgHExtended; i += 2)
  {
    const size_t offset = (i + posY - 3) * imgStride + posX - 3;

    const Pel *imgY0 = &srcExt[offset];
    const Pel *imgY1 = &srcExt[offset + imgStride];
    const Pel *imgY2 = &srcExt[offset + imgStride * 2];
    const Pel *imgY3 = &srcExt[offset + imgStride * 3];

    // pixel padding for gradient calculation
    int pos      = blkDst.pos().y - 2 + i;
    int posInCTU = pos & (vbCTUHeight - 1);
    if (pos > 0 && posInCTU == vbPos - 2)
    {
      imgY3 = imgY2;
    }
    else if (pos > 0 && posInCTU == vbPos)
    {
      imgY0 = imgY1;
    }

    __m128i prev = _mm_setzero_si128();

    for (int j = 0; j < imgWExtended; j += 8)
    {
      const __m128i x0 = _mm_loadu_si128((const __m128i *) (imgY0 + j));
      const __m128i x1 = _mm_loadu_si128((const __m128i *) (imgY1 + j));
      const __m128i x2 = _mm_loadu_si128((const __m128i *) (imgY2 + j));
      const __m128i x3 = _mm_loadu_si128((const __m128i *) (imgY3 + j));

      const __m128i x4 = _mm_loadu_si128((const __m128i *) (imgY0 + j + 2));
      const __m128i x5 = _mm_loadu_si128((const __m128i *) (imgY1 + j + 2));
      const __m128i x6 = _mm_loadu_si128((const __m128i *) (imgY2 + j + 2));
      const __m128i x7 = _mm_loadu_si128((const __m128i *) (imgY3 + j + 2));

      const __m128i nw = _mm_blend_epi16(x0, x1, 0xaa);
      const __m128i n  = _mm_blend_epi16(x0, x5, 0x55);
      const __m128i ne = _mm_blend_epi16(x4, x5, 0xaa);
      const __m128i w  = _mm_blend_epi16(x1, x2, 0xaa);
      const __m128i e  = _mm_blend_epi16(x5, x6, 0xaa);
      const __m128i sw = _mm_blend_epi16(x2, x3, 0xaa);
      const __m128i s  = _mm_blend_epi16(x2, x7, 0x55);
      const __m128i se = _mm_blend_epi16(x6, x7, 0xaa);

      __m128i c = _mm_blend_epi16(x1, x6, 0x55);
      c         = _mm_add_epi16(c, c);
      __m128i d = _mm_shuffle_epi8(c, _mm_setr_epi8(2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9, 14, 15, 12, 13));

      const __m128i ver = _mm_abs_epi16(_mm_sub_epi16(c, _mm_add_epi16(n, s)));
      const __m128i hor = _mm_abs_epi16(_mm_sub_epi16(d, _mm_add_epi16(w, e)));
      const __m128i di0 = _mm_abs_epi16(_mm_sub_epi16(d, _mm_add_epi16(nw, se)));
      const __m128i di1 = _mm_abs_epi16(_mm_sub_epi16(d, _mm_add_epi16(ne, sw)));

      const __m128i hv  = _mm_hadd_epi16(ver, hor);
      const __m128i di  = _mm_hadd_epi16(di0, di1);
      const __m128i all = _mm_hadd_epi16(hv, di);

      const __m128i t = _mm_blend_epi16(all, prev, 0xaa);
      _mm_storeu_si128((__m128i *) &colSums[i >> 1][j], _mm_hadd_epi16(t, all));
      prev = all;
    }
  }

  for (int i = 0; i < (blk.height >> 1); i += 4)
  {
    for (int j = 0; j < blk.width; j += 8)
    {
      __m128i x0, x1, x2, x3, x4, x5, x6, x7;

      const uint32_t z = (2 * i + blkDst.pos().y) & (vbCTUHeight - 1);
      const uint32_t z2 = (2 * i + 4 + blkDst.pos().y) & (vbCTUHeight - 1);

      x0 = (z == vbPos) ? _mm_setzero_si128() : _mm_loadu_si128((__m128i *) &colSums[i + 0][j + 4]);
      x1 = _mm_loadu_si128((__m128i *) &colSums[i + 1][j + 4]);
      x2 = _mm_loadu_si128((__m128i *) &colSums[i + 2][j + 4]);
      x3 = (z == vbPos - 4) ? _mm_setzero_si128() : _mm_loadu_si128((__m128i *) &colSums[i + 3][j + 4]);

      x4 = (z2 == vbPos) ? _mm_setzero_si128() : _mm_loadu_si128((__m128i *) &colSums[i + 2][j + 4]);
      x5 = _mm_loadu_si128((__m128i *) &colSums[i + 3][j + 4]);
      x6 = _mm_loadu_si128((__m128i *) &colSums[i + 4][j + 4]);
      x7 = (z2 == vbPos - 4) ? _mm_setzero_si128() : _mm_loadu_si128((__m128i *) &colSums[i + 5][j + 4]);

      __m128i x0l = _mm_cvtepu16_epi32(x0);
      __m128i x0h = _mm_unpackhi_epi16(x0, _mm_setzero_si128());
      __m128i x1l = _mm_cvtepu16_epi32(x1);
      __m128i x1h = _mm_unpackhi_epi16(x1, _mm_setzero_si128());
      __m128i x2l = _mm_cvtepu16_epi32(x2);
      __m128i x2h = _mm_unpackhi_epi16(x2, _mm_setzero_si128());
      __m128i x3l = _mm_cvtepu16_epi32(x3);
      __m128i x3h = _mm_unpackhi_epi16(x3, _mm_setzero_si128());
      __m128i x4l = _mm_cvtepu16_epi32(x4);
      __m128i x4h = _mm_unpackhi_epi16(x4, _mm_setzero_si128());
      __m128i x5l = _mm_cvtepu16_epi32(x5);
      __m128i x5h = _mm_unpackhi_epi16(x5, _mm_setzero_si128());
      __m128i x6l = _mm_cvtepu16_epi32(x6);
      __m128i x6h = _mm_unpackhi_epi16(x6, _mm_setzero_si128());
      __m128i x7l = _mm_cvtepu16_epi32(x7);
      __m128i x7h = _mm_unpackhi_epi16(x7, _mm_setzero_si128());

      x0l = _mm_add_epi32(x0l, x1l);
      x2l = _mm_add_epi32(x2l, x3l);
      x4l = _mm_add_epi32(x4l, x5l);
      x6l = _mm_add_epi32(x6l, x7l);
      x0h = _mm_add_epi32(x0h, x1h);
      x2h = _mm_add_epi32(x2h, x3h);
      x4h = _mm_add_epi32(x4h, x5h);
      x6h = _mm_add_epi32(x6h, x7h);

      x0l = _mm_add_epi32(x0l, x2l);
      x4l = _mm_add_epi32(x4l, x6l);
      x0h = _mm_add_epi32(x0h, x2h);
      x4h = _mm_add_epi32(x4h, x6h);

      x2l = _mm_unpacklo_epi32(x0l, x4l);
      x2h = _mm_unpackhi_epi32(x0l, x4l);
      x6l = _mm_unpacklo_epi32(x0h, x4h);
      x6h = _mm_unpackhi_epi32(x0h, x4h);

      __m128i sumV  = _mm_unpacklo_epi32(x2l, x6l);
      __m128i sumH  = _mm_unpackhi_epi32(x2l, x6l);
      __m128i sumD0 = _mm_unpacklo_epi32(x2h, x6h);
      __m128i sumD1 = _mm_unpackhi_epi32(x2h, x6h);

      //      uint32_t tempAct = sumV + sumH;
      __m128i tempAct = _mm_add_epi32(sumV, sumH);

      //      const uint32_t activity = std::min<uint32_t>(15, tempAct * scale >> shift);
      //      static const uint8_t th[16] = { 0, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4 };
      //      uint8_t classIdx = th[activity];
      const uint32_t scale  = (z == vbPos - 4 || z == vbPos) ? 96 : 64;
      const uint32_t scale2 = (z2 == vbPos - 4 || z2 == vbPos) ? 96 : 64;
      __m128i activity = _mm_mullo_epi32(tempAct, _mm_unpacklo_epi64(_mm_set1_epi32(scale), _mm_set1_epi32(scale2)));
      activity         = _mm_srl_epi32(activity, _mm_cvtsi32_si128(shift));
      activity         = _mm_min_epi32(activity, _mm_set1_epi32(15));
      __m128i classIdx = _mm_shuffle_epi8(_mm_setr_epi8(0, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4), activity);

      //      if (sumV > sumH)
      //      {
      //        hv1       = sumV;
      //        hv0       = sumH;
      //        dirTempHV = 0;
      //      }
      //      else
      //      {
      //        hv1       = sumH;
      //        hv0       = sumV;
      //        dirTempHV = 1;
      //      }
      __m128i dirTempHVMinus1 = _mm_cmpgt_epi32(sumV, sumH);
      __m128i hv1             = _mm_max_epi32(sumV, sumH);
      __m128i hv0             = _mm_min_epi32(sumV, sumH);

      //      if (sumD0 > sumD1)
      //      {
      //        d1       = sumD0;
      //        d0       = sumD1;
      //        dirTempD = 0;
      //      }
      //      else
      //      {
      //        d1       = sumD1;
      //        d0       = sumD0;
      //        dirTempD = 1;
      //      }
      __m128i dirTempDMinus1 = _mm_cmpgt_epi32(sumD0, sumD1);
      __m128i d1             = _mm_max_epi32(sumD0, sumD1);
      __m128i d0             = _mm_min_epi32(sumD0, sumD1);

      //      int dirIdx;
      //      if (d1 * hv0 > hv1 * d0)
      //      {
      //        hvd1   = d1;
      //        hvd0   = d0;
      //        dirIdx = 0;
      //      }
      //      else
      //      {
      //        hvd1   = hv1;
      //        hvd0   = hv0;
      //        dirIdx = 2;
      //      }
      __m128i a      = _mm_xor_si128(_mm_mullo_epi32(d1, hv0), _mm_set1_epi32(0x80000000));
      __m128i b      = _mm_xor_si128(_mm_mullo_epi32(hv1, d0), _mm_set1_epi32(0x80000000));
      __m128i dirIdx = _mm_cmpgt_epi32(a, b);
      __m128i hvd1   = _mm_blendv_epi8(hv1, d1, dirIdx);
      __m128i hvd0   = _mm_blendv_epi8(hv0, d0, dirIdx);

      //      if (hvd1 * 2 > 9 * hvd0)
      //      {
      //        classIdx += (dirIdx + 2) * 5;
      //      }
      //      else if (hvd1 > 2 * hvd0)
      //      {
      //        classIdx += (dirIdx + 1) * 5;
      //      }
      __m128i strength1 = _mm_cmpgt_epi32(hvd1, _mm_add_epi32(hvd0, hvd0));
      __m128i strength2 = _mm_cmpgt_epi32(_mm_add_epi32(hvd1, hvd1), _mm_add_epi32(hvd0, _mm_slli_epi32(hvd0, 3)));
      __m128i offset    = _mm_and_si128(strength1, _mm_set1_epi32(5));
      classIdx          = _mm_add_epi32(classIdx, offset);
      classIdx          = _mm_add_epi32(classIdx, _mm_and_si128(strength2, _mm_set1_epi32(5)));
      offset            = _mm_andnot_si128(dirIdx, offset);
      offset            = _mm_add_epi32(offset, offset);
      classIdx          = _mm_add_epi32(classIdx, offset);

      //      uint8_t transposeIdx = 2 * dirTempD + dirTempHV;
      __m128i transposeIdx = _mm_set1_epi32(3);
      transposeIdx         = _mm_add_epi32(transposeIdx, dirTempHVMinus1);
      transposeIdx         = _mm_add_epi32(transposeIdx, dirTempDMinus1);
      transposeIdx         = _mm_add_epi32(transposeIdx, dirTempDMinus1);

      int yOffset = 2 * i + blkDst.pos().y;
      int xOffset = j + blkDst.pos().x;

      static_assert(sizeof(AlfClassifier) == 2, "ALFClassifier type must be 16 bits wide");
      __m128i v;
      v = _mm_unpacklo_epi8(classIdx, transposeIdx);
      v = _mm_shuffle_epi8(v, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 8, 9, 8, 9, 8, 9, 8, 9));
      _mm_storeu_si128((__m128i *) (classifier[yOffset] + xOffset), v);
      _mm_storeu_si128((__m128i *) (classifier[yOffset + 1] + xOffset), v);
      _mm_storeu_si128((__m128i *) (classifier[yOffset + 2] + xOffset), v);
      _mm_storeu_si128((__m128i *) (classifier[yOffset + 3] + xOffset), v);
      v = _mm_unpackhi_epi8(classIdx, transposeIdx);
      v = _mm_shuffle_epi8(v, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 8, 9, 8, 9, 8, 9, 8, 9));
      _mm_storeu_si128((__m128i *) (classifier[yOffset + 4] + xOffset), v);
      _mm_storeu_si128((__m128i *) (classifier[yOffset + 5] + xOffset), v);
      _mm_storeu_si128((__m128i *) (classifier[yOffset + 6] + xOffset), v);
      _mm_storeu_si128((__m128i *) (classifier[yOffset + 7] + xOffset), v);
    }
  }
}
#endif

template<X86_VEXT vext>
static void simdFilter5x5Blk(AlfClassifier **classifier, const PelUnitBuf &recDst, const CPelUnitBuf &recSrc, const Area &blkDst, const Area &blk, const ComponentID compId, const short *filterSet,
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  const Pel *fClipSet, const ClpRng &clpRng, CodingStructure &cs
#else  
  const short *fClipSet, const ClpRng &clpRng, CodingStructure &cs
#endif
#if ALF_IMPROVEMENT
  , Pel ***fixedFilterResults, int fixedFilterSetIdx
#else
  , const int vbCTUHeight, int vbPos
#endif
)
{
#if !ALF_IMPROVEMENT
  CHECK((vbCTUHeight & (vbCTUHeight - 1)) != 0, "vbCTUHeight must be a power of 2");
  CHECK(!isChroma(compId), "ALF 5x5 filter is for chroma only");
#endif

  const CPelBuf srcBuffer = recSrc.get(compId);
  PelBuf        dstBuffer = recDst.get(compId);

  const size_t srcStride = srcBuffer.stride;
  const size_t dstStride = dstBuffer.stride;

  constexpr int SHIFT = AdaptiveLoopFilter::m_NUM_BITS - 1;
  constexpr int ROUND = 1 << (SHIFT - 1);
#if !ALF_IMPROVEMENT
  const __m128i mmOffset1 = _mm_set1_epi32((1 << ((SHIFT + 3) - 1)) - ROUND);
#endif

  const size_t width  = blk.width;
  const size_t height = blk.height;

  constexpr size_t STEP_X = 8;
#if ALF_IMPROVEMENT
  static_assert(sizeof(*filterSet) == 2, "ALF coeffs must be 16-bit wide");
  static_assert(sizeof(*fClipSet) == 2, "ALF clip values must be 16-bit wide");
  size_t STEP_Y = isChroma(compId) ? 4 : 2;
#else
  constexpr size_t STEP_Y = 4;
#endif

  CHECK(blk.y % STEP_Y, "Wrong startHeight in filtering");
  CHECK(blk.x % STEP_X, "Wrong startWidth in filtering");
  CHECK(height % STEP_Y, "Wrong endHeight in filtering");
  CHECK(width % 4, "Wrong endWidth in filtering");

  const Pel *src = srcBuffer.buf + blk.y * srcStride + blk.x;
  Pel *      dst = dstBuffer.buf + blkDst.y * dstStride + blkDst.x;



  const __m128i mmOffset = _mm_set1_epi32(ROUND);
  const __m128i mmMin = _mm_set1_epi16( clpRng.min );
  const __m128i mmMax = _mm_set1_epi16( clpRng.max );

#if !ALF_IMPROVEMENT
  __m128i params[2][3];
  __m128i fs   = _mm_loadu_si128((__m128i *) filterSet);
  params[0][0] = _mm_shuffle_epi32(fs, 0x00);
  params[0][1] = _mm_shuffle_epi32(fs, 0x55);
  params[0][2] = _mm_shuffle_epi32(fs, 0xaa);
  __m128i fc   = _mm_loadu_si128((__m128i *) fClipSet);
  params[1][0] = _mm_shuffle_epi32(fc, 0x00);
  params[1][1] = _mm_shuffle_epi32(fc, 0x55);
  params[1][2] = _mm_shuffle_epi32(fc, 0xaa);
#endif

  for (size_t i = 0; i < height; i += STEP_Y)
  {
#if ALF_IMPROVEMENT
    const AlfClassifier *pClass = isChroma(compId) ? nullptr : classifier[blkDst.y + i] + blkDst.x;
#endif
    for (size_t j = 0; j < width; j += STEP_X)
    {
#if ALF_IMPROVEMENT
      __m128i params[2][2][3];
      for (int k = 0; k < 2; ++k)
      {
        const int transposeIdx = pClass ? (pClass[j + 4 * k] & 0x3)  : 0;
        const int classIdx = pClass ? ( pClass[j + 4 * k] >> 2) : 0;
        const int transposeIdx1 = pClass ? ( pClass[j + 4 * k + 2] & 0x3) : 0;
        const int classIdx1 = pClass ? ( pClass[j + 4 * k + 2] >> 2 ) : 0;

        __m128i rawCoeff = _mm_loadu_si128((const __m128i *) (filterSet + classIdx * MAX_NUM_ALF_LUMA_COEFF));
        __m128i rawClip = _mm_loadu_si128((const __m128i *) (fClipSet + classIdx * MAX_NUM_ALF_LUMA_COEFF));

        __m128i s0 = _mm_loadu_si128((const __m128i *) shuffleTab5[transposeIdx]);
        __m128i s1 = _mm_xor_si128(s0, _mm_set1_epi8((char)0x80));

        rawCoeff = _mm_or_si128(_mm_shuffle_epi8(rawCoeff, s0), _mm_shuffle_epi8(rawCoeff, s1));
        rawClip = _mm_or_si128(_mm_shuffle_epi8(rawClip, s0), _mm_shuffle_epi8(rawClip, s1));

        __m128i rawCoeff1 = _mm_loadu_si128((const __m128i *) (filterSet + classIdx1 * MAX_NUM_ALF_LUMA_COEFF));
        __m128i rawClip1 = _mm_loadu_si128((const __m128i *) (fClipSet + classIdx1 * MAX_NUM_ALF_LUMA_COEFF));

        s0 = _mm_loadu_si128((const __m128i *) shuffleTab5[transposeIdx1]);
        s1 = _mm_xor_si128(s0, _mm_set1_epi8((char)0x80));

        rawCoeff1 = _mm_or_si128(_mm_shuffle_epi8(rawCoeff1, s0), _mm_shuffle_epi8(rawCoeff, s1));
        rawClip1 = _mm_or_si128(_mm_shuffle_epi8(rawClip1, s0), _mm_shuffle_epi8(rawClip, s1));

        params[k][0][0] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoeff, 0x00), _mm_shuffle_epi32(rawCoeff1, 0x00));
        params[k][0][1] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoeff, 0x55), _mm_shuffle_epi32(rawCoeff1, 0x55));
        params[k][0][2] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoeff, 0xaa), _mm_shuffle_epi32(rawCoeff1, 0xaa));
        params[k][1][0] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClip, 0x00), _mm_shuffle_epi32(rawClip1, 0x00));
        params[k][1][1] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClip, 0x55), _mm_shuffle_epi32(rawClip1, 0x55));
        params[k][1][2] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClip, 0xaa), _mm_shuffle_epi32(rawClip1, 0xaa));
      }
#endif
      for (size_t ii = 0; ii < STEP_Y; ii++)
      {
        const Pel *pImg0, *pImg1, *pImg2, *pImg3, *pImg4;

        pImg0 = src + j + ii * srcStride;
        pImg1 = pImg0 + srcStride;
        pImg2 = pImg0 - srcStride;
        pImg3 = pImg1 + srcStride;
        pImg4 = pImg2 - srcStride;
#if !ALF_IMPROVEMENT
        const int yVb = (blkDst.y + i + ii) & (vbCTUHeight - 1);
        if (yVb < vbPos && (yVb >= vbPos - 2))   // above
        {
          pImg1 = (yVb == vbPos - 1) ? pImg0 : pImg1;
          pImg3 = (yVb >= vbPos - 2) ? pImg1 : pImg3;

          pImg2 = (yVb == vbPos - 1) ? pImg0 : pImg2;
          pImg4 = (yVb >= vbPos - 2) ? pImg2 : pImg4;
        }
        else if (yVb >= vbPos && (yVb <= vbPos + 1))   // bottom
        {
          pImg2 = (yVb == vbPos) ? pImg0 : pImg2;
          pImg4 = (yVb <= vbPos + 1) ? pImg2 : pImg4;

          pImg1 = (yVb == vbPos) ? pImg0 : pImg1;
          pImg3 = (yVb <= vbPos + 1) ? pImg1 : pImg3;
        }
#endif
        __m128i cur = _mm_loadu_si128((const __m128i *) pImg0);
        __m128i accumA = mmOffset;
        __m128i accumB = mmOffset;

        auto process2coeffs = [&](const int i, const Pel *ptr0, const Pel *ptr1, const Pel *ptr2, const Pel *ptr3) {
          const __m128i val00 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) ptr0), cur);
          const __m128i val10 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) ptr2), cur);
          const __m128i val01 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) ptr1), cur);
          const __m128i val11 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) ptr3), cur);
          __m128i val01A = _mm_unpacklo_epi16(val00, val10);
          __m128i val01B = _mm_unpackhi_epi16(val00, val10);
          __m128i val01C = _mm_unpacklo_epi16(val01, val11);
          __m128i val01D = _mm_unpackhi_epi16(val01, val11);

#if ALF_IMPROVEMENT
          __m128i limit01A = params[0][1][i];
          __m128i limit01B = params[1][1][i];
          val01A = _mm_min_epi16(val01A, limit01A);
          val01B = _mm_min_epi16(val01B, limit01B);
          val01C = _mm_min_epi16(val01C, limit01A);
          val01D = _mm_min_epi16(val01D, limit01B);

          limit01A = _mm_sub_epi16(_mm_setzero_si128(), limit01A);
          limit01B = _mm_sub_epi16(_mm_setzero_si128(), limit01B);

          val01A = _mm_max_epi16(val01A, limit01A);
          val01B = _mm_max_epi16(val01B, limit01B);
          val01C = _mm_max_epi16(val01C, limit01A);
          val01D = _mm_max_epi16(val01D, limit01B);

          val01A = _mm_add_epi16(val01A, val01C);
          val01B = _mm_add_epi16(val01B, val01D);

          const __m128i coeff01A = params[0][0][i];
          const __m128i coeff01B = params[1][0][i];

          accumA = _mm_add_epi32(accumA, _mm_madd_epi16(val01A, coeff01A));
          accumB = _mm_add_epi32(accumB, _mm_madd_epi16(val01B, coeff01B));
#else
          __m128i limit01A = params[1][i];

          val01A = _mm_min_epi16(val01A, limit01A);
          val01B = _mm_min_epi16(val01B, limit01A);
          val01C = _mm_min_epi16(val01C, limit01A);
          val01D = _mm_min_epi16(val01D, limit01A);

          limit01A = _mm_sub_epi16(_mm_setzero_si128(), limit01A);

          val01A = _mm_max_epi16(val01A, limit01A);
          val01B = _mm_max_epi16(val01B, limit01A);
          val01C = _mm_max_epi16(val01C, limit01A);
          val01D = _mm_max_epi16(val01D, limit01A);

          val01A = _mm_add_epi16(val01A, val01C);
          val01B = _mm_add_epi16(val01B, val01D);

          __m128i coeff01A = params[0][i];

          accumA = _mm_add_epi32(accumA, _mm_madd_epi16(val01A, coeff01A));
          accumB = _mm_add_epi32(accumB, _mm_madd_epi16(val01B, coeff01A));
#endif
        };

        process2coeffs(0, pImg3 + 0, pImg4 + 0, pImg1 + 1, pImg2 - 1);
        process2coeffs(1, pImg1 + 0, pImg2 + 0, pImg1 - 1, pImg2 + 1);
        process2coeffs(2, pImg0 + 2, pImg0 - 2, pImg0 + 1, pImg0 - 1);
#if ALF_IMPROVEMENT
        accumA = _mm_srai_epi32( accumA, SHIFT );
        accumB = _mm_srai_epi32( accumB, SHIFT );
#else

        bool isNearVBabove = yVb < vbPos && (yVb >= vbPos - 1);
        bool isNearVBbelow = yVb >= vbPos && (yVb <= vbPos);
        if (!(isNearVBabove || isNearVBbelow))
        {
          accumA = _mm_srai_epi32(accumA, SHIFT);
          accumB = _mm_srai_epi32(accumB, SHIFT);
        }
        else
        {
          accumA = _mm_srai_epi32(_mm_add_epi32(accumA, mmOffset1), SHIFT + 3);
          accumB = _mm_srai_epi32(_mm_add_epi32(accumB, mmOffset1), SHIFT + 3);
        }
#endif
        accumA = _mm_packs_epi32(accumA, accumB);
        accumA = _mm_add_epi16(accumA, cur);
        accumA = _mm_min_epi16(mmMax, _mm_max_epi16(accumA, mmMin));

        if (j + STEP_X <= width)
        {
          _mm_storeu_si128((__m128i *) (dst + ii * dstStride + j), accumA);
        }
        else
        {
          _mm_storel_epi64((__m128i *) (dst + ii * dstStride + j), accumA);
        }
      }

    }

    src += srcStride * STEP_Y;
    dst += dstStride * STEP_Y;
  }
}

static const uint16_t shuffleTab[4][2][8] = {
  {
    { sh(0), sh(1), sh(2), sh(3), sh(4), sh(5), sh(6), sh(7) },
    { sh(8), sh(9), sh(10), sh(11), sh(12), sh(13), sh(14), sh(15) },
  },
  {
    { sh(9), sh(4), sh(10), sh(8), sh(1), sh(5), sh(11), sh(7) },
    { sh(3), sh(0), sh(2), sh(6), sh(12), sh(13), sh(14), sh(15) },
  },
  {
    { sh(0), sh(3), sh(2), sh(1), sh(8), sh(7), sh(6), sh(5) },
    { sh(4), sh(9), sh(10), sh(11), sh(12), sh(13), sh(14), sh(15) },
  },
  {
    { sh(9), sh(8), sh(10), sh(4), sh(3), sh(7), sh(11), sh(5) },
    { sh(1), sh(0), sh(2), sh(6), sh(12), sh(13), sh(14), sh(15) },
  },
};

template<X86_VEXT vext>
static void simdFilter7x7Blk(AlfClassifier **classifier, const PelUnitBuf &recDst, const CPelUnitBuf &recSrc, const Area &blkDst, const Area &blk, const ComponentID compId, const short *filterSet,
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  const Pel *fClipSet, const ClpRng &clpRng, CodingStructure &cs
#else  
  const short *fClipSet, const ClpRng &clpRng, CodingStructure &cs
#endif
#if ALF_IMPROVEMENT
  , Pel ***fixedFilterResults, int fixedFilterSetIdx
#else
  , const int vbCTUHeight, int vbPos
#endif
)
{
#if !ALF_IMPROVEMENT
  CHECK((vbCTUHeight & (vbCTUHeight - 1)) != 0, "vbCTUHeight must be a power of 2");
  CHECK(isChroma(compId), "7x7 ALF filter is meant for luma only");
#endif

  const CPelBuf srcBuffer = recSrc.get(compId);
  PelBuf        dstBuffer = recDst.get(compId);

  const size_t srcStride = srcBuffer.stride;
  const size_t dstStride = dstBuffer.stride;

  constexpr int SHIFT = AdaptiveLoopFilter::m_NUM_BITS - 1;
  constexpr int ROUND = 1 << (SHIFT - 1);

  const size_t width  = blk.width;
  const size_t height = blk.height;

  constexpr size_t STEP_X = 8;
#if ALF_IMPROVEMENT
  size_t STEP_Y = isChroma(compId) ? 4 : 2;
#else
  constexpr size_t STEP_Y = 4;
#endif

  CHECK(blk.y % STEP_Y, "Wrong startHeight in filtering");
  CHECK(blk.x % STEP_X, "Wrong startWidth in filtering");
  CHECK(height % STEP_Y, "Wrong endHeight in filtering");
  CHECK(width % STEP_X, "Wrong endWidth in filtering");

  const Pel *src = srcBuffer.buf + blk.y * srcStride + blk.x;
  Pel *      dst = dstBuffer.buf + blkDst.y * dstStride + blkDst.x;

  const __m128i mmOffset = _mm_set1_epi32(ROUND);
#if !ALF_IMPROVEMENT
  const __m128i mmOffset1 = _mm_set1_epi32((1 << ((SHIFT + 3) - 1)) - ROUND);
#endif
  const __m128i mmMin = _mm_set1_epi16( clpRng.min );
  const __m128i mmMax = _mm_set1_epi16( clpRng.max );


  for (size_t i = 0; i < height; i += STEP_Y)
  {
#if ALF_IMPROVEMENT
    const AlfClassifier *pClass = isChroma(compId) ? nullptr : classifier[blkDst.y + i] + blkDst.x;
#else
    const AlfClassifier *pClass = classifier[blkDst.y + i] + blkDst.x;
#endif

    for (size_t j = 0; j < width; j += STEP_X)
    {
      __m128i params[2][2][6];

      for (int k = 0; k < 2; ++k)
      {
#if ALF_IMPROVEMENT
        const int transposeIdx = pClass ? ( pClass[j + 4 * k] & 0x3 ) : 0;
        const int classIdx = pClass ? ( pClass[j + 4 * k] >> 2 ) : 0;
        const int transposeIdx1 = pClass ? ( pClass[j + 4 * k + 2] & 0x3 ) : 0;
        const int classIdx1 = pClass ? ( pClass[j + 4 * k + 2] >> 2 ) : 0;
#else
        const AlfClassifier &cl = pClass[j + 4 * k];

        const int transposeIdx = cl.transposeIdx;
        const int classIdx     = cl.classIdx;
#endif

        static_assert(sizeof(*filterSet) == 2, "ALF coeffs must be 16-bit wide");
        static_assert(sizeof(*fClipSet) == 2, "ALF clip values must be 16-bit wide");

        __m128i rawCoeff0, rawCoeff1;
        __m128i rawClip0, rawClip1;

          rawCoeff0 = _mm_loadu_si128((const __m128i *) (filterSet + classIdx * MAX_NUM_ALF_LUMA_COEFF));
          rawCoeff1 = _mm_loadl_epi64((const __m128i *) (filterSet + classIdx * MAX_NUM_ALF_LUMA_COEFF + 8));

          rawClip0 = _mm_loadu_si128((const __m128i *) (fClipSet + classIdx * MAX_NUM_ALF_LUMA_COEFF));
          rawClip1 = _mm_loadl_epi64((const __m128i *) (fClipSet + classIdx * MAX_NUM_ALF_LUMA_COEFF + 8));

        const __m128i s0 = _mm_loadu_si128((const __m128i *) shuffleTab[transposeIdx][0]);
        const __m128i s1 = _mm_xor_si128(s0, _mm_set1_epi8((char) 0x80));
        const __m128i s2 = _mm_loadu_si128((const __m128i *) shuffleTab[transposeIdx][1]);
        const __m128i s3 = _mm_xor_si128(s2, _mm_set1_epi8((char) 0x80));

        const __m128i rawCoeffLo = _mm_or_si128(_mm_shuffle_epi8(rawCoeff0, s0), _mm_shuffle_epi8(rawCoeff1, s1));
        const __m128i rawCoeffHi = _mm_or_si128(_mm_shuffle_epi8(rawCoeff0, s2), _mm_shuffle_epi8(rawCoeff1, s3));
        const __m128i rawClipLo  = _mm_or_si128(_mm_shuffle_epi8(rawClip0, s0), _mm_shuffle_epi8(rawClip1, s1));
        const __m128i rawClipHi  = _mm_or_si128(_mm_shuffle_epi8(rawClip0, s2), _mm_shuffle_epi8(rawClip1, s3));

#if ALF_IMPROVEMENT
        rawCoeff0 = _mm_loadu_si128((const __m128i *) (filterSet + classIdx1 * MAX_NUM_ALF_LUMA_COEFF));
        rawCoeff1 = _mm_loadl_epi64((const __m128i *) (filterSet + classIdx1 * MAX_NUM_ALF_LUMA_COEFF + 8));

        rawClip0 = _mm_loadu_si128((const __m128i *) (fClipSet + classIdx1 * MAX_NUM_ALF_LUMA_COEFF));
        rawClip1 = _mm_loadl_epi64((const __m128i *) (fClipSet + classIdx1 * MAX_NUM_ALF_LUMA_COEFF + 8));

        const __m128i s01 = _mm_loadu_si128((const __m128i *) shuffleTab[transposeIdx1][0]);
        const __m128i s11 = _mm_xor_si128(s01, _mm_set1_epi8((char)0x80));
        const __m128i s21 = _mm_loadu_si128((const __m128i *) shuffleTab[transposeIdx1][1]);
        const __m128i s31 = _mm_xor_si128(s21, _mm_set1_epi8((char)0x80));

        const __m128i rawCoeffLo1 = _mm_or_si128(_mm_shuffle_epi8(rawCoeff0, s01), _mm_shuffle_epi8(rawCoeff1, s11));
        const __m128i rawCoeffHi1 = _mm_or_si128(_mm_shuffle_epi8(rawCoeff0, s21), _mm_shuffle_epi8(rawCoeff1, s31));
        const __m128i rawClipLo1 = _mm_or_si128(_mm_shuffle_epi8(rawClip0, s01), _mm_shuffle_epi8(rawClip1, s11));
        const __m128i rawClipHi1 = _mm_or_si128(_mm_shuffle_epi8(rawClip0, s21), _mm_shuffle_epi8(rawClip1, s31));

        params[k][0][0] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoeffLo, 0x00), _mm_shuffle_epi32(rawCoeffLo1, 0x00));
        params[k][0][1] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoeffLo, 0x55), _mm_shuffle_epi32(rawCoeffLo1, 0x55));
        params[k][0][2] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoeffLo, 0xaa), _mm_shuffle_epi32(rawCoeffLo1, 0xaa));
        params[k][0][3] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoeffLo, 0xff), _mm_shuffle_epi32(rawCoeffLo1, 0xff));
        params[k][0][4] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoeffHi, 0x00), _mm_shuffle_epi32(rawCoeffHi1, 0x00));
        params[k][0][5] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoeffHi, 0x55), _mm_shuffle_epi32(rawCoeffHi1, 0x55));
        params[k][1][0] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClipLo, 0x00), _mm_shuffle_epi32(rawClipLo1, 0x00));
        params[k][1][1] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClipLo, 0x55), _mm_shuffle_epi32(rawClipLo1, 0x55));
        params[k][1][2] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClipLo, 0xaa), _mm_shuffle_epi32(rawClipLo1, 0xaa));
        params[k][1][3] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClipLo, 0xff), _mm_shuffle_epi32(rawClipLo1, 0xff));
        params[k][1][4] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClipHi, 0x00), _mm_shuffle_epi32(rawClipHi1, 0x00));
        params[k][1][5] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClipHi, 0x55), _mm_shuffle_epi32(rawClipHi1, 0x55));
#else

        params[k][0][0] = _mm_shuffle_epi32(rawCoeffLo, 0x00);
        params[k][0][1] = _mm_shuffle_epi32(rawCoeffLo, 0x55);
        params[k][0][2] = _mm_shuffle_epi32(rawCoeffLo, 0xaa);
        params[k][0][3] = _mm_shuffle_epi32(rawCoeffLo, 0xff);
        params[k][0][4] = _mm_shuffle_epi32(rawCoeffHi, 0x00);
        params[k][0][5] = _mm_shuffle_epi32(rawCoeffHi, 0x55);
        params[k][1][0] = _mm_shuffle_epi32(rawClipLo, 0x00);
        params[k][1][1] = _mm_shuffle_epi32(rawClipLo, 0x55);
        params[k][1][2] = _mm_shuffle_epi32(rawClipLo, 0xaa);
        params[k][1][3] = _mm_shuffle_epi32(rawClipLo, 0xff);
        params[k][1][4] = _mm_shuffle_epi32(rawClipHi, 0x00);
        params[k][1][5] = _mm_shuffle_epi32(rawClipHi, 0x55);
#endif
      }

      for (size_t ii = 0; ii < STEP_Y; ii++)
      {
        const Pel *pImg0, *pImg1, *pImg2, *pImg3, *pImg4, *pImg5, *pImg6;

        pImg0 = src + j + ii * srcStride;
        pImg1 = pImg0 + srcStride;
        pImg2 = pImg0 - srcStride;
        pImg3 = pImg1 + srcStride;
        pImg4 = pImg2 - srcStride;
        pImg5 = pImg3 + srcStride;
        pImg6 = pImg4 - srcStride;
#if !ALF_IMPROVEMENT
        const int yVb = (blkDst.y + i + ii) & (vbCTUHeight - 1);
        if (yVb < vbPos && (yVb >= vbPos - 4))   // above
        {
          pImg1 = (yVb == vbPos - 1) ? pImg0 : pImg1;
          pImg3 = (yVb >= vbPos - 2) ? pImg1 : pImg3;
          pImg5 = (yVb >= vbPos - 3) ? pImg3 : pImg5;

          pImg2 = (yVb == vbPos - 1) ? pImg0 : pImg2;
          pImg4 = (yVb >= vbPos - 2) ? pImg2 : pImg4;
          pImg6 = (yVb >= vbPos - 3) ? pImg4 : pImg6;
        }
        else if (yVb >= vbPos && (yVb <= vbPos + 3))   // bottom
        {
          pImg2 = (yVb == vbPos) ? pImg0 : pImg2;
          pImg4 = (yVb <= vbPos + 1) ? pImg2 : pImg4;
          pImg6 = (yVb <= vbPos + 2) ? pImg4 : pImg6;

          pImg1 = (yVb == vbPos) ? pImg0 : pImg1;
          pImg3 = (yVb <= vbPos + 1) ? pImg1 : pImg3;
          pImg5 = (yVb <= vbPos + 2) ? pImg3 : pImg5;
        }
#endif
        __m128i cur = _mm_loadu_si128((const __m128i *) pImg0);

        __m128i accumA = mmOffset;
        __m128i accumB = mmOffset;

        auto process2coeffs = [&](const int i, const Pel *ptr0, const Pel *ptr1, const Pel *ptr2, const Pel *ptr3) {
          const __m128i val00 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) ptr0), cur);
          const __m128i val10 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) ptr2), cur);
          const __m128i val01 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) ptr1), cur);
          const __m128i val11 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) ptr3), cur);

          __m128i val01A = _mm_unpacklo_epi16(val00, val10);
          __m128i val01B = _mm_unpackhi_epi16(val00, val10);
          __m128i val01C = _mm_unpacklo_epi16(val01, val11);
          __m128i val01D = _mm_unpackhi_epi16(val01, val11);

          __m128i limit01A = params[0][1][i];
          __m128i limit01B = params[1][1][i];

          val01A = _mm_min_epi16(val01A, limit01A);
          val01B = _mm_min_epi16(val01B, limit01B);
          val01C = _mm_min_epi16(val01C, limit01A);
          val01D = _mm_min_epi16(val01D, limit01B);

          limit01A = _mm_sub_epi16(_mm_setzero_si128(), limit01A);
          limit01B = _mm_sub_epi16(_mm_setzero_si128(), limit01B);

          val01A = _mm_max_epi16(val01A, limit01A);
          val01B = _mm_max_epi16(val01B, limit01B);
          val01C = _mm_max_epi16(val01C, limit01A);
          val01D = _mm_max_epi16(val01D, limit01B);

          val01A = _mm_add_epi16(val01A, val01C);
          val01B = _mm_add_epi16(val01B, val01D);

          const __m128i coeff01A = params[0][0][i];
          const __m128i coeff01B = params[1][0][i];

          accumA = _mm_add_epi32(accumA, _mm_madd_epi16(val01A, coeff01A));
          accumB = _mm_add_epi32(accumB, _mm_madd_epi16(val01B, coeff01B));
        };


        process2coeffs(0, pImg5 + 0, pImg6 + 0, pImg3 + 1, pImg4 - 1);
        process2coeffs(1, pImg3 + 0, pImg4 + 0, pImg3 - 1, pImg4 + 1);
        process2coeffs(2, pImg1 + 2, pImg2 - 2, pImg1 + 1, pImg2 - 1);
        process2coeffs(3, pImg1 + 0, pImg2 + 0, pImg1 - 1, pImg2 + 1);
        process2coeffs(4, pImg1 - 2, pImg2 + 2, pImg0 + 3, pImg0 - 3);
        process2coeffs(5, pImg0 + 2, pImg0 - 2, pImg0 + 1, pImg0 - 1);

#if ALF_IMPROVEMENT
        accumA = _mm_srai_epi32( accumA, SHIFT );
        accumB = _mm_srai_epi32( accumB, SHIFT );
#else
        bool isNearVBabove = yVb < vbPos && (yVb >= vbPos - 1);
        bool isNearVBbelow = yVb >= vbPos && (yVb <= vbPos);
        if (!(isNearVBabove || isNearVBbelow))
        {
          accumA = _mm_srai_epi32(accumA, SHIFT);
          accumB = _mm_srai_epi32(accumB, SHIFT);
        }
        else
        {
          accumA = _mm_srai_epi32(_mm_add_epi32(accumA, mmOffset1), SHIFT + 3);
          accumB = _mm_srai_epi32(_mm_add_epi32(accumB, mmOffset1), SHIFT + 3);
        }
#endif
        accumA = _mm_packs_epi32(accumA, accumB);
        accumA = _mm_add_epi16(accumA, cur);
        accumA = _mm_min_epi16(mmMax, _mm_max_epi16(accumA, mmMin));

        _mm_storeu_si128((__m128i *) (dst + ii * dstStride + j), accumA);
      }
    }

    src += srcStride * STEP_Y;
    dst += dstStride * STEP_Y;
  }
}

#if ALF_IMPROVEMENT
static const uint16_t shuffleTime9[4] = { 0, 3, 1, 3 };
static const uint16_t shuffleOp9[4][3][2] =
{
  {
    { 0, 1 },
    { 0, 1 },
    { 0, 1 },
  }, //0
  {
    { 0, 1 },
    { 0, 2 },
    { 1, 2 },
  }, //1
  {
    { 0, 1 },
    { 0, 1 },
    { 0, 1 },
  }, //2
  {
    { 0, 1 },
    { 0, 2 },
    { 1, 2 },
  }, //3
};

static const uint16_t shuffleTab9[4][3][2][8] =
{
  {
    {
      { sh(0) , sh(1) , sh(2) , sh(3) , sh(4) , sh(5) , sh(6) , sh(7) },
      { sh(8) , sh(9) , sh(10), sh(11), sh(12), sh(13), sh(14), sh(15) },
    },
    {
      { sh(0) , sh(1) , sh(2) , sh(3) , sh(4) , sh(5) , sh(6) , sh(7) },
      { sh(8) , sh(9) , sh(10), sh(11), sh(12), sh(13), sh(14), sh(15) },
    },
    {
      { sh(0) , sh(1) , sh(2) , sh(3) , sh(4) , sh(5) , sh(6) , sh(7) },
      { sh(8) , sh(9) , sh(10), sh(11), sh(12), sh(13), sh(14), sh(15) },
    },
  }, //0
  {
    {
      { sh(0) , sh(9) , sh(2) , sh(15) , sh(4) , sh(10) , sh(6) , sh(14) },
      { sh(8) , sh(1) , sh(5), sh(11), sh(12), sh(13), sh(7), sh(3) },
    },
    {
      { sh(8) , sh(1) , sh(9), sh(3) , sh(4) , sh(5) , sh(10) , sh(7) },
      { sh(0) , sh(2) , sh(6), sh(11), sh(12), sh(13), sh(14), sh(15) },
    },
    {
      { sh(0) , sh(1) , sh(2) , sh(3), sh(11), sh(5) , sh(6) , sh(7) },
      { sh(8) , sh(9) , sh(10), sh(4), sh(12), sh(13), sh(14), sh(15) },
    },
  }, //1
  {
    {
      { sh(0) , sh(3) , sh(2) , sh(1) , sh(8) , sh(7) , sh(6) , sh(5) },
      { sh(4) , sh(15), sh(14), sh(13), sh(12), sh(11), sh(10), sh(9) },
    },
    {
      { sh(0) , sh(1) , sh(2) , sh(3) , sh(4) , sh(5) , sh(6) , sh(7) },
      { sh(8) , sh(9) , sh(10), sh(11), sh(12), sh(13), sh(14), sh(15) },
    },
    {
      { sh(0) , sh(1) , sh(2) , sh(3) , sh(4) , sh(5) , sh(6) , sh(7) },
      { sh(8) , sh(9) , sh(10), sh(11), sh(12), sh(13), sh(14), sh(15) },
    },
  }, //2
  {
    {
      { sh(0) , sh(15), sh(2), sh(9), sh(8) , sh(14), sh(6), sh(10) },
      { sh(4) , sh(3) , sh(7), sh(13), sh(1), sh(11), sh(5), sh(12) },
    },
    {
      { sh(8) , sh(1) , sh(9), sh(3) , sh(4) , sh(5) , sh(10), sh(7) },
      { sh(0) , sh(2) , sh(6), sh(11), sh(12), sh(13), sh(14), sh(15) },
    },
    {
      { sh(0) , sh(1) , sh(2) , sh(3), sh(11), sh(5) , sh(6) , sh(4) },
      { sh(8) , sh(9) , sh(10), sh(7), sh(12), sh(13), sh(14), sh(15) },
    },
  }, //3
};


template<X86_VEXT vext>
static void simdFilter9x9Blk(AlfClassifier **classifier, const PelUnitBuf &recDst, const CPelUnitBuf &recSrc, const Area &blkDst, const Area &blk, const ComponentID compId, const short *filterSet, 
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  const Pel *fClipSet
#else
  const short *fClipSet
#endif
  , const ClpRng &clpRng, CodingStructure &cs, Pel ***fixedFilterResults, int fixedFilterSetIdx)
{
  const CPelBuf srcBuffer = recSrc.get(compId);
  PelBuf        dstBuffer = recDst.get(compId);

  const size_t srcStride = srcBuffer.stride;
  const size_t dstStride = dstBuffer.stride;

  constexpr int SHIFT = AdaptiveLoopFilter::m_NUM_BITS - 1;
  constexpr int ROUND = 1 << (SHIFT - 1);

  const size_t width = blk.width;
  const size_t height = blk.height;

  constexpr size_t STEP_X = 8;
  size_t STEP_Y = isChroma(compId) ? 4 : 2;

  static_assert(sizeof(*filterSet) == 2, "ALF coeffs must be 16-bit wide");
  static_assert(sizeof(*fClipSet) == 2, "ALF clip values must be 16-bit wide");

  CHECK(blk.y % STEP_Y, "Wrong startHeight in filtering");
  CHECK(blk.x % STEP_X, "Wrong startWidth in filtering");
  CHECK(height % STEP_Y, "Wrong endHeight in filtering");
  CHECK(width % STEP_X, "Wrong endWidth in filtering");

  const Pel *src = srcBuffer.buf + blk.y * srcStride + blk.x;
  Pel *      dst = dstBuffer.buf + blkDst.y * dstStride + blkDst.x;

  const __m128i mmOffset = _mm_set1_epi32(ROUND);
  const __m128i mmMin = _mm_set1_epi16(clpRng.min);
  const __m128i mmMax = _mm_set1_epi16(clpRng.max);

  for (size_t i = 0; i < height; i += STEP_Y)
  {
    const AlfClassifier *pClass = isChroma(compId) ? nullptr : classifier[blkDst.y + i] + blkDst.x;
    for (size_t j = 0; j < width; j += STEP_X)
    {
      __m128i params[2][2][10];
      for (int k = 0; k < 2; k++)
      {
        __m128i rawCoeff[2][3], rawClip[2][3], s0, s1, s2, s3, rawTmp0, rawTmp1;

        for (int l = 0; l < 2; l++)
        {
          const int transposeIdx = pClass ? (pClass[j + 4 * k + 2 * l] & 0x3 ) : 0;
          const int classIdx = pClass ? (pClass[j + 4 * k + 2 * l] >> 2): 0;

          rawCoeff[l][0] = _mm_loadu_si128((const __m128i *) (filterSet + classIdx * MAX_NUM_ALF_LUMA_COEFF));
          rawCoeff[l][1] = _mm_loadu_si128((const __m128i *) (filterSet + classIdx * MAX_NUM_ALF_LUMA_COEFF + 8));
          rawCoeff[l][2] = _mm_loadl_epi64((const __m128i *) (filterSet + classIdx * MAX_NUM_ALF_LUMA_COEFF + 16));
          rawClip[l][0] = _mm_loadu_si128((const __m128i *) (fClipSet + classIdx * MAX_NUM_ALF_LUMA_COEFF));
          rawClip[l][1] = _mm_loadu_si128((const __m128i *) (fClipSet + classIdx * MAX_NUM_ALF_LUMA_COEFF + 8));
          rawClip[l][2] = _mm_loadl_epi64((const __m128i *) (fClipSet + classIdx * MAX_NUM_ALF_LUMA_COEFF + 16));

          for (int m = 0; m < shuffleTime9[transposeIdx]; m++)
          {
            int op0 = shuffleOp9[transposeIdx][m][0];
            int op1 = shuffleOp9[transposeIdx][m][1];

            s0 = _mm_loadu_si128((const __m128i *) shuffleTab9[transposeIdx][m][0]);
            s1 = _mm_xor_si128(s0, _mm_set1_epi8((char)0x80));
            s2 = _mm_loadu_si128((const __m128i *) shuffleTab9[transposeIdx][m][1]);
            s3 = _mm_xor_si128(s2, _mm_set1_epi8((char)0x80));

            rawTmp0 = _mm_or_si128(_mm_shuffle_epi8(rawCoeff[l][op0], s0), _mm_shuffle_epi8(rawCoeff[l][op1], s1));
            rawTmp1 = _mm_or_si128(_mm_shuffle_epi8(rawCoeff[l][op0], s2), _mm_shuffle_epi8(rawCoeff[l][op1], s3));
            rawCoeff[l][op0] = rawTmp0;
            rawCoeff[l][op1] = rawTmp1;

            rawTmp0 = _mm_or_si128(_mm_shuffle_epi8(rawClip[l][op0], s0), _mm_shuffle_epi8(rawClip[l][op1], s1));
            rawTmp1 = _mm_or_si128(_mm_shuffle_epi8(rawClip[l][op0], s2), _mm_shuffle_epi8(rawClip[l][op1], s3));
            rawClip[l][op0] = rawTmp0;
            rawClip[l][op1] = rawTmp1;
          }
        }

        params[k][0][0] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoeff[0][0], 0x00), _mm_shuffle_epi32(rawCoeff[1][0], 0x00));
        params[k][0][1] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoeff[0][0], 0x55), _mm_shuffle_epi32(rawCoeff[1][0], 0x55));
        params[k][0][2] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoeff[0][0], 0xaa), _mm_shuffle_epi32(rawCoeff[1][0], 0xaa));
        params[k][0][3] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoeff[0][0], 0xff), _mm_shuffle_epi32(rawCoeff[1][0], 0xff));
        params[k][0][4] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoeff[0][1], 0x00), _mm_shuffle_epi32(rawCoeff[1][1], 0x00));
        params[k][0][5] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoeff[0][1], 0x55), _mm_shuffle_epi32(rawCoeff[1][1], 0x55));
        params[k][0][6] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoeff[0][1], 0xaa), _mm_shuffle_epi32(rawCoeff[1][1], 0xaa));
        params[k][0][7] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoeff[0][1], 0xff), _mm_shuffle_epi32(rawCoeff[1][1], 0xff));
        params[k][0][8] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoeff[0][2], 0x00), _mm_shuffle_epi32(rawCoeff[1][2], 0x00));
        params[k][0][9] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoeff[0][2], 0x55), _mm_shuffle_epi32(rawCoeff[1][2], 0x55));

        params[k][1][0] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClip[0][0], 0x00), _mm_shuffle_epi32(rawClip[1][0], 0x00));
        params[k][1][1] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClip[0][0], 0x55), _mm_shuffle_epi32(rawClip[1][0], 0x55));
        params[k][1][2] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClip[0][0], 0xaa), _mm_shuffle_epi32(rawClip[1][0], 0xaa));
        params[k][1][3] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClip[0][0], 0xff), _mm_shuffle_epi32(rawClip[1][0], 0xff));
        params[k][1][4] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClip[0][1], 0x00), _mm_shuffle_epi32(rawClip[1][1], 0x00));
        params[k][1][5] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClip[0][1], 0x55), _mm_shuffle_epi32(rawClip[1][1], 0x55));
        params[k][1][6] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClip[0][1], 0xaa), _mm_shuffle_epi32(rawClip[1][1], 0xaa));
        params[k][1][7] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClip[0][1], 0xff), _mm_shuffle_epi32(rawClip[1][1], 0xff));
        params[k][1][8] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClip[0][2], 0x00), _mm_shuffle_epi32(rawClip[1][2], 0x00));
        params[k][1][9] = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClip[0][2], 0x55), _mm_shuffle_epi32(rawClip[1][2], 0x55));
      }

      for (size_t ii = 0; ii < STEP_Y; ii++)
      {
        const Pel *pImg0, *pImg1, *pImg2, *pImg3, *pImg4, *pImg5, *pImg6, *pImg7, *pImg8;

        pImg0 = src + j + ii * srcStride;
        pImg1 = pImg0 + srcStride;
        pImg2 = pImg0 - srcStride;
        pImg3 = pImg1 + srcStride;
        pImg4 = pImg2 - srcStride;
        pImg5 = pImg3 + srcStride;
        pImg6 = pImg4 - srcStride;
        pImg7 = pImg5 + srcStride;
        pImg8 = pImg6 - srcStride;

        __m128i cur = _mm_loadu_si128((const __m128i *) pImg0);
        __m128i accumA = mmOffset;
        __m128i accumB = mmOffset;

        auto process2coeffs = [&](const int i, const Pel *ptr0, const Pel *ptr1, const Pel *ptr2, const Pel *ptr3) {
          const __m128i val00 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) ptr0), cur);
          const __m128i val10 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) ptr2), cur);
          const __m128i val01 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) ptr1), cur);
          const __m128i val11 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) ptr3), cur);

          __m128i val01A = _mm_unpacklo_epi16(val00, val10);
          __m128i val01B = _mm_unpackhi_epi16(val00, val10);
          __m128i val01C = _mm_unpacklo_epi16(val01, val11);
          __m128i val01D = _mm_unpackhi_epi16(val01, val11);

          __m128i limit01A = params[0][1][i];
          __m128i limit01B = params[1][1][i];

          val01A = _mm_min_epi16(val01A, limit01A);
          val01B = _mm_min_epi16(val01B, limit01B);
          val01C = _mm_min_epi16(val01C, limit01A);
          val01D = _mm_min_epi16(val01D, limit01B);

          limit01A = _mm_sub_epi16(_mm_setzero_si128(), limit01A);
          limit01B = _mm_sub_epi16(_mm_setzero_si128(), limit01B);

          val01A = _mm_max_epi16(val01A, limit01A);
          val01B = _mm_max_epi16(val01B, limit01B);
          val01C = _mm_max_epi16(val01C, limit01A);
          val01D = _mm_max_epi16(val01D, limit01B);

          val01A = _mm_add_epi16(val01A, val01C);
          val01B = _mm_add_epi16(val01B, val01D);

          const __m128i coeff01A = params[0][0][i];
          const __m128i coeff01B = params[1][0][i];

          accumA = _mm_add_epi32(accumA, _mm_madd_epi16(val01A, coeff01A));
          accumB = _mm_add_epi32(accumB, _mm_madd_epi16(val01B, coeff01B));
        };

        process2coeffs(0, pImg7 + 0, pImg8 + 0, pImg5 + 1, pImg6 - 1);
        process2coeffs(1, pImg5 + 0, pImg6 + 0, pImg5 - 1, pImg6 + 1);
        process2coeffs(2, pImg3 + 2, pImg4 - 2, pImg3 + 1, pImg4 - 1);
        process2coeffs(3, pImg3 + 0, pImg4 + 0, pImg3 - 1, pImg4 + 1);
        process2coeffs(4, pImg3 - 2, pImg4 + 2, pImg1 + 3, pImg2 - 3);
        process2coeffs(5, pImg1 + 2, pImg2 - 2, pImg1 + 1, pImg2 - 1);
        process2coeffs(6, pImg1 + 0, pImg2 + 0, pImg1 - 1, pImg2 + 1);
        process2coeffs(7, pImg1 - 2, pImg2 + 2, pImg1 - 3, pImg2 + 3);
        process2coeffs(8, pImg0 + 4, pImg0 - 4, pImg0 + 3, pImg0 - 3);
        process2coeffs(9, pImg0 + 2, pImg0 - 2, pImg0 + 1, pImg0 - 1);

        accumA = _mm_srai_epi32(accumA, SHIFT);
        accumB = _mm_srai_epi32(accumB, SHIFT);

        accumA = _mm_packs_epi32(accumA, accumB);
        accumA = _mm_add_epi16(accumA, cur);
        accumA = _mm_min_epi16(mmMax, _mm_max_epi16(accumA, mmMin));

        _mm_storeu_si128((__m128i *) (dst + ii * dstStride + j), accumA);
      } //for (size_t ii = 0; ii < STEP_Y; ii++)
    } //for (size_t j = 0; j < width; j += STEP_X)
    src += srcStride * STEP_Y;
    dst += dstStride * STEP_Y;
  }
}
#endif

#if ALF_IMPROVEMENT
template<X86_VEXT vext>
static void simdFilter9x9BlkExt(AlfClassifier **classifier, const PelUnitBuf &recDst, const CPelUnitBuf &recSrc, const Area &blkDst, const Area &blk, const ComponentID compId, const short *filterSet, 
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  const Pel *fClipSet
#else
  const short *fClipSet
#endif
  , const ClpRng &clpRng, CodingStructure &cs, Pel ***fixedFilterResults, int fixedFilterSetIdx)
{
  const CPelBuf srcBuffer = recSrc.get(compId);
  PelBuf        dstBuffer = recDst.get(compId);

  const size_t srcStride = srcBuffer.stride;
  const size_t dstStride = dstBuffer.stride;

  constexpr int SHIFT = AdaptiveLoopFilter::m_NUM_BITS - 1;
  constexpr int ROUND = 1 << (SHIFT - 1);

  const size_t width = blk.width;
  const size_t height = blk.height;

  constexpr size_t STEP_X = 8;
  size_t STEP_Y = 1;

  const __m128i mmOffset = _mm_set1_epi32(ROUND);
  const __m128i mmMin = _mm_set1_epi16(clpRng.min);
  const __m128i mmMax = _mm_set1_epi16(clpRng.max);

  static_assert(sizeof(*filterSet) == 2, "ALF coeffs must be 16-bit wide");
  static_assert(sizeof(*fClipSet) == 2, "ALF clip values must be 16-bit wide");

  const Pel *src = srcBuffer.buf + blk.y * srcStride + blk.x;
  Pel *      dst = dstBuffer.buf + blkDst.y * dstStride + blkDst.x;

  for (size_t i = 0; i < height; i += STEP_Y)
  {
    const AlfClassifier *pClass = classifier[blkDst.y + i] + blkDst.x;
    for (size_t j = 0; j < width; j += STEP_X)
    {
      __m128i params[2][2][11];
      for (int k = 0; k < 2; k++)
      {
        __m128i rawCoef[4][3], rawClip[4][3], s0, s1, s2, s3, rawTmp0, rawTmp1;
        for (int l = 0; l < 4; l++)
        {
          const int transposeIdx = pClass[j + 4 * k + l] & 0x3;
          const int classIdx = pClass[j + 4 * k + l] >> 2;

          rawCoef[l][0] = _mm_loadu_si128((const __m128i *) (filterSet + classIdx * MAX_NUM_ALF_LUMA_COEFF));
          rawCoef[l][1] = _mm_loadu_si128((const __m128i *) (filterSet + classIdx * MAX_NUM_ALF_LUMA_COEFF + 8));
          rawCoef[l][2] = _mm_loadu_si128((const __m128i *) (filterSet + classIdx * MAX_NUM_ALF_LUMA_COEFF + 16));
          rawClip[l][0] = _mm_loadu_si128((const __m128i *) (fClipSet + classIdx * MAX_NUM_ALF_LUMA_COEFF));
          rawClip[l][1] = _mm_loadu_si128((const __m128i *) (fClipSet + classIdx * MAX_NUM_ALF_LUMA_COEFF + 8));
          rawClip[l][2] = _mm_loadu_si128((const __m128i *) (fClipSet + classIdx * MAX_NUM_ALF_LUMA_COEFF + 16));

          for (int m = 0; m < shuffleTime9[transposeIdx]; m++)
          {
            int op0 = shuffleOp9[transposeIdx][m][0];
            int op1 = shuffleOp9[transposeIdx][m][1];

            s0 = _mm_loadu_si128((const __m128i *) shuffleTab9[transposeIdx][m][0]);
            s1 = _mm_xor_si128(s0, _mm_set1_epi8((char)0x80));
            s2 = _mm_loadu_si128((const __m128i *) shuffleTab9[transposeIdx][m][1]);
            s3 = _mm_xor_si128(s2, _mm_set1_epi8((char)0x80));

            rawTmp0 = _mm_or_si128(_mm_shuffle_epi8(rawCoef[l][op0], s0), _mm_shuffle_epi8(rawCoef[l][op1], s1));
            rawTmp1 = _mm_or_si128(_mm_shuffle_epi8(rawCoef[l][op0], s2), _mm_shuffle_epi8(rawCoef[l][op1], s3));
            rawCoef[l][op0] = rawTmp0;
            rawCoef[l][op1] = rawTmp1;

            rawTmp0 = _mm_or_si128(_mm_shuffle_epi8(rawClip[l][op0], s0), _mm_shuffle_epi8(rawClip[l][op1], s1));
            rawTmp1 = _mm_or_si128(_mm_shuffle_epi8(rawClip[l][op0], s2), _mm_shuffle_epi8(rawClip[l][op1], s3));
            rawClip[l][op0] = rawTmp0;
            rawClip[l][op1] = rawTmp1;
          }            
        }//for l

        for (unsigned char l = 0; l < 3; l++)
        {
          int m = l << 2;

          s0 = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoef[0][l], 0x00), _mm_shuffle_epi32(rawCoef[1][l], 0x00));
          s1 = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoef[2][l], 0x00), _mm_shuffle_epi32(rawCoef[3][l], 0x00));
          params[k][0][0 + m] = _mm_blend_epi16(_mm_shuffle_epi32(s0, 0x88), _mm_shuffle_epi32(s1, 0x88), 0xf0);
          s0 = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClip[0][l], 0x00), _mm_shuffle_epi32(rawClip[1][l], 0x00));
          s1 = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClip[2][l], 0x00), _mm_shuffle_epi32(rawClip[3][l], 0x00));
          params[k][1][0 + m] = _mm_blend_epi16(_mm_shuffle_epi32(s0, 0x88), _mm_shuffle_epi32(s1, 0x88), 0xf0);

          s0 = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoef[0][l], 0x55), _mm_shuffle_epi32(rawCoef[1][l], 0x55));
          s1 = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoef[2][l], 0x55), _mm_shuffle_epi32(rawCoef[3][l], 0x55));
          params[k][0][1 + m] = _mm_blend_epi16(_mm_shuffle_epi32(s0, 0x88), _mm_shuffle_epi32(s1, 0x88), 0xf0);
          s0 = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClip[0][l], 0x55), _mm_shuffle_epi32(rawClip[1][l], 0x55));
          s1 = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClip[2][l], 0x55), _mm_shuffle_epi32(rawClip[3][l], 0x55));
          params[k][1][1 + m] = _mm_blend_epi16(_mm_shuffle_epi32(s0, 0x88), _mm_shuffle_epi32(s1, 0x88), 0xf0);

          s0 = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoef[0][l], 0xaa), _mm_shuffle_epi32(rawCoef[1][l], 0xaa));
          s1 = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoef[2][l], 0xaa), _mm_shuffle_epi32(rawCoef[3][l], 0xaa));
          params[k][0][2 + m] = _mm_blend_epi16(_mm_shuffle_epi32(s0, 0x88), _mm_shuffle_epi32(s1, 0x88), 0xf0);
          s0 = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClip[0][l], 0xaa), _mm_shuffle_epi32(rawClip[1][l], 0xaa));
          s1 = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClip[2][l], 0xaa), _mm_shuffle_epi32(rawClip[3][l], 0xaa));
          params[k][1][2 + m] = _mm_blend_epi16(_mm_shuffle_epi32(s0, 0x88), _mm_shuffle_epi32(s1, 0x88), 0xf0);

          if (l < 2)
          {
            s0 = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoef[0][l], 0xff), _mm_shuffle_epi32(rawCoef[1][l], 0xff));
            s1 = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawCoef[2][l], 0xff), _mm_shuffle_epi32(rawCoef[3][l], 0xff));
            params[k][0][3 + m] = _mm_blend_epi16(_mm_shuffle_epi32(s0, 0x88), _mm_shuffle_epi32(s1, 0x88), 0xf0);
            s0 = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClip[0][l], 0xff), _mm_shuffle_epi32(rawClip[1][l], 0xff));
            s1 = _mm_unpacklo_epi64(_mm_shuffle_epi32(rawClip[2][l], 0xff), _mm_shuffle_epi32(rawClip[3][l], 0xff));
            params[k][1][3 + m] = _mm_blend_epi16(_mm_shuffle_epi32(s0, 0x88), _mm_shuffle_epi32(s1, 0x88), 0xf0);
          }
        }//for l
      }//for k

      const Pel *pImg0, *pImg1, *pImg2, *pImg3, *pImg4, *pImg5, *pImg6, *pImg7, *pImg8;
      pImg0 = src + j;
      pImg1 = pImg0 + srcStride;
      pImg2 = pImg0 - srcStride;
      pImg3 = pImg1 + srcStride;
      pImg4 = pImg2 - srcStride;
      pImg5 = pImg3 + srcStride;
      pImg6 = pImg4 - srcStride;
      pImg7 = pImg5 + srcStride;
      pImg8 = pImg6 - srcStride;

      __m128i cur = _mm_loadu_si128((const __m128i *) pImg0);
      __m128i accumA = mmOffset;
      __m128i accumB = mmOffset;

      auto process2coeffs = [&](const int i, const Pel *ptr0, const Pel *ptr1, const Pel *ptr2, const Pel *ptr3) {
        const __m128i val00 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) ptr0), cur);
        const __m128i val10 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) ptr2), cur);
        const __m128i val01 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) ptr1), cur);
        const __m128i val11 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) ptr3), cur);

        __m128i val01A = _mm_unpacklo_epi16(val00, val10);
        __m128i val01B = _mm_unpackhi_epi16(val00, val10);
        __m128i val01C = _mm_unpacklo_epi16(val01, val11);
        __m128i val01D = _mm_unpackhi_epi16(val01, val11);

        __m128i limit01A = params[0][1][i];
        __m128i limit01B = params[1][1][i];

        val01A = _mm_min_epi16(val01A, limit01A);
        val01B = _mm_min_epi16(val01B, limit01B);
        val01C = _mm_min_epi16(val01C, limit01A);
        val01D = _mm_min_epi16(val01D, limit01B);

        limit01A = _mm_sub_epi16(_mm_setzero_si128(), limit01A);
        limit01B = _mm_sub_epi16(_mm_setzero_si128(), limit01B);

        val01A = _mm_max_epi16(val01A, limit01A);
        val01B = _mm_max_epi16(val01B, limit01B);
        val01C = _mm_max_epi16(val01C, limit01A);
        val01D = _mm_max_epi16(val01D, limit01B);

        val01A = _mm_add_epi16(val01A, val01C);
        val01B = _mm_add_epi16(val01B, val01D);

        const __m128i coeff01A = params[0][0][i];
        const __m128i coeff01B = params[1][0][i];

        accumA = _mm_add_epi32(accumA, _mm_madd_epi16(val01A, coeff01A));
        accumB = _mm_add_epi32(accumB, _mm_madd_epi16(val01B, coeff01B));
      };

      process2coeffs(0, pImg7 + 0, pImg8 + 0, pImg5 + 1, pImg6 - 1);
      process2coeffs(1, pImg5 + 0, pImg6 + 0, pImg5 - 1, pImg6 + 1);
      process2coeffs(2, pImg3 + 2, pImg4 - 2, pImg3 + 1, pImg4 - 1);
      process2coeffs(3, pImg3 + 0, pImg4 + 0, pImg3 - 1, pImg4 + 1);
      process2coeffs(4, pImg3 - 2, pImg4 + 2, pImg1 + 3, pImg2 - 3);
      process2coeffs(5, pImg1 + 2, pImg2 - 2, pImg1 + 1, pImg2 - 1);
      process2coeffs(6, pImg1 + 0, pImg2 + 0, pImg1 - 1, pImg2 + 1);
      process2coeffs(7, pImg1 - 2, pImg2 + 2, pImg1 - 3, pImg2 + 3);
      process2coeffs(8, pImg0 + 4, pImg0 - 4, pImg0 + 3, pImg0 - 3);
      process2coeffs(9, pImg0 + 2, pImg0 - 2, pImg0 + 1, pImg0 - 1);
      const __m128i val00 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) (fixedFilterResults[fixedFilterSetIdx][blkDst.y + i] + blkDst.x + j)), cur);
      const __m128i val10 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) (fixedFilterResults[EXT_LENGTH + fixedFilterSetIdx][blkDst.y + i] + blkDst.x + j)), cur);
      __m128i val01A = _mm_unpacklo_epi16(val00, val10);
      __m128i val01B = _mm_unpackhi_epi16(val00, val10);
      __m128i limit01A = params[0][1][10];
      __m128i limit01B = params[1][1][10];
      val01A = _mm_min_epi16(val01A, limit01A);
      val01B = _mm_min_epi16(val01B, limit01B);
      limit01A = _mm_sub_epi16(_mm_setzero_si128(), limit01A);
      limit01B = _mm_sub_epi16(_mm_setzero_si128(), limit01B);
      val01A = _mm_max_epi16(val01A, limit01A);
      val01B = _mm_max_epi16(val01B, limit01B);
      const __m128i coeff01A = params[0][0][10];
      const __m128i coeff01B = params[1][0][10];
      accumA = _mm_add_epi32(accumA, _mm_madd_epi16(val01A, coeff01A));
      accumB = _mm_add_epi32(accumB, _mm_madd_epi16(val01B, coeff01B));

      accumA = _mm_srai_epi32(accumA, SHIFT);
      accumB = _mm_srai_epi32(accumB, SHIFT);

      accumA = _mm_packs_epi32(accumA, accumB);
      accumA = _mm_add_epi16(accumA, cur);
      accumA = _mm_min_epi16(mmMax, _mm_max_epi16(accumA, mmMin));

      _mm_storeu_si128((__m128i *) (dst + j), accumA);
    }//for j
    src += srcStride * STEP_Y;
    dst += dstStride * STEP_Y;
  }//for i
}

static const int8_t shTab[4][6][16] = {
  {
    {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
    {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
    {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15 },
    {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15 },
    {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15 },
    {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15 }
  },
  {
    { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
    {  8,  9,  6,  7,  4,  5,  2,  3,  0,  1, 14, 15, 12, 13, 10, 11 },
    {  8,  9,  6,  7,  4,  5,  2,  3,  0,  1, 14, 15, 12, 13, 10, 11 }, 
    {  6,  7,  4,  5,  2,  3,  0,  1, 10, 11,  8,  9, 12, 13, 14, 15 }, 
    {  0,  1,  2,  3, 10, 11,  8,  9,  6,  7,  4,  5, 14, 15, 12, 13 } 
  },
  {
    {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
    { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15 },
    {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15 }, 
    { 14, 15, 12, 13,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11 }, 
    {  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,  2,  3,  0,  1 } 
  },
  {
    { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
    {  8,  9,  6,  7,  4,  5,  2,  3,  0,  1, 14, 15, 12, 13, 10, 11 },
    {  8,  9,  6,  7,  4,  5,  2,  3,  0,  1, 14, 15, 12, 13, 10, 11 }, 
    { 14, 15, 12, 13,  6,  7,  4,  5,  2,  3,  0,  1, 10, 11,  8,  9 },
    { 10, 11,  8,  9,  6,  7,  4,  5, 14, 15, 12, 13,  2,  3,  0,  1 }
  }
};

template<X86_VEXT vext> 
static void simdFilter13x13Blk(
  AlfClassifier **classifier,
  const CPelBuf &srcLuma,
  const Area& curBlk,
  Pel ***fixedFilterResults,
  int picWidth,
  const int fixedFiltInd,
  const short classIndFixed[NUM_CLASSES_FIX],
  int fixedFiltQpInd, int dirWindSize,
  const ClpRng &clpRng,
  const Pel clippingValues[4]
)
{
  const int srcStride = srcLuma.stride;

  constexpr int SHIFT = AdaptiveLoopFilter::m_NUM_BITS_FIXED_FILTER - 1;
  constexpr int ROUND = 1 << (SHIFT - 1);

  const int width = curBlk.width;
  const int height = curBlk.height;

  constexpr int STEP_X = 8;
  constexpr int STEP_Y = 2;

  const Pel *src = srcLuma.buf + curBlk.y * srcStride + curBlk.x;
  const __m128i mmOffset = _mm_set1_epi32(ROUND);
  const __m128i mmMin = _mm_set1_epi16(clpRng.min);
  const __m128i mmMax = _mm_set1_epi16(clpRng.max);
  const __m128i mmClippingValues = _mm_loadl_epi64((const __m128i *)clippingValues);
  const __m128i mm11 = _mm_set1_epi8(1);
  const __m128i mm3 = _mm_set1_epi16(3);

  const std::array<std::array<short, 42>, NUM_FIXED_FILTERS>& filterCoeffFixed = packedDataFixedFilters[fixedFiltQpInd][dirWindSize];

  for (int i = 0; i < height; i += STEP_Y)
  {
    const AlfClassifier *pClass = classifier[curBlk.y + i] + curBlk.x;
    for (int j = 0; j < width; j += STEP_X)
    {
      __m128i params[21];
      __m128i rawCoef[4][6];
      for (int m = 0; m < 4; m++)
      {
        int transposeIdx = pClass[j + 2*m] & 0x3;
        const int filterIdx = classIndFixed[pClass[j + 2*m] >> 2];

        rawCoef[m][0] = _mm_loadu_si128((const __m128i*)(filterCoeffFixed[filterIdx].data()     ));
        rawCoef[m][1] = _mm_loadu_si128((const __m128i*)(filterCoeffFixed[filterIdx].data() +  6));
        rawCoef[m][2] = _mm_loadu_si128((const __m128i*)(filterCoeffFixed[filterIdx].data() + 12));
        rawCoef[m][3] = _mm_loadu_si128((const __m128i*)(filterCoeffFixed[filterIdx].data() + 20));
        rawCoef[m][4] = _mm_loadu_si128((const __m128i*)(filterCoeffFixed[filterIdx].data() + 28));
        rawCoef[m][5] = _mm_loadu_si128((const __m128i*)(filterCoeffFixed[filterIdx].data() + 34));

        //transpose
        {
          const __m128i s0 = _mm_loadu_si128((const __m128i*)shTab[transposeIdx][0]);
          const __m128i s1 = _mm_loadu_si128((const __m128i*)shTab[transposeIdx][1]);
          const __m128i s2 = _mm_loadu_si128((const __m128i*)shTab[transposeIdx][2]);
          const __m128i s3 = _mm_loadu_si128((const __m128i*)shTab[transposeIdx][3]);
          const __m128i s4 = _mm_loadu_si128((const __m128i*)shTab[transposeIdx][4]);
          const __m128i s5 = _mm_loadu_si128((const __m128i*)shTab[transposeIdx][5]);

          __m128i rawTmp[6];
          rawTmp[0] = rawCoef[m][0];
          rawTmp[1] = rawCoef[m][1];
          rawTmp[2] = _mm_shuffle_epi8(rawCoef[m][2], s2);
          rawTmp[3] = _mm_shuffle_epi8(rawCoef[m][3], s3);
          rawTmp[4] = _mm_shuffle_epi8(rawCoef[m][4], s4);
          rawTmp[5] = _mm_shuffle_epi8(rawCoef[m][5], s5);

          rawCoef[m][0] = _mm_add_epi16(rawTmp[0], _mm_and_si128(s0, _mm_sub_epi16(rawTmp[1], rawTmp[0])));
          rawCoef[m][1] = _mm_sub_epi16(rawTmp[1], _mm_and_si128(s0, _mm_sub_epi16(rawTmp[1], rawTmp[0])));
          rawCoef[m][2] = _mm_add_epi16(rawTmp[2], _mm_and_si128(s1, _mm_sub_epi16(rawTmp[3], rawTmp[2])));
          rawCoef[m][3] = _mm_sub_epi16(rawTmp[3], _mm_and_si128(s1, _mm_sub_epi16(rawTmp[3], rawTmp[2])));
          rawCoef[m][4] = _mm_add_epi16(rawTmp[4], _mm_and_si128(s1, _mm_sub_epi16(rawTmp[5], rawTmp[4])));
          rawCoef[m][5] = _mm_sub_epi16(rawTmp[5], _mm_and_si128(s1, _mm_sub_epi16(rawTmp[5], rawTmp[4])));
        }
      }//for(m)

      params[ 0] = _mm_unpacklo_epi64(_mm_unpacklo_epi32(rawCoef[0][0], rawCoef[1][0]), _mm_unpacklo_epi32(rawCoef[2][0], rawCoef[3][0]));
      params[ 1] = _mm_unpackhi_epi64(_mm_unpacklo_epi32(rawCoef[0][0], rawCoef[1][0]), _mm_unpacklo_epi32(rawCoef[2][0], rawCoef[3][0]));
      params[ 2] = _mm_unpacklo_epi64(_mm_unpackhi_epi32(rawCoef[0][0], rawCoef[1][0]), _mm_unpackhi_epi32(rawCoef[2][0], rawCoef[3][0]));

      params[ 3] = _mm_unpacklo_epi64(_mm_unpacklo_epi32(rawCoef[0][1], rawCoef[1][1]), _mm_unpacklo_epi32(rawCoef[2][1], rawCoef[3][1]));
      params[ 4] = _mm_unpackhi_epi64(_mm_unpacklo_epi32(rawCoef[0][1], rawCoef[1][1]), _mm_unpacklo_epi32(rawCoef[2][1], rawCoef[3][1]));
      params[ 5] = _mm_unpacklo_epi64(_mm_unpackhi_epi32(rawCoef[0][1], rawCoef[1][1]), _mm_unpackhi_epi32(rawCoef[2][1], rawCoef[3][1]));

      params[ 6] = _mm_unpacklo_epi64(_mm_unpacklo_epi32(rawCoef[0][2], rawCoef[1][2]), _mm_unpacklo_epi32(rawCoef[2][2], rawCoef[3][2]));
      params[ 7] = _mm_unpackhi_epi64(_mm_unpacklo_epi32(rawCoef[0][2], rawCoef[1][2]), _mm_unpacklo_epi32(rawCoef[2][2], rawCoef[3][2]));
      params[ 8] = _mm_unpacklo_epi64(_mm_unpackhi_epi32(rawCoef[0][2], rawCoef[1][2]), _mm_unpackhi_epi32(rawCoef[2][2], rawCoef[3][2]));
      params[ 9] = _mm_unpackhi_epi64(_mm_unpackhi_epi32(rawCoef[0][2], rawCoef[1][2]), _mm_unpackhi_epi32(rawCoef[2][2], rawCoef[3][2]));

      params[10] = _mm_unpacklo_epi64(_mm_unpacklo_epi32(rawCoef[0][3], rawCoef[1][3]), _mm_unpacklo_epi32(rawCoef[2][3], rawCoef[3][3]));
      params[11] = _mm_unpackhi_epi64(_mm_unpacklo_epi32(rawCoef[0][3], rawCoef[1][3]), _mm_unpacklo_epi32(rawCoef[2][3], rawCoef[3][3]));
      params[12] = _mm_unpacklo_epi64(_mm_unpackhi_epi32(rawCoef[0][3], rawCoef[1][3]), _mm_unpackhi_epi32(rawCoef[2][3], rawCoef[3][3]));
      params[13] = _mm_unpackhi_epi64(_mm_unpackhi_epi32(rawCoef[0][3], rawCoef[1][3]), _mm_unpackhi_epi32(rawCoef[2][3], rawCoef[3][3]));

      params[14] = _mm_unpacklo_epi64(_mm_unpacklo_epi32(rawCoef[0][4], rawCoef[1][4]), _mm_unpacklo_epi32(rawCoef[2][4], rawCoef[3][4]));
      params[15] = _mm_unpackhi_epi64(_mm_unpacklo_epi32(rawCoef[0][4], rawCoef[1][4]), _mm_unpacklo_epi32(rawCoef[2][4], rawCoef[3][4]));
      params[16] = _mm_unpacklo_epi64(_mm_unpackhi_epi32(rawCoef[0][4], rawCoef[1][4]), _mm_unpackhi_epi32(rawCoef[2][4], rawCoef[3][4]));

      params[17] = _mm_unpackhi_epi64(_mm_unpackhi_epi32(rawCoef[0][4], rawCoef[1][4]), _mm_unpackhi_epi32(rawCoef[2][4], rawCoef[3][4]));

      params[18] = _mm_unpackhi_epi64(_mm_unpacklo_epi32(rawCoef[0][5], rawCoef[1][5]), _mm_unpacklo_epi32(rawCoef[2][5], rawCoef[3][5]));
      params[19] = _mm_unpacklo_epi64(_mm_unpackhi_epi32(rawCoef[0][5], rawCoef[1][5]), _mm_unpackhi_epi32(rawCoef[2][5], rawCoef[3][5]));
      params[20] = _mm_unpackhi_epi64(_mm_unpackhi_epi32(rawCoef[0][5], rawCoef[1][5]), _mm_unpackhi_epi32(rawCoef[2][5], rawCoef[3][5]));

      for (int ii = 0; ii < STEP_Y; ii++)
      {
        const Pel *pImg0, *pImg1, *pImg2, *pImg3, *pImg4, *pImg5, *pImg6, *pImg7, *pImg8, *pImg9, *pImg10, *pImg11, *pImg12;
        pImg0 = src + j + ii * srcStride;
        pImg1 = pImg0 + srcStride;
        pImg2 = pImg0 - srcStride;
        pImg3 = pImg1 + srcStride;
        pImg4 = pImg2 - srcStride;
        pImg5 = pImg3 + srcStride;
        pImg6 = pImg4 - srcStride;
        pImg7 = pImg5 + srcStride;
        pImg8 = pImg6 - srcStride;
        pImg9 = pImg7 + srcStride;
        pImg10 = pImg8 - srcStride;
        pImg11 = pImg9 + srcStride;
        pImg12 = pImg10 - srcStride;

        __m128i cur = _mm_loadu_si128((const __m128i *) pImg0);
        __m128i accumA = mmOffset;
        __m128i accumB = mmOffset;

        auto process2coeffs = [&](const int i, const Pel *ptr0, const Pel *ptr1, const Pel *ptr2, const Pel *ptr3) {
          const __m128i val00 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) ptr0), cur);
          const __m128i val10 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) ptr2), cur);
          const __m128i val01 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) ptr1), cur);
          const __m128i val11 = _mm_sub_epi16(_mm_loadu_si128((const __m128i *) ptr3), cur);

          __m128i val01A = _mm_blend_epi16(val00, _mm_slli_si128(val10, 2), 0xAA);
          __m128i val01B = _mm_blend_epi16(_mm_srli_si128(val00, 2), val10, 0xAA);
          __m128i val01C = _mm_blend_epi16(val01, _mm_slli_si128(val11, 2), 0xAA);
          __m128i val01D = _mm_blend_epi16(_mm_srli_si128(val01, 2), val11, 0xAA);

          __m128i mmClippingFixed = _mm_and_si128(params[i], mm3);
          
          __m128i mmClippingFixed2 = _mm_packs_epi16(mmClippingFixed, mmClippingFixed);
          mmClippingFixed2 = _mm_add_epi8(mmClippingFixed2, mmClippingFixed2);
          __m128i xmm2 = _mm_add_epi8(mmClippingFixed2, mm11);
          __m128i xmmA = _mm_unpacklo_epi8(mmClippingFixed2, xmm2);
          __m128i limit = _mm_shuffle_epi8(mmClippingValues, xmmA);

          val01A = _mm_min_epi16(val01A, limit);
          val01B = _mm_min_epi16(val01B, limit);
          val01C = _mm_min_epi16(val01C, limit);
          val01D = _mm_min_epi16(val01D, limit);

          limit = _mm_sub_epi16(_mm_setzero_si128(), limit);

          val01A = _mm_max_epi16(val01A, limit);
          val01B = _mm_max_epi16(val01B, limit);
          val01C = _mm_max_epi16(val01C, limit);
          val01D = _mm_max_epi16(val01D, limit);

          val01A = _mm_add_epi16(val01A, val01C);
          val01B = _mm_add_epi16(val01B, val01D);

          const __m128i coeff = _mm_srai_epi16(params[i], 2);

          accumA = _mm_add_epi32(accumA, _mm_madd_epi16(val01A, coeff));
          accumB = _mm_add_epi32(accumB, _mm_madd_epi16(val01B, coeff));
        };

        process2coeffs(0, pImg11 + 0, pImg12 + 0, pImg9 + 0, pImg10 - 0);
        process2coeffs(1, pImg7 + 0, pImg8 + 0, pImg5 - 0, pImg6 + 0);
        process2coeffs(2, pImg3 + 0, pImg4 - 0, pImg1 + 0, pImg2 - 0);

        process2coeffs(3, pImg0 + 6, pImg0 - 6, pImg0 + 5, pImg0 - 5);
        process2coeffs(4, pImg0 + 4, pImg0 - 4, pImg0 + 3, pImg0 - 3);
        process2coeffs(5, pImg0 + 2, pImg0 - 2, pImg0 + 1, pImg0 - 1);

        process2coeffs(6, pImg9 + 1, pImg10 - 1, pImg7 + 2, pImg8 - 2);
        process2coeffs(7, pImg5 + 3, pImg6 - 3, pImg3 + 4, pImg4 - 4);
        process2coeffs(8, pImg1 + 5, pImg2 - 5, pImg5 + 1, pImg6 - 1);
        process2coeffs(9, pImg3 + 2, pImg4 - 2, pImg1 + 3, pImg2 - 3);
        
        process2coeffs(10, pImg9 - 1, pImg10 + 1, pImg7 - 2, pImg8 + 2);
        process2coeffs(11, pImg5 - 3, pImg6 + 3, pImg3 - 4, pImg4 + 4);
        process2coeffs(12, pImg1 - 5, pImg2 + 5, pImg5 - 1, pImg6 + 1);
        process2coeffs(13, pImg3 - 2, pImg4 + 2, pImg1 - 3, pImg2 + 3);
        
        process2coeffs(14, pImg7 + 1, pImg8 - 1, pImg5 + 2, pImg6 - 2);
        process2coeffs(15, pImg3 + 3, pImg4 - 3, pImg1 + 4, pImg2 - 4);
        process2coeffs(16, pImg3 + 1, pImg4 - 1, pImg1 + 2, pImg2 - 2);

        process2coeffs(17, pImg1 + 1, pImg2 - 1, pImg1 - 1, pImg2 + 1);
        
        process2coeffs(18, pImg7 - 1, pImg8 + 1, pImg5 - 2, pImg6 + 2);
        process2coeffs(19, pImg3 - 3, pImg4 + 3, pImg1 - 4, pImg2 + 4);
        process2coeffs(20, pImg3 - 1, pImg4 + 1, pImg1 - 2, pImg2 + 2);

        accumA = _mm_srai_epi32(accumA, SHIFT);
        accumB = _mm_srai_epi32(accumB, SHIFT);

        accumA = _mm_blend_epi16(accumA, _mm_slli_si128(accumB, 2), 0xAA);
        accumA = _mm_add_epi16(accumA, cur);
        accumA = _mm_min_epi16(mmMax, _mm_max_epi16(accumA, mmMin));

        _mm_storeu_si128((__m128i *)(&(fixedFilterResults[fixedFiltInd][curBlk.y + i + ii][curBlk.x + j])), accumA);
      } //for (size_t ii = 0; ii < STEP_Y; ii++)
    }//for (size_t j = 0; j < width; j += STEP_X)
    src += srcStride * STEP_Y;
  }
}

static void simdDeriveClassificationLaplacian(const CPelBuf &srcLuma, const Area &blkDst, const Area &blk, uint32_t **laplacian[NUM_DIRECTIONS])
{
  const size_t imgStride = srcLuma.stride;
  const Pel *  srcExt = srcLuma.buf;
  const size_t flP1 = ALF_CLASSIFIER_FL + 1;
  const size_t fl2 = ALF_CLASSIFIER_FL << 1;
  const int imgHExtended = blk.height + fl2;
  const int imgWExtended = blk.width + fl2;

  const int posX = blk.pos().x;
  const int posY = blk.pos().y;

  const size_t stride2 = imgStride * 2;
  const size_t stride3 = imgStride * 3;

  uint16_t colSums2x2[NUM_DIRECTIONS][(AdaptiveLoopFilter::m_CLASSIFICATION_BLK_SIZE + 10) >> 1]
    [((AdaptiveLoopFilter::m_CLASSIFICATION_BLK_SIZE + 16) >> 1) + 4];

  for (int i = 0; i < imgHExtended; i += 2)
  {
    const size_t offset = (i + posY - flP1) * imgStride + posX - flP1;

    const Pel *imgY0 = &srcExt[offset];
    const Pel *imgY1 = &srcExt[offset + imgStride];
    const Pel *imgY2 = &srcExt[offset + stride2];
    const Pel *imgY3 = &srcExt[offset + stride3];

    for (int j = 0; j < imgWExtended; j += 8)
    {
      const __m128i x0 = _mm_loadu_si128((const __m128i *) (imgY0 + j));
      const __m128i x1 = _mm_loadu_si128((const __m128i *) (imgY1 + j));
      const __m128i x2 = _mm_loadu_si128((const __m128i *) (imgY2 + j));
      const __m128i x3 = _mm_loadu_si128((const __m128i *) (imgY3 + j));
      const __m128i x0next = _mm_loadu_si128((const __m128i *) (imgY0 + j + 8));
      const __m128i x1next = _mm_loadu_si128((const __m128i *) (imgY1 + j + 8));
      const __m128i x2next = _mm_loadu_si128((const __m128i *) (imgY2 + j + 8));
      const __m128i x3next = _mm_loadu_si128((const __m128i *) (imgY3 + j + 8));

      const __m128i pixel1 = _mm_slli_epi16(_mm_alignr_epi8(x1next, x1, 2), 1);
      const __m128i pixel2 = _mm_slli_epi16(_mm_alignr_epi8(x2next, x2, 2), 1);

      //ver
      __m128i ver1 = _mm_add_epi16(_mm_alignr_epi8(x0next, x0, 2), _mm_alignr_epi8(x2next, x2, 2));
      ver1 = _mm_sub_epi16(pixel1, ver1);
      ver1 = _mm_abs_epi16(ver1);
      __m128i ver2 = _mm_add_epi16(_mm_alignr_epi8(x1next, x1, 2), _mm_alignr_epi8(x3next, x3, 2));
      ver2 = _mm_sub_epi16(pixel2, ver2);
      ver2 = _mm_abs_epi16(ver2);
      ver1 = _mm_add_epi16(ver1, ver2);  //8 ver (each is 2 values in a col)

      //hor
      __m128i hor1 = _mm_add_epi16(_mm_alignr_epi8(x1next, x1, 4), x1);
      hor1 = _mm_sub_epi16(pixel1, hor1);
      hor1 = _mm_abs_epi16(hor1);
      __m128i hor2 = _mm_add_epi16(_mm_alignr_epi8(x2next, x2, 4), x2);
      hor2 = _mm_sub_epi16(pixel2, hor2);
      hor2 = _mm_abs_epi16(hor2);
      hor1 = _mm_add_epi16(hor1, hor2);

      //dig0
      __m128i di01 = _mm_add_epi16(_mm_alignr_epi8(x2next, x2, 4), x0);
      di01 = _mm_sub_epi16(pixel1, di01);
      di01 = _mm_abs_epi16(di01);
      __m128i di02 = _mm_add_epi16(_mm_alignr_epi8(x3next, x3, 4), x1);
      di02 = _mm_sub_epi16(pixel2, di02);
      di02 = _mm_abs_epi16(di02);
      di01 = _mm_add_epi16(di01, di02);

      //dig1
      __m128i di11 = _mm_add_epi16(_mm_alignr_epi8(x0next, x0, 4), x2);
      di11 = _mm_sub_epi16(pixel1, di11);
      di11 = _mm_abs_epi16(di11);
      __m128i di12 = _mm_add_epi16(_mm_alignr_epi8(x1next, x1, 4), x3);
      di12 = _mm_sub_epi16(pixel2, di12);
      di12 = _mm_abs_epi16(di12);
      di11 = _mm_add_epi16(di11, di12);

      __m128i vh = _mm_hadd_epi16(ver1, hor1);   //ver: 2x2, 2x2, 2x2, 2x2; hor: 2x2, 2x2, 2x2, 2x2
      __m128i di = _mm_hadd_epi16(di01, di11);

      _mm_storel_epi64((__m128i *) &colSums2x2[VER][i >> 1][j >> 1], vh);
      _mm_storel_epi64((__m128i *) &colSums2x2[HOR][i >> 1][j >> 1], _mm_srli_si128(vh, 8));
      _mm_storel_epi64((__m128i *) &colSums2x2[DIAG0][i >> 1][j >> 1], di);
      _mm_storel_epi64((__m128i *) &colSums2x2[DIAG1][i >> 1][j >> 1], _mm_srli_si128(di, 8));
    }//(int j = 0; j < imgWExtended; j += 8)
  }//for (int i = 0; i < imgHExtended; i += 2)

  //get 4x4 sums
  for (int i = 0; i < (imgHExtended >> 1); i++)
  {
    for (int j = 0; j < (imgWExtended >> 1); j+=4)
    {   
      __m128i x0v = _mm_cvtepu16_epi32(_mm_loadl_epi64((__m128i *)&colSums2x2[VER][i][j]));
      __m128i x1v = _mm_cvtepu16_epi32(_mm_loadl_epi64((__m128i *)&colSums2x2[VER][i][j + 1]));
      __m128i x0h = _mm_cvtepu16_epi32(_mm_loadl_epi64((__m128i *)&colSums2x2[HOR][i][j]));
      __m128i x1h = _mm_cvtepu16_epi32(_mm_loadl_epi64((__m128i *)&colSums2x2[HOR][i][j + 1]));
      __m128i x0d0 = _mm_cvtepu16_epi32(_mm_loadl_epi64((__m128i *)&colSums2x2[DIAG0][i][j]));
      __m128i x1d0 = _mm_cvtepu16_epi32(_mm_loadl_epi64((__m128i *)&colSums2x2[DIAG0][i][j + 1]));
      __m128i x0d1 = _mm_cvtepu16_epi32(_mm_loadl_epi64((__m128i *)&colSums2x2[DIAG1][i][j]));
      __m128i x1d1 = _mm_cvtepu16_epi32(_mm_loadl_epi64((__m128i *)&colSums2x2[DIAG1][i][j + 1]));

      x0v = _mm_add_epi32(x0v, x1v);
      x0h = _mm_add_epi32(x0h, x1h);
      x0d0 = _mm_add_epi32(x0d0, x1d0);
      x0d1 = _mm_add_epi32(x0d1, x1d1);
      _mm_storeu_si128((__m128i *) &laplacian[VER][i][j], x0v);
      _mm_storeu_si128((__m128i *) &laplacian[HOR][i][j], x0h);
      _mm_storeu_si128((__m128i *) &laplacian[DIAG0][i][j], x0d0);
      _mm_storeu_si128((__m128i *) &laplacian[DIAG1][i][j], x0d1); //2x4

      if (i > 0) //4x4
      {
        x1v = _mm_loadu_si128((const __m128i *) &laplacian[VER][i-1][j]);
        x1h = _mm_loadu_si128((const __m128i *) &laplacian[HOR][i-1][j]);
        x1d0 = _mm_loadu_si128((const __m128i *) &laplacian[DIAG0][i-1][j]);
        x1d1 = _mm_loadu_si128((const __m128i *) &laplacian[DIAG1][i-1][j]);
        x0v = _mm_add_epi32(x0v, x1v);
        x0h = _mm_add_epi32(x0h, x1h);
        x0d0 = _mm_add_epi32(x0d0, x1d0);
        x0d1 = _mm_add_epi32(x0d1, x1d1);
        _mm_storeu_si128((__m128i *) &laplacian[VER][i - 1][j], x0v);
        _mm_storeu_si128((__m128i *) &laplacian[HOR][i - 1][j], x0h);
        _mm_storeu_si128((__m128i *) &laplacian[DIAG0][i - 1][j], x0d0);
        _mm_storeu_si128((__m128i *) &laplacian[DIAG1][i - 1][j], x0d1);
      }
    }//for (int j = 0; j < (imgWExtended >> 1); j+=4)
  }// for (int i = 0; i < (imgHExtended >> 1); i++)
}

static void simdDeriveClassificationLaplacianBig(const Area &curBlk, uint32_t **laplacian[NUM_DIRECTIONS])
{
  //get 12x12 sums, laplacian stores 4x4 sums
  int fl2 = ALF_CLASSIFIER_FL << 1;
  const int imgHExtended = curBlk.height + fl2;
  const int imgWExtended = curBlk.width + fl2;

  //4x12 sums
  for (int i = 0; i < (imgHExtended >> 1); i++)
  {
    __m128i x0v = _mm_loadu_si128((const __m128i *) &laplacian[VER][i][0]);
    __m128i x0h = _mm_loadu_si128((const __m128i *) &laplacian[HOR][i][0]);
    __m128i x0d0 = _mm_loadu_si128((const __m128i *) &laplacian[DIAG0][i][0]);
    __m128i x0d1 = _mm_loadu_si128((const __m128i *) &laplacian[DIAG1][i][0]);
    for (int j = 4; j < (imgWExtended >> 1); j += 4)
    {      
      __m128i x2v = _mm_loadu_si128((const __m128i *) &laplacian[VER][i][j]);
      __m128i x1v = _mm_shuffle_epi32(_mm_blend_epi16(x0v, x2v, 0x0f), 0x4e);
      x0v = _mm_add_epi32(x0v, x1v);
      x0v = _mm_add_epi32(x0v, x2v);
      
      __m128i x2h = _mm_loadu_si128((const __m128i *) &laplacian[HOR][i][j]);
      __m128i x1h = _mm_shuffle_epi32(_mm_blend_epi16(x0h, x2h, 0x0f), 0x4e);
      x0h = _mm_add_epi32(x0h, x1h);
      x0h = _mm_add_epi32(x0h, x2h);

      __m128i x2d0 = _mm_loadu_si128((const __m128i *) &laplacian[DIAG0][i][j]);
      __m128i x1d0 = _mm_shuffle_epi32(_mm_blend_epi16(x0d0, x2d0, 0x0f), 0x4e);
      x0d0 = _mm_add_epi32(x0d0, x1d0);
      x0d0 = _mm_add_epi32(x0d0, x2d0);

      __m128i x2d1 = _mm_loadu_si128((const __m128i *) &laplacian[DIAG1][i][j]);
      __m128i x1d1 = _mm_shuffle_epi32(_mm_blend_epi16(x0d1, x2d1, 0x0f), 0x4e);
      x0d1 = _mm_add_epi32(x0d1, x1d1);
      x0d1 = _mm_add_epi32(x0d1, x2d1);   

      _mm_storeu_si128((__m128i *) &laplacian[VER][i][j-4], x0v);
      _mm_storeu_si128((__m128i *) &laplacian[HOR][i][j-4], x0h);
      _mm_storeu_si128((__m128i *) &laplacian[DIAG0][i][j-4], x0d0);
      _mm_storeu_si128((__m128i *) &laplacian[DIAG1][i][j-4], x0d1);
      x0v = x2v;
      x0h = x2h;
      x0d0 = x2d0;
      x0d1 = x2d1;
    }
  }

  for (int i = 0; i < (imgHExtended >> 1) - 5; i++)
  {
    for (int j = 0; j < (imgWExtended >> 1) - 5; j += 4)
    {
      __m128i x0v = _mm_loadu_si128((const __m128i *) &laplacian[VER][i][j]);
      __m128i x1v = _mm_loadu_si128((const __m128i *) &laplacian[VER][i + 2][j]);
      __m128i x2v = _mm_loadu_si128((const __m128i *) &laplacian[VER][i + 4][j]);
      x0v = _mm_add_epi32(x0v, x1v);
      x0v = _mm_add_epi32(x0v, x2v);

      __m128i x0h = _mm_loadu_si128((const __m128i *) &laplacian[HOR][i][j]);
      __m128i x1h = _mm_loadu_si128((const __m128i *) &laplacian[HOR][i + 2][j]);
      __m128i x2h = _mm_loadu_si128((const __m128i *) &laplacian[HOR][i + 4][j]);
      x0h = _mm_add_epi32(x0h, x1h);
      x0h = _mm_add_epi32(x0h, x2h);

      __m128i x0d0 = _mm_loadu_si128((const __m128i *) &laplacian[DIAG0][i][j]);
      __m128i x1d0 = _mm_loadu_si128((const __m128i *) &laplacian[DIAG0][i + 2][j]);
      __m128i x2d0 = _mm_loadu_si128((const __m128i *) &laplacian[DIAG0][i + 4][j]);
      x0d0 = _mm_add_epi32(x0d0, x1d0);
      x0d0 = _mm_add_epi32(x0d0, x2d0);

      __m128i x0d1 = _mm_loadu_si128((const __m128i *) &laplacian[DIAG1][i][j]);
      __m128i x1d1 = _mm_loadu_si128((const __m128i *) &laplacian[DIAG1][i + 2][j]);
      __m128i x2d1 = _mm_loadu_si128((const __m128i *) &laplacian[DIAG1][i + 4][j]);
      x0d1 = _mm_add_epi32(x0d1, x1d1);
      x0d1 = _mm_add_epi32(x0d1, x2d1);

      _mm_storeu_si128((__m128i *) &laplacian[VER][i][j], x0v);
      _mm_storeu_si128((__m128i *) &laplacian[HOR][i][j], x0h);
      _mm_storeu_si128((__m128i *) &laplacian[DIAG0][i][j], x0d0);
      _mm_storeu_si128((__m128i *) &laplacian[DIAG1][i][j], x0d1);
    }
  }
}

static void simdCalcClass0(AlfClassifier **classifier, const Area &blkDst, const Area &curBlk, int dirWindSize, int classDir, int noDir, int noAct, int bitDepth, int subBlkSize, int mappingDir[NUM_DIR_FIX][NUM_DIR_FIX], uint32_t **laplacian[NUM_DIRECTIONS])
{
  const __m128i shift = _mm_cvtsi32_si128(9 + bitDepth - 4);
  const int multTab[] = { 5628, 1407, 624, 351, 225, 156 };
  const __m128i mult = _mm_set1_epi32(multTab[dirWindSize]);
  const __m128i dirOff = _mm_set1_epi32(noDir * (noDir + 1));
  const __m128i ones = _mm_set1_epi32(1);
  const __m128i zeros = _mm_setzero_si128();
  const __m128i scale = _mm_set1_epi32(192);

  int lapOffset = (dirWindSize == 1) ? 2 : 0;
  for (int i = 0; i < curBlk.height; i += 2)
  {
    int iOffset = (i >> 1) + lapOffset;
    for (int j = 0; j < curBlk.width; j += 8)
    {
      int jOffset = (j >> 1) + lapOffset;
      __m128i sumV = _mm_loadu_si128((const __m128i *) &laplacian[VER][iOffset][jOffset]);  //4 32-bit values
      __m128i sumH = _mm_loadu_si128((const __m128i *) &laplacian[HOR][iOffset][jOffset]);
      __m128i sumD0 = _mm_loadu_si128((const __m128i *) &laplacian[DIAG0][iOffset][jOffset]);
      __m128i sumD1 = _mm_loadu_si128((const __m128i *) &laplacian[DIAG1][iOffset][jOffset]);

      //sum += sumV + sumH;
      __m128i tempAct = _mm_add_epi32(sumV, sumH);
      __m128i activity = _mm_mullo_epi32(tempAct, mult);
      activity = _mm_srl_epi32(activity, shift);
      activity = _mm_min_epi32(activity, scale);

      __m128i xmm2 = activity;
      __m128i xmm0 = _mm_setzero_si128();
      __m128i xmm15 = _mm_cmpeq_epi32(xmm0, xmm0);
      __m128i xmm1 = _mm_srli_epi32(xmm15, 31);
      __m128i xmm7 = _mm_srli_epi32(xmm15, 29);
      __m128i xmm9 = _mm_add_epi32(_mm_slli_epi32(xmm7, 2), xmm1);

      __m128i LUT192 = _mm_set_epi32(0x0C020A00, 0x0E040608, 0x0E040608, 0x0C020A00);

      xmm2 = _mm_or_si128(xmm2, _mm_srli_epi32(xmm2, 1));
      xmm2 = _mm_or_si128(xmm2, _mm_srli_epi32(xmm2, 2));
      xmm2 = _mm_or_si128(xmm2, _mm_srli_epi32(xmm2, 4));
      xmm2 = _mm_mullo_epi16(xmm2, xmm9);
      xmm2 = _mm_and_si128(_mm_srli_epi32(xmm2, 5), xmm7);
      xmm2 = _mm_shuffle_epi8(LUT192, xmm2);

       __m128i xmm4 = _mm_xor_si128(activity, _mm_srli_epi32(activity, 1));
      xmm4 = _mm_cmplt_epi32(xmm4, activity);
      xmm4 = _mm_or_si128(_mm_cmpeq_epi32(activity, xmm1), xmm4);
      xmm4 = _mm_and_si128(xmm4, xmm1);

      activity = _mm_or_si128(xmm2, xmm4);
      
      __m128i hv1 = _mm_max_epi32(sumV, sumH);
      __m128i hv0 = _mm_min_epi32(sumV, sumH);

      __m128i d1 = _mm_max_epi32(sumD0, sumD1);
      __m128i d0 = _mm_min_epi32(sumD0, sumD1);

      //edgeStrengthHV, to optimize
      __m128i hv0Two = _mm_slli_epi32(hv0, 1);
      __m128i hv0Eight = _mm_slli_epi32(hv0, 3);
      __m128i hv1Two = _mm_slli_epi32(hv1, 1);
      __m128i strength = _mm_cmpgt_epi32(_mm_slli_epi32(hv1, 2), _mm_add_epi32(hv0, _mm_slli_epi32(hv0, 2)));  //4, 5
      __m128i edgeStrengthHV = _mm_and_si128(strength, ones);

      strength = _mm_cmpgt_epi32(hv1Two, _mm_add_epi32(hv0, hv0Two)); //2, 3
      edgeStrengthHV = _mm_add_epi32(edgeStrengthHV, _mm_and_si128(strength, ones));

      strength = _mm_cmpgt_epi32(hv1, hv0Two); //1, 2
      edgeStrengthHV = _mm_add_epi32(edgeStrengthHV, _mm_and_si128(strength, ones));

      strength = _mm_cmpgt_epi32(hv1, _mm_add_epi32(hv0, hv0Two)); //1, 3
      edgeStrengthHV = _mm_add_epi32(edgeStrengthHV, _mm_and_si128(strength, ones));

      strength = _mm_cmpgt_epi32(hv1Two, _mm_add_epi32(hv0, hv0Eight)); //2, 9
      edgeStrengthHV = _mm_add_epi32(edgeStrengthHV, _mm_and_si128(strength, ones));

      strength = _mm_cmpgt_epi32(hv1, hv0Eight); //1, 8
      edgeStrengthHV = _mm_add_epi32(edgeStrengthHV, _mm_and_si128(strength, ones));

      //edgeStrengthD, to optimize
      __m128i d0Two = _mm_slli_epi32(d0, 1);
      __m128i d0Eight = _mm_slli_epi32(d0, 3);
      __m128i d1Two = _mm_slli_epi32(d1, 1);
      strength = _mm_cmpgt_epi32(_mm_slli_epi32(d1, 2), _mm_add_epi32(d0, _mm_slli_epi32(d0, 2))); //4, 5
      __m128i edgeStrengthD = _mm_and_si128(strength, ones);

      strength = _mm_cmpgt_epi32(d1Two, _mm_add_epi32(d0, d0Two)); //2, 3
      edgeStrengthD = _mm_add_epi32(edgeStrengthD, _mm_and_si128(strength, ones));

      strength = _mm_cmpgt_epi32(d1, d0Two); //1, 2
      edgeStrengthD = _mm_add_epi32(edgeStrengthD, _mm_and_si128(strength, ones)); 

      strength = _mm_cmpgt_epi32(d1, _mm_add_epi32(d0, d0Two)); //1, 3
      edgeStrengthD = _mm_add_epi32(edgeStrengthD, _mm_and_si128(strength, ones));

      strength = _mm_cmpgt_epi32(d1Two, _mm_add_epi32(d0, d0Eight));//2, 9
      edgeStrengthD = _mm_add_epi32(edgeStrengthD, _mm_and_si128(strength, ones));

      strength = _mm_cmpgt_epi32(d1, d0Eight);//1, 8
      edgeStrengthD = _mm_add_epi32(edgeStrengthD, _mm_and_si128(strength, ones));

      const __m128i hv1Xd0e = _mm_mul_epi32(hv1, d0); 
      const __m128i hv0Xd1e = _mm_mul_epi32(hv0, d1);
      const __m128i hv1Xd0o = _mm_mul_epi32(_mm_srli_si128(hv1, 4), _mm_srli_si128(d0, 4));
      const __m128i hv0Xd1o = _mm_mul_epi32(_mm_srli_si128(hv0, 4), _mm_srli_si128(d1, 4));

      const __m128i xmme = _mm_sub_epi64(hv0Xd1e, hv1Xd0e);
      const __m128i xmmo = _mm_sub_epi64(hv0Xd1o, hv1Xd0o);

      __m128i dirCondition = _mm_srai_epi32(_mm_blend_epi16(_mm_srli_si128(xmme, 4), xmmo, 0xCC), 31);

      __m128i cx = _mm_blendv_epi8(edgeStrengthHV, edgeStrengthD, dirCondition);  //x
      __m128i cy = _mm_blendv_epi8(edgeStrengthD, edgeStrengthHV, dirCondition);  //y
      __m128i dirOffset = _mm_blendv_epi8(_mm_set1_epi32(28), zeros, dirCondition);
      //direction = (y*(y+1))/2 + x
      __m128i direction = _mm_mullo_epi32(cy, cy);
      direction = _mm_add_epi32(direction, cy);
      direction = _mm_srli_epi32(direction, 1);
      direction = _mm_add_epi32(direction, cx);
      direction = _mm_andnot_si128(_mm_cmpgt_epi32(cx, cy), direction);
      direction = _mm_add_epi32(direction, dirOffset);

      __m128i classIdx = _mm_mullo_epi32(dirOff, activity);
      classIdx = _mm_add_epi32(classIdx, direction);     

      //transpose
      __m128i dirTempHVMinus1 = _mm_cmpgt_epi32(sumV, sumH);
      __m128i dirTempDMinus1 = _mm_cmpgt_epi32(sumD0, sumD1);
      __m128i transposeIdx = _mm_set1_epi32(3);
      transposeIdx = _mm_add_epi32(transposeIdx, dirTempHVMinus1);
      transposeIdx = _mm_add_epi32(transposeIdx, dirTempDMinus1);
      transposeIdx = _mm_add_epi32(transposeIdx, dirTempDMinus1);

      classIdx = _mm_slli_epi16(classIdx, 2);
      classIdx = _mm_add_epi16(classIdx, transposeIdx);
      classIdx = _mm_shuffle_epi8(classIdx, _mm_setr_epi8(0, 1, 0, 1, 4, 5, 4, 5, 8, 9, 8, 9, 12, 13, 12, 13));
      _mm_storeu_si128((__m128i *) &classifier[curBlk.pos().y + i][curBlk.pos().x + j], classIdx);
      _mm_storeu_si128((__m128i *) &classifier[curBlk.pos().y + i + 1][curBlk.pos().x + j], classIdx);
    }//for (int j = 0; j < curBlk.width; j += 8)
  }//for (int i = 0; i < curBlk.height; i += 2)
}

static void simdCalcClass1(AlfClassifier **classifier, const Area &blkDst, const Area &curBlk, int dirWindSize, int classDir, int noDir, int noAct, int bitDepth, int subBlkSize, int mappingDir[NUM_DIR_FIX][NUM_DIR_FIX], uint32_t **laplacian[NUM_DIRECTIONS])
{
  const __m128i shift = _mm_cvtsi32_si128(9 + bitDepth);
  const int multTab[] = { 5628, 1407, 624, 351, 225, 156 };
  const __m128i mult = _mm_set1_epi32(multTab[dirWindSize]);
  const __m128i scale = _mm_set1_epi32(15);

  for (int i = 0; i < curBlk.height; i += 2)
  {
    int iOffset = i >> 1;
    for (int j = 0; j < curBlk.width; j += 8)
    {
      int jOffset = j >> 1;
      __m128i sumV = _mm_loadu_si128((const __m128i *) &laplacian[VER][iOffset][jOffset]);  //4 32-bit values
      __m128i sumH = _mm_loadu_si128((const __m128i *) &laplacian[HOR][iOffset][jOffset]);
      __m128i sumD0 = _mm_loadu_si128((const __m128i *) &laplacian[DIAG0][iOffset][jOffset]);
      __m128i sumD1 = _mm_loadu_si128((const __m128i *) &laplacian[DIAG1][iOffset][jOffset]);

      //sum += sumV + sumH;
      __m128i tempAct = _mm_add_epi32(sumV, sumH);
      __m128i activity = _mm_mullo_epi32(tempAct, mult);
      activity = _mm_srl_epi32(activity, shift);
      activity = _mm_min_epi32(activity, scale);
      __m128i classIdx = _mm_shuffle_epi8(_mm_setr_epi8(0, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4), activity);
      classIdx = _mm_add_epi32(classIdx, _mm_slli_epi32(classIdx, 2));  // activity * 5

      __m128i hv1 = _mm_max_epi32(sumV, sumH);
      __m128i hv0 = _mm_min_epi32(sumV, sumH);

      __m128i d1 = _mm_max_epi32(sumD0, sumD1);
      __m128i d0 = _mm_min_epi32(sumD0, sumD1);

      const __m128i hv1Xd0e = _mm_mul_epi32(hv1, d0);
      const __m128i hv0Xd1e = _mm_mul_epi32(hv0, d1);
      const __m128i hv1Xd0o = _mm_mul_epi32(_mm_srli_si128(hv1, 4), _mm_srli_si128(d0, 4));
      const __m128i hv0Xd1o = _mm_mul_epi32(_mm_srli_si128(hv0, 4), _mm_srli_si128(d1, 4));

      const __m128i xmme = _mm_sub_epi64(hv1Xd0e, hv0Xd1e);
      const __m128i xmmo = _mm_sub_epi64(hv1Xd0o, hv0Xd1o);

      __m128i dirCondition = _mm_srai_epi32(_mm_blend_epi16(_mm_srli_si128(xmme, 4), xmmo, 0xCC), 31);

      __m128i hvd1 = _mm_blendv_epi8(hv1, d1, dirCondition);
      __m128i hvd0 = _mm_blendv_epi8(hv0, d0, dirCondition);
      __m128i strength1 = _mm_cmpgt_epi32(hvd1, _mm_add_epi32(hvd0, hvd0));
      __m128i strength2 = _mm_cmpgt_epi32(_mm_add_epi32(hvd1, hvd1), _mm_add_epi32(hvd0, _mm_slli_epi32(hvd0, 3)));
      __m128i offset = _mm_and_si128(strength1, _mm_set1_epi32(1));
      __m128i direction = _mm_add_epi32(offset, _mm_and_si128(strength2, _mm_set1_epi32(1)));
      direction = _mm_add_epi32(direction, _mm_andnot_si128(dirCondition, _mm_set1_epi32(2)));      
      direction = _mm_and_si128(direction, strength1);      
      classIdx = _mm_add_epi32(direction, classIdx);

      //transpose
      __m128i dirTempHVMinus1 = _mm_cmpgt_epi32(sumV, sumH);
      __m128i dirTempDMinus1 = _mm_cmpgt_epi32(sumD0, sumD1);
      __m128i transposeIdx = _mm_set1_epi32(3);
      transposeIdx = _mm_add_epi32(transposeIdx, dirTempHVMinus1);
      transposeIdx = _mm_add_epi32(transposeIdx, dirTempDMinus1);
      transposeIdx = _mm_add_epi32(transposeIdx, dirTempDMinus1);
      classIdx = _mm_slli_epi16(classIdx, 2);
      classIdx = _mm_add_epi16(classIdx, transposeIdx);
      classIdx = _mm_shuffle_epi8(classIdx, _mm_setr_epi8(0, 1, 0, 1, 4, 5, 4, 5, 8, 9, 8, 9, 12, 13, 12, 13));
      _mm_storeu_si128((__m128i *) &classifier[curBlk.pos().y + i][curBlk.pos().x + j], classIdx);
      _mm_storeu_si128((__m128i *) &classifier[curBlk.pos().y + i + 1][curBlk.pos().x + j], classIdx);
    }//for (int j = 0; j < curBlk.width; j += 8)
  }//for (int i = 0; i < curBlk.height; i += 2)
}
#endif


template <X86_VEXT vext>
void AdaptiveLoopFilter::_initAdaptiveLoopFilterX86()
{
  m_filter5x5Blk = simdFilter5x5Blk<vext>;
  m_filter7x7Blk = simdFilter7x7Blk<vext>;
#if ALF_IMPROVEMENT
  m_filter9x9Blk = simdFilter9x9Blk<vext>;
  m_filter9x9BlkExt = simdFilter9x9BlkExt<vext>;
  m_filter13x13Blk = simdFilter13x13Blk<vext>;
  m_deriveClassificationLaplacian = simdDeriveClassificationLaplacian;
  m_deriveClassificationLaplacianBig = simdDeriveClassificationLaplacianBig;
  m_calcClass0 = simdCalcClass0;
  m_calcClass1 = simdCalcClass1;

  for( int i = 0; i < NUM_SETS_FIXED_FILTERS; i++ )
  {
    for( int j = 0; j < NUM_CLASSIFIER; j++ )
    {
      for( int k = 0; k < NUM_FIXED_FILTERS; k++ )
      {
        for( int m = 0; m < 42; m++ )
        {
          packedDataFixedFilters[i][j][k][m] = ( m_filterCoeffFixed[i][j][k][p0[m]] << 2 ) | m_clippingFixed[i][j][k][p0[m]];
        }
      }
    }
  }
#else
  m_deriveClassificationBlk = simdDeriveClassificationBlk<vext>;
#endif
}

template void AdaptiveLoopFilter::_initAdaptiveLoopFilterX86<SIMDX86>();
#endif   // TARGET_SIMD_X86
