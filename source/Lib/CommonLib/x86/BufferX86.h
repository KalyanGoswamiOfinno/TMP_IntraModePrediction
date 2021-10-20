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

/** \file     YuvX86.cpp
    \brief    SIMD averaging.
*/

//! \ingroup CommonLib
//! \{


#include "CommonLib/CommonDef.h"
#include "CommonDefX86.h"
#include "CommonLib/Unit.h"
#include "CommonLib/Buffer.h"
#include "CommonLib/InterpolationFilter.h"

#if ENABLE_SIMD_OPT_BUFFER
#ifdef TARGET_SIMD_X86
#if JVET_W0097_GPM_MMVD_TM
template< X86_VEXT vext >
void roundBD_SSE(const Pel* srcp, const int srcStride, Pel* dest, const int destStride, int width, int height, const ClpRng& clpRng)
{
  const int32_t clipbd = clpRng.bd;
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  const int32_t shiftDefault = IF_INTERNAL_FRAC_BITS(clipbd);
#else
  const int32_t shiftDefault = std::max<int>(2, (IF_INTERNAL_PREC - clipbd));
#endif
  const int32_t offsetDefault = (1 << (shiftDefault - 1)) + IF_INTERNAL_OFFS;

#if USE_AVX2
  if (vext >= AVX2 && (width & 0x0f) == 0)
  {
    __m256i voffset = _mm256_set1_epi16((short)offsetDefault);
    __m256i vibdmin = _mm256_set1_epi16((short)clpRng.min);
    __m256i vibdmax = _mm256_set1_epi16((short)clpRng.max);
    __m256i vsrc;
    for (int row = 0; row < height; row++)
    {
      for (int col = 0; col < width; col += 16)
      {
        vsrc = _mm256_lddqu_si256((__m256i *)&srcp[col]);
        vsrc = _mm256_adds_epi16(vsrc, voffset);
        vsrc = _mm256_srai_epi16(vsrc, shiftDefault);
        vsrc = _mm256_min_epi16(vibdmax, _mm256_max_epi16(vibdmin, vsrc));
        _mm256_storeu_si256((__m256i *)&dest[col], vsrc);
      }
      srcp += srcStride;
      dest += destStride;
    }
  }
  else
  {
#endif
    __m128i voffset = _mm_set1_epi16((short)offsetDefault);
    __m128i vibdmin = _mm_set1_epi16((short)clpRng.min);
    __m128i vibdmax = _mm_set1_epi16((short)clpRng.max);
    __m128i vsrc;
    for (int row = 0; row < height; row++)
    {
      int col = 0;
      for (; col < ((width >> 3) << 3); col += 8)
      {
        vsrc = _mm_lddqu_si128((__m128i *)&srcp[col]);
        vsrc = _mm_adds_epi16(vsrc, voffset);
        vsrc = _mm_srai_epi16(vsrc, shiftDefault);
        vsrc = _mm_min_epi16(vibdmax, _mm_max_epi16(vibdmin, vsrc));
        _mm_storeu_si128((__m128i *)&dest[col], vsrc);
      }
      for (; col < ((width >> 2) << 2); col += 4)
      {
        vsrc = _mm_loadl_epi64((__m128i *)&srcp[col]);
        vsrc = _mm_adds_epi16(vsrc, voffset);
        vsrc = _mm_srai_epi16(vsrc, shiftDefault);
        vsrc = _mm_min_epi16(vibdmax, _mm_max_epi16(vibdmin, vsrc));
        _mm_storel_epi64((__m128i *)&dest[col], vsrc);
      }
      for (; col < width; col++)
      {
        dest[col] = ClipPel(rightShift(srcp[col] + offsetDefault, shiftDefault), clpRng);
      }
      srcp += srcStride;
      dest += destStride;
    }
#if USE_AVX2
  }
#endif
}

template< X86_VEXT vext >
void weightedAvg_SSE(const Pel* src0, const unsigned src0Stride, const Pel* src1, const unsigned src1Stride, Pel* dest, const unsigned destStride, const int8_t w0, const int8_t w1, int width, int height, const ClpRng& clpRng)
{
  const int8_t log2WeightBase = g_BcwLog2WeightBase;
  const int    clipbd = clpRng.bd;
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  const int shiftNum = IF_INTERNAL_FRAC_BITS(clipbd) + log2WeightBase;
#else
  const int shiftNum = std::max<int>(2, (IF_INTERNAL_PREC - clipbd)) + log2WeightBase;
#endif
  const int offset = (1 << (shiftNum - 1)) + (IF_INTERNAL_OFFS << log2WeightBase);

#if USE_AVX2
  if ((vext >= AVX2) && (width & 0x7) == 0)
  {
    __m256i mw = _mm256_unpacklo_epi16(_mm256_set1_epi16(w0), _mm256_set1_epi16(w1));
    __m256i voffset = _mm256_set1_epi32(offset);
    __m256i vibdmin = _mm256_set1_epi16((short)clpRng.min);
    __m256i vibdmax = _mm256_set1_epi16((short)clpRng.max);
    __m256i msrc0, msrc1, msum0, msum1;

    for (int row = 0; row < height; row++)
    {
      for (int col = 0; col < width; col += 8)
      {
        msrc0 = _mm256_castsi128_si256(_mm_lddqu_si128((__m128i*)(&src0[col])));
        msrc1 = _mm256_castsi128_si256(_mm_lddqu_si128((__m128i*)(&src1[col])));
        msum0 = _mm256_unpacklo_epi16(msrc0, msrc1);
        msum1 = _mm256_unpackhi_epi16(msrc0, msrc1);
        msum0 = _mm256_madd_epi16(msum0, mw);
        msum1 = _mm256_madd_epi16(msum1, mw);
        msum0 = _mm256_add_epi32(msum0, voffset);
        msum1 = _mm256_add_epi32(msum1, voffset);
        msum0 = _mm256_srai_epi32(msum0, shiftNum);
        msum1 = _mm256_srai_epi32(msum1, shiftNum);
        msum0 = _mm256_packs_epi32(msum0, msum1);
        msum0 = _mm256_min_epi16(vibdmax, _mm256_max_epi16(vibdmin, msum0));
        _mm_storeu_si128((__m128i *)&dest[col], _mm256_castsi256_si128(msum0));
      }
      src0 += src0Stride;
      src1 += src1Stride;
      dest += destStride;
    }
  }
  else
  {
#endif
    __m128i mw = _mm_unpacklo_epi16(_mm_set1_epi16(w0), _mm_set1_epi16(w1));
    __m128i voffset = _mm_set1_epi32(offset);
    __m128i vibdmin = _mm_set1_epi16((short)clpRng.min);
    __m128i vibdmax = _mm_set1_epi16((short)clpRng.max);

    for (int row = 0; row < height; row++)
    {
      int col = 0;
      for (; col < ((width >> 2) << 2); col += 4)
      {
        __m128i msrc = _mm_unpacklo_epi16(_mm_loadl_epi64((__m128i *)&src0[col]), _mm_loadl_epi64((__m128i *)&src1[col]));
        msrc = _mm_madd_epi16(msrc, mw);
        msrc = _mm_add_epi32(msrc, voffset);
        msrc = _mm_srai_epi32(msrc, shiftNum);
        msrc = _mm_packs_epi32(msrc, msrc);
        msrc = _mm_min_epi16(vibdmax, _mm_max_epi16(vibdmin, msrc));
        _mm_storel_epi64((__m128i *)&dest[col], msrc);
      }
      for (; col < width; col++)
      {
        dest[col] = ClipPel(rightShift(src0[col] * w0 + src1[col] * w1 + offset, shiftNum), clpRng);
      }
      src0 += src0Stride;
      src1 += src1Stride;
      dest += destStride;
    }
#if USE_AVX2
  }
#endif
}

template< X86_VEXT vext >
void copyClip_SSE(const Pel* srcp, const unsigned srcStride, Pel* dest, const unsigned destStride, int width, int height, const ClpRng& clpRng)
{
#if USE_AVX2
  if (vext >= AVX2 && (width & 0x0f) == 0)
  {
    __m256i vibdmin = _mm256_set1_epi16((short)clpRng.min);
    __m256i vibdmax = _mm256_set1_epi16((short)clpRng.max);
    __m256i vsrc;
    for (int row = 0; row < height; row++)
    {
      for (int col = 0; col < width; col += 16)
      {
        vsrc = _mm256_lddqu_si256((__m256i *)&srcp[col]);
        vsrc = _mm256_min_epi16(vibdmax, _mm256_max_epi16(vibdmin, vsrc));
        _mm256_storeu_si256((__m256i *)&dest[col], vsrc);
      }
      srcp += srcStride;
      dest += destStride;
    }
  }
  else
  {
#endif
    __m128i vibdmin = _mm_set1_epi16((short)clpRng.min);
    __m128i vibdmax = _mm_set1_epi16((short)clpRng.max);
    __m128i vsrc;
    for (int row = 0; row < height; row++)
    {
      int col = 0;
      for (; col < ((width >> 3) << 3); col += 8)
      {
        vsrc = _mm_lddqu_si128((__m128i *)&srcp[col]);
        vsrc = _mm_min_epi16(vibdmax, _mm_max_epi16(vibdmin, vsrc));
        _mm_storeu_si128((__m128i *)&dest[col], vsrc);
      }
      for (; col < ((width >> 2) << 2); col += 4)
      {
        vsrc = _mm_loadl_epi64((__m128i *)&srcp[col]);
        vsrc = _mm_min_epi16(vibdmax, _mm_max_epi16(vibdmin, vsrc));
        _mm_storel_epi64((__m128i *)&dest[col], vsrc);
      }
      for (; col < width; col++)
      {
        dest[col] = ClipPel(srcp[col], clpRng);
      }
      srcp += srcStride;
      dest += destStride;
    }
#if USE_AVX2
  }
#endif
}
#endif

template< X86_VEXT vext, int W >
void addAvg_SSE( const int16_t* src0, int src0Stride, const int16_t* src1, int src1Stride, int16_t *dst, int dstStride, int width, int height, int shift, int offset, const ClpRng& clpRng )
{
  if( W == 8 )
  {
    CHECK(offset & 1, "offset must be even");
    CHECK(offset < -32768 || offset > 32767, "offset must be a 16-bit value");

    __m128i vibdimin = _mm_set1_epi16(clpRng.min);
    __m128i vibdimax = _mm_set1_epi16(clpRng.max);

    for (int row = 0; row < height; row++)
    {
      for (int col = 0; col < width; col += 8)
      {
        __m128i vsrc0 = _mm_loadu_si128((const __m128i *) &src0[col]);
        __m128i vsrc1 = _mm_loadu_si128((const __m128i *) &src1[col]);

        vsrc0 = _mm_xor_si128(vsrc0, _mm_set1_epi16(0x7fff));
        vsrc1 = _mm_xor_si128(vsrc1, _mm_set1_epi16(0x7fff));
        vsrc0 = _mm_avg_epu16(vsrc0, vsrc1);
        vsrc0 = _mm_xor_si128(vsrc0, _mm_set1_epi16(0x7fff));
        vsrc0 = _mm_adds_epi16(vsrc0, _mm_set1_epi16(offset >> 1));
        vsrc0 = _mm_sra_epi16(vsrc0, _mm_cvtsi32_si128(shift - 1));
        vsrc0 = _mm_max_epi16(vsrc0, vibdimin);
        vsrc0 = _mm_min_epi16(vsrc0, vibdimax);
        _mm_storeu_si128((__m128i *) &dst[col], vsrc0);
      }

      src0 += src0Stride;
      src1 += src1Stride;
      dst += dstStride;
    }
  }
  else if( W == 4 )
  {
    __m128i vzero     = _mm_setzero_si128();
    __m128i voffset   = _mm_set1_epi32( offset );
    __m128i vibdimin  = _mm_set1_epi16( clpRng.min );
    __m128i vibdimax  = _mm_set1_epi16( clpRng.max );

    for( int row = 0; row < height; row++ )
    {
      for( int col = 0; col < width; col += 4 )
      {
        __m128i vsum = _mm_loadl_epi64  ( ( const __m128i * )&src0[col] );
        __m128i vdst = _mm_loadl_epi64  ( ( const __m128i * )&src1[col] );
        vsum = _mm_cvtepi16_epi32       ( vsum );
        vdst = _mm_cvtepi16_epi32       ( vdst );
        vsum = _mm_add_epi32            ( vsum, vdst );
        vsum = _mm_add_epi32            ( vsum, voffset );
        vsum = _mm_srai_epi32           ( vsum, shift );
        vsum = _mm_packs_epi32          ( vsum, vzero );

        vsum = _mm_min_epi16( vibdimax, _mm_max_epi16( vibdimin, vsum ) );
        _mm_storel_epi64( ( __m128i * )&dst[col], vsum );
      }

      src0 += src0Stride;
      src1 += src1Stride;
      dst  +=  dstStride;
    }
  }
  else
  {
    THROW( "Unsupported size" );
  }
}

template<X86_VEXT vext>
void copyBufferSimd(Pel *src, int srcStride, Pel *dst, int dstStride, int width, int height)
{
  if (width < 8)
  {
    CHECK(width < 4, "width must be at least 4");

    for (size_t x = 0; x < width; x += 4)
    {
      if (x > width - 4)
        x = width - 4;
      for (size_t y = 0; y < height; y++)
      {
        __m128i val = _mm_loadl_epi64((const __m128i *) (src + y * srcStride + x));
        _mm_storel_epi64((__m128i *) (dst + y * dstStride + x), val);
      }
    }
  }
  else
  {
    for (size_t x = 0; x < width; x += 8)
    {
      if (x > width - 8)
        x = width - 8;
      for (size_t y = 0; y < height; y++)
      {
        __m128i val = _mm_loadu_si128((const __m128i *) (src + y * srcStride + x));
        _mm_storeu_si128((__m128i *) (dst + y * dstStride + x), val);
      }
    }
  }
}

template<X86_VEXT vext>
void paddingSimd(Pel *dst, int stride, int width, int height, int padSize)
{
  size_t extWidth = width + 2 * padSize;
  CHECK(extWidth < 8, "width plus 2 times padding size must be at least 8");

  if (padSize == 1)
  {
    for (size_t i = 0; i < height; i++)
    {
      Pel left                = dst[i * stride];
      Pel right               = dst[i * stride + width - 1];
      dst[i * stride - 1]     = left;
      dst[i * stride + width] = right;
    }

    dst -= 1;

    for (size_t i = 0; i < extWidth - 8; i++)
    {
      __m128i top = _mm_loadu_si128((const __m128i *) (dst + i));
      _mm_storeu_si128((__m128i *) (dst - stride + i), top);
    }
    __m128i top = _mm_loadu_si128((const __m128i *) (dst + extWidth - 8));
    _mm_storeu_si128((__m128i *) (dst - stride + extWidth - 8), top);

    dst += height * stride;

    for (size_t i = 0; i < extWidth - 8; i++)
    {
      __m128i bottom = _mm_loadu_si128((const __m128i *) (dst - stride + i));
      _mm_storeu_si128((__m128i *) (dst + i), bottom);
    }
    __m128i bottom = _mm_loadu_si128((const __m128i *) (dst - stride + extWidth - 8));
    _mm_storeu_si128((__m128i *) (dst + extWidth - 8), bottom);
  }
  else if (padSize == 2)
  {
    for (size_t i = 0; i < height; i++)
    {
      Pel left                    = dst[i * stride];
      Pel right                   = dst[i * stride + width - 1];
      dst[i * stride - 2]         = left;
      dst[i * stride - 1]         = left;
      dst[i * stride + width]     = right;
      dst[i * stride + width + 1] = right;
    }

    dst -= 2;

    for (size_t i = 0; i < extWidth - 8; i++)
    {
      __m128i top = _mm_loadu_si128((const __m128i *) (dst + i));
      _mm_storeu_si128((__m128i *) (dst - 2 * stride + i), top);
      _mm_storeu_si128((__m128i *) (dst - stride + i), top);
    }
    __m128i top = _mm_loadu_si128((const __m128i *) (dst + extWidth - 8));
    _mm_storeu_si128((__m128i *) (dst - 2 * stride + extWidth - 8), top);
    _mm_storeu_si128((__m128i *) (dst - stride + extWidth - 8), top);

    dst += height * stride;

    for (size_t i = 0; i < extWidth - 8; i++)
    {
      __m128i bottom = _mm_loadu_si128((const __m128i *) (dst - stride + i));
      _mm_storeu_si128((__m128i *) (dst + i), bottom);
      _mm_storeu_si128((__m128i *) (dst + stride + i), bottom);
    }
    __m128i bottom = _mm_loadu_si128((const __m128i *) (dst - stride + extWidth - 8));
    _mm_storeu_si128((__m128i *) (dst + extWidth - 8), bottom);
    _mm_storeu_si128((__m128i *) (dst + stride + extWidth - 8), bottom);
  }
  else
  {
    CHECK(false, "padding size must be 1 or 2");
  }
}

template< X86_VEXT vext >
void addBIOAvg4_SSE(const Pel* src0, int src0Stride, const Pel* src1, int src1Stride, Pel *dst, int dstStride, const Pel *gradX0, const Pel *gradX1, const Pel *gradY0, const Pel*gradY1, int gradStride, int width, int height, int tmpx, int tmpy, int shift, int offset, const ClpRng& clpRng)
{
  __m128i c        = _mm_unpacklo_epi16(_mm_set1_epi16(tmpx), _mm_set1_epi16(tmpy));
  __m128i vibdimin = _mm_set1_epi16(clpRng.min);
  __m128i vibdimax = _mm_set1_epi16(clpRng.max);

  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x += 4)
    {
      __m128i a   = _mm_unpacklo_epi16(_mm_loadl_epi64((const __m128i *) (gradX0 + x)),
                                     _mm_loadl_epi64((const __m128i *) (gradY0 + x)));
      __m128i b   = _mm_unpacklo_epi16(_mm_loadl_epi64((const __m128i *) (gradX1 + x)),
                                     _mm_loadl_epi64((const __m128i *) (gradY1 + x)));
      a           = _mm_sub_epi16(a, b);
      __m128i sum = _mm_madd_epi16(a, c);

      a   = _mm_unpacklo_epi16(_mm_loadl_epi64((const __m128i *) (src0 + x)),
                             _mm_loadl_epi64((const __m128i *) (src1 + x)));
      sum = _mm_add_epi32(sum, _mm_madd_epi16(a, _mm_set1_epi16(1)));
      sum = _mm_add_epi32(sum, _mm_set1_epi32(offset));
      sum = _mm_sra_epi32(sum, _mm_cvtsi32_si128(shift));
      sum = _mm_packs_epi32(sum, sum);
      sum = _mm_max_epi16(sum, vibdimin);
      sum = _mm_min_epi16(sum, vibdimax);
      _mm_storel_epi64((__m128i *) (dst + x), sum);
    }
    dst += dstStride;       src0 += src0Stride;     src1 += src1Stride;
    gradX0 += gradStride; gradX1 += gradStride; gradY0 += gradStride; gradY1 += gradStride;
  }
}

#if MULTI_PASS_DMVR || SAMPLE_BASED_BDOF
template< X86_VEXT vext >
void calcBIOParameter_SSE(const Pel* srcY0Tmp, const Pel* srcY1Tmp, Pel* gradX0, Pel* gradX1, Pel* gradY0, Pel* gradY1, int width, int height, const int src0Stride, const int src1Stride, const int widthG, const int bitDepth, Pel* absGX, Pel* absGY, Pel* dIX, Pel* dIY, Pel* signGY_GX, Pel* dI)
{
  width -= 2;
  height -= 2;
  const int bioParamOffset = widthG + 1;
  srcY0Tmp += src0Stride + 1;
  srcY1Tmp += src1Stride + 1;
  gradX0 += bioParamOffset;  gradX1 += bioParamOffset;
  gradY0 += bioParamOffset;  gradY1 += bioParamOffset;
  absGX  += bioParamOffset;  absGY  += bioParamOffset;
  dIX    += bioParamOffset;  dIY    += bioParamOffset;
  signGY_GX += bioParamOffset;
  if (dI)
  {
    dI += bioParamOffset;
  }
  int shift4 = 4;
  int shift5 = 1;
#ifdef USE_AVX2
  if (width == 8)
  {
    for (int y = 0; y < height; y++)
    {
      __m128i shiftSrcY0Tmp = _mm_srai_epi16(_mm_loadu_si128((const __m128i*)srcY0Tmp), shift4);
      __m128i shiftSrcY1Tmp = _mm_srai_epi16(_mm_loadu_si128((const __m128i*)srcY1Tmp), shift4);
      __m128i sumGradX = _mm_add_epi16(_mm_loadu_si128((const __m128i*)gradX0), _mm_loadu_si128((const __m128i*)gradX1));
      __m128i sumGradY = _mm_add_epi16(_mm_loadu_si128((const __m128i*)gradY0), _mm_loadu_si128((const __m128i*)gradY1));

      __m128i subTemp1 = _mm_sub_epi16(shiftSrcY1Tmp, shiftSrcY0Tmp);
      __m128i packTempX = _mm_srai_epi16(sumGradX, shift5);
      __m128i packTempY = _mm_srai_epi16(sumGradY, shift5);

      __m128i gX_tmp = _mm_abs_epi16(packTempX);
      __m128i gY_tmp = _mm_abs_epi16(packTempY);
      if (dI)
      {
        _mm_storeu_si128((__m128i *) dI, subTemp1);
      }
      __m128i dIX_tmp       = _mm_sign_epi16(subTemp1,  packTempX );
      __m128i dIY_tmp       = _mm_sign_epi16(subTemp1,  packTempY );
      __m128i signGY_GX_tmp = _mm_sign_epi16(packTempX, packTempY );

      _mm_storeu_si128( ( __m128i * )absGX, gX_tmp );
      _mm_storeu_si128( ( __m128i * )absGY, gY_tmp );
      _mm_storeu_si128( ( __m128i * )dIX, dIX_tmp );
      _mm_storeu_si128( ( __m128i * )dIY, dIY_tmp );
      _mm_storeu_si128( ( __m128i * )signGY_GX, signGY_GX_tmp );
      srcY0Tmp += src0Stride;
      srcY1Tmp += src1Stride;
      gradX0 += widthG;
      gradX1 += widthG;
      gradY0 += widthG;
      gradY1 += widthG;
      absGX += widthG;
      absGY += widthG;
      if (dI)
      {
        dI += widthG;
      }
      dIX += widthG;
      dIY += widthG;
      signGY_GX += widthG;
    }
  }
  else // width = 12, 20, 36, 68, 132, 260
  {
    for (int y = 0; y < height; y++)
    {
      for (int x = 0; x < width; x += 16)
      {
        __m256i shiftSrcY0Tmp = _mm256_srai_epi16(_mm256_loadu_si256((const __m256i*)(srcY0Tmp + x)), shift4);
        __m256i shiftSrcY1Tmp = _mm256_srai_epi16(_mm256_loadu_si256((const __m256i*)(srcY1Tmp + x)), shift4);
        __m256i sumGradX = _mm256_add_epi16(_mm256_loadu_si256((const __m256i*)(gradX0 + x)), _mm256_loadu_si256((const __m256i*)(gradX1 + x)));
        __m256i sumGradY = _mm256_add_epi16(_mm256_loadu_si256((const __m256i*)(gradY0 + x)), _mm256_loadu_si256((const __m256i*)(gradY1 + x)));

        __m256i subTemp1 = _mm256_sub_epi16(shiftSrcY1Tmp, shiftSrcY0Tmp);
        __m256i packTempX = _mm256_srai_epi16(sumGradX, shift5);
        __m256i packTempY = _mm256_srai_epi16(sumGradY, shift5);

        __m256i gX_tmp = _mm256_abs_epi16(packTempX);
        __m256i gY_tmp = _mm256_abs_epi16(packTempY);
        if (dI)
        {
          _mm256_storeu_si256((__m256i *) (dI + x), subTemp1);
        }
        __m256i dIX_tmp       = _mm256_sign_epi16(subTemp1,  packTempX );
        __m256i dIY_tmp       = _mm256_sign_epi16(subTemp1,  packTempY );
        __m256i signGY_GX_tmp = _mm256_sign_epi16(packTempX, packTempY );

        _mm256_storeu_si256( ( __m256i * ) ( absGX + x ), gX_tmp );
        _mm256_storeu_si256( ( __m256i * ) ( absGY + x ), gY_tmp );
        _mm256_storeu_si256( ( __m256i * ) ( dIX + x ), dIX_tmp );
        _mm256_storeu_si256( ( __m256i * ) ( dIY + x ), dIY_tmp );
        _mm256_storeu_si256( ( __m256i * ) ( signGY_GX + x ), signGY_GX_tmp );
      }
      srcY0Tmp += src0Stride;
      srcY1Tmp += src1Stride;
      gradX0 += widthG;
      gradX1 += widthG;
      gradY0 += widthG;
      gradY1 += widthG;
      absGX += widthG;
      absGY += widthG;
      if (dI)
      {
        dI += widthG;
      }
      dIX += widthG;
      dIY += widthG;
      signGY_GX += widthG;
    }
  }
#else
  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x += 8)
    {
      __m128i shiftSrcY0Tmp = _mm_srai_epi16(_mm_loadu_si128((const __m128i*)(srcY0Tmp + x)), shift4);
      __m128i shiftSrcY1Tmp = _mm_srai_epi16(_mm_loadu_si128((const __m128i*)(srcY1Tmp + x)), shift4);
      __m128i sumGradX = _mm_add_epi16(_mm_loadu_si128((const __m128i*)(gradX0 + x)), _mm_loadu_si128((const __m128i*)(gradX1 + x)));
      __m128i sumGradY = _mm_add_epi16(_mm_loadu_si128((const __m128i*)(gradY0 + x)), _mm_loadu_si128((const __m128i*)(gradY1 + x)));

      __m128i subTemp1 = _mm_sub_epi16(shiftSrcY1Tmp, shiftSrcY0Tmp);
      __m128i packTempX = _mm_srai_epi16(sumGradX, shift5);
      __m128i packTempY = _mm_srai_epi16(sumGradY, shift5);

      __m128i gX_tmp = _mm_abs_epi16(packTempX);
      __m128i gY_tmp = _mm_abs_epi16(packTempY);
      if (dI)
      {
        _mm_storeu_si128((__m128i *) (dI + x), subTemp1);
      }
      __m128i dIX_tmp       = _mm_sign_epi16(subTemp1,  packTempX );
      __m128i dIY_tmp       = _mm_sign_epi16(subTemp1,  packTempY );
      __m128i signGY_GX_tmp = _mm_sign_epi16(packTempX, packTempY );

      _mm_storeu_si128( ( __m128i * ) ( absGX + x ), gX_tmp );
      _mm_storeu_si128( ( __m128i * ) ( absGY + x ), gY_tmp );
      _mm_storeu_si128( ( __m128i * ) ( dIX + x ), dIX_tmp );
      _mm_storeu_si128( ( __m128i * ) ( dIY + x ), dIY_tmp );
      _mm_storeu_si128( ( __m128i * ) ( signGY_GX + x ), signGY_GX_tmp );
    }
    srcY0Tmp += src0Stride;
    srcY1Tmp += src1Stride;
    gradX0 += widthG;
    gradX1 += widthG;
    gradY0 += widthG;
    gradY1 += widthG;
    absGX += widthG;
    absGY += widthG;
    dIX += widthG;
    dIY += widthG;
    if (dI)
    {
      dI += widthG;
    }
    signGY_GX += widthG;
  }
#endif
}
template< X86_VEXT vext >
void calcBIOParamSum5_SSE(Pel* absGX, Pel* absGY, Pel* dIX, Pel* dIY, Pel* signGY_GX, const int widthG, const int width, const int height, int* sumAbsGX, int* sumAbsGY, int* sumDIX, int* sumDIY, int* sumSignGY_GX)
{
  __m128i vzero = _mm_setzero_si128();
  __m128i vmask = _mm_setr_epi16(1, 1, 1, 1, 1, 0, 0, 0);

  __m128i sumAbsGXTmp16 = vzero;
  __m128i sumDIXTmp16 = vzero;
  __m128i sumAbsGYTmp16 = vzero;
  __m128i sumDIYTmp16 = vzero;
  __m128i sumSignGyGxTmp16 = vzero;
  __m128i sumAbsGXTmp32 = vzero;
  __m128i sumDIXTmp32 = vzero;
  __m128i sumAbsGYTmp32 = vzero;
  __m128i sumDIYTmp32 = vzero;
  __m128i sumSignGyGxTmp32 = vzero;
  __m128i a12 = vzero;
  __m128i a3 = vzero;
  __m128i b12 = vzero;
  __m128i b3 = vzero;
  __m128i c1 = vzero;

  const int width_height = (width * height);
  const int width_N = (1 - width_height);
  const int widthG_height = (widthG * height);
  const int widthG_2 = (widthG << 1);
  const int widthG_3 = widthG_2 + widthG;
  const int widthG_4 = widthG_3 + widthG;
  const int widthG_N = (1 - widthG_height);
  for (int x = 0; x < width; x++)
  {
    sumAbsGXTmp16  = _mm_add_epi16(_mm_loadu_si128((const __m128i*)absGX), _mm_loadu_si128((const __m128i*)(absGX + widthG)));
    sumAbsGXTmp16  = _mm_add_epi16(sumAbsGXTmp16, _mm_loadu_si128((const __m128i*)(absGX + widthG_2)));
    sumAbsGXTmp16  = _mm_add_epi16(sumAbsGXTmp16, _mm_loadu_si128((const __m128i*)(absGX + widthG_3)));
    sumDIXTmp16  = _mm_add_epi16(_mm_loadu_si128((const __m128i*)dIX), _mm_loadu_si128((const __m128i*)(dIX + widthG)));
    sumDIXTmp16  = _mm_add_epi16(sumDIXTmp16, _mm_loadu_si128((const __m128i*)(dIX + widthG_2)));
    sumDIXTmp16  = _mm_add_epi16(sumDIXTmp16, _mm_loadu_si128((const __m128i*)(dIX + widthG_3)));
    sumAbsGYTmp16  = _mm_add_epi16(_mm_loadu_si128((const __m128i*)absGY), _mm_loadu_si128((const __m128i*)(absGY + widthG)));
    sumAbsGYTmp16  = _mm_add_epi16(sumAbsGYTmp16, _mm_loadu_si128((const __m128i*)(absGY + widthG_2)));
    sumAbsGYTmp16  = _mm_add_epi16(sumAbsGYTmp16, _mm_loadu_si128((const __m128i*)(absGY + widthG_3)));
    sumDIYTmp16  = _mm_add_epi16(_mm_loadu_si128((const __m128i*)dIY), _mm_loadu_si128((const __m128i*)(dIY + widthG)));
    sumDIYTmp16  = _mm_add_epi16(sumDIYTmp16, _mm_loadu_si128((const __m128i*)(dIY + widthG_2)));
    sumDIYTmp16  = _mm_add_epi16(sumDIYTmp16, _mm_loadu_si128((const __m128i*)(dIY + widthG_3)));
    sumSignGyGxTmp16  = _mm_add_epi16(_mm_loadu_si128((const __m128i*)signGY_GX), _mm_loadu_si128((const __m128i*)(signGY_GX + widthG)));
    sumSignGyGxTmp16  = _mm_add_epi16(sumSignGyGxTmp16, _mm_loadu_si128((const __m128i*)(signGY_GX + widthG_2)));
    sumSignGyGxTmp16  = _mm_add_epi16(sumSignGyGxTmp16, _mm_loadu_si128((const __m128i*)(signGY_GX + widthG_3)));

    __m128i absGXOneRow = vzero;
    __m128i dIXOneRow = vzero;
    __m128i absGYOneRow = vzero;
    __m128i dIYOneRow = vzero;
    __m128i signGyGxOneRow = vzero;
    for (int y = 0; y < height; y++)
    {
      sumAbsGXTmp16 = _mm_sub_epi16(sumAbsGXTmp16, absGXOneRow);
      absGXOneRow = _mm_loadu_si128((const __m128i*)absGX);
      sumAbsGXTmp16 = _mm_add_epi16(sumAbsGXTmp16, _mm_loadu_si128((const __m128i*)(absGX + widthG_4)));
      sumDIXTmp16  = _mm_sub_epi16(sumDIXTmp16, dIXOneRow);
      dIXOneRow = _mm_loadu_si128((const __m128i*)dIX);
      sumDIXTmp16  = _mm_add_epi16(sumDIXTmp16, _mm_loadu_si128((const __m128i*)(dIX + widthG_4)));
      sumAbsGYTmp16  = _mm_sub_epi16(sumAbsGYTmp16, absGYOneRow);
      absGYOneRow = _mm_loadu_si128((const __m128i*)absGY);
      sumAbsGYTmp16  = _mm_add_epi16(sumAbsGYTmp16, _mm_loadu_si128((const __m128i*)(absGY + widthG_4)));
      sumDIYTmp16  = _mm_sub_epi16(sumDIYTmp16, dIYOneRow);
      dIYOneRow = _mm_loadu_si128((const __m128i*)dIY);
      sumDIYTmp16  = _mm_add_epi16(sumDIYTmp16, _mm_loadu_si128((const __m128i*)(dIY + widthG_4)));
      sumSignGyGxTmp16  = _mm_sub_epi16(sumSignGyGxTmp16, signGyGxOneRow);
      signGyGxOneRow = _mm_loadu_si128((const __m128i*)signGY_GX);
      sumSignGyGxTmp16  = _mm_add_epi16(sumSignGyGxTmp16, _mm_loadu_si128((const __m128i*)(signGY_GX+ widthG_4)));

      sumAbsGXTmp32 = _mm_madd_epi16(sumAbsGXTmp16, vmask);
      sumAbsGYTmp32 = _mm_madd_epi16(sumAbsGYTmp16, vmask);
      sumDIXTmp32 = _mm_madd_epi16(sumDIXTmp16, vmask);
      sumDIYTmp32 = _mm_madd_epi16(sumDIYTmp16, vmask);
      a12 = _mm_unpacklo_epi32(sumAbsGXTmp32, sumAbsGYTmp32);
      a3  = _mm_unpackhi_epi32(sumAbsGXTmp32, sumAbsGYTmp32);
      b12 = _mm_unpacklo_epi32(sumDIXTmp32, sumDIYTmp32);
      b3  = _mm_unpackhi_epi32(sumDIXTmp32, sumDIYTmp32);

      c1  = _mm_unpacklo_epi64(a12, b12);
      c1 = _mm_add_epi32(c1, _mm_unpackhi_epi64(a12, b12));
      c1 = _mm_add_epi32(c1, _mm_unpacklo_epi64(a3, b3));

      *sumAbsGX = _mm_cvtsi128_si32(c1);
      *sumAbsGY = _mm_cvtsi128_si32(_mm_shuffle_epi32(c1, 0x55));
      *sumDIX   = _mm_cvtsi128_si32(_mm_shuffle_epi32(c1, 0xaa));
      *sumDIY   = _mm_cvtsi128_si32(_mm_shuffle_epi32(c1, 0xff));

      sumSignGyGxTmp32 = _mm_madd_epi16(sumSignGyGxTmp16, vmask);
      sumSignGyGxTmp32 = _mm_add_epi32(sumSignGyGxTmp32, _mm_shuffle_epi32(sumSignGyGxTmp32, 0x4e));   // 01001110
      sumSignGyGxTmp32 = _mm_add_epi32(sumSignGyGxTmp32, _mm_shuffle_epi32(sumSignGyGxTmp32, 0xb1));   // 10110001
      *sumSignGY_GX = _mm_cvtsi128_si32(sumSignGyGxTmp32);

      // bio parameter increment
      absGX += widthG;
      absGY += widthG;
      dIX += widthG;
      dIY += widthG;
      signGY_GX += widthG;
      // sum parameter increment
      sumAbsGX += width;
      sumAbsGY += width;
      sumDIX += width;
      sumDIY += width;
      sumSignGY_GX += width;
    }
    // bio parameter back to first row
    absGX += widthG_N;
    absGY += widthG_N;
    dIX += widthG_N;
    dIY += widthG_N;
    signGY_GX += widthG_N;
    // sum parameter back to first row
    sumAbsGX += width_N;
    sumAbsGY += width_N;
    sumDIX += width_N;
    sumDIY += width_N;
    sumSignGY_GX += width_N;
  }
  sumDIX -= width;
  sumDIY -= width;
#ifdef USE_AVX2
  for (int idx = 0; idx < width_height; idx += 8)
  {
    __m256i sumDIXTmp256 = _mm256_loadu_si256((const __m256i*)sumDIX);
    __m256i sumDIYTmp256 = _mm256_loadu_si256((const __m256i*)sumDIY);
    _mm256_storeu_si256( ( __m256i * )sumDIX, _mm256_slli_epi32(sumDIXTmp256, 2) );
    _mm256_storeu_si256( ( __m256i * )sumDIY, _mm256_slli_epi32(sumDIYTmp256, 2) );
    sumDIX += 8;
    sumDIY += 8;
  }
#else
  for (int idx = 0; idx < width_height; idx += 4)
  {
    sumDIXTmp32 = _mm_loadu_si128((const __m128i*)sumDIX);
    sumDIYTmp32 = _mm_loadu_si128((const __m128i*)sumDIY);
    _mm_storeu_si128( ( __m128i * )sumDIX, _mm_slli_epi32(sumDIXTmp32, 2) );
    _mm_storeu_si128( ( __m128i * )sumDIY, _mm_slli_epi32(sumDIYTmp32, 2) );
    sumDIX += 4;
    sumDIY += 4;
  }
#endif
}
template< X86_VEXT vext >
void calcBIOParamSum4_SSE(Pel* absGX, Pel* absGY, Pel* dIX, Pel* dIY, Pel* signGY_GX, int width, int height, const int widthG, int* sumAbsGX, int* sumAbsGY, int* sumDIX, int* sumDIY, int* sumSignGY_GX)
{
  __m128i vzero = _mm_setzero_si128();

  __m128i sumAbsGXTmp16 = vzero;
  __m128i sumDIXTmp16 = vzero;
  __m128i sumAbsGYTmp16 = vzero;
  __m128i sumDIYTmp16 = vzero;
  __m128i sumSignGyGxTmp16 = vzero;
  if (width == 8)
  {
    for (int y = 0; y < height; y++)
    {
      sumAbsGXTmp16  = _mm_add_epi16(sumAbsGXTmp16, _mm_loadu_si128((const __m128i*)absGX));
      sumDIXTmp16  = _mm_add_epi16(sumDIXTmp16, _mm_loadu_si128((const __m128i*)dIX));
      sumAbsGYTmp16  = _mm_add_epi16(sumAbsGYTmp16, _mm_loadu_si128((const __m128i*)absGY));
      sumDIYTmp16  = _mm_add_epi16(sumDIYTmp16, _mm_loadu_si128((const __m128i*)dIY));
      sumSignGyGxTmp16  = _mm_add_epi16(sumSignGyGxTmp16, _mm_loadu_si128((const __m128i*)signGY_GX));
      // bio parameter increment
      absGX += widthG;
      absGY += widthG;
      dIX += widthG;
      dIY += widthG;
      signGY_GX += widthG;
    }
  }
  else // (width == 12)
  {
    for (int y = 0; y < height; y++)
    {
      sumAbsGXTmp16  = _mm_add_epi16(sumAbsGXTmp16, _mm_loadu_si128((const __m128i*)absGX));
      sumAbsGXTmp16  = _mm_add_epi16(sumAbsGXTmp16, _mm_loadl_epi64((const __m128i*)(absGX + 8)));
      sumDIXTmp16  = _mm_add_epi16(sumDIXTmp16, _mm_loadu_si128((const __m128i*)dIX));
      sumDIXTmp16  = _mm_add_epi16(sumDIXTmp16, _mm_loadl_epi64((const __m128i*)(dIX + 8)));
      sumAbsGYTmp16  = _mm_add_epi16(sumAbsGYTmp16, _mm_loadu_si128((const __m128i*)absGY));
      sumAbsGYTmp16  = _mm_add_epi16(sumAbsGYTmp16, _mm_loadl_epi64((const __m128i*)(absGY + 8)));
      sumDIYTmp16  = _mm_add_epi16(sumDIYTmp16, _mm_loadu_si128((const __m128i*)dIY));
      sumDIYTmp16  = _mm_add_epi16(sumDIYTmp16, _mm_loadl_epi64((const __m128i*)(dIY + 8)));
      sumSignGyGxTmp16  = _mm_add_epi16(sumSignGyGxTmp16, _mm_loadu_si128((const __m128i*)signGY_GX));
      sumSignGyGxTmp16  = _mm_add_epi16(sumSignGyGxTmp16, _mm_loadl_epi64((const __m128i*)(signGY_GX + 8)));
      // bio parameter increment
      absGX += widthG;
      absGY += widthG;
      dIX += widthG;
      dIY += widthG;
      signGY_GX += widthG;
    }
  }

  __m128i sumAbsGXTmp32 = _mm_add_epi32(_mm_unpacklo_epi16(sumAbsGXTmp16, vzero), _mm_unpackhi_epi16(sumAbsGXTmp16, vzero));
  __m128i sumAbsGYTmp32 = _mm_add_epi32(_mm_unpacklo_epi16(sumAbsGYTmp16, vzero), _mm_unpackhi_epi16(sumAbsGYTmp16, vzero));
  __m128i sumDIXTmp32 = _mm_madd_epi16(sumDIXTmp16, _mm_set1_epi16(1));
  __m128i sumDIYTmp32 = _mm_madd_epi16(sumDIYTmp16, _mm_set1_epi16(1));
  __m128i a12 = _mm_unpacklo_epi32(sumAbsGXTmp32, sumAbsGYTmp32);
  __m128i a3  = _mm_unpackhi_epi32(sumAbsGXTmp32, sumAbsGYTmp32);
  __m128i b12 = _mm_unpacklo_epi32(sumDIXTmp32, sumDIYTmp32);
  __m128i b3  = _mm_unpackhi_epi32(sumDIXTmp32, sumDIYTmp32);

  __m128i c1  = _mm_unpacklo_epi64(a12, b12);
  c1 = _mm_add_epi32(c1, _mm_unpackhi_epi64(a12, b12));
  c1 = _mm_add_epi32(c1, _mm_unpacklo_epi64(a3, b3));
  c1 = _mm_add_epi32(c1, _mm_unpackhi_epi64(a3, b3));

  *sumAbsGX = _mm_cvtsi128_si32(c1);
  *sumAbsGY = _mm_cvtsi128_si32(_mm_shuffle_epi32(c1, 0x55));
  *sumDIX   = _mm_cvtsi128_si32(_mm_shuffle_epi32(c1, 0xaa));
  *sumDIY   = _mm_cvtsi128_si32(_mm_shuffle_epi32(c1, 0xff));

  __m128i sumSignGyGxTmp32 = _mm_madd_epi16(sumSignGyGxTmp16, _mm_set1_epi16(1));
  sumSignGyGxTmp32 = _mm_add_epi32(sumSignGyGxTmp32, _mm_shuffle_epi32(sumSignGyGxTmp32, 0x4e));   // 01001110
  sumSignGyGxTmp32 = _mm_add_epi32(sumSignGyGxTmp32, _mm_shuffle_epi32(sumSignGyGxTmp32, 0xb1));   // 10110001
  *sumSignGY_GX = _mm_cvtsi128_si32(sumSignGyGxTmp32);
}

template< X86_VEXT vext >
void calcBIOClippedVxVy_SSE(int* sumDIX_pixel_32bit, int* sumAbsGX_pixel_32bit, int* sumDIY_pixel_32bit, int* sumAbsGY_pixel_32bit, int* sumSignGY_GX_pixel_32bit, const int limit, const int bioSubblockSize, int* tmpx_pixel_32bit, int* tmpy_pixel_32bit)
{
#ifdef USE_AVX2
  __m256i vibdimin = _mm256_set1_epi32(-limit);
  __m256i vibdimax = _mm256_set1_epi32(limit);
  __m256i tmp256  = _mm256_setzero_si256();

  for (int idx = 0; idx < bioSubblockSize; idx += 8)
  {
    tmp256 = _mm256_loadu_si256((const __m256i*)sumDIX_pixel_32bit);
    tmp256 = _mm256_srav_epi32(tmp256, _mm256_loadu_si256((const __m256i*)sumAbsGX_pixel_32bit));
    tmp256 = _mm256_max_epi32(tmp256, vibdimin);
    tmp256 = _mm256_min_epi32(tmp256, vibdimax);
    _mm256_storeu_si256( ( __m256i * )tmpx_pixel_32bit, tmp256);
    tmp256 = _mm256_mullo_epi32(tmp256, _mm256_loadu_si256((const __m256i*)sumSignGY_GX_pixel_32bit));
    tmp256 = _mm256_srai_epi32(tmp256, 1);
    tmp256 = _mm256_sub_epi32(_mm256_loadu_si256((const __m256i*)sumDIY_pixel_32bit), tmp256);
    tmp256 = _mm256_srav_epi32(tmp256, _mm256_loadu_si256((const __m256i*)sumAbsGY_pixel_32bit));
    tmp256 = _mm256_max_epi32(tmp256, vibdimin);
    tmp256 = _mm256_min_epi32(tmp256, vibdimax);
    _mm256_storeu_si256( ( __m256i * )tmpy_pixel_32bit, tmp256);
    sumDIX_pixel_32bit += 8;
    sumAbsGX_pixel_32bit += 8;
    sumDIY_pixel_32bit += 8;
    sumAbsGY_pixel_32bit += 8;
    sumSignGY_GX_pixel_32bit += 8;
    tmpx_pixel_32bit += 8;
    tmpy_pixel_32bit += 8;
  }
#else
  __m128i vibdimin = _mm_set1_epi32(-limit);
  __m128i vibdimax = _mm_set1_epi32(limit);
  __m128i tmp128  = _mm_setzero_si128();

  for (int idx = 0; idx < bioSubblockSize; idx += 4)
  {
    *sumDIX_pixel_32bit = (*sumDIX_pixel_32bit) >> (*sumAbsGX_pixel_32bit);
    *(sumDIX_pixel_32bit + 1) = (*(sumDIX_pixel_32bit + 1)) >> (*(sumAbsGX_pixel_32bit + 1));
    *(sumDIX_pixel_32bit + 2) = (*(sumDIX_pixel_32bit + 2)) >> (*(sumAbsGX_pixel_32bit + 2));
    *(sumDIX_pixel_32bit + 3) = (*(sumDIX_pixel_32bit + 3)) >> (*(sumAbsGX_pixel_32bit + 3));
    tmp128 = _mm_loadu_si128((const __m128i*)sumDIX_pixel_32bit);
    tmp128 = _mm_max_epi32(tmp128, vibdimin);
    tmp128 = _mm_min_epi32(tmp128, vibdimax);
    _mm_storeu_si128( ( __m128i * )tmpx_pixel_32bit, tmp128);
    tmp128 = _mm_mullo_epi32(tmp128, _mm_loadu_si128((const __m128i*)sumSignGY_GX_pixel_32bit));
    tmp128 = _mm_srai_epi32(tmp128, 1);
    tmp128 = _mm_sub_epi32(_mm_loadu_si128((const __m128i*)sumDIY_pixel_32bit), tmp128);
    _mm_storeu_si128( ( __m128i * )sumDIY_pixel_32bit, tmp128);
    *sumDIY_pixel_32bit = (*sumDIY_pixel_32bit) >> (*sumAbsGY_pixel_32bit);
    *(sumDIY_pixel_32bit + 1) = (*(sumDIY_pixel_32bit + 1)) >> (*(sumAbsGY_pixel_32bit + 1));
    *(sumDIY_pixel_32bit + 2) = (*(sumDIY_pixel_32bit + 2)) >> (*(sumAbsGY_pixel_32bit + 2));
    *(sumDIY_pixel_32bit + 3) = (*(sumDIY_pixel_32bit + 3)) >> (*(sumAbsGY_pixel_32bit + 3));
    tmp128 = _mm_loadu_si128((const __m128i*)sumDIY_pixel_32bit);
    tmp128 = _mm_max_epi32(tmp128, vibdimin);
    tmp128 = _mm_min_epi32(tmp128, vibdimax);
    _mm_storeu_si128( ( __m128i * )tmpy_pixel_32bit, tmp128);
    sumDIX_pixel_32bit += 4;
    sumAbsGX_pixel_32bit += 4;
    sumDIY_pixel_32bit += 4;
    sumAbsGY_pixel_32bit += 4;
    sumSignGY_GX_pixel_32bit += 4;
    tmpx_pixel_32bit += 4;
    tmpy_pixel_32bit += 4;
  }
#endif
}

template< X86_VEXT vext >
void addBIOAvgN_SSE(const Pel* src0, int src0Stride, const Pel* src1, int src1Stride, Pel *dst, int dstStride, const Pel *gradX0, const Pel *gradX1, const Pel *gradY0, const Pel *gradY1, int gradStride, int width, int height, int* tmpx, int* tmpy, int shift, int offset, const ClpRng& clpRng)
{
#ifdef USE_AVX2
  __m128i vibdimin = _mm_set1_epi16(clpRng.min);
  __m128i vibdimax = _mm_set1_epi16(clpRng.max);
  __m128i var1 = _mm_setzero_si128();
  __m128i var2 = _mm_setzero_si128();

  if (width == 8)
  {
    for (int y = 0; y < height; y++)
    {
      var1        = _mm_loadu_si128((const __m128i *)gradX0);
      var1        = _mm_sub_epi16(var1,  _mm_loadu_si128((const __m128i *)gradX1));
      var2        = _mm_loadu_si128((const __m128i *)gradY0);
      var2        = _mm_sub_epi16(var2,  _mm_loadu_si128((const __m128i *)gradY1));
      __m256i a   = _mm256_inserti128_si256(_mm256_castsi128_si256(_mm_unpacklo_epi16(var1, var2)),
                                            _mm_unpackhi_epi16(var1, var2), 1);

      __m256i vClipTmpx256 = _mm256_loadu_si256((const __m256i *)tmpx);
      __m256i vClipTmpy256 = _mm256_loadu_si256((const __m256i *)tmpy);
      __m256i c   = _mm256_packs_epi32 (_mm256_unpacklo_epi32(vClipTmpx256, vClipTmpy256), _mm256_unpackhi_epi32(vClipTmpx256, vClipTmpy256));

      __m256i sum32 = _mm256_madd_epi16(a, c);  // 8 32bits

      var1        = _mm_loadu_si128((const __m128i *)src0);
      var2        = _mm_loadu_si128((const __m128i *)src1);
      a           = _mm256_inserti128_si256(_mm256_castsi128_si256(_mm_unpacklo_epi16(var1, var2)),
                                            _mm_unpackhi_epi16(var1, var2), 1);

      sum32 = _mm256_add_epi32(sum32, _mm256_madd_epi16(a, _mm256_set1_epi16(1)));
      sum32 = _mm256_add_epi32(sum32, _mm256_set1_epi32(offset));
      sum32 = _mm256_sra_epi32(sum32, _mm_cvtsi32_si128(shift));
      __m128i sum16 = _mm_packs_epi32(_mm256_castsi256_si128(sum32), _mm256_extractf128_si256(sum32, 1));
      sum16 = _mm_max_epi16(sum16, vibdimin);
      sum16 = _mm_min_epi16(sum16, vibdimax);
      _mm_storeu_si128((__m128i *)dst, sum16);
      dst += dstStride;       src0 += src0Stride;     src1 += src1Stride;
      gradX0 += gradStride; gradX1 += gradStride; gradY0 += gradStride; gradY1 += gradStride;
      tmpx += width; tmpy += width;
    }
  }
  else // (width == 4)
  {
    for (int y = 0; y < height; y++)
    {
      __m128i a   = _mm_unpacklo_epi16(_mm_loadl_epi64((const __m128i *)gradX0),
                                     _mm_loadl_epi64((const __m128i *)gradY0));
      __m128i b   = _mm_unpacklo_epi16(_mm_loadl_epi64((const __m128i *)gradX1),
                                     _mm_loadl_epi64((const __m128i *)gradY1));
      a           = _mm_sub_epi16(a, b);

      __m128i vClipTmpx128   = _mm_loadu_si128((const __m128i*)tmpx);
      __m128i vClipTmpy128   = _mm_loadu_si128((const __m128i*)tmpy);
      __m128i c   = _mm_packs_epi32(_mm_unpacklo_epi32(vClipTmpx128, vClipTmpy128), _mm_unpackhi_epi32(vClipTmpx128, vClipTmpy128));

      __m128i sum = _mm_madd_epi16(a, c);

      a   = _mm_unpacklo_epi16(_mm_loadl_epi64((const __m128i *)src0),
                             _mm_loadl_epi64((const __m128i *)src1));
      sum = _mm_add_epi32(sum, _mm_madd_epi16(a, _mm_set1_epi16(1)));
      sum = _mm_add_epi32(sum, _mm_set1_epi32(offset));
      sum = _mm_sra_epi32(sum, _mm_cvtsi32_si128(shift));
      sum = _mm_packs_epi32(sum, sum);
      sum = _mm_max_epi16(sum, vibdimin);
      sum = _mm_min_epi16(sum, vibdimax);
      _mm_storel_epi64((__m128i *)dst, sum);
      dst += dstStride;       src0 += src0Stride;     src1 += src1Stride;
      gradX0 += gradStride; gradX1 += gradStride; gradY0 += gradStride; gradY1 += gradStride;
      tmpx += width; tmpy += width;
    }
  }
#else
  __m128i vibdimin = _mm_set1_epi16(clpRng.min);
  __m128i vibdimax = _mm_set1_epi16(clpRng.max);
  __m128i vClipTmpx = _mm_setzero_si128();
  __m128i vClipTmpy = _mm_setzero_si128();

  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x += 4)
    {
      __m128i a   = _mm_unpacklo_epi16(_mm_loadl_epi64((const __m128i *) (gradX0 + x)),
                                     _mm_loadl_epi64((const __m128i *) (gradY0 + x)));
      __m128i b   = _mm_unpacklo_epi16(_mm_loadl_epi64((const __m128i *) (gradX1 + x)),
                                     _mm_loadl_epi64((const __m128i *) (gradY1 + x)));
      a           = _mm_sub_epi16(a, b);

      vClipTmpx   = _mm_loadu_si128((const __m128i*)(tmpx + x));
      vClipTmpy   = _mm_loadu_si128((const __m128i*)(tmpy + x));
      __m128i c   = _mm_packs_epi32(_mm_unpacklo_epi32(vClipTmpx, vClipTmpy), _mm_unpackhi_epi32(vClipTmpx, vClipTmpy));
      __m128i sum = _mm_madd_epi16(a, c);

      a   = _mm_unpacklo_epi16(_mm_loadl_epi64((const __m128i *) (src0 + x)),
                             _mm_loadl_epi64((const __m128i *) (src1 + x)));
      sum = _mm_add_epi32(sum, _mm_madd_epi16(a, _mm_set1_epi16(1)));
      sum = _mm_add_epi32(sum, _mm_set1_epi32(offset));
      sum = _mm_sra_epi32(sum, _mm_cvtsi32_si128(shift));
      sum = _mm_packs_epi32(sum, sum);
      sum = _mm_max_epi16(sum, vibdimin);
      sum = _mm_min_epi16(sum, vibdimax);
      _mm_storel_epi64((__m128i *) (dst + x), sum);
    }
    dst += dstStride;       src0 += src0Stride;     src1 += src1Stride;
    gradX0 += gradStride; gradX1 += gradStride; gradY0 += gradStride; gradY1 += gradStride;
    tmpx += width; tmpy += width;
  }
#endif
}

template< X86_VEXT vext >
void calAbsSum_SSE(const Pel* diff, int stride, int width, int height, int* absSum)
{
  __m128i vzero = _mm_setzero_si128();
  __m128i vsum32 = vzero;
  __m128i vdiff = vzero;

  if (width == 8)
  {
    for (int y = 0; y < height; y++)
    {
      vdiff = _mm_abs_epi16(_mm_loadu_si128((const __m128i*)(diff)));
      vsum32 = _mm_add_epi32(vsum32, _mm_unpacklo_epi16(vdiff, vzero));
      vsum32 = _mm_add_epi32(vsum32, _mm_unpackhi_epi16(vdiff, vzero));
      diff += stride;
    }
  }
  else // (width == 4)
  {
    const int strideDouble = (stride << 1);
    for (int y = 0; y < height; y += 2)
    {
      vdiff = _mm_abs_epi16(_mm_unpacklo_epi16(_mm_loadl_epi64((const __m128i*)diff), _mm_loadl_epi64((const __m128i*)(diff + stride))));
      vsum32 = _mm_add_epi32(vsum32, _mm_unpacklo_epi16(vdiff, vzero));
      vsum32 = _mm_add_epi32(vsum32, _mm_unpackhi_epi16(vdiff, vzero));
      diff += strideDouble;
    }
  }

  vsum32 = _mm_add_epi32(vsum32, _mm_shuffle_epi32(vsum32, 0x4e));   // 01001110
  vsum32 = _mm_add_epi32(vsum32, _mm_shuffle_epi32(vsum32, 0xb1));   // 10110001
  *absSum = _mm_cvtsi128_si32(vsum32);
}
#endif

template< X86_VEXT vext >
void calcBIOSums_SSE(const Pel* srcY0Tmp, const Pel* srcY1Tmp, Pel* gradX0, Pel* gradX1, Pel* gradY0, Pel* gradY1, int xu, int yu, const int src0Stride, const int src1Stride, const int widthG, const int bitDepth, int* sumAbsGX, int* sumAbsGY, int* sumDIX, int* sumDIY, int* sumSignGY_GX)

{
  int shift4 = 4;
  int shift5 = 1;

  __m128i sumAbsGXTmp = _mm_setzero_si128();
  __m128i sumDIXTmp = _mm_setzero_si128();
  __m128i sumAbsGYTmp = _mm_setzero_si128();
  __m128i sumDIYTmp = _mm_setzero_si128();
  __m128i sumSignGyGxTmp = _mm_setzero_si128();

  for (int y = 0; y < 6; y++)
  {
    // Note: loading 8 values also works, but valgrind doesn't like it
    auto load6values = [](const Pel *ptr) {
      __m128i a = _mm_loadl_epi64((const __m128i *) ptr);
      __m128i b = _mm_cvtsi32_si128(*(uint32_t *) (ptr + 4));
      return _mm_unpacklo_epi64(a, b);
    };

    __m128i shiftSrcY0Tmp = _mm_srai_epi16(load6values(srcY0Tmp), shift4);
    __m128i shiftSrcY1Tmp = _mm_srai_epi16(load6values(srcY1Tmp), shift4);
    __m128i loadGradX0    = load6values(gradX0);
    __m128i loadGradX1    = load6values(gradX1);
    __m128i loadGradY0    = load6values(gradY0);
    __m128i loadGradY1    = load6values(gradY1);

    __m128i subTemp1 = _mm_sub_epi16(shiftSrcY1Tmp, shiftSrcY0Tmp);
    __m128i packTempX = _mm_srai_epi16(_mm_add_epi16(loadGradX0, loadGradX1), shift5);
    __m128i packTempY = _mm_srai_epi16(_mm_add_epi16(loadGradY0, loadGradY1), shift5);
    __m128i gX = _mm_abs_epi16(packTempX);
    __m128i gY = _mm_abs_epi16(packTempY);
    __m128i dIX       = _mm_sign_epi16(subTemp1,  packTempX );
    __m128i dIY       = _mm_sign_epi16(subTemp1,  packTempY );
    __m128i signGY_GX = _mm_sign_epi16(packTempX, packTempY );

    sumAbsGXTmp = _mm_add_epi16(sumAbsGXTmp, gX);
    sumDIXTmp = _mm_add_epi16(sumDIXTmp, dIX);
    sumAbsGYTmp = _mm_add_epi16(sumAbsGYTmp, gY);
    sumDIYTmp = _mm_add_epi16(sumDIYTmp, dIY);
    sumSignGyGxTmp = _mm_add_epi16(sumSignGyGxTmp, signGY_GX);
    srcY0Tmp += src0Stride;
    srcY1Tmp += src1Stride;
    gradX0 += widthG;
    gradX1 += widthG;
    gradY0 += widthG;
    gradY1 += widthG;
  }

  sumAbsGXTmp    = _mm_madd_epi16(sumAbsGXTmp, _mm_setr_epi16(1, 1, 1, 1, 1, 1, 0, 0));
  sumDIXTmp      = _mm_madd_epi16(sumDIXTmp, _mm_setr_epi16(1, 1, 1, 1, 1, 1, 0, 0));
  sumAbsGYTmp    = _mm_madd_epi16(sumAbsGYTmp, _mm_setr_epi16(1, 1, 1, 1, 1, 1, 0, 0));
  sumDIYTmp      = _mm_madd_epi16(sumDIYTmp, _mm_setr_epi16(1, 1, 1, 1, 1, 1, 0, 0));
  sumSignGyGxTmp = _mm_madd_epi16(sumSignGyGxTmp, _mm_setr_epi16(1, 1, 1, 1, 1, 1, 0, 0));

  __m128i a12 = _mm_unpacklo_epi32(sumAbsGXTmp, sumAbsGYTmp);
  __m128i a3  = _mm_unpackhi_epi32(sumAbsGXTmp, sumAbsGYTmp);
  __m128i b12 = _mm_unpacklo_epi32(sumDIXTmp, sumDIYTmp);
  __m128i b3  = _mm_unpackhi_epi32(sumDIXTmp, sumDIYTmp);
  __m128i c1  = _mm_unpacklo_epi64(a12, b12);
  __m128i c2  = _mm_unpackhi_epi64(a12, b12);
  __m128i c3  = _mm_unpacklo_epi64(a3, b3);

  c1 = _mm_add_epi32(c1, c2);
  c1 = _mm_add_epi32(c1, c3);

  *sumAbsGX = _mm_cvtsi128_si32(c1);
  *sumAbsGY = _mm_cvtsi128_si32(_mm_shuffle_epi32(c1, 0x55));
  *sumDIX   = _mm_cvtsi128_si32(_mm_shuffle_epi32(c1, 0xaa));
  *sumDIY   = _mm_cvtsi128_si32(_mm_shuffle_epi32(c1, 0xff));

  sumSignGyGxTmp = _mm_add_epi32(sumSignGyGxTmp, _mm_shuffle_epi32(sumSignGyGxTmp, 0x4e));   // 01001110
  sumSignGyGxTmp = _mm_add_epi32(sumSignGyGxTmp, _mm_shuffle_epi32(sumSignGyGxTmp, 0xb1));   // 10110001
  *sumSignGY_GX  = _mm_cvtsi128_si32(sumSignGyGxTmp);
}

template< X86_VEXT vext >
void applyPROF_SSE(Pel* dstPel, int dstStride, const Pel* srcPel, int srcStride, int width, int height, const Pel* gradX, const Pel* gradY, int gradStride, const int* dMvX, const int* dMvY, int dMvStride, const bool& bi, int shiftNum, Pel offset, const ClpRng& clpRng)
{
  CHECKD((width & 3), "block width error!");

  const int dILimit = 1 << std::max<int>(clpRng.bd + 1, 13);

#ifdef USE_AVX2
  __m256i mm_dmvx, mm_dmvy, mm_gradx, mm_grady, mm_dI, mm_dI0, mm_src;
  __m256i mm_offset = _mm256_set1_epi16(offset);
  __m256i vibdimin = _mm256_set1_epi16(clpRng.min);
  __m256i vibdimax = _mm256_set1_epi16(clpRng.max);
  __m256i mm_dimin = _mm256_set1_epi32(-dILimit);
  __m256i mm_dimax = _mm256_set1_epi32(dILimit - 1);
#else
  __m128i mm_dmvx, mm_dmvy, mm_gradx, mm_grady, mm_dI, mm_dI0;
  __m128i mm_offset = _mm_set1_epi16(offset);
  __m128i vibdimin = _mm_set1_epi16(clpRng.min);
  __m128i vibdimax = _mm_set1_epi16(clpRng.max);
  __m128i mm_dimin = _mm_set1_epi32(-dILimit);
  __m128i mm_dimax = _mm_set1_epi32(dILimit - 1);
#endif

#if USE_AVX2
  for (int h = 0; h < height; h += 4)
#else
  for (int h = 0; h < height; h += 2)
#endif
  {
    const int* vX = dMvX;
    const int* vY = dMvY;
    const Pel* gX = gradX;
    const Pel* gY = gradY;
    const Pel* src = srcPel;
    Pel*       dst = dstPel;

    for (int w = 0; w < width; w += 4)
    {
#if USE_AVX2
      const int *vX0 = vX, *vY0 = vY;
      const Pel *gX0 = gX, *gY0 = gY;

      // first two rows
      mm_dmvx = _mm256_inserti128_si256(_mm256_castsi128_si256(_mm_loadu_si128((const __m128i *)vX0)), _mm_loadu_si128((const __m128i *)(vX0 + dMvStride)), 1);
      mm_dmvy = _mm256_inserti128_si256(_mm256_castsi128_si256(_mm_loadu_si128((const __m128i *)vY0)), _mm_loadu_si128((const __m128i *)(vY0 + dMvStride)), 1);
      mm_gradx = _mm256_inserti128_si256(
        _mm256_castsi128_si256(_mm_cvtepi16_epi32(_mm_loadl_epi64((__m128i*)gX0))),
        _mm_cvtepi16_epi32(_mm_loadl_epi64((__m128i*)(gX0 + gradStride))), 1);
      mm_grady = _mm256_inserti128_si256(
        _mm256_castsi128_si256(_mm_cvtepi16_epi32(_mm_loadl_epi64((__m128i*)gY0))),
        _mm_cvtepi16_epi32(_mm_loadl_epi64((__m128i*)(gY0 + gradStride))), 1);
      mm_dI0 = _mm256_add_epi32(_mm256_mullo_epi32(mm_dmvx, mm_gradx), _mm256_mullo_epi32(mm_dmvy, mm_grady));
      mm_dI0 = _mm256_min_epi32(mm_dimax, _mm256_max_epi32(mm_dimin, mm_dI0));

      // next two rows
      vX0 += (dMvStride << 1); vY0 += (dMvStride << 1); gX0 += (gradStride << 1); gY0 += (gradStride << 1);
      mm_dmvx = _mm256_inserti128_si256(_mm256_castsi128_si256(_mm_loadu_si128((const __m128i *)vX0)), _mm_loadu_si128((const __m128i *)(vX0 + dMvStride)), 1);
      mm_dmvy = _mm256_inserti128_si256(_mm256_castsi128_si256(_mm_loadu_si128((const __m128i *)vY0)), _mm_loadu_si128((const __m128i *)(vY0 + dMvStride)), 1);
      mm_gradx = _mm256_inserti128_si256(
        _mm256_castsi128_si256(_mm_cvtepi16_epi32(_mm_loadl_epi64((__m128i*)gX0))),
        _mm_cvtepi16_epi32(_mm_loadl_epi64((__m128i*)(gX0 + gradStride))), 1);
      mm_grady = _mm256_inserti128_si256(
        _mm256_castsi128_si256(_mm_cvtepi16_epi32(_mm_loadl_epi64((__m128i*)gY0))),
        _mm_cvtepi16_epi32(_mm_loadl_epi64((__m128i*)(gY0 + gradStride))), 1);
      mm_dI = _mm256_add_epi32(_mm256_mullo_epi32(mm_dmvx, mm_gradx), _mm256_mullo_epi32(mm_dmvy, mm_grady));
      mm_dI = _mm256_min_epi32(mm_dimax, _mm256_max_epi32(mm_dimin, mm_dI));

      // combine four rows
      mm_dI = _mm256_packs_epi32(mm_dI0, mm_dI);
      const Pel* src0 = src + srcStride;
      mm_src = _mm256_inserti128_si256(
        _mm256_castsi128_si256(_mm_unpacklo_epi64(_mm_loadl_epi64((const __m128i *)src), _mm_loadl_epi64((const __m128i *)(src + (srcStride << 1))))),
        _mm_unpacklo_epi64(_mm_loadl_epi64((const __m128i *)src0), _mm_loadl_epi64((const __m128i *)(src0 + (srcStride << 1)))),
        1
      );
      mm_dI = _mm256_add_epi16(mm_dI, mm_src);
      if (!bi)
      {
        mm_dI = _mm256_srai_epi16(_mm256_adds_epi16(mm_dI, mm_offset), shiftNum);
        mm_dI = _mm256_min_epi16(vibdimax, _mm256_max_epi16(vibdimin, mm_dI));
      }

      // store final results
      __m128i dITmp = _mm256_extractf128_si256(mm_dI, 1);
      Pel* dst0 = dst;
      _mm_storel_epi64((__m128i *)dst0, _mm256_castsi256_si128(mm_dI));
      dst0 += dstStride; _mm_storel_epi64((__m128i *)dst0, dITmp);
      dst0 += dstStride; _mm_storel_epi64((__m128i *)dst0, _mm_unpackhi_epi64(_mm256_castsi256_si128(mm_dI), _mm256_castsi256_si128(mm_dI)));
      dst0 += dstStride; _mm_storel_epi64((__m128i *)dst0, _mm_unpackhi_epi64(dITmp, dITmp));
#else
      // first row
      mm_dmvx = _mm_loadu_si128((const __m128i *)vX);
      mm_dmvy = _mm_loadu_si128((const __m128i *)vY);
      mm_gradx = _mm_cvtepi16_epi32(_mm_loadl_epi64((__m128i*)gX));
      mm_grady = _mm_cvtepi16_epi32(_mm_loadl_epi64((__m128i*)gY));
      mm_dI0 = _mm_add_epi32(_mm_mullo_epi32(mm_dmvx, mm_gradx), _mm_mullo_epi32(mm_dmvy, mm_grady));
      mm_dI0 = _mm_min_epi32(mm_dimax, _mm_max_epi32(mm_dimin, mm_dI0));

      // second row
      mm_dmvx = _mm_loadu_si128((const __m128i *)(vX + dMvStride));
      mm_dmvy = _mm_loadu_si128((const __m128i *)(vY + dMvStride));
      mm_gradx = _mm_cvtepi16_epi32(_mm_loadl_epi64((__m128i*)(gX + gradStride)));
      mm_grady = _mm_cvtepi16_epi32(_mm_loadl_epi64((__m128i*)(gY + gradStride)));
      mm_dI = _mm_add_epi32(_mm_mullo_epi32(mm_dmvx, mm_gradx), _mm_mullo_epi32(mm_dmvy, mm_grady));
      mm_dI = _mm_min_epi32(mm_dimax, _mm_max_epi32(mm_dimin, mm_dI));

      // combine both rows
      mm_dI = _mm_packs_epi32(mm_dI0, mm_dI);
      mm_dI = _mm_add_epi16(_mm_unpacklo_epi64(_mm_loadl_epi64((const __m128i *)src), _mm_loadl_epi64((const __m128i *)(src + srcStride))), mm_dI);
      if (!bi)
      {
        mm_dI = _mm_srai_epi16(_mm_adds_epi16(mm_dI, mm_offset), shiftNum);
        mm_dI = _mm_min_epi16(vibdimax, _mm_max_epi16(vibdimin, mm_dI));
      }

      _mm_storel_epi64((__m128i *)dst, mm_dI);
      _mm_storel_epi64((__m128i *)(dst + dstStride), _mm_unpackhi_epi64(mm_dI, mm_dI));
#endif
      vX += 4; vY += 4; gX += 4; gY += 4; src += 4; dst += 4;
    }

#if USE_AVX2
    dMvX += (dMvStride << 2);
    dMvY += (dMvStride << 2);
    gradX += (gradStride << 2);
    gradY += (gradStride << 2);
    srcPel += (srcStride << 2);
    dstPel += (dstStride << 2);
#else
    dMvX += (dMvStride << 1);
    dMvY += (dMvStride << 1);
    gradX += (gradStride << 1);
    gradY += (gradStride << 1);
    srcPel += (srcStride << 1);
    dstPel += (dstStride << 1);
#endif
  }
}


template< X86_VEXT vext >
void roundIntVector_SIMD(int* v, int size, unsigned int nShift, const int dmvLimit)
{
  CHECKD(size % 16 != 0, "Size must be multiple of 16!");
#ifdef USE_AVX512
  if (vext >= AVX512 && size >= 16)
  {
    __m512i dMvMin = _mm256_set1_epi32(-dmvLimit);
    __m512i dMvMax = _mm256_set1_epi32( dmvLimit );
    __m512i nOffset = _mm512_set1_epi32((1 << (nShift - 1)));
    __m512i vones = _mm512_set1_epi32(1);
    __m512i vzero = _mm512_setzero_si512();
    for (int i = 0; i < size; i += 16, v += 16)
    {
      __m512i src = _mm512_loadu_si512(v);
      __mmask16 mask = _mm512_cmpge_epi32_mask(src, vzero);
      src = __mm512_add_epi32(src, nOffset);
      __mm512i dst = _mm512_srai_epi32(_mm512_mask_sub_epi32(src, mask, src, vones), nShift);
      dst = _mm512_min_epi32(dMvMax, _mm512_max_epi32(dMvMin, dst));
      _mm512_storeu_si512(v, dst);
    }
  }
  else
#endif
#ifdef USE_AVX2
  if (vext >= AVX2 && size >= 8)
  {
    __m256i dMvMin = _mm256_set1_epi32(-dmvLimit);
    __m256i dMvMax = _mm256_set1_epi32( dmvLimit );
    __m256i nOffset = _mm256_set1_epi32(1 << (nShift - 1));
    __m256i vzero = _mm256_setzero_si256();
    for (int i = 0; i < size; i += 8, v += 8)
    {
      __m256i src = _mm256_lddqu_si256((__m256i*)v);
      __m256i of  = _mm256_cmpgt_epi32(src, vzero);
      __m256i dst = _mm256_srai_epi32(_mm256_add_epi32(_mm256_add_epi32(src, nOffset), of), nShift);
      dst = _mm256_min_epi32(dMvMax, _mm256_max_epi32(dMvMin, dst));
      _mm256_storeu_si256((__m256i*)v, dst);
    }
  }
  else
#endif
  {
    __m128i dMvMin = _mm_set1_epi32(-dmvLimit);
    __m128i dMvMax = _mm_set1_epi32( dmvLimit );
    __m128i nOffset = _mm_set1_epi32((1 << (nShift - 1)));
    __m128i vzero = _mm_setzero_si128();
    for (int i = 0; i < size; i += 4, v += 4)
    {
      __m128i src = _mm_loadu_si128((__m128i*)v);
      __m128i of  = _mm_cmpgt_epi32(src, vzero);
      __m128i dst = _mm_srai_epi32(_mm_add_epi32(_mm_add_epi32(src, nOffset), of), nShift);
      dst = _mm_min_epi32(dMvMax, _mm_max_epi32(dMvMin, dst));
      _mm_storeu_si128((__m128i*)v, dst);
    }
  }
}

template< X86_VEXT vext, bool PAD = true>
void gradFilter_SSE(Pel* src, int srcStride, int width, int height, int gradStride, Pel* gradX, Pel* gradY, const int bitDepth)
{
  Pel* srcTmp = src + srcStride + 1;
  Pel* gradXTmp = gradX + gradStride + 1;
  Pel* gradYTmp = gradY + gradStride + 1;

#if MULTI_PASS_DMVR || SAMPLE_BASED_BDOF
  int widthInside = width - 2;
  int heightInside = height - 2;
#else
  int widthInside = width - 2 * BIO_EXTEND_SIZE;
  int heightInside = height - 2 * BIO_EXTEND_SIZE;
#endif
  int shift1 = 6;
  __m128i mmShift1 = _mm_cvtsi32_si128( shift1 );
  assert((widthInside & 3) == 0);

#if MULTI_PASS_DMVR || SAMPLE_BASED_BDOF
  if ( widthInside > 4 )
#else
  if ( ( widthInside & 7 ) == 0 )
#endif
  {
    for (int y = 0; y < heightInside; y++)
    {
      int x = 0;
      for ( ; x < widthInside; x += 8 )
      {
        __m128i mmPixTop    = _mm_sra_epi16( _mm_loadu_si128( ( __m128i* ) ( srcTmp + x - srcStride ) ), mmShift1 );
        __m128i mmPixBottom = _mm_sra_epi16( _mm_loadu_si128( ( __m128i* ) ( srcTmp + x + srcStride ) ), mmShift1 );
        __m128i mmPixLeft   = _mm_sra_epi16( _mm_loadu_si128( ( __m128i* ) ( srcTmp + x - 1 ) ), mmShift1 );
        __m128i mmPixRight  = _mm_sra_epi16( _mm_loadu_si128( ( __m128i* ) ( srcTmp + x + 1 ) ), mmShift1 );

        __m128i mmGradVer = _mm_sub_epi16( mmPixBottom, mmPixTop );
        __m128i mmGradHor = _mm_sub_epi16( mmPixRight, mmPixLeft );

        _mm_storeu_si128( ( __m128i * ) ( gradYTmp + x ), mmGradVer );
        _mm_storeu_si128( ( __m128i * ) ( gradXTmp + x ), mmGradHor );
      }
      gradXTmp += gradStride;
      gradYTmp += gradStride;
      srcTmp += srcStride;
    }
  }
  else
  {
    __m128i mmPixTop = _mm_sra_epi16( _mm_unpacklo_epi64( _mm_loadl_epi64( (__m128i*) ( srcTmp - srcStride ) ), _mm_loadl_epi64( (__m128i*) ( srcTmp ) ) ), mmShift1 );
    for ( int y = 0; y < heightInside; y += 2 )
    {
      __m128i mmPixBottom = _mm_sra_epi16( _mm_unpacklo_epi64( _mm_loadl_epi64( (__m128i*) ( srcTmp + srcStride ) ), _mm_loadl_epi64( (__m128i*) ( srcTmp + ( srcStride << 1 ) ) ) ), mmShift1 );
      __m128i mmPixLeft   = _mm_sra_epi16( _mm_unpacklo_epi64( _mm_loadl_epi64( (__m128i*) ( srcTmp - 1 ) ), _mm_loadl_epi64( (__m128i*) ( srcTmp - 1 + srcStride ) ) ), mmShift1 );
      __m128i mmPixRight  = _mm_sra_epi16( _mm_unpacklo_epi64( _mm_loadl_epi64( (__m128i*) ( srcTmp + 1 ) ), _mm_loadl_epi64( (__m128i*) ( srcTmp + 1 + srcStride ) ) ), mmShift1 );

      __m128i mmGradVer = _mm_sub_epi16( mmPixBottom, mmPixTop );
      __m128i mmGradHor = _mm_sub_epi16( mmPixRight, mmPixLeft );

      _mm_storel_epi64( (__m128i *) gradYTmp, mmGradVer );
      _mm_storel_epi64( (__m128i *) ( gradYTmp + gradStride ), _mm_unpackhi_epi64( mmGradVer, mmGradHor ) );
      _mm_storel_epi64( (__m128i *) gradXTmp, mmGradHor );
      _mm_storel_epi64( (__m128i *) ( gradXTmp + gradStride ), _mm_unpackhi_epi64( mmGradHor, mmGradVer ) );

      mmPixTop = mmPixBottom;
      gradXTmp += gradStride << 1;
      gradYTmp += gradStride << 1;
      srcTmp   += srcStride << 1;
    }
  }

#if !MULTI_PASS_DMVR && !SAMPLE_BASED_BDOF
  if (PAD)
  {
  gradXTmp = gradX + gradStride + 1;
  gradYTmp = gradY + gradStride + 1;
  for (int y = 0; y < heightInside; y++)
  {
    gradXTmp[-1] = gradXTmp[0];
    gradXTmp[widthInside] = gradXTmp[widthInside - 1];
    gradXTmp += gradStride;

    gradYTmp[-1] = gradYTmp[0];
    gradYTmp[widthInside] = gradYTmp[widthInside - 1];
    gradYTmp += gradStride;
  }

  gradXTmp = gradX + gradStride;
  gradYTmp = gradY + gradStride;
  ::memcpy(gradXTmp - gradStride, gradXTmp, sizeof(Pel)*(width));
  ::memcpy(gradXTmp + heightInside*gradStride, gradXTmp + (heightInside - 1)*gradStride, sizeof(Pel)*(width));
  ::memcpy(gradYTmp - gradStride, gradYTmp, sizeof(Pel)*(width));
  ::memcpy(gradYTmp + heightInside*gradStride, gradYTmp + (heightInside - 1)*gradStride, sizeof(Pel)*(width));
  }
#endif
}


template< X86_VEXT vext >
void calcBlkGradient_SSE(int sx, int sy, int     *arraysGx2, int     *arraysGxGy, int     *arraysGxdI, int     *arraysGy2, int     *arraysGydI, int     &sGx2, int     &sGy2, int     &sGxGy, int     &sGxdI, int     &sGydI, int width, int height, int unitSize)
{
  int     *Gx2 = arraysGx2;
  int     *Gy2 = arraysGy2;
  int     *GxGy = arraysGxGy;
  int     *GxdI = arraysGxdI;
  int     *GydI = arraysGydI;

  // set to the above row due to JVET_K0485_BIO_EXTEND_SIZE
  Gx2 -= (BIO_EXTEND_SIZE*width);
  Gy2 -= (BIO_EXTEND_SIZE*width);
  GxGy -= (BIO_EXTEND_SIZE*width);
  GxdI -= (BIO_EXTEND_SIZE*width);
  GydI -= (BIO_EXTEND_SIZE*width);

  __m128i vzero = _mm_setzero_si128();
  __m128i mmGx2Total = _mm_setzero_si128();
  __m128i mmGy2Total = _mm_setzero_si128();
  __m128i mmGxGyTotal = _mm_setzero_si128();
  __m128i mmGxdITotal = _mm_setzero_si128();
  __m128i mmGydITotal = _mm_setzero_si128();

  for (int y = -BIO_EXTEND_SIZE; y < unitSize + BIO_EXTEND_SIZE; y++)
  {
    __m128i mmsGx2 = _mm_loadu_si128((__m128i*)(Gx2 - 1));   __m128i mmsGx2Sec = _mm_loadl_epi64((__m128i*)(Gx2 + 3));
    __m128i mmsGy2 = _mm_loadu_si128((__m128i*)(Gy2 - 1));   __m128i mmsGy2Sec = _mm_loadl_epi64((__m128i*)(Gy2 + 3));
    __m128i mmsGxGy = _mm_loadu_si128((__m128i*)(GxGy - 1));  __m128i mmsGxGySec = _mm_loadl_epi64((__m128i*)(GxGy + 3));
    __m128i mmsGxdI = _mm_loadu_si128((__m128i*)(GxdI - 1));  __m128i mmsGxdISec = _mm_loadl_epi64((__m128i*)(GxdI + 3));
    __m128i mmsGydI = _mm_loadu_si128((__m128i*)(GydI - 1));  __m128i mmsGydISec = _mm_loadl_epi64((__m128i*)(GydI + 3));

    mmsGx2 = _mm_add_epi32(mmsGx2, mmsGx2Sec);
    mmsGy2 = _mm_add_epi32(mmsGy2, mmsGy2Sec);
    mmsGxGy = _mm_add_epi32(mmsGxGy, mmsGxGySec);
    mmsGxdI = _mm_add_epi32(mmsGxdI, mmsGxdISec);
    mmsGydI = _mm_add_epi32(mmsGydI, mmsGydISec);


    mmGx2Total = _mm_add_epi32(mmGx2Total, mmsGx2);
    mmGy2Total = _mm_add_epi32(mmGy2Total, mmsGy2);
    mmGxGyTotal = _mm_add_epi32(mmGxGyTotal, mmsGxGy);
    mmGxdITotal = _mm_add_epi32(mmGxdITotal, mmsGxdI);
    mmGydITotal = _mm_add_epi32(mmGydITotal, mmsGydI);

    Gx2 += width;
    Gy2 += width;
    GxGy += width;
    GxdI += width;
    GydI += width;
  }

  mmGx2Total = _mm_hadd_epi32(_mm_hadd_epi32(mmGx2Total, vzero), vzero);
  mmGy2Total = _mm_hadd_epi32(_mm_hadd_epi32(mmGy2Total, vzero), vzero);
  mmGxGyTotal = _mm_hadd_epi32(_mm_hadd_epi32(mmGxGyTotal, vzero), vzero);
  mmGxdITotal = _mm_hadd_epi32(_mm_hadd_epi32(mmGxdITotal, vzero), vzero);
  mmGydITotal = _mm_hadd_epi32(_mm_hadd_epi32(mmGydITotal, vzero), vzero);

  sGx2 = _mm_cvtsi128_si32(mmGx2Total);
  sGy2 = _mm_cvtsi128_si32(mmGy2Total);
  sGxGy = _mm_cvtsi128_si32(mmGxGyTotal);
  sGxdI = _mm_cvtsi128_si32(mmGxdITotal);
  sGydI = _mm_cvtsi128_si32(mmGydITotal);
}

template< X86_VEXT vext, int W >
void reco_SSE( const int16_t* src0, int src0Stride, const int16_t* src1, int src1Stride, int16_t *dst, int dstStride, int width, int height, const ClpRng& clpRng )
{
  if( W == 8 )
  {
    if( vext >= AVX2 && ( width & 15 ) == 0 )
    {
#if USE_AVX2
      __m256i vbdmin = _mm256_set1_epi16( clpRng.min );
      __m256i vbdmax = _mm256_set1_epi16( clpRng.max );

      for( int row = 0; row < height; row++ )
      {
        for( int col = 0; col < width; col += 16 )
        {
          __m256i vdest = _mm256_lddqu_si256( ( const __m256i * )&src0[col] );
          __m256i vsrc1 = _mm256_lddqu_si256( ( const __m256i * )&src1[col] );

          vdest = _mm256_adds_epi16( vdest, vsrc1 );
          vdest = _mm256_min_epi16( vbdmax, _mm256_max_epi16( vbdmin, vdest ) );

          _mm256_storeu_si256( ( __m256i * )&dst[col], vdest );
        }

        src0 += src0Stride;
        src1 += src1Stride;
        dst  += dstStride;
      }
#endif
    }
    else
    {
      __m128i vbdmin = _mm_set1_epi16( clpRng.min );
      __m128i vbdmax = _mm_set1_epi16( clpRng.max );

      for( int row = 0; row < height; row++ )
      {
        for( int col = 0; col < width; col += 8 )
        {
          __m128i vdest = _mm_loadu_si128( ( const __m128i * )&src0[col] );
          __m128i vsrc1 = _mm_loadu_si128( ( const __m128i * )&src1[col] );

          vdest = _mm_adds_epi16( vdest, vsrc1 );
          vdest = _mm_min_epi16( vbdmax, _mm_max_epi16( vbdmin, vdest ) );

          _mm_storeu_si128( ( __m128i * )&dst[col], vdest );
        }

        src0 += src0Stride;
        src1 += src1Stride;
        dst  += dstStride;
      }
    }
  }
  else if( W == 4 )
  {
    __m128i vbdmin = _mm_set1_epi16( clpRng.min );
    __m128i vbdmax = _mm_set1_epi16( clpRng.max );

    for( int row = 0; row < height; row++ )
    {
      for( int col = 0; col < width; col += 4 )
      {
        __m128i vsrc = _mm_loadl_epi64( ( const __m128i * )&src0[col] );
        __m128i vdst = _mm_loadl_epi64( ( const __m128i * )&src1[col] );

        vdst = _mm_adds_epi16( vdst, vsrc );
        vdst = _mm_min_epi16( vbdmax, _mm_max_epi16( vbdmin, vdst ) );

        _mm_storel_epi64( ( __m128i * )&dst[col], vdst );
      }

      src0 += src0Stride;
      src1 += src1Stride;
      dst  +=  dstStride;
    }
  }
  else
  {
    THROW( "Unsupported size" );
  }
}

#if ENABLE_SIMD_OPT_BCW
template< X86_VEXT vext, int W >
void removeWeightHighFreq_SSE(int16_t* src0, int src0Stride, const int16_t* src1, int src1Stride, int width, int height, int shift, int bcwWeight)
{
  int normalizer = ((1 << 16) + (bcwWeight>0 ? (bcwWeight >> 1) : -(bcwWeight >> 1))) / bcwWeight;
  int weight0 = normalizer << g_BcwLog2WeightBase;
  int weight1 = (g_BcwWeightBase - bcwWeight)*normalizer;
  int offset = 1 << (shift - 1);
  if (W == 8)
  {
    __m128i vzero = _mm_setzero_si128();
    __m128i voffset = _mm_set1_epi32(offset);
    __m128i vw0 = _mm_set1_epi32(weight0);
    __m128i vw1 = _mm_set1_epi32(weight1);

    for (int row = 0; row < height; row++)
    {
      for (int col = 0; col < width; col += 8)
      {
        __m128i vsrc0 = _mm_loadu_si128( (const __m128i *)&src0[col] );
        __m128i vsrc1 = _mm_loadu_si128( (const __m128i *)&src1[col] );

        __m128i vtmp, vdst, vsrc;
        vdst = _mm_cvtepi16_epi32(vsrc0);
        vsrc = _mm_cvtepi16_epi32(vsrc1);
        vdst = _mm_mullo_epi32(vdst, vw0);
        vsrc = _mm_mullo_epi32(vsrc, vw1);
        vtmp = _mm_add_epi32(_mm_sub_epi32(vdst, vsrc), voffset);
        vtmp = _mm_srai_epi32(vtmp, shift);

        vsrc0 = _mm_unpackhi_epi64(vsrc0, vzero);
        vsrc1 = _mm_unpackhi_epi64(vsrc1, vzero);
        vdst = _mm_cvtepi16_epi32(vsrc0);
        vsrc = _mm_cvtepi16_epi32(vsrc1);
        vdst = _mm_mullo_epi32(vdst, vw0);
        vsrc = _mm_mullo_epi32(vsrc, vw1);
        vdst = _mm_add_epi32(_mm_sub_epi32(vdst, vsrc), voffset);
        vdst = _mm_srai_epi32(vdst, shift);
        vdst = _mm_packs_epi32(vtmp, vdst);

        _mm_store_si128((__m128i *)&src0[col], vdst);
      }
      src0 += src0Stride;
      src1 += src1Stride;
    }
  }
  else if (W == 4)
  {
    __m128i vzero = _mm_setzero_si128();
    __m128i voffset = _mm_set1_epi32(offset);
    __m128i vw0 = _mm_set1_epi32(weight0);
    __m128i vw1 = _mm_set1_epi32(weight1);

    for (int row = 0; row < height; row++)
    {
      __m128i vsum = _mm_loadl_epi64((const __m128i *)src0);
      __m128i vdst = _mm_loadl_epi64((const __m128i *)src1);

      vsum = _mm_cvtepi16_epi32(vsum);
      vdst = _mm_cvtepi16_epi32(vdst);
      vsum = _mm_mullo_epi32(vsum, vw0);
      vdst = _mm_mullo_epi32(vdst, vw1);
      vsum = _mm_add_epi32(_mm_sub_epi32(vsum, vdst), voffset);
      vsum = _mm_srai_epi32(vsum, shift);
      vsum = _mm_packs_epi32(vsum, vzero);

      _mm_storel_epi64((__m128i *)src0, vsum);

      src0 += src0Stride;
      src1 += src1Stride;
    }
  }
  else
  {
    THROW("Unsupported size");
  }
}

template< X86_VEXT vext, int W >
void removeHighFreq_SSE(int16_t* src0, int src0Stride, const int16_t* src1, int src1Stride, int width, int height)
{
  if (W == 8)
  {
    // TODO: AVX2 impl
    {
      for (int row = 0; row < height; row++)
      {
        for (int col = 0; col < width; col += 8)
        {
          __m128i vsrc0 = _mm_loadu_si128( (const __m128i *)&src0[col] );
          __m128i vsrc1 = _mm_loadu_si128( (const __m128i *)&src1[col] );

          vsrc0 = _mm_sub_epi16(_mm_slli_epi16(vsrc0, 1), vsrc1);
          _mm_store_si128((__m128i *)&src0[col], vsrc0);
        }

        src0 += src0Stride;
        src1 += src1Stride;
      }
    }
  }
  else if (W == 4)
  {
    for (int row = 0; row < height; row += 2)
    {
      __m128i vsrc0 = _mm_loadl_epi64((const __m128i *)src0);
      __m128i vsrc1 = _mm_loadl_epi64((const __m128i *)src1);
      __m128i vsrc0_2 = _mm_loadl_epi64((const __m128i *)(src0 + src0Stride));
      __m128i vsrc1_2 = _mm_loadl_epi64((const __m128i *)(src1 + src1Stride));

      vsrc0 = _mm_unpacklo_epi64(vsrc0, vsrc0_2);
      vsrc1 = _mm_unpacklo_epi64(vsrc1, vsrc1_2);

      vsrc0 = _mm_sub_epi16(_mm_slli_epi16(vsrc0, 1), vsrc1);
      _mm_storel_epi64((__m128i *)src0, vsrc0);
      _mm_storel_epi64((__m128i *)(src0 + src0Stride), _mm_unpackhi_epi64(vsrc0, vsrc0));

      src0 += (src0Stride << 1);
      src1 += (src1Stride << 1);
    }
  }
  else
  {
    THROW("Unsupported size");
  }
}
#endif

template<bool doShift, bool shiftR, typename T> static inline void do_shift( T &vreg, int num );
#if USE_AVX2
template<> inline void do_shift<true,  true , __m256i>( __m256i &vreg, int num ) { vreg = _mm256_srai_epi32( vreg, num ); }
template<> inline void do_shift<true,  false, __m256i>( __m256i &vreg, int num ) { vreg = _mm256_slli_epi32( vreg, num ); }
template<> inline void do_shift<false, true , __m256i>( __m256i &vreg, int num ) { }
template<> inline void do_shift<false, false, __m256i>( __m256i &vreg, int num ) { }
#endif
template<> inline void do_shift<true,  true , __m128i>( __m128i &vreg, int num ) { vreg = _mm_srai_epi32( vreg, num ); }
template<> inline void do_shift<true,  false, __m128i>( __m128i &vreg, int num ) { vreg = _mm_slli_epi32( vreg, num ); }
template<> inline void do_shift<false, true , __m128i>( __m128i &vreg, int num ) { }
template<> inline void do_shift<false, false, __m128i>( __m128i &vreg, int num ) { }

template<bool mult, typename T> static inline void do_mult( T& vreg, T& vmult );
template<> inline void do_mult<false, __m128i>( __m128i&, __m128i& ) { }
#if USE_AVX2
template<> inline void do_mult<false, __m256i>( __m256i&, __m256i& ) { }
#endif
template<> inline void do_mult<true,   __m128i>( __m128i& vreg, __m128i& vmult ) { vreg = _mm_mullo_epi32   ( vreg, vmult ); }
#if USE_AVX2
template<> inline void do_mult<true,   __m256i>( __m256i& vreg, __m256i& vmult ) { vreg = _mm256_mullo_epi32( vreg, vmult ); }
#endif

template<bool add, typename T> static inline void do_add( T& vreg, T& vadd );
template<> inline void do_add<false, __m128i>( __m128i&, __m128i& ) { }
#if USE_AVX2
template<> inline void do_add<false, __m256i>( __m256i&, __m256i& ) { }
#endif
template<> inline void do_add<true,  __m128i>( __m128i& vreg, __m128i& vadd ) { vreg = _mm_add_epi32( vreg, vadd ); }
#if USE_AVX2
template<> inline void do_add<true,  __m256i>( __m256i& vreg, __m256i& vadd ) { vreg = _mm256_add_epi32( vreg, vadd ); }
#endif

template<bool clip, typename T> static inline void do_clip( T& vreg, T& vbdmin, T& vbdmax );
template<> inline void do_clip<false, __m128i>( __m128i&, __m128i&, __m128i& ) { }
#if USE_AVX2
template<> inline void do_clip<false, __m256i>( __m256i&, __m256i&, __m256i& ) { }
#endif
template<> inline void do_clip<true,  __m128i>( __m128i& vreg, __m128i& vbdmin, __m128i& vbdmax ) { vreg = _mm_min_epi16   ( vbdmax, _mm_max_epi16   ( vbdmin, vreg ) ); }
#if USE_AVX2
template<> inline void do_clip<true,  __m256i>( __m256i& vreg, __m256i& vbdmin, __m256i& vbdmax ) { vreg = _mm256_min_epi16( vbdmax, _mm256_max_epi16( vbdmin, vreg ) ); }
#endif


template<X86_VEXT vext, int W, bool doAdd, bool mult, bool doShift, bool shiftR, bool clip>
void linTf_SSE( const Pel* src, int srcStride, Pel *dst, int dstStride, int width, int height, int scale, int shift, int offset, const ClpRng& clpRng )
{
  if( vext >= AVX2 && ( width & 7 ) == 0 && W == 8 )
  {
#if USE_AVX2
    __m256i vzero    = _mm256_setzero_si256();
    __m256i vbdmin   = _mm256_set1_epi16( clpRng.min );
    __m256i vbdmax   = _mm256_set1_epi16( clpRng.max );
    __m256i voffset  = _mm256_set1_epi32( offset );
    __m256i vscale   = _mm256_set1_epi32( scale );

    for( int row = 0; row < height; row++ )
    {
      for( int col = 0; col < width; col += 8 )
      {
        __m256i val;
        val = _mm256_cvtepi16_epi32       (  _mm_loadu_si128( ( const __m128i * )&src[col] ) );
        do_mult<mult, __m256i>            ( val, vscale );
        do_shift<doShift, shiftR, __m256i>( val, shift );
        do_add<doAdd, __m256i>            ( val, voffset );
        val = _mm256_packs_epi32          ( val, vzero );
        do_clip<clip, __m256i>            ( val, vbdmin, vbdmax );
        val = _mm256_permute4x64_epi64    ( val, ( 0 << 0 ) + ( 2 << 2 ) + ( 1 << 4 ) + ( 1 << 6 ) );

        _mm_storeu_si128                  ( ( __m128i * )&dst[col], _mm256_castsi256_si128( val ) );
      }

      src += srcStride;
      dst += dstStride;
    }
#endif
  }
  else
  {
    __m128i vzero   = _mm_setzero_si128();
    __m128i vbdmin  = _mm_set1_epi16   ( clpRng.min );
    __m128i vbdmax  = _mm_set1_epi16   ( clpRng.max );
    __m128i voffset = _mm_set1_epi32   ( offset );
    __m128i vscale  = _mm_set1_epi32   ( scale );

    for( int row = 0; row < height; row++ )
    {
      for( int col = 0; col < width; col += 4 )
      {
        __m128i val;
        val = _mm_loadl_epi64             ( ( const __m128i * )&src[col] );
        val = _mm_cvtepi16_epi32          ( val );
        do_mult<mult, __m128i>            ( val, vscale );
        do_shift<doShift, shiftR, __m128i>( val, shift );
        do_add<doAdd, __m128i>            ( val, voffset );
        val = _mm_packs_epi32             ( val, vzero );
        do_clip<clip, __m128i>            ( val, vbdmin, vbdmax );

        _mm_storel_epi64                  ( ( __m128i * )&dst[col], val );
      }

      src += srcStride;
      dst += dstStride;
    }
  }
}

template<X86_VEXT vext, int W>
void linTf_SSE_entry( const Pel* src, int srcStride, Pel *dst, int dstStride, int width, int height, int scale, int shift, int offset, const ClpRng& clpRng, bool clip )
{
  int fn = ( offset == 0 ? 16 : 0 ) + ( scale == 1 ? 8 : 0 ) + ( shift == 0 ? 4 : 0 ) + ( shift < 0 ? 2 : 0 ) + ( !clip ? 1 : 0 );

  switch( fn )
  {
  case  0: linTf_SSE<vext, W, true,  true,  true,  true,  true >( src, srcStride, dst, dstStride, width, height, scale,  shift, offset, clpRng ); break;
  case  1: linTf_SSE<vext, W, true,  true,  true,  true,  false>( src, srcStride, dst, dstStride, width, height, scale,  shift, offset, clpRng ); break;
  case  2: linTf_SSE<vext, W, true,  true,  true,  false, true >( src, srcStride, dst, dstStride, width, height, scale, -shift, offset, clpRng ); break;
  case  3: linTf_SSE<vext, W, true,  true,  true,  false, false>( src, srcStride, dst, dstStride, width, height, scale, -shift, offset, clpRng ); break;
  case  4: linTf_SSE<vext, W, true,  true,  false, true,  true >( src, srcStride, dst, dstStride, width, height, scale,  shift, offset, clpRng ); break;
  case  5: linTf_SSE<vext, W, true,  true,  false, true,  false>( src, srcStride, dst, dstStride, width, height, scale,  shift, offset, clpRng ); break;
  case  6: linTf_SSE<vext, W, true,  true,  false, false, true >( src, srcStride, dst, dstStride, width, height, scale, -shift, offset, clpRng ); break;
  case  7: linTf_SSE<vext, W, true,  true,  false, false, false>( src, srcStride, dst, dstStride, width, height, scale, -shift, offset, clpRng ); break;
  case  8: linTf_SSE<vext, W, true,  false, true,  true,  true >( src, srcStride, dst, dstStride, width, height, scale,  shift, offset, clpRng ); break;
  case  9: linTf_SSE<vext, W, true,  false, true,  true,  false>( src, srcStride, dst, dstStride, width, height, scale,  shift, offset, clpRng ); break;
  case 10: linTf_SSE<vext, W, true,  false, true,  false, true >( src, srcStride, dst, dstStride, width, height, scale, -shift, offset, clpRng ); break;
  case 11: linTf_SSE<vext, W, true,  false, true,  false, false>( src, srcStride, dst, dstStride, width, height, scale, -shift, offset, clpRng ); break;
  case 12: linTf_SSE<vext, W, true,  false, false, true,  true >( src, srcStride, dst, dstStride, width, height, scale,  shift, offset, clpRng ); break;
  case 13: linTf_SSE<vext, W, true,  false, false, true,  false>( src, srcStride, dst, dstStride, width, height, scale,  shift, offset, clpRng ); break;
  case 14: linTf_SSE<vext, W, true,  false, false, false, true >( src, srcStride, dst, dstStride, width, height, scale, -shift, offset, clpRng ); break;
  case 15: linTf_SSE<vext, W, true,  false, false, false, false>( src, srcStride, dst, dstStride, width, height, scale, -shift, offset, clpRng ); break;
  case 16: linTf_SSE<vext, W, false, true,  true,  true,  true >( src, srcStride, dst, dstStride, width, height, scale,  shift, offset, clpRng ); break;
  case 17: linTf_SSE<vext, W, false, true,  true,  true,  false>( src, srcStride, dst, dstStride, width, height, scale,  shift, offset, clpRng ); break;
  case 18: linTf_SSE<vext, W, false, true,  true,  false, true >( src, srcStride, dst, dstStride, width, height, scale, -shift, offset, clpRng ); break;
  case 19: linTf_SSE<vext, W, false, true,  true,  false, false>( src, srcStride, dst, dstStride, width, height, scale, -shift, offset, clpRng ); break;
  case 20: linTf_SSE<vext, W, false, true,  false, true,  true >( src, srcStride, dst, dstStride, width, height, scale,  shift, offset, clpRng ); break;
  case 21: linTf_SSE<vext, W, false, true,  false, true,  false>( src, srcStride, dst, dstStride, width, height, scale,  shift, offset, clpRng ); break;
  case 22: linTf_SSE<vext, W, false, true,  false, false, true >( src, srcStride, dst, dstStride, width, height, scale, -shift, offset, clpRng ); break;
  case 23: linTf_SSE<vext, W, false, true,  false, false, false>( src, srcStride, dst, dstStride, width, height, scale, -shift, offset, clpRng ); break;
  case 24: linTf_SSE<vext, W, false, false, true,  true,  true >( src, srcStride, dst, dstStride, width, height, scale,  shift, offset, clpRng ); break;
  case 25: linTf_SSE<vext, W, false, false, true,  true,  false>( src, srcStride, dst, dstStride, width, height, scale,  shift, offset, clpRng ); break;
  case 26: linTf_SSE<vext, W, false, false, true,  false, true >( src, srcStride, dst, dstStride, width, height, scale, -shift, offset, clpRng ); break;
  case 27: linTf_SSE<vext, W, false, false, true,  false, false>( src, srcStride, dst, dstStride, width, height, scale, -shift, offset, clpRng ); break;
  case 28: linTf_SSE<vext, W, false, false, false, true,  true >( src, srcStride, dst, dstStride, width, height, scale,  shift, offset, clpRng ); break;
  case 29: linTf_SSE<vext, W, false, false, false, true,  false>( src, srcStride, dst, dstStride, width, height, scale,  shift, offset, clpRng ); break;
  case 30: linTf_SSE<vext, W, false, false, false, false, true >( src, srcStride, dst, dstStride, width, height, scale, -shift, offset, clpRng ); break;
  case 31: linTf_SSE<vext, W, false, false, false, false, false>( src, srcStride, dst, dstStride, width, height, scale, -shift, offset, clpRng ); break;
  default:
    THROW( "Unknown parametrization of the linear transformation" );
    break;
  }
}

#if TM_AMVP || TM_MRG
template<X86_VEXT vext>
int64_t getSumOfDifference_SSE(const Pel* src0, int src0Stride, const Pel* src1, int src1Stride, int width, int height, int rowSubShift, int bitDepth)
{
  if (width < 4 || (bitDepth - rowSubShift) > 10) // note: it seems 12b also works
  {
    return getSumOfDifferenceCore(src0, src0Stride, src1, src1Stride, width, height, rowSubShift, bitDepth);
  }

  const short* pOrg = src0;
  const short* pCur = src1;

  const int subShift  = rowSubShift;
  const int subStep   = 1 << subShift;
  const int strideOrg = src0Stride * subStep;
  const int strideCur = src1Stride * subStep;
  int deltaAvg = 0;

  // internal bit-depth must be 12-bit or lower

  if (width & 7) // multiple of 4
  {
    __m128i vzero = _mm_setzero_si128();
    __m128i vsum32 = vzero;

    for (; height != 0; height -= subStep)
    {
      __m128i vsum16 = vzero;

      for (int n = 0; n < width; n += 4)
      {
        __m128i org = _mm_loadl_epi64((__m128i*)(pOrg + n));
        __m128i cur = _mm_loadl_epi64((__m128i*)(pCur + n));
        vsum16 = _mm_adds_epi16(vsum16, _mm_sub_epi16(org, cur));
      }

      __m128i vsign = _mm_cmpgt_epi16(vzero, vsum16);
      vsum32 = _mm_add_epi32(vsum32, _mm_unpacklo_epi16(vsum16, vsign));

      pOrg += strideOrg;
      pCur += strideCur;
    }

    vsum32 = _mm_add_epi32(vsum32, _mm_shuffle_epi32(vsum32, 0x4e));   // 01001110
    vsum32 = _mm_add_epi32(vsum32, _mm_shuffle_epi32(vsum32, 0xb1));   // 10110001
    deltaAvg = _mm_cvtsi128_si32(vsum32);
  }
#ifdef USE_AVX2
  else if (vext >= AVX2 && width >= 16) // multiple of 16
  {
    __m256i vzero = _mm256_setzero_si256();
    __m256i vsum32 = vzero;

    for (; height != 0; height -= subStep)
    {
      __m256i vsum16 = vzero;

      for (int n = 0; n < width; n += 16)
      {
        __m256i org = _mm256_lddqu_si256((__m256i*)(pOrg + n));
        __m256i cur = _mm256_lddqu_si256((__m256i*)(pCur + n));
        vsum16 = _mm256_adds_epi16(vsum16, _mm256_sub_epi16(org, cur));
      }

      __m256i vsign = _mm256_cmpgt_epi16(vzero, vsum16);
      __m256i vsumtemp = _mm256_add_epi32(_mm256_unpacklo_epi16(vsum16, vsign), _mm256_unpackhi_epi16(vsum16, vsign));
      vsum32 = _mm256_add_epi32(vsum32, vsumtemp);

      pOrg += strideOrg;
      pCur += strideCur;
    }

    vsum32 = _mm256_hadd_epi32(vsum32, vzero);
    vsum32 = _mm256_hadd_epi32(vsum32, vzero);
    deltaAvg = _mm_cvtsi128_si32(_mm256_castsi256_si128(vsum32)) + _mm_cvtsi128_si32(_mm256_castsi256_si128(_mm256_permute2x128_si256(vsum32, vsum32, 0x11)));
  }
#endif
  else // multiple of 8
  {
    __m128i vzero = _mm_setzero_si128();
    __m128i vsum32 = vzero;

    for (; height != 0; height -= subStep)
    {
      __m128i vsum16 = vzero;

      for (int n = 0; n < width; n += 8)
      {
        __m128i org = _mm_lddqu_si128((__m128i*)(pOrg + n));
        __m128i cur = _mm_lddqu_si128((__m128i*)(pCur + n));
        vsum16 = _mm_adds_epi16(vsum16, _mm_sub_epi16(org, cur));
      }

      __m128i vsign = _mm_cmpgt_epi16(vzero, vsum16);
      __m128i vsumtemp = _mm_add_epi32(_mm_unpacklo_epi16(vsum16, vsign), _mm_unpackhi_epi16(vsum16, vsign));
      vsum32 = _mm_add_epi32(vsum32, vsumtemp);

      pOrg += strideOrg;
      pCur += strideCur;
    }

    vsum32 = _mm_add_epi32(vsum32, _mm_shuffle_epi32(vsum32, 0x4e));   // 01001110
    vsum32 = _mm_add_epi32(vsum32, _mm_shuffle_epi32(vsum32, 0xb1));   // 10110001
    deltaAvg = _mm_cvtsi128_si32(vsum32);
  }

  deltaAvg <<= subShift;
  return deltaAvg;
}
#endif

template<X86_VEXT vext>
void PelBufferOps::_initPelBufOpsX86()
{
#if JVET_W0097_GPM_MMVD_TM
  roundBD = roundBD_SSE<vext>;
  weightedAvg = weightedAvg_SSE<vext>;
  copyClip = copyClip_SSE<vext>;
#endif
  addAvg8 = addAvg_SSE<vext, 8>;
  addAvg4 = addAvg_SSE<vext, 4>;

  addBIOAvg4      = addBIOAvg4_SSE<vext>;
#if MULTI_PASS_DMVR || SAMPLE_BASED_BDOF
  calcBIOParameter   = calcBIOParameter_SSE<vext>;
  calcBIOParamSum5   = calcBIOParamSum5_SSE<vext>;
  calcBIOParamSum4   = calcBIOParamSum4_SSE<vext>;
  calcBIOClippedVxVy = calcBIOClippedVxVy_SSE<vext>;
  addBIOAvgN         = addBIOAvgN_SSE<vext>;
  calAbsSum          = calAbsSum_SSE<vext>;
  bioGradFilter      = gradFilter_SSE<vext, false>;
#else
  bioGradFilter   = gradFilter_SSE<vext>;
#endif
  calcBIOSums = calcBIOSums_SSE<vext>;

  copyBuffer = copyBufferSimd<vext>;
#if !MULTI_PASS_DMVR
  padding    = paddingSimd<vext>;
#endif
  reco8 = reco_SSE<vext, 8>;
  reco4 = reco_SSE<vext, 4>;

  linTf8 = linTf_SSE_entry<vext, 8>;
  linTf4 = linTf_SSE_entry<vext, 4>;
#if ENABLE_SIMD_OPT_BCW
  removeWeightHighFreq8 = removeWeightHighFreq_SSE<vext, 8>;
  removeWeightHighFreq4 = removeWeightHighFreq_SSE<vext, 4>;
  removeHighFreq8 = removeHighFreq_SSE<vext, 8>;
  removeHighFreq4 = removeHighFreq_SSE<vext, 4>;
#endif
  profGradFilter = gradFilter_SSE<vext, false>;
  applyPROF      = applyPROF_SSE<vext>;
  roundIntVector = roundIntVector_SIMD<vext>;
#if TM_AMVP || TM_MRG
  getSumOfDifference = getSumOfDifference_SSE<vext>;
#endif
}

template void PelBufferOps::_initPelBufOpsX86<SIMDX86>();

#endif // TARGET_SIMD_X86
#endif
//! \}
