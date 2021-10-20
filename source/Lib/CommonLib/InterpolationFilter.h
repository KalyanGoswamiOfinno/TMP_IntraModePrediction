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
 * \file
 * \brief Declaration of InterpolationFilter class
 */

#ifndef __INTERPOLATIONFILTER__
#define __INTERPOLATIONFILTER__

#include "CommonDef.h"
#include "CacheModel.h"

//! \ingroup CommonLib
//! \{

#define IF_INTERNAL_PREC 14 ///< Number of bits for internal precision
#if IF_12TAP
#define IF_FILTER_PREC    8 ///< Log2 of sum of filter taps
#else
#define IF_FILTER_PREC    6 ///< Log2 of sum of filter taps
#endif
#define IF_INTERNAL_OFFS (1<<(IF_INTERNAL_PREC-1)) ///< Offset used internally
#define IF_INTERNAL_PREC_BILINEAR 10 ///< Number of bits for internal precision
#define IF_FILTER_PREC_BILINEAR   4  ///< Bilinear filter coeff precision so that intermediate value will not exceed 16 bit for SIMD - bit exact
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
#define IF_INTERNAL_FRAC_BITS(bd) std::max(2, IF_INTERNAL_PREC - int(bd))
#endif
/**
 * \brief Interpolation filter class
 */
class InterpolationFilter
{
#if IF_12TAP
  static const TFilterCoeff m_lumaFilter4x4[LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][8];
#else
  static const TFilterCoeff m_lumaFilter4x4[LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][NTAPS_LUMA];
#endif
public:
#if IF_12TAP
  static const TFilterCoeff m_lumaFilter12[LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS + 1][12];     ///< Luma filter taps //+1 added by //kolya
  static const TFilterCoeff m_lumaFilter[LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][8]; ///< Luma filter taps, for affine 
  static const TFilterCoeff m_chromaFilter[CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][NTAPS_CHROMA]; ///< Chroma filter taps
#if INTRA_6TAP
  static const TFilterCoeff m_lumaIntraFilter[CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][6]; ///< Chroma filter 6 taps
  static const TFilterCoeff m_weak4TapFilter[CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][NTAPS_CHROMA]; ///< Weak filter 4 taps
#if JVET_W0123_TIMD_FUSION
  static const TFilterCoeff m_lumaIntraFilterExt[CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS << 1][6]; ///< Chroma filter 6 taps
#endif
#endif
  static const TFilterCoeff m_lumaFilterRPR1[LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][8]; ///< Luma filter taps 1.5x
  static const TFilterCoeff m_lumaFilterRPR2[LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][8]; ///< Luma filter taps 2x
  static const TFilterCoeff m_chromaFilterRPR1[CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][NTAPS_CHROMA]; ///< Chroma filter taps 1.5x
  static const TFilterCoeff m_chromaFilterRPR2[CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][NTAPS_CHROMA]; ///< Chroma filter taps 2x

  static const TFilterCoeff m_affineLumaFilterRPR1[LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][8]; ///< Luma filter taps 1.5x
  static const TFilterCoeff m_affineLumaFilterRPR2[LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][8]; ///< Luma filter taps 2x

private:
  static const TFilterCoeff m_lumaAltHpelIFilter[8]; ///< Luma filter taps
#else
  static const TFilterCoeff m_lumaFilter[LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][NTAPS_LUMA]; ///< Luma filter taps
  static const TFilterCoeff m_chromaFilter[CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][NTAPS_CHROMA]; ///< Chroma filter taps
#if INTRA_6TAP
  static const TFilterCoeff m_lumaIntraFilter[CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][6]; ///< Chroma filter 6 taps
  static const TFilterCoeff m_weak4TapFilter[CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][NTAPS_CHROMA]; ///< Weak filter 4 taps
#endif
  static const TFilterCoeff m_lumaFilterRPR1[LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][NTAPS_LUMA]; ///< Luma filter taps 1.5x
  static const TFilterCoeff m_lumaFilterRPR2[LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][NTAPS_LUMA]; ///< Luma filter taps 2x
  static const TFilterCoeff m_chromaFilterRPR1[CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][NTAPS_CHROMA]; ///< Chroma filter taps 1.5x
  static const TFilterCoeff m_chromaFilterRPR2[CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][NTAPS_CHROMA]; ///< Chroma filter taps 2x
  static const TFilterCoeff m_affineLumaFilterRPR1[LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][NTAPS_LUMA]; ///< Luma filter taps 1.5x
  static const TFilterCoeff m_affineLumaFilterRPR2[LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][NTAPS_LUMA]; ///< Luma filter taps 2x
private:
  static const TFilterCoeff m_lumaAltHpelIFilter[NTAPS_LUMA]; ///< Luma filter taps
#endif
#if JVET_W0123_TIMD_FUSION
  static const TFilterCoeff g_aiExtIntraCubicFilter[CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS<<1][NTAPS_CHROMA]; ///< Chroma filter taps
  static const TFilterCoeff g_aiExtIntraGaussFilter[CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS<<1][NTAPS_CHROMA]; ///< Chroma filter taps
#endif
  static const TFilterCoeff m_bilinearFilter[LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][NTAPS_BILINEAR]; ///< bilinear filter taps
  static const TFilterCoeff m_bilinearFilterPrec4[LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS][NTAPS_BILINEAR]; ///< bilinear filter taps
public:
  template<bool isFirst, bool isLast>
  static void filterCopy( const ClpRng& clpRng, const Pel *src, int srcStride, Pel *dst, int dstStride, int width, int height, bool biMCForDMVR);

  template<int N, bool isVertical, bool isFirst, bool isLast>
  static void filter(const ClpRng& clpRng, Pel const *src, int srcStride, Pel *dst, int dstStride, int width, int height, TFilterCoeff const *coeff, bool biMCForDMVR);

  void filter4x4( const ClpRng& clpRng, Pel const *src, int srcStride, Pel *dst, int dstStride, int xFrac, int yFrac, bool isLast);

  template<int N>
  void filterHor(const ClpRng& clpRng, Pel const* src, int srcStride, Pel *dst, int dstStride, int width, int height, bool isLast, TFilterCoeff const *coeff, bool biMCForDMVR);

  template<int N>
  void filterVer(const ClpRng& clpRng, Pel const* src, int srcStride, Pel *dst, int dstStride, int width, int height, bool isFirst, bool isLast, TFilterCoeff const *coeff, bool biMCForDMVR);

  static void xWeightedGeoBlk(const PredictionUnit &pu, const uint32_t width, const uint32_t height, const ComponentID compIdx, const uint8_t splitDir, PelUnitBuf& predDst, PelUnitBuf& predSrc0, PelUnitBuf& predSrc1);
  void weightedGeoBlk(const PredictionUnit &pu, const uint32_t width, const uint32_t height, const ComponentID compIdx, const uint8_t splitDir, PelUnitBuf& predDst, PelUnitBuf& predSrc0, PelUnitBuf& predSrc1);
protected:
#if JVET_J0090_MEMORY_BANDWITH_MEASURE
  static CacheModel* m_cacheModel;
#endif
public:
  InterpolationFilter();
  ~InterpolationFilter() {}
#if IF_12TAP
#if SIMD_4x4_12
  void(*m_filter4x4[2])(const int16_t* src, int srcStride, int16_t *dst, int dstStride, int shiftH, int offsetH, int shiftV, int offsetV, int16_t const *coeffH, int16_t const *coeffV, int ibdimin, int ibdimax); //kolya
#endif

  void(*m_filterHor[4][2][2])(const ClpRng& clpRng, Pel const *src, int srcStride, Pel *dst, int dstStride, int width, int height, TFilterCoeff const *coeff, bool biMCForDMVR);
  void(*m_filterVer[4][2][2])(const ClpRng& clpRng, Pel const *src, int srcStride, Pel *dst, int dstStride, int width, int height, TFilterCoeff const *coeff, bool biMCForDMVR);
#else
  void(*m_filterHor[3][2][2])(const ClpRng& clpRng, Pel const *src, int srcStride, Pel *dst, int dstStride, int width, int height, TFilterCoeff const *coeff, bool biMCForDMVR);
  void(*m_filterVer[3][2][2])(const ClpRng& clpRng, Pel const *src, int srcStride, Pel *dst, int dstStride, int width, int height, TFilterCoeff const *coeff, bool biMCForDMVR);
#endif
  void( *m_filterCopy[2][2] )  ( const ClpRng& clpRng, Pel const *src, int srcStride, Pel *dst, int dstStride, int width, int height, bool biMCForDMVR);
  void( *m_weightedGeoBlk )(const PredictionUnit &pu, const uint32_t width, const uint32_t height, const ComponentID compIdx, const uint8_t splitDir, PelUnitBuf& predDst, PelUnitBuf& predSrc0, PelUnitBuf& predSrc1);

  void initInterpolationFilter( bool enable );
#ifdef TARGET_SIMD_X86
  void initInterpolationFilterX86();
  template <X86_VEXT vext>
  void _initInterpolationFilterX86();
#endif
  void filterHor(const ComponentID compID, Pel const* src, int srcStride, Pel *dst, int dstStride, int width, int height, int frac,               bool isLast, const ChromaFormat fmt, const ClpRng& clpRng, int nFilterIdx = 0, bool biMCForDMVR = false, bool useAltHpelIf = false);
  void filterVer(const ComponentID compID, Pel const* src, int srcStride, Pel *dst, int dstStride, int width, int height, int frac, bool isFirst, bool isLast, const ChromaFormat fmt, const ClpRng& clpRng, int nFilterIdx = 0, bool biMCForDMVR = false, bool useAltHpelIf = false);
#if JVET_J0090_MEMORY_BANDWITH_MEASURE
  void cacheAssign( CacheModel *cache ) { m_cacheModel = cache; }
#endif

#if INTRA_6TAP
  static TFilterCoeff const * const getIntraLumaFilterTable(const int deltaFract) { return m_lumaIntraFilter[deltaFract]; };
  static TFilterCoeff const * const getWeak4TapFilterTable(const int deltaFract) { return m_weak4TapFilter[deltaFract]; };
#if JVET_W0123_TIMD_FUSION
  static TFilterCoeff const * const getIntraLumaFilterTableExt(const int deltaFract) { return m_lumaIntraFilterExt[deltaFract]; };
#endif
#endif
  static TFilterCoeff const * const getChromaFilterTable(const int deltaFract) { return m_chromaFilter[deltaFract]; };
#if JVET_W0123_TIMD_FUSION
  static TFilterCoeff const * const getExtIntraCubicFilter(const int deltaFract) { return g_aiExtIntraCubicFilter[deltaFract]; };
  static TFilterCoeff const * const getExtIntraGaussFilter(const int deltaFract) { return g_aiExtIntraGaussFilter[deltaFract]; };
#endif
};

//! \}

#endif
