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

#pragma once

#include "CommonDefX86.h"
#include "../BilateralFilter.h"

#ifdef TARGET_SIMD_X86
#if defined _MSC_VER
#include <tmmintrin.h>
#else
#include <x86intrin.h>
#endif

#if ENABLE_SIMD_BILATERAL_FILTER
#if JVET_W0066_CCSAO
template<X86_VEXT vext>
void BilateralFilter::simdFilterDiamond5x5NoClip(uint32_t uiWidth, uint32_t uiHeight, int16_t block[], int16_t blkFilt[], const ClpRng& clpRng, Pel* recPtr, int recStride, int iWidthExtSIMD, int bfac, int bif_round_add, int bif_round_shift, bool isRDO, const char* LUTrowPtr)
{
  //if( uiWidth < 4 || ( uiWidth < 8 && isRDO ) )
  if (uiWidth < 4)
  {
    return blockBilateralFilterDiamond5x5NoClip(uiWidth, uiHeight, block, blkFilt, clpRng, recPtr, recStride, iWidthExtSIMD, bfac, bif_round_add, bif_round_shift, isRDO, LUTrowPtr);
  }

  int pad = 2;
  int padwidth = iWidthExtSIMD;

  __m128i center, left, right, up, down, lu, ld, ru, rd, diffabs, four, fifteen, lut, acc, temp, round_add, clipmin, clipmax, inputVals;
  __m128i ll, rr, uu, dd;

  four = _mm_set1_epi16(4);
  fifteen = _mm_set1_epi16(15);
  round_add = _mm_set1_epi16(bif_round_add);
  clipmin = _mm_set1_epi16(clpRng.min);
  clipmax = _mm_set1_epi16(clpRng.max);

  lut = _mm_loadu_si128((__m128i*)(LUTrowPtr));
  acc = _mm_set1_epi32(0);

  // Copy back parameters
  Pel* tempBlockPtr = (short*)blkFilt + (((padwidth + 4) << 1) + 2);
  int tempBlockStride = padwidth + 4;

  for (int col = 0; col < uiWidth; col += 8)
  {
    for (int row = 0; row < uiHeight; row++)
    {
      acc = _mm_set1_epi32(0);
      int16_t* point = &block[(row + pad) * padwidth + pad + col];

      center = _mm_loadu_si128((__m128i*)(point));

      //load neighbours
      left = _mm_loadu_si128((__m128i*)(point - 1));
      right = _mm_loadu_si128((__m128i*)(point + 1));
      up = _mm_loadu_si128((__m128i*)(point - padwidth));
      down = _mm_loadu_si128((__m128i*)(point + padwidth));

      lu = _mm_loadu_si128((__m128i*)(point - 1 - padwidth));
      ld = _mm_loadu_si128((__m128i*)(point - 1 + padwidth));
      ru = _mm_loadu_si128((__m128i*)(point + 1 - padwidth));
      rd = _mm_loadu_si128((__m128i*)(point + 1 + padwidth));

      ll = _mm_loadu_si128((__m128i*)(point - 2));
      rr = _mm_loadu_si128((__m128i*)(point + 2));
      uu = _mm_loadu_si128((__m128i*)(point - 2 * padwidth));
      dd = _mm_loadu_si128((__m128i*)(point + 2 * padwidth));

      //calculate diffs
      left = _mm_sub_epi16(left, center);
      right = _mm_sub_epi16(right, center);
      up = _mm_sub_epi16(up, center);
      down = _mm_sub_epi16(down, center);

      lu = _mm_sub_epi16(lu, center);
      ld = _mm_sub_epi16(ld, center);
      ru = _mm_sub_epi16(ru, center);
      rd = _mm_sub_epi16(rd, center);

      ll = _mm_sub_epi16(ll, center);
      rr = _mm_sub_epi16(rr, center);
      uu = _mm_sub_epi16(uu, center);
      dd = _mm_sub_epi16(dd, center);

      //LEFT!
      //calculate abs
      diffabs = _mm_abs_epi16(left); //abs
      diffabs = _mm_add_epi16(diffabs, four); //+4
      diffabs = _mm_srai_epi16(diffabs, 3); //>>3
      diffabs = _mm_min_epi16(diffabs, fifteen); //min(x,15)
      diffabs = _mm_packus_epi16(diffabs, diffabs); //convert to 8
      diffabs = _mm_shuffle_epi8(lut, diffabs);//lut
      diffabs = _mm_cvtepi8_epi16(diffabs);//back to 16-bit
      diffabs = _mm_sign_epi16(diffabs, left);//fix sign!
      acc = _mm_add_epi16(diffabs, acc); //add to acc
      //RIGHT!
      //calculate abs
      diffabs = _mm_abs_epi16(right); //abs
      diffabs = _mm_add_epi16(diffabs, four); //+4
      diffabs = _mm_srai_epi16(diffabs, 3); //>>3
      diffabs = _mm_min_epi16(diffabs, fifteen); //min(x,15)
      diffabs = _mm_packus_epi16(diffabs, diffabs); //convert to 8
      diffabs = _mm_shuffle_epi8(lut, diffabs);//lut
      diffabs = _mm_cvtepi8_epi16(diffabs);//back to 16-bit
      diffabs = _mm_sign_epi16(diffabs, right);//fix sign!
      acc = _mm_add_epi16(diffabs, acc); //add to acc
      //UP!
      //calculate abs
      diffabs = _mm_abs_epi16(up); //abs
      diffabs = _mm_add_epi16(diffabs, four); //+4
      diffabs = _mm_srai_epi16(diffabs, 3); //>>3
      diffabs = _mm_min_epi16(diffabs, fifteen); //min(x,15)
      diffabs = _mm_packus_epi16(diffabs, diffabs); //convert to 8
      diffabs = _mm_shuffle_epi8(lut, diffabs);//lut
      diffabs = _mm_cvtepi8_epi16(diffabs);//back to 16-bit
      diffabs = _mm_sign_epi16(diffabs, up);//fix sign!
      acc = _mm_add_epi16(diffabs, acc); //add to acc

      //DOWN!
      //calculate abs
      diffabs = _mm_abs_epi16(down); //abs
      diffabs = _mm_add_epi16(diffabs, four); //+4
      diffabs = _mm_srai_epi16(diffabs, 3); //>>3
      diffabs = _mm_min_epi16(diffabs, fifteen); //min(x,15)
      diffabs = _mm_packus_epi16(diffabs, diffabs); //convert to 8
      diffabs = _mm_shuffle_epi8(lut, diffabs);//lut
      diffabs = _mm_cvtepi8_epi16(diffabs);//back to 16-bit
      diffabs = _mm_sign_epi16(diffabs, down);//fix sign!
      acc = _mm_add_epi16(diffabs, acc); //add to acc

      //lu!
      //calculate abs
      diffabs = _mm_abs_epi16(lu); //abs
      diffabs = _mm_add_epi16(diffabs, four); //+4
      diffabs = _mm_srai_epi16(diffabs, 3); //>>3
      diffabs = _mm_min_epi16(diffabs, fifteen); //min(x,15)
      diffabs = _mm_packus_epi16(diffabs, diffabs); //convert to 8
      diffabs = _mm_shuffle_epi8(lut, diffabs);//lut
      diffabs = _mm_cvtepi8_epi16(diffabs);//back to 16-bit
      diffabs = _mm_srai_epi16(diffabs, 1);//diagonal shift!
      diffabs = _mm_sign_epi16(diffabs, lu);//fix sign!
      acc = _mm_add_epi16(diffabs, acc); //add to acc
      //ld!
      //calculate abs
      diffabs = _mm_abs_epi16(ld); //abs
      diffabs = _mm_add_epi16(diffabs, four); //+4
      diffabs = _mm_srai_epi16(diffabs, 3); //>>3
      diffabs = _mm_min_epi16(diffabs, fifteen); //min(x,15)
      diffabs = _mm_packus_epi16(diffabs, diffabs); //convert to 8
      diffabs = _mm_shuffle_epi8(lut, diffabs);//lut
      diffabs = _mm_cvtepi8_epi16(diffabs);//back to 16-bit
      diffabs = _mm_srai_epi16(diffabs, 1);//diagonal shift!
      diffabs = _mm_sign_epi16(diffabs, ld);//fix sign!
      acc = _mm_add_epi16(diffabs, acc); //add to acc
      //ru!
      //calculate abs
      diffabs = _mm_abs_epi16(ru); //abs
      diffabs = _mm_add_epi16(diffabs, four); //+4
      diffabs = _mm_srai_epi16(diffabs, 3); //>>3
      diffabs = _mm_min_epi16(diffabs, fifteen); //min(x,15)
      diffabs = _mm_packus_epi16(diffabs, diffabs); //convert to 8
      diffabs = _mm_shuffle_epi8(lut, diffabs);//lut
      diffabs = _mm_cvtepi8_epi16(diffabs);//back to 16-bit
      diffabs = _mm_srai_epi16(diffabs, 1);//diagonal shift!
      diffabs = _mm_sign_epi16(diffabs, ru);//fix sign!
      acc = _mm_add_epi16(diffabs, acc); //add to acc
      //rd!
      //calculate abs
      diffabs = _mm_abs_epi16(rd); //abs
      diffabs = _mm_add_epi16(diffabs, four); //+4
      diffabs = _mm_srai_epi16(diffabs, 3); //>>3
      diffabs = _mm_min_epi16(diffabs, fifteen); //min(x,15)
      diffabs = _mm_packus_epi16(diffabs, diffabs); //convert to 8
      diffabs = _mm_shuffle_epi8(lut, diffabs);//lut
      diffabs = _mm_cvtepi8_epi16(diffabs);//back to 16-bit
      diffabs = _mm_srai_epi16(diffabs, 1);//diagonal shift!
      diffabs = _mm_sign_epi16(diffabs, rd);//fix sign!
      acc = _mm_add_epi16(diffabs, acc); //add to acc

      //ll!
      //calculate abs
      diffabs = _mm_abs_epi16(ll); //abs
      diffabs = _mm_add_epi16(diffabs, four); //+4
      diffabs = _mm_srai_epi16(diffabs, 3); //>>3
      diffabs = _mm_min_epi16(diffabs, fifteen); //min(x,15)
      diffabs = _mm_packus_epi16(diffabs, diffabs); //convert to 8
      diffabs = _mm_shuffle_epi8(lut, diffabs);//lut
      diffabs = _mm_cvtepi8_epi16(diffabs);//back to 16-bit
      diffabs = _mm_srai_epi16(diffabs, 1);//diagonal shift!
      diffabs = _mm_sign_epi16(diffabs, ll);//fix sign!
      acc = _mm_add_epi16(diffabs, acc); //add to acc
      //rr!
      //calculate abs
      diffabs = _mm_abs_epi16(rr); //abs
      diffabs = _mm_add_epi16(diffabs, four); //+4
      diffabs = _mm_srai_epi16(diffabs, 3); //>>3
      diffabs = _mm_min_epi16(diffabs, fifteen); //min(x,15)
      diffabs = _mm_packus_epi16(diffabs, diffabs); //convert to 8
      diffabs = _mm_shuffle_epi8(lut, diffabs);//lut
      diffabs = _mm_cvtepi8_epi16(diffabs);//back to 16-bit
      diffabs = _mm_srai_epi16(diffabs, 1);//diagonal shift!
      diffabs = _mm_sign_epi16(diffabs, rr);//fix sign!
      acc = _mm_add_epi16(diffabs, acc); //add to acc
      //uu!
      //calculate abs
      diffabs = _mm_abs_epi16(uu); //abs
      diffabs = _mm_add_epi16(diffabs, four); //+4
      diffabs = _mm_srai_epi16(diffabs, 3); //>>3
      diffabs = _mm_min_epi16(diffabs, fifteen); //min(x,15)
      diffabs = _mm_packus_epi16(diffabs, diffabs); //convert to 8
      diffabs = _mm_shuffle_epi8(lut, diffabs);//lut
      diffabs = _mm_cvtepi8_epi16(diffabs);//back to 16-bit
      diffabs = _mm_srai_epi16(diffabs, 1);//diagonal shift!
      diffabs = _mm_sign_epi16(diffabs, uu);//fix sign!
      acc = _mm_add_epi16(diffabs, acc); //add to acc
      //dd!
      //calculate abs
      diffabs = _mm_abs_epi16(dd); //abs
      diffabs = _mm_add_epi16(diffabs, four); //+4
      diffabs = _mm_srai_epi16(diffabs, 3); //>>3
      diffabs = _mm_min_epi16(diffabs, fifteen); //min(x,15)
      diffabs = _mm_packus_epi16(diffabs, diffabs); //convert to 8
      diffabs = _mm_shuffle_epi8(lut, diffabs);//lut
      diffabs = _mm_cvtepi8_epi16(diffabs);//back to 16-bit
      diffabs = _mm_srai_epi16(diffabs, 1);//diagonal shift!
      diffabs = _mm_sign_epi16(diffabs, dd);//fix sign!
      acc = _mm_add_epi16(diffabs, acc); //add to acc

      if (bfac == 2)
      {
        acc = _mm_slli_epi16(acc, 1);   // Shift left to get 2*
      }
      else if (bfac == 3)
      {
        temp = _mm_slli_epi16(acc, 1);  // Multiply by two by shifting left
        acc = _mm_add_epi16(acc, temp); // Add original value to get 3*
      }

      // Add 16 and shift 5
      acc = _mm_add_epi16(acc, round_add);
      acc = _mm_srai_epi16(acc, bif_round_shift);

      // Instead we add our input values to the delta
      if (isRDO)
      {
        acc = _mm_add_epi16(acc, center);
      }
      else
      {
        int16_t* recpoint = &recPtr[row * recStride + col];
        inputVals = _mm_loadu_si128((__m128i*)(recpoint));
        acc = _mm_add_epi16(acc, inputVals);
      }

      // Clip
      if (isRDO)
      {
        acc = _mm_max_epi16(acc, clipmin); // No clipping applied in this function, will be clipped later on in CCSAO
        acc = _mm_min_epi16(acc, clipmax);
      }

      _mm_store_si128((__m128i*)(blkFilt + (row + pad) * (padwidth + 4) + col + pad), acc);
    }
  }

  // Copy back from tempbufFilter to recBuf
  int onerow = uiWidth * sizeof(Pel);
  for (uint32_t yy = 0; yy < uiHeight; yy++)
  {
    std::memcpy(recPtr, tempBlockPtr, onerow);
    recPtr += recStride;
    tempBlockPtr += tempBlockStride;
  }
}
#endif
template<X86_VEXT vext>
void BilateralFilter::simdFilterDiamond5x5(uint32_t uiWidth, uint32_t uiHeight, int16_t block[], int16_t blkFilt[], const ClpRng& clpRng, Pel* recPtr, int recStride, int iWidthExtSIMD, int bfac, int bif_round_add, int bif_round_shift, bool isRDO, const char* LUTrowPtr )
{
  //if( uiWidth < 4 || ( uiWidth < 8 && isRDO ) )
  if( uiWidth < 4 )
  {
    return blockBilateralFilterDiamond5x5(uiWidth, uiHeight, block, blkFilt, clpRng, recPtr, recStride, iWidthExtSIMD, bfac, bif_round_add, bif_round_shift, isRDO, LUTrowPtr );
  }

  int pad = 2;
  int padwidth = iWidthExtSIMD;

  __m128i center, left, right, up, down, lu, ld, ru, rd, diffabs, four, fifteen, lut, acc, temp, round_add, clipmin, clipmax, inputVals;
  __m128i ll, rr, uu, dd;

  four = _mm_set1_epi16(4);
  fifteen = _mm_set1_epi16(15);
  round_add = _mm_set1_epi16(bif_round_add);
  clipmin = _mm_set1_epi16(clpRng.min);
  clipmax = _mm_set1_epi16(clpRng.max);

  lut = _mm_loadu_si128((__m128i*)(LUTrowPtr));
  acc = _mm_set1_epi32(0);
  
  // Copy back parameters
  Pel *tempBlockPtr = (short*)blkFilt + (((padwidth+4) << 1) + 2);
  int tempBlockStride = padwidth+4;  
  
  for (int col = 0; col < uiWidth; col += 8)
  {
    for (int row = 0; row < uiHeight; row++)
    {
      acc = _mm_set1_epi32(0);
      int16_t *point = &block[(row + pad)*padwidth + pad + col];
      
      center = _mm_loadu_si128((__m128i*)(point));
      
      //load neighbours
      left = _mm_loadu_si128((__m128i*)(point - 1));
      right = _mm_loadu_si128((__m128i*)(point + 1));
      up = _mm_loadu_si128((__m128i*)(point - padwidth));
      down = _mm_loadu_si128((__m128i*)(point + padwidth));
      
      lu = _mm_loadu_si128((__m128i*)(point - 1 - padwidth));
      ld = _mm_loadu_si128((__m128i*)(point - 1 + padwidth));
      ru = _mm_loadu_si128((__m128i*)(point + 1 - padwidth));
      rd = _mm_loadu_si128((__m128i*)(point + 1 + padwidth));

      ll = _mm_loadu_si128((__m128i*)(point - 2));
      rr = _mm_loadu_si128((__m128i*)(point + 2));
      uu = _mm_loadu_si128((__m128i*)(point - 2*padwidth));
      dd = _mm_loadu_si128((__m128i*)(point + 2*padwidth));
      
      //calculate diffs
      left = _mm_sub_epi16(left, center);
      right = _mm_sub_epi16(right, center);
      up = _mm_sub_epi16(up, center);
      down = _mm_sub_epi16(down, center);
      
      lu = _mm_sub_epi16(lu, center);
      ld = _mm_sub_epi16(ld, center);
      ru = _mm_sub_epi16(ru, center);
      rd = _mm_sub_epi16(rd, center);

      ll = _mm_sub_epi16(ll, center);
      rr = _mm_sub_epi16(rr, center);
      uu = _mm_sub_epi16(uu, center);
      dd = _mm_sub_epi16(dd, center);
      
      //LEFT!
      //calculate abs
      diffabs = _mm_abs_epi16(left); //abs
      diffabs = _mm_add_epi16(diffabs, four); //+4
      diffabs = _mm_srai_epi16(diffabs, 3); //>>3
      diffabs = _mm_min_epi16(diffabs, fifteen); //min(x,15)
      diffabs = _mm_packus_epi16(diffabs, diffabs); //convert to 8
      diffabs = _mm_shuffle_epi8(lut, diffabs);//lut
      diffabs = _mm_cvtepi8_epi16(diffabs);//back to 16-bit
      diffabs = _mm_sign_epi16(diffabs, left);//fix sign!
      acc = _mm_add_epi16(diffabs, acc); //add to acc
      //RIGHT!
      //calculate abs
      diffabs = _mm_abs_epi16(right); //abs
      diffabs = _mm_add_epi16(diffabs, four); //+4
      diffabs = _mm_srai_epi16(diffabs, 3); //>>3
      diffabs = _mm_min_epi16(diffabs, fifteen); //min(x,15)
      diffabs = _mm_packus_epi16(diffabs, diffabs); //convert to 8
      diffabs = _mm_shuffle_epi8(lut, diffabs);//lut
      diffabs = _mm_cvtepi8_epi16(diffabs);//back to 16-bit
      diffabs = _mm_sign_epi16(diffabs, right);//fix sign!
      acc = _mm_add_epi16(diffabs, acc); //add to acc
      //UP!
      //calculate abs
      diffabs = _mm_abs_epi16(up); //abs
      diffabs = _mm_add_epi16(diffabs, four); //+4
      diffabs = _mm_srai_epi16(diffabs, 3); //>>3
      diffabs = _mm_min_epi16(diffabs, fifteen); //min(x,15)
      diffabs = _mm_packus_epi16(diffabs, diffabs); //convert to 8
      diffabs = _mm_shuffle_epi8(lut, diffabs);//lut
      diffabs = _mm_cvtepi8_epi16(diffabs);//back to 16-bit
      diffabs = _mm_sign_epi16(diffabs, up);//fix sign!
      acc = _mm_add_epi16(diffabs, acc); //add to acc
      
      //DOWN!
      //calculate abs
      diffabs = _mm_abs_epi16(down); //abs
      diffabs = _mm_add_epi16(diffabs, four); //+4
      diffabs = _mm_srai_epi16(diffabs, 3); //>>3
      diffabs = _mm_min_epi16(diffabs, fifteen); //min(x,15)
      diffabs = _mm_packus_epi16(diffabs, diffabs); //convert to 8
      diffabs = _mm_shuffle_epi8(lut, diffabs);//lut
      diffabs = _mm_cvtepi8_epi16(diffabs);//back to 16-bit
      diffabs = _mm_sign_epi16(diffabs, down);//fix sign!
      acc = _mm_add_epi16(diffabs, acc); //add to acc
      
      //lu!
      //calculate abs
      diffabs = _mm_abs_epi16(lu); //abs
      diffabs = _mm_add_epi16(diffabs, four); //+4
      diffabs = _mm_srai_epi16(diffabs, 3); //>>3
      diffabs = _mm_min_epi16(diffabs, fifteen); //min(x,15)
      diffabs = _mm_packus_epi16(diffabs, diffabs); //convert to 8
      diffabs = _mm_shuffle_epi8(lut, diffabs);//lut
      diffabs = _mm_cvtepi8_epi16(diffabs);//back to 16-bit
      diffabs = _mm_srai_epi16(diffabs, 1);//diagonal shift!
      diffabs = _mm_sign_epi16(diffabs, lu);//fix sign!
      acc = _mm_add_epi16(diffabs, acc); //add to acc
      //ld!
      //calculate abs
      diffabs = _mm_abs_epi16(ld); //abs
      diffabs = _mm_add_epi16(diffabs, four); //+4
      diffabs = _mm_srai_epi16(diffabs, 3); //>>3
      diffabs = _mm_min_epi16(diffabs, fifteen); //min(x,15)
      diffabs = _mm_packus_epi16(diffabs, diffabs); //convert to 8
      diffabs = _mm_shuffle_epi8(lut, diffabs);//lut
      diffabs = _mm_cvtepi8_epi16(diffabs);//back to 16-bit
      diffabs = _mm_srai_epi16(diffabs, 1);//diagonal shift!
      diffabs = _mm_sign_epi16(diffabs, ld);//fix sign!
      acc = _mm_add_epi16(diffabs, acc); //add to acc
      //ru!
      //calculate abs
      diffabs = _mm_abs_epi16(ru); //abs
      diffabs = _mm_add_epi16(diffabs, four); //+4
      diffabs = _mm_srai_epi16(diffabs, 3); //>>3
      diffabs = _mm_min_epi16(diffabs, fifteen); //min(x,15)
      diffabs = _mm_packus_epi16(diffabs, diffabs); //convert to 8
      diffabs = _mm_shuffle_epi8(lut, diffabs);//lut
      diffabs = _mm_cvtepi8_epi16(diffabs);//back to 16-bit
      diffabs = _mm_srai_epi16(diffabs, 1);//diagonal shift!
      diffabs = _mm_sign_epi16(diffabs, ru);//fix sign!
      acc = _mm_add_epi16(diffabs, acc); //add to acc
      //rd!
      //calculate abs
      diffabs = _mm_abs_epi16(rd); //abs
      diffabs = _mm_add_epi16(diffabs, four); //+4
      diffabs = _mm_srai_epi16(diffabs, 3); //>>3
      diffabs = _mm_min_epi16(diffabs, fifteen); //min(x,15)
      diffabs = _mm_packus_epi16(diffabs, diffabs); //convert to 8
      diffabs = _mm_shuffle_epi8(lut, diffabs);//lut
      diffabs = _mm_cvtepi8_epi16(diffabs);//back to 16-bit
      diffabs = _mm_srai_epi16(diffabs, 1);//diagonal shift!
      diffabs = _mm_sign_epi16(diffabs, rd);//fix sign!
      acc = _mm_add_epi16(diffabs, acc); //add to acc

      //ll!
      //calculate abs
      diffabs = _mm_abs_epi16(ll); //abs
      diffabs = _mm_add_epi16(diffabs, four); //+4
      diffabs = _mm_srai_epi16(diffabs, 3); //>>3
      diffabs = _mm_min_epi16(diffabs, fifteen); //min(x,15)
      diffabs = _mm_packus_epi16(diffabs, diffabs); //convert to 8
      diffabs = _mm_shuffle_epi8(lut, diffabs);//lut
      diffabs = _mm_cvtepi8_epi16(diffabs);//back to 16-bit
      diffabs = _mm_srai_epi16(diffabs, 1);//diagonal shift!
      diffabs = _mm_sign_epi16(diffabs, ll);//fix sign!
      acc = _mm_add_epi16(diffabs, acc); //add to acc
      //rr!
      //calculate abs
      diffabs = _mm_abs_epi16(rr); //abs
      diffabs = _mm_add_epi16(diffabs, four); //+4
      diffabs = _mm_srai_epi16(diffabs, 3); //>>3
      diffabs = _mm_min_epi16(diffabs, fifteen); //min(x,15)
      diffabs = _mm_packus_epi16(diffabs, diffabs); //convert to 8
      diffabs = _mm_shuffle_epi8(lut, diffabs);//lut
      diffabs = _mm_cvtepi8_epi16(diffabs);//back to 16-bit
      diffabs = _mm_srai_epi16(diffabs, 1);//diagonal shift!
      diffabs = _mm_sign_epi16(diffabs, rr);//fix sign!
      acc = _mm_add_epi16(diffabs, acc); //add to acc
      //uu!
      //calculate abs
      diffabs = _mm_abs_epi16(uu); //abs
      diffabs = _mm_add_epi16(diffabs, four); //+4
      diffabs = _mm_srai_epi16(diffabs, 3); //>>3
      diffabs = _mm_min_epi16(diffabs, fifteen); //min(x,15)
      diffabs = _mm_packus_epi16(diffabs, diffabs); //convert to 8
      diffabs = _mm_shuffle_epi8(lut, diffabs);//lut
      diffabs = _mm_cvtepi8_epi16(diffabs);//back to 16-bit
      diffabs = _mm_srai_epi16(diffabs, 1);//diagonal shift!
      diffabs = _mm_sign_epi16(diffabs, uu);//fix sign!
      acc = _mm_add_epi16(diffabs, acc); //add to acc
      //dd!
      //calculate abs
      diffabs = _mm_abs_epi16(dd); //abs
      diffabs = _mm_add_epi16(diffabs, four); //+4
      diffabs = _mm_srai_epi16(diffabs, 3); //>>3
      diffabs = _mm_min_epi16(diffabs, fifteen); //min(x,15)
      diffabs = _mm_packus_epi16(diffabs, diffabs); //convert to 8
      diffabs = _mm_shuffle_epi8(lut, diffabs);//lut
      diffabs = _mm_cvtepi8_epi16(diffabs);//back to 16-bit
      diffabs = _mm_srai_epi16(diffabs, 1);//diagonal shift!
      diffabs = _mm_sign_epi16(diffabs, dd);//fix sign!
      acc = _mm_add_epi16(diffabs, acc); //add to acc
      
      if (bfac == 2)
      {
        acc = _mm_slli_epi16(acc, 1);   // Shift left to get 2*
      }
      else if (bfac == 3)
      {
        temp = _mm_slli_epi16(acc, 1);  // Multiply by two by shifting left
        acc = _mm_add_epi16(acc, temp); // Add original value to get 3*
      }
      
      // Add 16 and shift 5
      acc = _mm_add_epi16(acc, round_add);
      acc = _mm_srai_epi16(acc, bif_round_shift);
      
      // Instead we add our input values to the delta
      if(isRDO)
      {
        acc = _mm_add_epi16(acc, center);
      }
      else
      {
        int16_t *recpoint = &recPtr[row * recStride + col];
        inputVals = _mm_loadu_si128((__m128i*)(recpoint));
        acc = _mm_add_epi16(acc, inputVals);
      }
      
      // Clip
      acc = _mm_max_epi16(acc, clipmin);
      acc = _mm_min_epi16(acc, clipmax);

      _mm_store_si128((__m128i*)(blkFilt + (row + pad) * (padwidth + 4) + col + pad), acc);
    }
  }
  
  // Copy back from tempbufFilter to recBuf
  int onerow = uiWidth * sizeof(Pel);
  for(uint32_t yy = 0; yy < uiHeight; yy++)
  {
    std::memcpy(recPtr, tempBlockPtr, onerow);
    recPtr += recStride;
    tempBlockPtr += tempBlockStride;
  }
}

template <X86_VEXT vext>
void BilateralFilter::_initBilateralFilterX86()
{
#if JVET_W0066_CCSAO   
  m_bilateralFilterDiamond5x5NoClip = simdFilterDiamond5x5NoClip<vext>;
#endif
  m_bilateralFilterDiamond5x5 = simdFilterDiamond5x5<vext>;  
}

template void BilateralFilter::_initBilateralFilterX86<SIMDX86>();
#endif
#endif   // TARGET_SIMD_X86
