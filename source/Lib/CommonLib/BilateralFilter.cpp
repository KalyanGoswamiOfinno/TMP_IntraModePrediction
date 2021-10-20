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

#include "BilateralFilter.h"
#if JVET_V0094_BILATERAL_FILTER
#include "Unit.h"
#include "UnitTools.h"
#endif

#if JVET_V0094_BILATERAL_FILTER
#include <tmmintrin.h>
#include <smmintrin.h>
#include <immintrin.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <algorithm>

#include "CommonLib/UnitTools.h"

BilateralFilter::BilateralFilter()
{
  m_bilateralFilterDiamond5x5 = blockBilateralFilterDiamond5x5;
#if JVET_W0066_CCSAO
  m_bilateralFilterDiamond5x5NoClip = blockBilateralFilterDiamond5x5NoClip;
#endif

#if ENABLE_SIMD_BILATERAL_FILTER
#ifdef TARGET_SIMD_X86
  initBilateralFilterX86();
#endif
#endif
}

BilateralFilter::~BilateralFilter()
{
}

void BilateralFilter::create()
{
}

void BilateralFilter::destroy()
{
}

const char* BilateralFilter::getFilterLutParameters( const int size, const PredMode predMode, const int32_t qp, int& bfac )
{
  if( size <= 4 )
  {
    bfac = 3;
  }
  else if( size >= 16 )
  {
    bfac = 1;
  }
  else
  {
    bfac = 2;
  }

  if( predMode == MODE_INTER )
  {
    if( size <= 4 )
    {
      bfac = 2;
    }
    else if( size >= 16 )
    {
      bfac = 1;
    }
    else
    {
      bfac = 2;
    }
  }

  int sqp = qp;

  if( sqp < 17 )
  {
    sqp = 17;
  }

  if( sqp > 42 )
  {
    sqp = 42;
  }

  return m_wBIF[sqp - 17];
}

#if JVET_W0066_CCSAO
void BilateralFilter::blockBilateralFilterDiamond5x5NoClip(uint32_t uiWidth, uint32_t uiHeight, int16_t block[], int16_t blkFilt[], const ClpRng& clpRng, Pel* recPtr, int recStride, int iWidthExtSIMD, int bfac, int bif_round_add, int bif_round_shift, bool isRDO, const char* LUTrowPtr)
{
  int pad = 2;

  int padwidth = iWidthExtSIMD;
  int downbuffer[64];
  int downleftbuffer[65];
  int downrightbuffer[2][65];
  int Shift, sg0, v0, idx, w0;
  Shift = sizeof(int) * 8 - 1;
  downbuffer[0] = 0;

  for (int x = 0; x < uiWidth; x++)
  {
    int pixel = block[(-1 + pad) * padwidth + x + pad];
    int below = block[(-1 + pad + 1) * padwidth + x + pad];
    int diff = below - pixel;
    sg0 = diff >> Shift;
    v0 = (diff + sg0) ^ sg0;
    v0 = (v0 + 4) >> 3;
    idx = 15 + ((v0 - 15) & ((v0 - 15) >> Shift));
    w0 = LUTrowPtr[idx];
    int mod = (w0 + sg0) ^ sg0;
    downbuffer[x] = mod;

    int belowright = block[(-1 + pad + 1) * padwidth + x + pad + 1];
    diff = belowright - pixel;
    sg0 = diff >> Shift;
    v0 = (diff + sg0) ^ sg0;
    v0 = (v0 + 4) >> 3;
    idx = 15 + ((v0 - 15) & ((v0 - 15) >> Shift));
    w0 = LUTrowPtr[idx] >> 1;
    mod = (w0 + sg0) ^ sg0;
    downrightbuffer[1][x + 1] = mod;

    int belowleft = block[(-1 + pad + 1) * padwidth + x + pad - 1];
    diff = belowleft - pixel;
    sg0 = diff >> Shift;
    v0 = (diff + sg0) ^ sg0;
    v0 = (v0 + 4) >> 3;
    idx = 15 + ((v0 - 15) & ((v0 - 15) >> Shift));
    w0 = LUTrowPtr[idx] >> 1;
    mod = (w0 + sg0) ^ sg0;
    downleftbuffer[x] = mod;
  }
  int width = uiWidth;
  for (int y = 0; y < uiHeight; y++)
  {
    int diff;

    int16_t* rowStart = &block[(y + pad) * padwidth + pad];

    int pixel = rowStart[-1];

    int right = rowStart[0];
    diff = right - pixel;
    sg0 = diff >> Shift;
    v0 = (diff + sg0) ^ sg0;
    v0 = (v0 + 4) >> 3;
    idx = 15 + ((v0 - 15) & ((v0 - 15) >> Shift));
    w0 = LUTrowPtr[idx];
    int mod = (w0 + sg0) ^ sg0;
    int rightmod = mod;

    pixel = rowStart[-padwidth - 1];
    int belowright = right;
    diff = belowright - pixel;
    sg0 = diff >> Shift;
    v0 = (diff + sg0) ^ sg0;
    v0 = (v0 + 4) >> 3;
    idx = 15 + ((v0 - 15) & ((v0 - 15) >> Shift));
    w0 = LUTrowPtr[idx] >> 1;
    mod = (w0 + sg0) ^ sg0;
    downrightbuffer[(y + 1) % 2][0] = mod;

    pixel = rowStart[-padwidth + width];
    int belowleft = rowStart[width - 1];
    diff = belowleft - pixel;
    sg0 = diff >> Shift;
    v0 = (diff + sg0) ^ sg0;
    v0 = (v0 + 4) >> 3;
    idx = 15 + ((v0 - 15) & ((v0 - 15) >> Shift));
    w0 = LUTrowPtr[idx] >> 1;
    mod = (w0 + sg0) ^ sg0;
    downleftbuffer[width] = mod;

    for (int x = 0; x < uiWidth; x++)
    {
      pixel = rowStart[x];
      int modsum = 0;

      int abovemod = -downbuffer[x];
      modsum += abovemod;

      int leftmod = -rightmod;
      modsum += leftmod;

      right = rowStart[x + 1];
      diff = right - pixel;
      sg0 = diff >> Shift;
      v0 = (diff + sg0) ^ sg0;
      v0 = (v0 + 4) >> 3;
      idx = 15 + ((v0 - 15) & ((v0 - 15) >> Shift));
      w0 = LUTrowPtr[idx];
      mod = (w0 + sg0) ^ sg0;

      modsum += mod;
      rightmod = mod;

      int below = rowStart[x + padwidth];
      diff = below - pixel;
      sg0 = diff >> Shift;
      v0 = (diff + sg0) ^ sg0;
      v0 = (v0 + 4) >> 3;
      idx = 15 + ((v0 - 15) & ((v0 - 15) >> Shift));
      w0 = LUTrowPtr[idx];
      mod = (w0 + sg0) ^ sg0;
      modsum += mod;
      downbuffer[x] = mod;

      int aboverightmod = -downleftbuffer[x + 1];
      // modsum += ((int16_t)((uint16_t)((aboverightmod) >> 1)));
      modsum += aboverightmod;

      int aboveleftmod = -downrightbuffer[(y + 1) % 2][x];
      // modsum += ((int16_t)((uint16_t)((aboveleftmod) >> 1)));
      modsum += aboveleftmod;

      int belowleft = rowStart[x + padwidth - 1];
      diff = belowleft - pixel;
      sg0 = diff >> Shift;
      v0 = (diff + sg0) ^ sg0;
      v0 = (v0 + 4) >> 3;
      idx = 15 + ((v0 - 15) & ((v0 - 15) >> Shift));
      w0 = LUTrowPtr[idx] >> 1;
      mod = (w0 + sg0) ^ sg0;
      // modsum += ((int16_t)((uint16_t)((mod) >> 1)));
      modsum += mod;
      downleftbuffer[x] = mod;

      int belowright = rowStart[x + padwidth + 1];
      diff = belowright - pixel;
      sg0 = diff >> Shift;
      v0 = (diff + sg0) ^ sg0;
      v0 = (v0 + 4) >> 3;
      idx = 15 + ((v0 - 15) & ((v0 - 15) >> Shift));
      w0 = LUTrowPtr[idx] >> 1;
      mod = (w0 + sg0) ^ sg0;
      //modsum += ((int16_t)((uint16_t)((mod) >> 1)));
      modsum += mod;
      downrightbuffer[y % 2][x + 1] = mod;

      // For samples two pixels out, we do not reuse previously calculated
      // values even though that is possible. Doing so would likely increase
      // speed when SIMD is turned off.

      int above = rowStart[x - 2 * padwidth];
      diff = above - pixel;
      sg0 = diff >> Shift;
      v0 = (diff + sg0) ^ sg0;
      v0 = (v0 + 4) >> 3;
      idx = 15 + ((v0 - 15) & ((v0 - 15) >> Shift));
      w0 = LUTrowPtr[idx] >> 1;
      mod = (w0 + sg0) ^ sg0;
      modsum += mod;

      below = rowStart[x + 2 * padwidth];
      diff = below - pixel;
      sg0 = diff >> Shift;
      v0 = (diff + sg0) ^ sg0;
      v0 = (v0 + 4) >> 3;
      idx = 15 + ((v0 - 15) & ((v0 - 15) >> Shift));
      w0 = LUTrowPtr[idx] >> 1;
      mod = (w0 + sg0) ^ sg0;
      modsum += mod;

      int left = rowStart[x - 2];
      diff = left - pixel;
      sg0 = diff >> Shift;
      v0 = (diff + sg0) ^ sg0;
      v0 = (v0 + 4) >> 3;
      idx = 15 + ((v0 - 15) & ((v0 - 15) >> Shift));
      w0 = LUTrowPtr[idx] >> 1;
      mod = (w0 + sg0) ^ sg0;
      modsum += mod;

      right = rowStart[x + 2];
      diff = right - pixel;
      sg0 = diff >> Shift;
      v0 = (diff + sg0) ^ sg0;
      v0 = (v0 + 4) >> 3;
      idx = 15 + ((v0 - 15) & ((v0 - 15) >> Shift));
      w0 = LUTrowPtr[idx] >> 1;
      mod = (w0 + sg0) ^ sg0;
      modsum += mod;

      blkFilt[(y + pad) * (padwidth + 4) + x + pad] = ((int16_t)((uint16_t)((modsum * bfac + bif_round_add) >> bif_round_shift)));
    }
  }

  // Copy back
  Pel* tempBlockPtr = (short*)blkFilt + (((padwidth + 4) << 1) + 2);
  int tempBlockStride = padwidth + 4;
  if (isRDO)
  {
    Pel* srcBlockPtr = (short*)block + (((padwidth) << 1) + 2);
    int srcBlockStride = padwidth;
    for (uint32_t yy = 0; yy < uiHeight; yy++)
    {
      for (uint32_t xx = 0; xx < uiWidth; xx++)
      {
        recPtr[xx] = ClipPel(srcBlockPtr[xx] + tempBlockPtr[xx], clpRng);
      }
      recPtr += recStride;
      tempBlockPtr += tempBlockStride;
      srcBlockPtr += srcBlockStride;
    }
  }
  else
  {
    for (uint32_t yy = 0; yy < uiHeight; yy++)
    {
      for (uint32_t xx = 0; xx < uiWidth; xx++)
      {
        // new result = old result (which is SAO-treated already) + diff due to bilateral filtering
        //recPtr[xx] = ClipPel<int>(recPtr[xx] + tempBlockPtr[xx], clpRng);
        recPtr[xx] = recPtr[xx] + tempBlockPtr[xx]; // clipping is done jointly for SAO/BIF/CCSAO
      }
      recPtr += recStride;
      tempBlockPtr += tempBlockStride;
    }
  }
}
#endif
void BilateralFilter::blockBilateralFilterDiamond5x5( uint32_t uiWidth, uint32_t uiHeight, int16_t block[], int16_t blkFilt[], const ClpRng& clpRng, Pel* recPtr, int recStride, int iWidthExtSIMD, int bfac, int bif_round_add, int bif_round_shift, bool isRDO, const char* LUTrowPtr )
{
  int pad = 2;

  int padwidth = iWidthExtSIMD;
  int downbuffer[64];
  int downleftbuffer[65];
  int downrightbuffer[2][65];
  int Shift, sg0, v0, idx, w0;
  Shift = sizeof( int ) * 8 - 1;
  downbuffer[0] = 0;

  for( int x = 0; x < uiWidth; x++ )
  {
    int pixel = block[(-1 + pad)*padwidth + x + pad];
    int below = block[(-1 + pad + 1)*padwidth + x + pad];
    int diff = below - pixel;
    sg0 = diff >> Shift;
    v0 = (diff + sg0) ^ sg0;
    v0 = (v0 + 4) >> 3;
    idx = 15 + ((v0 - 15)&((v0 - 15) >> Shift));
    w0 = LUTrowPtr[idx];
    int mod = (w0 + sg0) ^ sg0;
    downbuffer[x] = mod;

    int belowright = block[(-1 + pad + 1)*padwidth + x + pad + 1];
    diff = belowright - pixel;
    sg0 = diff >> Shift;
    v0 = (diff + sg0) ^ sg0;
    v0 = (v0 + 4) >> 3;
    idx = 15 + ((v0 - 15)&((v0 - 15) >> Shift));
    w0 = LUTrowPtr[idx] >> 1;
    mod = (w0 + sg0) ^ sg0;
    downrightbuffer[1][x + 1] = mod;

    int belowleft = block[(-1 + pad + 1)*padwidth + x + pad - 1];
    diff = belowleft - pixel;
    sg0 = diff >> Shift;
    v0 = (diff + sg0) ^ sg0;
    v0 = (v0 + 4) >> 3;
    idx = 15 + ((v0 - 15)&((v0 - 15) >> Shift));
    w0 = LUTrowPtr[idx] >> 1;
    mod = (w0 + sg0) ^ sg0;
    downleftbuffer[x] = mod;
  }
  int width = uiWidth;
  for( int y = 0; y < uiHeight; y++ )
  {
    int diff;

    int16_t *rowStart = &block[(y + pad)*padwidth + pad];

    int pixel = rowStart[-1];

    int right = rowStart[0];
    diff = right - pixel;
    sg0 = diff >> Shift;
    v0 = (diff + sg0) ^ sg0;
    v0 = (v0 + 4) >> 3;
    idx = 15 + ((v0 - 15)&((v0 - 15) >> Shift));
    w0 = LUTrowPtr[idx];
    int mod = (w0 + sg0) ^ sg0;
    int rightmod = mod;

    pixel = rowStart[-padwidth - 1];
    int belowright = right;
    diff = belowright - pixel;
    sg0 = diff >> Shift;
    v0 = (diff + sg0) ^ sg0;
    v0 = (v0 + 4) >> 3;
    idx = 15 + ((v0 - 15)&((v0 - 15) >> Shift));
    w0 = LUTrowPtr[idx] >> 1;
    mod = (w0 + sg0) ^ sg0;
    downrightbuffer[(y + 1) % 2][0] = mod;

    pixel = rowStart[-padwidth + width];
    int belowleft = rowStart[width - 1];
    diff = belowleft - pixel;
    sg0 = diff >> Shift;
    v0 = (diff + sg0) ^ sg0;
    v0 = (v0 + 4) >> 3;
    idx = 15 + ((v0 - 15)&((v0 - 15) >> Shift));
    w0 = LUTrowPtr[idx] >> 1;
    mod = (w0 + sg0) ^ sg0;
    downleftbuffer[width] = mod;

    for( int x = 0; x < uiWidth; x++ )
    {
      pixel = rowStart[x];
      int modsum = 0;

      int abovemod = -downbuffer[x];
      modsum += abovemod;

      int leftmod = -rightmod;
      modsum += leftmod;

      right = rowStart[x + 1];
      diff = right - pixel;
      sg0 = diff >> Shift;
      v0 = (diff + sg0) ^ sg0;
      v0 = (v0 + 4) >> 3;
      idx = 15 + ((v0 - 15)&((v0 - 15) >> Shift));
      w0 = LUTrowPtr[idx];
      mod = (w0 + sg0) ^ sg0;

      modsum += mod;
      rightmod = mod;

      int below = rowStart[x + padwidth];
      diff = below - pixel;
      sg0 = diff >> Shift;
      v0 = (diff + sg0) ^ sg0;
      v0 = (v0 + 4) >> 3;
      idx = 15 + ((v0 - 15)&((v0 - 15) >> Shift));
      w0 = LUTrowPtr[idx];
      mod = (w0 + sg0) ^ sg0;
      modsum += mod;
      downbuffer[x] = mod;

      int aboverightmod = -downleftbuffer[x + 1];
      // modsum += ((int16_t)((uint16_t)((aboverightmod) >> 1)));
      modsum += aboverightmod;

      int aboveleftmod = -downrightbuffer[(y + 1) % 2][x];
      // modsum += ((int16_t)((uint16_t)((aboveleftmod) >> 1)));
      modsum += aboveleftmod;

      int belowleft = rowStart[x + padwidth - 1];
      diff = belowleft - pixel;
      sg0 = diff >> Shift;
      v0 = (diff + sg0) ^ sg0;
      v0 = (v0 + 4) >> 3;
      idx = 15 + ((v0 - 15)&((v0 - 15) >> Shift));
      w0 = LUTrowPtr[idx] >> 1;
      mod = (w0 + sg0) ^ sg0;
      // modsum += ((int16_t)((uint16_t)((mod) >> 1)));
      modsum += mod;
      downleftbuffer[x] = mod;

      int belowright = rowStart[x + padwidth + 1];
      diff = belowright - pixel;
      sg0 = diff >> Shift;
      v0 = (diff + sg0) ^ sg0;
      v0 = (v0 + 4) >> 3;
      idx = 15 + ((v0 - 15)&((v0 - 15) >> Shift));
      w0 = LUTrowPtr[idx] >> 1;
      mod = (w0 + sg0) ^ sg0;
      //modsum += ((int16_t)((uint16_t)((mod) >> 1)));
      modsum += mod;
      downrightbuffer[y % 2][x + 1] = mod;

      // For samples two pixels out, we do not reuse previously calculated
      // values even though that is possible. Doing so would likely increase
      // speed when SIMD is turned off.

      int above = rowStart[x - 2 * padwidth];
      diff = above - pixel;
      sg0 = diff >> Shift;
      v0 = (diff + sg0) ^ sg0;
      v0 = (v0 + 4) >> 3;
      idx = 15 + ((v0 - 15)&((v0 - 15) >> Shift));
      w0 = LUTrowPtr[idx] >> 1;
      mod = (w0 + sg0) ^ sg0;
      modsum += mod;

      below = rowStart[x + 2 * padwidth];
      diff = below - pixel;
      sg0 = diff >> Shift;
      v0 = (diff + sg0) ^ sg0;
      v0 = (v0 + 4) >> 3;
      idx = 15 + ((v0 - 15)&((v0 - 15) >> Shift));
      w0 = LUTrowPtr[idx] >> 1;
      mod = (w0 + sg0) ^ sg0;
      modsum += mod;

      int left = rowStart[x - 2];
      diff = left - pixel;
      sg0 = diff >> Shift;
      v0 = (diff + sg0) ^ sg0;
      v0 = (v0 + 4) >> 3;
      idx = 15 + ((v0 - 15)&((v0 - 15) >> Shift));
      w0 = LUTrowPtr[idx] >> 1;
      mod = (w0 + sg0) ^ sg0;
      modsum += mod;

      right = rowStart[x + 2];
      diff = right - pixel;
      sg0 = diff >> Shift;
      v0 = (diff + sg0) ^ sg0;
      v0 = (v0 + 4) >> 3;
      idx = 15 + ((v0 - 15)&((v0 - 15) >> Shift));
      w0 = LUTrowPtr[idx] >> 1;
      mod = (w0 + sg0) ^ sg0;
      modsum += mod;

      blkFilt[(y + pad)*(padwidth + 4) + x + pad] = (( int16_t ) (( uint16_t ) ((modsum*bfac + bif_round_add) >> bif_round_shift)));
    }
  }

  // Copy back
  Pel *tempBlockPtr = ( short* ) blkFilt + (((padwidth + 4) << 1) + 2);
  int tempBlockStride = padwidth + 4;
  if( isRDO )
  {
    Pel *srcBlockPtr = ( short* ) block + (((padwidth) << 1) + 2);
    int srcBlockStride = padwidth;
    for( uint32_t yy = 0; yy < uiHeight; yy++ )
    {
      for( uint32_t xx = 0; xx < uiWidth; xx++ )
      {
        recPtr[xx] = ClipPel( srcBlockPtr[xx] + tempBlockPtr[xx], clpRng );
      }
      recPtr += recStride;
      tempBlockPtr += tempBlockStride;
      srcBlockPtr += srcBlockStride;
    }
  }
  else
  {
    for( uint32_t yy = 0; yy < uiHeight; yy++ )
    {
      for( uint32_t xx = 0; xx < uiWidth; xx++ )
      {
        // new result = old result (which is SAO-treated already) + diff due to bilateral filtering
        recPtr[xx] = ClipPel<int>( recPtr[xx] + tempBlockPtr[xx], clpRng );
      }
      recPtr += recStride;
      tempBlockPtr += tempBlockStride;
    }
  }
}

void BilateralFilter::bilateralFilterRDOdiamond5x5(PelBuf& resiBuf, const CPelBuf& predBuf, PelBuf& recoBuf, int32_t qp, const CPelBuf& recIPredBuf, const ClpRng& clpRng, TransformUnit & currTU, bool useReco, bool doReshape, std::vector<Pel>& pLUT)
{
  const unsigned uiWidth = predBuf.width;
  const unsigned uiHeight = predBuf.height;
  
  int bfac = 1;
  int bif_round_add = (BIF_ROUND_ADD) >> (currTU.cs->pps->getBIFStrength());
  int bif_round_shift = ( BIF_ROUND_SHIFT ) -(currTU.cs->pps->getBIFStrength());

  const char* LUTrowPtr = getFilterLutParameters( std::min( uiWidth, uiHeight ), currTU.cu->predMode, qp + currTU.cs->pps->getBIFQPOffset(), bfac );

  const unsigned uiPredStride = predBuf.stride;
  const unsigned uiStrideRes = resiBuf.stride;
  const unsigned uiRecStride = recoBuf.stride;
  const Pel *piPred = predBuf.buf;
  Pel *piResi = resiBuf.buf;
  Pel *piReco = recoBuf.buf;

  const Pel *piRecIPred = recIPredBuf.buf;
  const unsigned uiRecIPredStride = recIPredBuf.stride;

  const Pel *piPredTemp = piPred;
  Pel *piResiTemp = piResi;
  Pel *piRecoTemp = piReco;
  // Reco = Pred + Resi
  
  Pel *tempBlockPtr;
  
  uint32_t   uiWidthExt = uiWidth + (NUMBER_PADDED_SAMPLES << 1);
  uint32_t   uiHeightExt = uiHeight + (NUMBER_PADDED_SAMPLES << 1);
  
  int iWidthExtSIMD = uiWidthExt;
  if( uiWidth < 8 )
  {
    iWidthExtSIMD = 8 + (NUMBER_PADDED_SAMPLES << 1);
  }
  
  memset(tempblock, 0, iWidthExtSIMD*uiHeightExt * sizeof(short));
  tempBlockPtr = tempblock + (NUMBER_PADDED_SAMPLES)* iWidthExtSIMD + NUMBER_PADDED_SAMPLES;
  
  //// Clip and move block to temporary block
  if (useReco)
  {
    for (uint32_t uiY = 0; uiY < uiHeight; ++uiY)
    {
      std::memcpy(tempBlockPtr, piReco, uiWidth * sizeof(Pel));
      piReco += uiRecStride;
      tempBlockPtr += iWidthExtSIMD;
    }
    piReco = piRecoTemp;
  }
  else
  {
    for (uint32_t uiY = 0; uiY < uiHeight; ++uiY)
    {
      for (uint32_t uiX = 0; uiX < uiWidth; ++uiX)
      {
        tempBlockPtr[uiX] = ClipPel(piPred[uiX] + piResi[uiX], clpRng);
      }
      piPred += uiPredStride;
      piResi += uiStrideRes;
      piReco += uiRecStride;
      tempBlockPtr += iWidthExtSIMD;
    }
  }
  
  piPred = piPredTemp;
  piResi = piResiTemp;
  piReco = piRecoTemp;
  
  // Now do non-local filtering
  //
  // If a surrounding block is available, use samples from that block.
  // Otherwise, use block padded using odd padding:  . . a b c -> c b a b c
  
  // Pad entire block first and then overwrite with samples from surrounding blocks
  // if they pass the test.
  for (int yy = 1; yy< uiHeightExt -1 ; yy++)
  {
    tempblock[yy*iWidthExtSIMD + 1] = tempblock[yy*iWidthExtSIMD + 2];
    tempblock[yy*iWidthExtSIMD + uiWidthExt - 2] = tempblock[yy*iWidthExtSIMD + uiWidthExt - 3];
  }
  for (int xx = 1; xx< uiWidthExt - 1; xx++)
  {
    tempblock[1 * iWidthExtSIMD + xx] = tempblock[2 * iWidthExtSIMD + xx];
    tempblock[(uiHeightExt - 2)*iWidthExtSIMD + xx] = tempblock[(uiHeightExt - 3)*iWidthExtSIMD + xx];
  }
  
  bool subTuVer = currTU.lx() > currTU.cu->lx() ? true : false;
  bool subTuHor = currTU.ly() > currTU.cu->ly() ? true : false;
  
  bool isCTUboundary = currTU.ly() % currTU.cs->slice->getSPS()->getCTUSize() == 0;
  
  bool topAvailable = (currTU.ly() - 2 >= 0) && (currTU.ly() == currTU.cu->ly());
  topAvailable &= !isCTUboundary;
  bool leftAvailable = (currTU.lx() - 2 >= 0) && (currTU.lx() == currTU.cu->lx());
  
  // Fill in samples from blocks that pass the test
  if (topAvailable || leftAvailable || subTuVer || subTuHor)
  {
    const CompArea &area = currTU.blocks[COMPONENT_Y];
    CodingStructure& cs = *currTU.cs;
    
    if (topAvailable && leftAvailable)
    {
      // top left pixels
      tempblock[iWidthExtSIMD + 1] = *(piRecIPred - (uiRecIPredStride)-1);
      // Reshape copied pixels if necessary.
      if(doReshape)
      {
        tempblock[iWidthExtSIMD + 1] = pLUT[tempblock[iWidthExtSIMD + 1]];
      }
    }
    
    // top row
    if (topAvailable)
    {
      for (int blockx = 0; blockx < area.width; blockx += 1)
      {
        // copy 4 pixels one line above block from block to blockx + 3
        std::copy(piRecIPred - (uiRecIPredStride)+blockx, piRecIPred - (uiRecIPredStride)+blockx + 1, tempblock + 2 + iWidthExtSIMD + blockx);
        if( doReshape )
        {
          for( int xx = 0; xx < 1; xx++ )
          {
            tempblock[2 + iWidthExtSIMD + blockx + xx] = pLUT[tempblock[2 + iWidthExtSIMD + blockx + xx]];
          }
        }
      }
    }
    else if (subTuHor)
    {
      const CompArea &prevHalfArea = currTU.prev->blocks[COMPONENT_Y];
      CPelBuf earlierHalfBuf = cs.getPredBuf(prevHalfArea);
      earlierHalfBuf = cs.getRecoBuf(prevHalfArea);
      const unsigned earlierStride = earlierHalfBuf.stride;
      const Pel *earlierPel = earlierHalfBuf.buf + (currTU.prev->lheight() - 1)*earlierStride;
      
      std::copy(earlierPel, earlierPel + area.width, tempblock + 2 + iWidthExtSIMD);
      if( doReshape )
      {
        for( int xx = 0; xx < area.width; xx++ )
        {
          tempblock[2 + iWidthExtSIMD + xx] = pLUT[tempblock[2 + iWidthExtSIMD + xx]];
        }
      }
    }
    // left column
    if (leftAvailable)
    {
      for (int blocky = 0; blocky < area.height; blocky += 1)
      {
        tempblock[(iWidthExtSIMD << 1) + (blocky + 0) * iWidthExtSIMD + 1] = *(piRecIPred + (blocky + 0)*uiRecIPredStride - 1); // 1 pel out
        if(doReshape)
        {
          tempblock[(iWidthExtSIMD << 1) + (blocky + 0) * iWidthExtSIMD + 1] = pLUT[tempblock[(iWidthExtSIMD << 1) + (blocky + 0) * iWidthExtSIMD + 1]];
        }
      }
    }
    else if (subTuVer)
    {
      const CompArea &prevHalfArea = currTU.prev->blocks[COMPONENT_Y];
      CPelBuf earlierHalfBuf = cs.getPredBuf(prevHalfArea);
      earlierHalfBuf = cs.getRecoBuf(prevHalfArea);
      const unsigned earlierStride = earlierHalfBuf.stride;
      const Pel *earlierPel = earlierHalfBuf.buf + (currTU.prev->lwidth() - 1); // second to last pixel of first row of previous block
      
      for (int yy = 0; yy < currTU.lheight(); yy++)
      {
        tempblock[(iWidthExtSIMD << 1) + yy * iWidthExtSIMD + 1] = *(earlierPel + yy*earlierStride + 0);
      }
      if(doReshape)
      {
        for (int yy = 0; yy < currTU.lheight(); yy++)
        {
          tempblock[(iWidthExtSIMD << 1) + yy * iWidthExtSIMD + 1] = pLUT[tempblock[(iWidthExtSIMD << 1) + yy * iWidthExtSIMD + 1]];
        }
      }
    }
  }
  
  // Sloppy copying of outer layer
  for(int yy = 0; yy < uiHeight+2; yy++)
  {
    tempblock[iWidthExtSIMD + yy*iWidthExtSIMD] = tempblock[iWidthExtSIMD + yy*iWidthExtSIMD + 1];
    tempblock[iWidthExtSIMD + uiWidthExt - 1 + yy*iWidthExtSIMD] = tempblock[iWidthExtSIMD + uiWidthExt - 2 + yy*iWidthExtSIMD];
  }
  std::copy(tempblock + iWidthExtSIMD, tempblock + iWidthExtSIMD + uiWidthExt, tempblock);
  std::copy(tempblock  + iWidthExtSIMD*(uiHeightExt-2), tempblock  + iWidthExtSIMD*(uiHeightExt-2) + uiWidthExt, tempblock + iWidthExtSIMD*(uiHeightExt-1));

  m_bilateralFilterDiamond5x5(uiWidth, uiHeight, tempblock, tempblockFiltered, clpRng, piReco, uiRecStride, iWidthExtSIMD, bfac, bif_round_add, bif_round_shift, true, LUTrowPtr );

  if (!useReco)
  {
    // need to be performed if residual  is used
    // Resi' = Reco' - Pred
    for (uint32_t uiY = 0; uiY < uiHeight; ++uiY)
    {
      for (uint32_t uiX = 0; uiX < uiWidth; ++uiX)
      {
        piResi[uiX] = piReco[uiX] - piPred[uiX];
      }
      piPred += uiPredStride;
      piResi += uiStrideRes;
      piReco += uiRecStride;
    }
  }
}

#if JVET_W0066_CCSAO
void BilateralFilter::bilateralFilterDiamond5x5NoClip(const CPelUnitBuf& src, PelUnitBuf& rec, int32_t qp, const ClpRng& clpRng, TransformUnit& currTU)
{
  CompArea& compArea = currTU.block(COMPONENT_Y);

  const unsigned uiWidth = compArea.width;
  const unsigned uiHeight = compArea.height;

  bool topAltAvailable;
  bool leftAltAvailable;

  int srcStride = src.get(COMPONENT_Y).stride;
  const Pel* srcPtr = src.get(COMPONENT_Y).bufAt(compArea);
  const Pel* srcPtrTemp = srcPtr;

  int recStride = rec.get(COMPONENT_Y).stride;
  Pel* recPtr = rec.get(COMPONENT_Y).bufAt(compArea);

  int bfac = 1;
  const char* LUTrowPtr = getFilterLutParameters(std::min(uiWidth, uiHeight), currTU.cu->predMode, qp + currTU.cs->pps->getBIFQPOffset(), bfac);

  int bif_round_add = (BIF_ROUND_ADD) >> (currTU.cs->pps->getBIFStrength());
  int bif_round_shift = (BIF_ROUND_SHIFT)-(currTU.cs->pps->getBIFStrength());

  const CompArea& myArea = currTU.blocks[COMPONENT_Y];
  topAltAvailable = myArea.y - 2 >= 0;
  leftAltAvailable = myArea.x - 2 >= 0;
  bool bottomAltAvailable = myArea.y + myArea.height + 1 < currTU.cu->slice->getSPS()->getMaxPicHeightInLumaSamples();
  bool rightAltAvailable = myArea.x + myArea.width + 1 < currTU.cu->slice->getSPS()->getMaxPicWidthInLumaSamples();

  uint32_t   uiWidthExt = uiWidth + (NUMBER_PADDED_SAMPLES << 1);
  uint32_t   uiHeightExt = uiHeight + (NUMBER_PADDED_SAMPLES << 1);

  int iWidthExtSIMD = uiWidthExt;
  if (uiWidth < 8)
  {
    iWidthExtSIMD = 8 + (NUMBER_PADDED_SAMPLES << 1);
  }

  Pel* tempBlockPtr;

  bool allAvail = topAltAvailable && bottomAltAvailable && leftAltAvailable && rightAltAvailable;

  memset(tempblock, 0, iWidthExtSIMD * uiHeightExt * sizeof(short));

  if (allAvail)
  {
    // set pointer two rows up and two pixels to the left from the start of the block
    tempBlockPtr = tempblock;

    // same with image data
    srcPtr = srcPtr - 2 * srcStride - 2;

    //// Move block to temporary block

    // Check if the block a the top block of a CTU.
    bool isCTUboundary = myArea.y % currTU.cs->slice->getSPS()->getCTUSize() == 0;
    if (isCTUboundary)
    {
      // The samples two lines up are out of bounds. (One line above the CTU is OK, since SAO uses that line.)
      // Hence the top line of tempblock is unavailable if the block is the top block of a CTU.
      // Therefore, copy samples from one line up instead of from two lines up by updating srcPtr *before* copy.
      srcPtr += srcStride;
      std::memcpy(tempBlockPtr, srcPtr, (uiWidthExt) * sizeof(Pel));
    }
    else
    {
      std::memcpy(tempBlockPtr, srcPtr, (uiWidthExt) * sizeof(Pel));
      srcPtr += srcStride;
    }
    tempBlockPtr += iWidthExtSIMD;
    // Copy samples that are not out of bounds.
    for (uint32_t uiY = 1; uiY < uiHeightExt - 1; ++uiY)
    {
      std::memcpy(tempBlockPtr, srcPtr, (uiWidthExt) * sizeof(Pel));
      srcPtr += srcStride;
      tempBlockPtr += iWidthExtSIMD;
    }
    // Check if the block is a bottom block of a CTU.
    isCTUboundary = (myArea.y + uiHeight) % currTU.cs->slice->getSPS()->getCTUSize() == 0;
    if (isCTUboundary)
    {
      // The samples two lines down are out of bounds. (One line below the CTU is OK, since SAO uses that line.)
      // Hence the bottom line of tempblock is unavailable if the block at the bottom of a CTU.
      // Therefore, copy samples from the second to last line instead of the last line by subtracting srcPtr before copy.
      srcPtr -= srcStride;
      std::memcpy(tempBlockPtr, srcPtr, (uiWidthExt) * sizeof(Pel));
    }
    else
    {
      std::memcpy(tempBlockPtr, srcPtr, (uiWidthExt) * sizeof(Pel));
    }
    return m_bilateralFilterDiamond5x5NoClip(uiWidth, uiHeight, tempblock, tempblockFiltered, clpRng, recPtr, recStride, iWidthExtSIMD, bfac, bif_round_add, bif_round_shift, false, LUTrowPtr);
  }
  else
  {
    tempBlockPtr = tempblock + (NUMBER_PADDED_SAMPLES)*iWidthExtSIMD + NUMBER_PADDED_SAMPLES;

    //// Move block to temporary block
    for (uint32_t uiY = 0; uiY < uiHeight; ++uiY)
    {
      std::memcpy(tempBlockPtr, srcPtr, uiWidth * sizeof(Pel));
      srcPtr += srcStride;
      tempBlockPtr += iWidthExtSIMD;
    }
    srcPtr = srcPtrTemp;

    if (topAltAvailable)
    {
      std::copy(srcPtr - 2 * srcStride, srcPtr - 2 * srcStride + uiWidth, tempblock + 2);
      std::copy(srcPtr - srcStride, srcPtr - srcStride + uiWidth, tempblock + iWidthExtSIMD + 2);
    }
    if (bottomAltAvailable)
    {
      std::copy(srcPtr + (uiHeight + 1) * srcStride, srcPtr + (uiHeight + 1) * srcStride + uiWidth, tempblock + (uiHeightExt - 1) * iWidthExtSIMD + 2);
      std::copy(srcPtr + uiHeight * srcStride, srcPtr + uiHeight * srcStride + uiWidth, tempblock + (uiHeightExt - 2) * iWidthExtSIMD + 2);
    }
    if (leftAltAvailable)
    {
      for (int yy = 0; yy < uiHeight; yy++)
      {
        tempblock[(iWidthExtSIMD << 1) + yy * iWidthExtSIMD + 0] = *(srcPtr + yy * srcStride - 2);
        tempblock[(iWidthExtSIMD << 1) + yy * iWidthExtSIMD + 1] = *(srcPtr + yy * srcStride - 1);
      }
    }
    if (rightAltAvailable)
    {
      for (int yy = 0; yy < uiHeight; yy++)
      {
        tempblock[(iWidthExtSIMD << 1) + uiWidthExt - 1 + yy * iWidthExtSIMD] = *(srcPtr + uiWidth + yy * srcStride + 1);
        tempblock[(iWidthExtSIMD << 1) + uiWidthExt - 2 + yy * iWidthExtSIMD] = *(srcPtr + uiWidth + yy * srcStride);
      }
    }

    // if not all available, copy from inside tempbuffer
    if (!topAltAvailable)
    {
      std::copy(tempblock + iWidthExtSIMD * 2 + 2, tempblock + iWidthExtSIMD * 2 + 2 + uiWidth, tempblock + 2);
      std::copy(tempblock + iWidthExtSIMD * 2 + 2, tempblock + iWidthExtSIMD * 2 + 2 + uiWidth, tempblock + iWidthExtSIMD + 2);
    }
    if (!bottomAltAvailable)
    {
      std::copy(tempblock + (uiHeightExt - 3) * iWidthExtSIMD + 2, tempblock + (uiHeightExt - 3) * iWidthExtSIMD + 2 + uiWidth, tempblock + (uiHeightExt - 2) * iWidthExtSIMD + 2);
      std::copy(tempblock + (uiHeightExt - 3) * iWidthExtSIMD + 2, tempblock + (uiHeightExt - 3) * iWidthExtSIMD + 2 + uiWidth, tempblock + (uiHeightExt - 1) * iWidthExtSIMD + 2);
    }
    if (!leftAltAvailable)
    {
      for (int yy = 0; yy < uiHeight; yy++)
      {
        tempblock[(iWidthExtSIMD << 1) + yy * iWidthExtSIMD + 0] = tempblock[(iWidthExtSIMD << 1) + yy * iWidthExtSIMD + 2];
        tempblock[(iWidthExtSIMD << 1) + yy * iWidthExtSIMD + 1] = tempblock[(iWidthExtSIMD << 1) + yy * iWidthExtSIMD + 2];
      }
    }
    if (!rightAltAvailable)
    {
      for (int yy = 0; yy < uiHeight; yy++)
      {
        tempblock[(iWidthExtSIMD << 1) + uiWidthExt - 2 + yy * iWidthExtSIMD] = tempblock[(iWidthExtSIMD << 1) + uiWidthExt - 2 + yy * iWidthExtSIMD - 1];
        tempblock[(iWidthExtSIMD << 1) + uiWidthExt - 1 + yy * iWidthExtSIMD] = tempblock[(iWidthExtSIMD << 1) + uiWidthExt - 2 + yy * iWidthExtSIMD - 1];
      }
    }

    // All sides are available, easy to just copy corners also.
    if (topAltAvailable && leftAltAvailable)
    {
      tempblock[0] = *(srcPtr - 2 * srcStride - 2);                                               // a     top left corner
      tempblock[1] = *(srcPtr - 2 * srcStride - 1);                                               // b     a b|x x
      tempblock[iWidthExtSIMD + 0] = *(srcPtr - srcStride - 2);                                   // c     c d|x x
      tempblock[iWidthExtSIMD + 1] = *(srcPtr - srcStride - 1);                                   // d     -------
    }
    else
    {
      tempblock[0] = tempblock[iWidthExtSIMD * 2 + 2];                                            // extend top left
      tempblock[1] = tempblock[iWidthExtSIMD * 2 + 2];                                            // extend top left
      tempblock[iWidthExtSIMD + 0] = tempblock[iWidthExtSIMD * 2 + 2];                            // extend top left
      tempblock[iWidthExtSIMD + 1] = tempblock[iWidthExtSIMD * 2 + 2];                            // extend top left
    }

    if (topAltAvailable && rightAltAvailable)
    {
      tempblock[iWidthExtSIMD - 2] = *(srcPtr - 2 * srcStride + uiWidth);                         // a
      tempblock[iWidthExtSIMD - 1] = *(srcPtr - 2 * srcStride + uiWidth + 1);                     // b
      tempblock[iWidthExtSIMD + uiWidthExt - 2] = *(srcPtr - srcStride + uiWidth);                // c
      tempblock[iWidthExtSIMD + uiWidthExt - 1] = *(srcPtr - srcStride + uiWidth + 1);            // d
    }
    else
    {
      tempblock[iWidthExtSIMD - 2] = tempblock[iWidthExtSIMD * 2 + uiWidthExt - 3];               // extend top right
      tempblock[iWidthExtSIMD - 1] = tempblock[iWidthExtSIMD * 2 + uiWidthExt - 3];               // extend top right
      tempblock[iWidthExtSIMD + uiWidthExt - 2] = tempblock[iWidthExtSIMD * 2 + uiWidthExt - 3];  // extend top right
      tempblock[iWidthExtSIMD + uiWidthExt - 1] = tempblock[iWidthExtSIMD * 2 + uiWidthExt - 3];  // extend top right
    }

    if (bottomAltAvailable && leftAltAvailable)
    {
      tempblock[iWidthExtSIMD * (uiHeightExt - 2) + 0] = *(srcPtr + uiHeight * srcStride - 2);          // a
      tempblock[iWidthExtSIMD * (uiHeightExt - 2) + 1] = *(srcPtr + uiHeight * srcStride - 1);          // b
      tempblock[iWidthExtSIMD * (uiHeightExt - 1) + 0] = *(srcPtr + (uiHeight + 1) * srcStride - 2);    // c
      tempblock[iWidthExtSIMD * (uiHeightExt - 1) + 1] = *(srcPtr + (uiHeight + 1) * srcStride - 1);    // d
    }
    else
    {
      tempblock[iWidthExtSIMD * (uiHeightExt - 2) + 0] = tempblock[iWidthExtSIMD * (uiHeightExt - 3) + 2];  // bot avail: mirror left/right
      tempblock[iWidthExtSIMD * (uiHeightExt - 2) + 1] = tempblock[iWidthExtSIMD * (uiHeightExt - 3) + 2];  // bot avail: mirror left/right
      tempblock[iWidthExtSIMD * (uiHeightExt - 1) + 0] = tempblock[iWidthExtSIMD * (uiHeightExt - 3) + 2];  // bot avail: mirror left/right
      tempblock[iWidthExtSIMD * (uiHeightExt - 1) + 1] = tempblock[iWidthExtSIMD * (uiHeightExt - 3) + 2];  // bot avail: mirror left/right
    }

    if (bottomAltAvailable && rightAltAvailable)
    {
      tempblock[iWidthExtSIMD * (uiHeightExt - 2) + uiWidthExt - 2] = *(srcPtr + uiHeight * srcStride + uiWidth);                // a
      tempblock[iWidthExtSIMD * (uiHeightExt - 2) + uiWidthExt - 1] = *(srcPtr + uiHeight * srcStride + uiWidth + 1);            // b
      tempblock[iWidthExtSIMD * (uiHeightExt - 1) + uiWidthExt - 2] = *(srcPtr + (uiHeight + 1) * srcStride + uiWidth);          // c
      tempblock[iWidthExtSIMD * (uiHeightExt - 1) + uiWidthExt - 1] = *(srcPtr + (uiHeight + 1) * srcStride + uiWidth + 1);      // d
    }
    else
    {
      tempblock[iWidthExtSIMD * (uiHeightExt - 2) + uiWidthExt - 2] = tempblock[iWidthExtSIMD * (uiHeightExt - 3) + uiWidthExt - 3];
      tempblock[iWidthExtSIMD * (uiHeightExt - 2) + uiWidthExt - 1] = tempblock[iWidthExtSIMD * (uiHeightExt - 3) + uiWidthExt - 3];
      tempblock[iWidthExtSIMD * (uiHeightExt - 1) + uiWidthExt - 2] = tempblock[iWidthExtSIMD * (uiHeightExt - 3) + uiWidthExt - 3];
      tempblock[iWidthExtSIMD * (uiHeightExt - 1) + uiWidthExt - 1] = tempblock[iWidthExtSIMD * (uiHeightExt - 3) + uiWidthExt - 3];
    }
  }

  m_bilateralFilterDiamond5x5NoClip(uiWidth, uiHeight, tempblock, tempblockFiltered, clpRng, recPtr, recStride, iWidthExtSIMD, bfac, bif_round_add, bif_round_shift, false, LUTrowPtr);
}
#endif

void BilateralFilter::bilateralFilterDiamond5x5(const CPelUnitBuf& src, PelUnitBuf& rec, int32_t qp, const ClpRng& clpRng, TransformUnit & currTU)
{
  CompArea &compArea = currTU.block(COMPONENT_Y);
  
  const unsigned uiWidth = compArea.width;
  const unsigned uiHeight = compArea.height;
  
  bool topAltAvailable;
  bool leftAltAvailable;
  
  int srcStride = src.get(COMPONENT_Y).stride;
  const Pel *srcPtr = src.get(COMPONENT_Y).bufAt(compArea);
  const Pel *srcPtrTemp = srcPtr;
  
  int recStride = rec.get(COMPONENT_Y).stride;
  Pel *recPtr = rec.get(COMPONENT_Y).bufAt(compArea);
  
  int bfac = 1;  
  const char* LUTrowPtr = getFilterLutParameters( std::min( uiWidth, uiHeight ), currTU.cu->predMode, qp + currTU.cs->pps->getBIFQPOffset(), bfac );
  
  int bif_round_add = (BIF_ROUND_ADD) >> (currTU.cs->pps->getBIFStrength());
  int bif_round_shift = (BIF_ROUND_SHIFT) - (currTU.cs->pps->getBIFStrength());
  
  const CompArea &myArea = currTU.blocks[COMPONENT_Y];
  topAltAvailable = myArea.y - 2 >= 0;
  leftAltAvailable = myArea.x - 2 >= 0;
  bool bottomAltAvailable = myArea.y + myArea.height + 1 < currTU.cu->slice->getSPS()->getMaxPicHeightInLumaSamples();
  bool rightAltAvailable = myArea.x + myArea.width + 1 < currTU.cu->slice->getSPS()->getMaxPicWidthInLumaSamples();
  
  uint32_t   uiWidthExt = uiWidth + (NUMBER_PADDED_SAMPLES << 1);
  uint32_t   uiHeightExt = uiHeight + (NUMBER_PADDED_SAMPLES << 1);
  
  int iWidthExtSIMD = uiWidthExt;
  if( uiWidth < 8 )
  {
    iWidthExtSIMD = 8 + (NUMBER_PADDED_SAMPLES << 1);
  }
  
  Pel *tempBlockPtr;
  
  bool allAvail = topAltAvailable && bottomAltAvailable && leftAltAvailable && rightAltAvailable;

  memset(tempblock, 0, iWidthExtSIMD*uiHeightExt * sizeof(short));

  if(allAvail)
  {
    // set pointer two rows up and two pixels to the left from the start of the block
    tempBlockPtr = tempblock;
    
    // same with image data
    srcPtr = srcPtr - 2*srcStride - 2;

    //// Move block to temporary block

    // Check if the block a the top block of a CTU.
    bool isCTUboundary = myArea.y % currTU.cs->slice->getSPS()->getCTUSize() == 0;
    if(isCTUboundary)
    {
      // The samples two lines up are out of bounds. (One line above the CTU is OK, since SAO uses that line.)
      // Hence the top line of tempblock is unavailable if the block is the top block of a CTU.
      // Therefore, copy samples from one line up instead of from two lines up by updating srcPtr *before* copy.
      srcPtr += srcStride;
      std::memcpy(tempBlockPtr, srcPtr, (uiWidthExt) * sizeof(Pel));
    }
    else
    {
      std::memcpy(tempBlockPtr, srcPtr, (uiWidthExt) * sizeof(Pel));
      srcPtr += srcStride;
    }
    tempBlockPtr += iWidthExtSIMD;
    // Copy samples that are not out of bounds.
    for (uint32_t uiY = 1; uiY < uiHeightExt-1; ++uiY)
    {
      std::memcpy(tempBlockPtr, srcPtr, (uiWidthExt) * sizeof(Pel));
      srcPtr += srcStride;
      tempBlockPtr += iWidthExtSIMD;
    }
    // Check if the block is a bottom block of a CTU.
    isCTUboundary = (myArea.y + uiHeight) % currTU.cs->slice->getSPS()->getCTUSize() == 0;
    if(isCTUboundary)
    {
      // The samples two lines down are out of bounds. (One line below the CTU is OK, since SAO uses that line.)
      // Hence the bottom line of tempblock is unavailable if the block at the bottom of a CTU.
      // Therefore, copy samples from the second to last line instead of the last line by subtracting srcPtr before copy.
      srcPtr -= srcStride;
      std::memcpy(tempBlockPtr, srcPtr, (uiWidthExt) * sizeof(Pel));
    }
    else
    {
      std::memcpy(tempBlockPtr, srcPtr, (uiWidthExt) * sizeof(Pel));
    }
    return m_bilateralFilterDiamond5x5(uiWidth, uiHeight, tempblock, tempblockFiltered, clpRng, recPtr, recStride, iWidthExtSIMD, bfac, bif_round_add, bif_round_shift, false, LUTrowPtr );
  }
  else
  {
    tempBlockPtr = tempblock + (NUMBER_PADDED_SAMPLES)* iWidthExtSIMD + NUMBER_PADDED_SAMPLES;

    //// Move block to temporary block
    for (uint32_t uiY = 0; uiY < uiHeight; ++uiY)
    {
      std::memcpy(tempBlockPtr, srcPtr, uiWidth * sizeof(Pel));
      srcPtr += srcStride;
      tempBlockPtr += iWidthExtSIMD;
    }
    srcPtr = srcPtrTemp;
    
    if(topAltAvailable)
    {
      std::copy(srcPtr - 2*srcStride, srcPtr - 2*srcStride + uiWidth, tempblock + 2);
      std::copy(srcPtr - srcStride, srcPtr - srcStride + uiWidth, tempblock + iWidthExtSIMD + 2);
    }
    if(bottomAltAvailable)
    {
      std::copy(srcPtr + (uiHeight+1)*srcStride, srcPtr +(uiHeight+1)*srcStride + uiWidth, tempblock + (uiHeightExt-1)*iWidthExtSIMD + 2);
      std::copy(srcPtr + uiHeight*srcStride, srcPtr +uiHeight*srcStride + uiWidth, tempblock + (uiHeightExt-2)*iWidthExtSIMD + 2);
    }
    if(leftAltAvailable)
    {
      for(int yy = 0; yy<uiHeight; yy++)
      {
        tempblock[(iWidthExtSIMD<<1) + yy*iWidthExtSIMD + 0] = *(srcPtr + yy*srcStride -2);
        tempblock[(iWidthExtSIMD<<1) + yy*iWidthExtSIMD + 1] = *(srcPtr + yy*srcStride -1);
      }
    }
    if(rightAltAvailable)
    {
      for(int yy = 0; yy<uiHeight; yy++)
      {
        tempblock[(iWidthExtSIMD<<1) + uiWidthExt-1 + yy*iWidthExtSIMD] = *(srcPtr + uiWidth + yy*srcStride + 1);
        tempblock[(iWidthExtSIMD<<1) + uiWidthExt-2 + yy*iWidthExtSIMD] = *(srcPtr + uiWidth + yy*srcStride);
      }
    }
    
    // if not all available, copy from inside tempbuffer
    if(!topAltAvailable)
    {
      std::copy(tempblock + iWidthExtSIMD*2 + 2, tempblock + iWidthExtSIMD*2 + 2 + uiWidth, tempblock + 2);
      std::copy(tempblock + iWidthExtSIMD*2 + 2, tempblock + iWidthExtSIMD*2 + 2 + uiWidth, tempblock + iWidthExtSIMD + 2);
    }
    if(!bottomAltAvailable)
    {
      std::copy(tempblock + (uiHeightExt-3)*iWidthExtSIMD + 2, tempblock + (uiHeightExt-3)*iWidthExtSIMD + 2 + uiWidth, tempblock + (uiHeightExt-2)*iWidthExtSIMD + 2);
      std::copy(tempblock + (uiHeightExt-3)*iWidthExtSIMD + 2, tempblock + (uiHeightExt-3)*iWidthExtSIMD + 2 + uiWidth, tempblock + (uiHeightExt-1)*iWidthExtSIMD + 2);
    }
    if(!leftAltAvailable)
    {
      for(int yy = 0; yy<uiHeight; yy++)
      {
        tempblock[(iWidthExtSIMD<<1) + yy*iWidthExtSIMD + 0] = tempblock[(iWidthExtSIMD<<1) + yy*iWidthExtSIMD + 2];
        tempblock[(iWidthExtSIMD<<1) + yy*iWidthExtSIMD + 1] = tempblock[(iWidthExtSIMD<<1) + yy*iWidthExtSIMD + 2];
      }
    }
    if(!rightAltAvailable)
    {
      for(int yy = 0; yy<uiHeight; yy++)
      {
        tempblock[(iWidthExtSIMD<<1) + uiWidthExt-2 + yy*iWidthExtSIMD] = tempblock[(iWidthExtSIMD<<1) + uiWidthExt-2 + yy*iWidthExtSIMD - 1];
        tempblock[(iWidthExtSIMD<<1) + uiWidthExt-1 + yy*iWidthExtSIMD] = tempblock[(iWidthExtSIMD<<1) + uiWidthExt-2 + yy*iWidthExtSIMD - 1];
      }
    }
    
    // All sides are available, easy to just copy corners also.
    if(topAltAvailable && leftAltAvailable)
    {
      tempblock[0] = *(srcPtr - 2*srcStride -2);                                                // a     top left corner
      tempblock[1] = *(srcPtr - 2*srcStride -1);                                                // b     a b|x x
      tempblock[iWidthExtSIMD + 0] = *(srcPtr - srcStride -2);                                  // c     c d|x x
      tempblock[iWidthExtSIMD + 1] = *(srcPtr - srcStride -1);                                  // d     -------
    }
    else
    {
      tempblock[0] = tempblock[iWidthExtSIMD*2 + 2];                                            // extend top left
      tempblock[1] = tempblock[iWidthExtSIMD*2 + 2];                                            // extend top left
      tempblock[iWidthExtSIMD + 0] = tempblock[iWidthExtSIMD*2 + 2];                            // extend top left
      tempblock[iWidthExtSIMD + 1] = tempblock[iWidthExtSIMD*2 + 2];                            // extend top left
    }
    
    if(topAltAvailable && rightAltAvailable)
    {
      tempblock[iWidthExtSIMD - 2] = *(srcPtr - 2*srcStride + uiWidth);                         // a
      tempblock[iWidthExtSIMD - 1] = *(srcPtr - 2*srcStride + uiWidth + 1);                     // b
      tempblock[iWidthExtSIMD + uiWidthExt - 2] = *(srcPtr - srcStride + uiWidth);              // c
      tempblock[iWidthExtSIMD + uiWidthExt - 1] = *(srcPtr - srcStride + uiWidth + 1);          // d
    }
    else
    {
      tempblock[iWidthExtSIMD - 2] = tempblock[iWidthExtSIMD*2 + uiWidthExt - 3];               // extend top right
      tempblock[iWidthExtSIMD - 1] = tempblock[iWidthExtSIMD*2 + uiWidthExt - 3];               // extend top right
      tempblock[iWidthExtSIMD + uiWidthExt - 2] = tempblock[iWidthExtSIMD*2 + uiWidthExt - 3];  // extend top right
      tempblock[iWidthExtSIMD + uiWidthExt - 1] = tempblock[iWidthExtSIMD*2 + uiWidthExt - 3];  // extend top right
    }
    
    if(bottomAltAvailable && leftAltAvailable)
    {
      tempblock[iWidthExtSIMD*(uiHeightExt-2) + 0] = *(srcPtr + uiHeight*srcStride -2);          // a
      tempblock[iWidthExtSIMD*(uiHeightExt-2) + 1] = *(srcPtr + uiHeight*srcStride -1);          // b
      tempblock[iWidthExtSIMD*(uiHeightExt-1) + 0] = *(srcPtr + (uiHeight+1)*srcStride -2);      // c
      tempblock[iWidthExtSIMD*(uiHeightExt-1) + 1] = *(srcPtr + (uiHeight+1)*srcStride -1);      // d
    }
    else
    {
      tempblock[iWidthExtSIMD*(uiHeightExt-2) + 0] = tempblock[iWidthExtSIMD*(uiHeightExt-3) + 2];  // bot avail: mirror left/right
      tempblock[iWidthExtSIMD*(uiHeightExt-2) + 1] = tempblock[iWidthExtSIMD*(uiHeightExt-3) + 2];  // bot avail: mirror left/right
      tempblock[iWidthExtSIMD*(uiHeightExt-1) + 0] = tempblock[iWidthExtSIMD*(uiHeightExt-3) + 2];  // bot avail: mirror left/right
      tempblock[iWidthExtSIMD*(uiHeightExt-1) + 1] = tempblock[iWidthExtSIMD*(uiHeightExt-3) + 2];  // bot avail: mirror left/right
    }
    
    if(bottomAltAvailable && rightAltAvailable)
    {
      tempblock[iWidthExtSIMD*(uiHeightExt-2) + uiWidthExt - 2] = *(srcPtr + uiHeight*srcStride + uiWidth);                // a
      tempblock[iWidthExtSIMD*(uiHeightExt-2) + uiWidthExt - 1] = *(srcPtr + uiHeight*srcStride + uiWidth + 1);            // b
      tempblock[iWidthExtSIMD*(uiHeightExt-1) + uiWidthExt - 2] = *(srcPtr + (uiHeight+1)*srcStride + uiWidth);            // c
      tempblock[iWidthExtSIMD*(uiHeightExt-1) + uiWidthExt - 1] = *(srcPtr + (uiHeight+1)*srcStride + uiWidth + 1);        // d
    }
    else
    {
      tempblock[iWidthExtSIMD*(uiHeightExt-2) + uiWidthExt - 2] = tempblock[iWidthExtSIMD*(uiHeightExt-3) + uiWidthExt - 3];
      tempblock[iWidthExtSIMD*(uiHeightExt-2) + uiWidthExt - 1] = tempblock[iWidthExtSIMD*(uiHeightExt-3) + uiWidthExt - 3];
      tempblock[iWidthExtSIMD*(uiHeightExt-1) + uiWidthExt - 2] = tempblock[iWidthExtSIMD*(uiHeightExt-3) + uiWidthExt - 3];
      tempblock[iWidthExtSIMD*(uiHeightExt-1) + uiWidthExt - 1] = tempblock[iWidthExtSIMD*(uiHeightExt-3) + uiWidthExt - 3];
    }
  }
  
  m_bilateralFilterDiamond5x5(uiWidth, uiHeight, tempblock, tempblockFiltered, clpRng, recPtr, recStride, iWidthExtSIMD, bfac, bif_round_add, bif_round_shift, false, LUTrowPtr );
}

void BilateralFilter::clipNotBilaterallyFilteredBlocks(const CPelUnitBuf& src, PelUnitBuf& rec, const ClpRng& clpRng, TransformUnit & currTU)
{
  PelUnitBuf myRecBuf = currTU.cs->getRecoBuf(currTU);
  if(myRecBuf.bufs[0].width > 1)
  {
    // new result = old result (which is SAO-treated already) + diff due to bilateral filtering
    myRecBuf.bufs[0].copyClip(myRecBuf.bufs[0], clpRng);
  }
  else
  {
    CompArea &compArea = currTU.block(COMPONENT_Y);
    const unsigned uiHeight = compArea.height;
    int recStride = rec.get(COMPONENT_Y).stride;
    Pel *recPtr = rec.get(COMPONENT_Y).bufAt(compArea);
    
    for(uint32_t yy = 0; yy < uiHeight; yy++)
    {
      // new result = old result (which is SAO-treated already) + diff due to bilateral filtering
      recPtr[0] = ClipPel<int>(recPtr[0], clpRng);
      recPtr += recStride;
    }
  }
}

static double getDist(PelBuf &RecBuf, PelBuf &OrigBuf)
{
  double Dist = 0;
  assert((RecBuf.height == OrigBuf.height));
  assert((RecBuf.width == OrigBuf.width));
  
  int recStride = RecBuf.stride;
  const Pel *recPtr = RecBuf.buf;
  
  int orgStride = OrigBuf.stride;
  const Pel *orgPtr = OrigBuf.buf;
  
  for(int yy = 0; yy<RecBuf.height; yy++)
  {
    for(int xx = 0; xx<RecBuf.width; xx++)
    {
      int diff = (*recPtr++)-(*orgPtr++);
      Dist = Dist + ((double)diff)*((double)diff);
    }
    recPtr += recStride-RecBuf.width;
    orgPtr += orgStride-RecBuf.width;
  }
  return Dist;
}

void copyBack(PelBuf &srcBuf, PelBuf &dstBuf)
{
  assert((srcBuf.height == dstBuf.height));
  assert((srcBuf.width == dstBuf.width));
  
  int srcStride = srcBuf.stride;
  const Pel *srcPtr = srcBuf.buf;
  
  int dstStride = dstBuf.stride;
  Pel *dstPtr = dstBuf.buf;
  
  for(int yy = 0; yy<srcBuf.height; yy++)
  {
    for(int xx = 0; xx<srcBuf.width; xx++)
    {
      (*dstPtr++) = (*srcPtr++);
    }
    srcPtr += srcStride-srcBuf.width;
    dstPtr += dstStride-srcBuf.width;
  }
}

void BilateralFilter::bilateralFilterPicRDOperCTU(CodingStructure& cs, PelUnitBuf& src, BIFCabacEst* BifCABACEstimator)
{
  // We must have already copied recobuf into src before running this
  // such as src.copyFrom(rec);
  
  const PreCalcValues& pcv = *cs.pcv;

  PelUnitBuf rec = cs.getRecoBuf();

  double frameMSEnoBIF = 0;
  double frameMSEallBIF = 0;
  double frameMSEswitchBIF = 0;
  BifParams& bifParams = cs.picture->getBifParam();
  int ctuIdx = 0;
  
  for (int y = 0; y < pcv.heightInCtus; y++)
  {
    for (int x = 0; x < pcv.widthInCtus; x++)
    {
      UnitArea ctuArea(pcv.chrFormat, Area(x << pcv.maxCUWidthLog2, y << pcv.maxCUHeightLog2, pcv.maxCUWidth, pcv.maxCUWidth));

      ctuArea = clipArea(ctuArea, *cs.slice->getPic());
      PelBuf piOrg = cs.getOrgBuf(ctuArea).Y();
      PelBuf piSrc = src.subBuf(ctuArea).Y();
      double MSEnoBIF = getDist(piSrc, piOrg);

      for (auto &currCU : cs.traverseCUs(CS::getArea(cs, ctuArea, CH_L), CH_L))
      {
        for (auto &currTU : CU::traverseTUs(currCU))
        {
          bool isInter = (currCU.predMode == MODE_INTER) ? true : false;
          // We should ideally also check the CTU-BIF-flag here. However, given that this function
          // is only called by the encoder, and the encoder always has CTU-BIF-flag on, there is no
          // need to check.
         
          if ((TU::getCbf(currTU, COMPONENT_Y) || isInter == false) && (currTU.cu->qp > 17) && (128 > std::max(currTU.lumaSize().width, currTU.lumaSize().height)) && ((isInter == false) || (32 > std::min(currTU.lumaSize().width, currTU.lumaSize().height))))
          {
            bilateralFilterDiamond5x5(src, rec, currTU.cu->qp, cs.slice->clpRng(COMPONENT_Y), currTU);
          }
        }
      }

      PelBuf piRec = rec.subBuf(ctuArea).Y();
      double MSEafterBIF = getDist(piRec, piOrg);
      frameMSEnoBIF += MSEnoBIF;
      frameMSEallBIF += MSEafterBIF;
      if(MSEnoBIF<MSEafterBIF)
      {
        frameMSEswitchBIF += MSEnoBIF;
        bifParams.ctuOn[ctuIdx] = 0;
      }
      else
      {
        frameMSEswitchBIF += MSEafterBIF;
        bifParams.ctuOn[ctuIdx] = 1;

      }
      ctuIdx++;

    }
  }

  double lambda = cs.picture->slices[0]->getLambdas()[0];
  double costAllCTUsBIF  = frameMSEallBIF + lambda * 1;      // To turn everything on, only slice_bif_all_ctb_enabled_flag = 1, so one bit.
  double costNoCTUsBIF = frameMSEnoBIF + lambda * 2;         // To turn everything off, slice_bif_all_ctb_enabled = 0 && slice_bif_enabled_flag = 0, so two bits.
  double costSwitchCTUsBIF;
  
  // Does CABAC estimation instead
  const double FracBitsScale = 1.0 / double(1 << SCALE_BITS);
  bifParams.frmOn = 1;
  bifParams.allCtuOn = 0;
  double ctuSwitchBits = FracBitsScale*BifCABACEstimator->getBits(*cs.slice, bifParams);
  costSwitchCTUsBIF  = frameMSEswitchBIF  + lambda * ctuSwitchBits;
 
  double bestCost = MAX_DOUBLE;
  if (costAllCTUsBIF < bestCost)
  {
    // If everything should be BIF-filtered, we do not need to change any of the samples,
    // since they are already filtered.
    bestCost = costAllCTUsBIF;
    bifParams.frmOn = 1;
    bifParams.allCtuOn = 1;
  }
  if (costSwitchCTUsBIF < bestCost)
  {
    bestCost = costSwitchCTUsBIF;
    bifParams.frmOn       = 1;
    bifParams.allCtuOn  = 0;

    // If only some CTUs should be BIF-filtered, we need to restore the ones
    // that should not be filtered. This test must be done before the above one
    // since it is partly destroying our filtered data.
    ctuIdx = 0;
    for (int y = 0; y < pcv.heightInCtus; y++)
    {
      for (int x = 0; x < pcv.widthInCtus; x++)
      {
        UnitArea ctuArea(pcv.chrFormat, Area(x << pcv.maxCUWidthLog2, y << pcv.maxCUHeightLog2, pcv.maxCUWidth, pcv.maxCUWidth));
        ctuArea = clipArea(ctuArea, *cs.slice->getPic());
        PelBuf piRec = rec.subBuf(ctuArea).Y();
        PelBuf piSrc = src.subBuf(ctuArea).Y();
        if(bifParams.ctuOn[ctuIdx++]==0)
          copyBack(piSrc, piRec); // Copy unfiltered samples back to rec
      }
    }
  }
  if (costNoCTUsBIF < bestCost)
  {
    bestCost = costNoCTUsBIF;
    bifParams.frmOn       = 0;
    bifParams.allCtuOn  = 0;
    // If no CTUs should be BIF-filtered, we need to restore all CTUs.
    // Note that this test must be done last since it is destroying all
    // of our filtered data.
    rec.copyFrom(src);
  }

  if( bifParams.frmOn == 0 )
  {
    std::fill( bifParams.ctuOn.begin(), bifParams.ctuOn.end(), 0 );
  }
  else if( bifParams.allCtuOn )
  {
    std::fill( bifParams.ctuOn.begin(), bifParams.ctuOn.end(), 1 );
  }
}

#endif
