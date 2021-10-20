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

/** \file     EncAdaptiveLoopFilter.cpp
 \brief    estimation part of adaptive loop filter class
 */
#include "EncAdaptiveLoopFilter.h"

#include "CommonLib/Picture.h"
#include "CommonLib/CodingStructure.h"

#define AlfCtx(c) SubCtx( Ctx::Alf, c)
std::vector<double> EncAdaptiveLoopFilter::m_lumaLevelToWeightPLUT;

#include <algorithm>

#if MAX_NUM_CC_ALF_FILTERS>1
struct FilterIdxCount
{
  uint64_t count;
  uint8_t filterIdx;
};

bool compareCounts(FilterIdxCount a, FilterIdxCount b) { return a.count > b.count; }
#endif

void AlfCovariance::getClipMax(const AlfFilterShape& alfShape, int *clip_max) const
{
  for( int k = 0; k < numCoeff-1; ++k )
  {
    clip_max[k] = 0;

    bool inc = true;
    while( inc && clip_max[k]+1 < numBins && y[clip_max[k]+1][k] == y[clip_max[k]][k] )
    {
      for( int l = 0; inc && l < numCoeff; ++l )
      {
        if( E[clip_max[k]][0][k][l] != E[clip_max[k]+1][0][k][l] )
        {
          inc = false;
        }
      }
      if( inc )
      {
        ++clip_max[k];
      }
    }
  }
  clip_max[numCoeff-1] = 0;
}

void AlfCovariance::reduceClipCost(const AlfFilterShape& alfShape, int *clip) const
{
  for( int k = 0; k < numCoeff-1; ++k )
  {
    bool dec = true;
    while( dec && clip[k] > 0 && y[clip[k]-1][k] == y[clip[k]][k] )
    {
      for( int l = 0; dec && l < numCoeff; ++l )
      {
        if( E[clip[k]][clip[l]][k][l] != E[clip[k]-1][clip[l]][k][l] )
        {
          dec = false;
        }
      }
      if( dec )
      {
        --clip[k];
      }
    }
  }
}

double AlfCovariance::optimizeFilter(const AlfFilterShape& alfShape, int* clip, double *f, bool optimize_clip) const
{
  const int size = alfShape.numCoeff;
  int clip_max[MAX_NUM_ALF_LUMA_COEFF];
  double fBest[MAX_NUM_ALF_LUMA_COEFF];
  int copySize = sizeof( double ) * size;
  double err_best, err_last;

  TE kE;
  Ty ky;

  if( optimize_clip )
  {
    // Start by looking for min clipping that has no impact => max_clipping
    getClipMax(alfShape, clip_max);
    for (int k=0; k<size; ++k)
    {
      clip[k] = std::max(clip_max[k], clip[k]);
      clip[k] = std::min(clip[k], numBins-1);
    }
  }

  setEyFromClip( clip, kE, ky, size );

  gnsSolveByChol( kE, ky, f, size );
  err_best = calculateError( clip, f, size );
  memcpy( fBest, f, copySize );
  int step = optimize_clip ? (numBins+1)/2 : 0;

  while( step > 0 )
  {
    double err_min = err_best;
    int idx_min = -1;
    int inc_min = 0;

    for( int k = 0; k < size-1; ++k )
    {
      if( clip[k] - step >= clip_max[k] )
      {
        clip[k] -= step;
        ky[k] = y[clip[k]][k];
        for( int l = 0; l < size; l++ )
        {
          kE[k][l] = E[clip[k]][clip[l]][k][l];
          kE[l][k] = E[clip[l]][clip[k]][l][k];
        }

        gnsSolveByChol( kE, ky, f, size );
        err_last = calculateError( clip, f, size );

        if( err_last < err_min )
        {
          err_min = err_last;
          idx_min = k;
          inc_min = -step;
          memcpy(fBest, f, copySize);
        }
        clip[k] += step;
      }
      if( clip[k] + step < numBins )
      {
        clip[k] += step;
        ky[k] = y[clip[k]][k];
        for( int l = 0; l < size; l++ )
        {
          kE[k][l] = E[clip[k]][clip[l]][k][l];
          kE[l][k] = E[clip[l]][clip[k]][l][k];
        }

        gnsSolveByChol( kE, ky, f, size );
        err_last = calculateError( clip, f, size );

        if( err_last < err_min )
        {
          err_min = err_last;
          idx_min = k;
          inc_min = step;
          memcpy(fBest, f, copySize);
        }
        clip[k] -= step;

      }
      ky[k] = y[clip[k]][k];
      for( int l = 0; l < size; l++ )
      {
        kE[k][l] = E[clip[k]][clip[l]][k][l];
        kE[l][k] = E[clip[l]][clip[k]][l][k];
      }
    }

    if( idx_min >= 0 )
    {
      err_best = err_min;
      clip[idx_min] += inc_min;
      ky[idx_min] = y[clip[idx_min]][idx_min];
      for( int l = 0; l < size; l++ )
      {
        kE[idx_min][l] = E[clip[idx_min]][clip[l]][idx_min][l];
        kE[l][idx_min] = E[clip[l]][clip[idx_min]][l][idx_min];
      }
    }
    else
    {
      --step;
    }
  }

  if (optimize_clip)
  {
    // test all max
    for( int k = 0; k < size-1; ++k )
    {
      clip_max[k] = 0;
    }
    TE kE_max;
    Ty ky_max;
    setEyFromClip( clip_max, kE_max, ky_max, size );

    gnsSolveByChol( kE_max, ky_max, f, size );
    err_last = calculateError( clip_max, f, size );
    if( err_last < err_best )
    {
      err_best = err_last;
      for (int k=0; k<size; ++k)
      {
        clip[k] = clip_max[k];
      }
    }
    else
    {
      // update clip to reduce coding cost
      reduceClipCost(alfShape, clip);
      memcpy( f, fBest, copySize );
    }
  }

  return err_best;
}

double AlfCovariance::calcErrorForCoeffs( const int *clip, const int *coeff, const int numCoeff, const int bitDepth ) const
{
  double factor = 1 << ( bitDepth - 1 );
  double error = 0;

  for( int i = 0; i < numCoeff; i++ )   //diagonal
  {
    double sum = 0;
    for( int j = i + 1; j < numCoeff; j++ )
    {
      // E[j][i] = E[i][j], sum will be multiplied by 2 later
      sum += E[clip[i]][clip[j]][i][j] * coeff[j];
    }
    error += ( ( E[clip[i]][clip[i]][i][i] * coeff[i] + sum * 2 ) / factor - 2 * y[clip[i]][i] ) * coeff[i];
  }

  return error / factor;
}

double AlfCovariance::calcErrorForCcAlfCoeffs(const int16_t *coeff, const int numCoeff, const int bitDepth) const
{
  double factor = 1 << (bitDepth - 1);
  double error = 0;

  for (int i = 0; i < numCoeff; i++)   // diagonal
  {
    double sum = 0;
    for (int j = i + 1; j < numCoeff; j++)
    {
      // E[j][i] = E[i][j], sum will be multiplied by 2 later
      sum += E[0][0][i][j] * coeff[j];
    }
    error += ((E[0][0][i][i] * coeff[i] + sum * 2) / factor - 2 * y[0][i]) * coeff[i];
  }

  return error / factor;
}

double AlfCovariance::calculateError( const int *clip, const double *coeff, const int numCoeff ) const
{
  double sum = 0;
  for( int i = 0; i < numCoeff; i++ )
  {
    sum += coeff[i] * y[clip[i]][i];
  }

  return pixAcc - sum;
}

double AlfCovariance::calculateError( const int *clip ) const
{
  Ty c;

  return optimizeFilter( clip, c, numCoeff );
}
//********************************
// Cholesky decomposition
//********************************

#define ROUND(a)  (((a) < 0)? (int)((a) - 0.5) : (int)((a) + 0.5))
#define REG              0.0001
#define REG_SQR          0.0000001

//Find filter coeff related
int AlfCovariance::gnsCholeskyDec( TE inpMatr, TE outMatr, int numEq ) const
{
  for( int i = 0; i < numEq; i++ )
  {
    double* inputM  = inpMatr[i];
    double* outputM = outMatr[i];
    double scale    = inputM[i];

    for( int k = i - 1; k >= 0; k-- )
    {
      scale -= outMatr[k][i] * outMatr[k][i];
    }

    if( scale <= REG_SQR ) // inpMatr is singular
    {
      return 0;
    }

    outputM[i] = sqrt( scale );
    double tmp = 1 / outputM[i];

    for( int j = i + 1; j < numEq; j++ )
    {
      scale = inputM[j];
      for( int k = i - 1; k >= 0; k-- )
      {
        scale -= outMatr[k][j] * outMatr[k][i];
      }

      outputM[j] = scale * tmp; // Upper triangular
      //outMatr[j][i] = 0.0;      // Lower triangular part
    }
  }

  return 1; // Signal that Cholesky factorization is successfully performed
}

void AlfCovariance::gnsTransposeBacksubstitution( TE U, double* rhs, double* x, int order ) const
{
  /* Backsubstitution starts */
  x[0] = rhs[0] / U[0][0];               /* First row of U'                   */
  for( int i = 1; i < order; i++ )
  {         /* For the rows 1..order-1           */

    double sum = 0; //Holds backsubstitution from already handled rows

    for( int j = 0; j < i; j++ ) /* Backsubst already solved unknowns */
    {
      sum += x[j] * U[j][i];
    }

    x[i] = ( rhs[i] - sum ) / U[i][i];       /* i'th component of solution vect.  */
  }
}

void AlfCovariance::gnsBacksubstitution( TE R, double* z, int size, double* A ) const
{
  size--;
  A[size] = z[size] / R[size][size];

  for( int i = size - 1; i >= 0; i-- )
  {
    double sum = 0;

    for( int j = i + 1; j <= size; j++ )
    {
      sum += R[i][j] * A[j];
    }

    A[i] = ( z[i] - sum ) / R[i][i];
  }
}

int AlfCovariance::gnsSolveByChol( const int *clip, double *x, int numEq ) const
{
  TE LHS;
  Ty rhs;

  setEyFromClip( clip, LHS, rhs, numEq );
  return gnsSolveByChol( LHS, rhs, x, numEq );
}

int AlfCovariance::gnsSolveByChol( TE LHS, double* rhs, double *x, int numEq ) const
{
  Ty aux;     /* Auxiliary vector */
  TE U;    /* Upper triangular Cholesky factor of LHS */

  int res = 1;  // Signal that Cholesky factorization is successfully performed

                /* The equation to be solved is LHSx = rhs */

                /* Compute upper triangular U such that U'*U = LHS */

  if( gnsCholeskyDec( LHS, U, numEq ) ) /* If Cholesky decomposition has been successful */
  {
    /* Now, the equation is  U'*U*x = rhs, where U is upper triangular
    * Solve U'*aux = rhs for aux
    */
    gnsTransposeBacksubstitution( U, rhs, aux, numEq );

    /* The equation is now U*x = aux, solve it for x (new motion coefficients) */
    gnsBacksubstitution( U, aux, numEq, x );

  }
  else /* LHS was singular */
  {
    res = 0;

    /* Regularize LHS */
    for( int i = 0; i < numEq; i++ )
    {
      LHS[i][i] += REG;
    }

    /* Compute upper triangular U such that U'*U = regularized LHS */
    res = gnsCholeskyDec( LHS, U, numEq );

    if( !res )
    {
      std::memset( x, 0, sizeof( double )*numEq );
      return 0;
    }

    /* Solve  U'*aux = rhs for aux */
    gnsTransposeBacksubstitution( U, rhs, aux, numEq );

    /* Solve U*x = aux for x */
    gnsBacksubstitution( U, aux, numEq, x );
  }
  return res;
}
//////////////////////////////////////////////////////////////////////////////////////////

EncAdaptiveLoopFilter::EncAdaptiveLoopFilter( int& apsIdStart )
  : m_CABACEstimator( nullptr )
  , m_apsIdStart( apsIdStart )
{
  for( int i = 0; i < MAX_NUM_COMPONENT; i++ )
  {
    m_alfCovariance[i] = nullptr;
  }
  for( int i = 0; i < MAX_NUM_CHANNEL_TYPE; i++ )
  {
    m_alfCovarianceFrame[i] = nullptr;
  }
  m_filterCoeffSet = nullptr;
  m_filterClippSet = nullptr;
  m_diffFilterCoeff = nullptr;

  m_alfWSSD = 0;

  m_alfCovarianceCcAlf[0] = nullptr;
  m_alfCovarianceCcAlf[1] = nullptr;
  m_alfCovarianceFrameCcAlf[0] = nullptr;
  m_alfCovarianceFrameCcAlf[1] = nullptr;

}

void EncAdaptiveLoopFilter::create( const EncCfg* encCfg, const int picWidth, const int picHeight, const ChromaFormat chromaFormatIDC, const int maxCUWidth, const int maxCUHeight, const int maxCUDepth, const int inputBitDepth[MAX_NUM_CHANNEL_TYPE], const int internalBitDepth[MAX_NUM_CHANNEL_TYPE] )
{
  AdaptiveLoopFilter::create( picWidth, picHeight, chromaFormatIDC, maxCUWidth, maxCUHeight, maxCUDepth, inputBitDepth );
  CHECK( encCfg == nullptr, "encCfg must not be null" );
  m_encCfg = encCfg;


  for( int channelIdx = 0; channelIdx < MAX_NUM_CHANNEL_TYPE; channelIdx++ )
  {
    ChannelType chType = (ChannelType)channelIdx;
    int numClasses = channelIdx ? 1 : MAX_NUM_ALF_CLASSES;
    m_alfCovarianceFrame[chType] = new AlfCovariance*[m_filterShapes[chType].size()];

    for( int i = 0; i != m_filterShapes[chType].size(); i++ )
    {
#if ALF_IMPROVEMENT
      m_alfCovarianceFrame[chType][i] = nullptr;
      if (m_filterTypeTest[chType][m_filterShapes[chType][i].filterType] == false)
      {
        continue;
      }
#endif
      m_alfCovarianceFrame[chType][i] = new AlfCovariance[numClasses];

      for( int k = 0; k < numClasses; k++ )
      {
        m_alfCovarianceFrame[chType][i][k].create( m_filterShapes[chType][i].numCoeff );
      }
    }
  }

  for( int compIdx = 0; compIdx < MAX_NUM_COMPONENT; compIdx++ )
  {
    m_ctuEnableFlagTmp[compIdx] = new uint8_t[m_numCTUsInPic];
    m_ctuEnableFlagTmp2[compIdx] = new uint8_t[m_numCTUsInPic];
#if !ALF_IMPROVEMENT
    if( isLuma( ComponentID(compIdx) ) )
    {
      m_ctuAlternativeTmp[compIdx] = nullptr;
    }
    else
#endif
    {
      m_ctuAlternativeTmp[compIdx] = new uint8_t[m_numCTUsInPic];
      std::fill_n( m_ctuAlternativeTmp[compIdx], m_numCTUsInPic, 0 );
    }
    ChannelType chType = toChannelType( ComponentID( compIdx ) );
    int numClasses = compIdx ? 1 : MAX_NUM_ALF_CLASSES;
#if ALF_IMPROVEMENT
    m_alfCovariance[compIdx] = new AlfCovariance***[m_filterShapes[chType].size()];
    for( int i = 0; i != m_filterShapes[chType].size(); i++ )
    {
      m_alfCovariance[compIdx][i] = nullptr;
      if( m_filterTypeTest[chType][m_filterShapes[chType][i].filterType] == false )
      {
        continue;
      }
      m_alfCovariance[compIdx][i] = new AlfCovariance**[m_numCTUsInPic];
      int numFixedFilterSets = ( m_filterShapes[chType][i].filterType == ALF_FILTER_EXT || m_filterShapes[chType][i].filterType == ALF_FILTER_9_EXT ) ? 2 : 1;
      for( int j = 0; j < m_numCTUsInPic; j++ )
      {
        m_alfCovariance[compIdx][i][j] = new AlfCovariance*[numFixedFilterSets];
        for( int fixedFilterSetIdx = 0; fixedFilterSetIdx < numFixedFilterSets; fixedFilterSetIdx++)
        {
          m_alfCovariance[compIdx][i][j][fixedFilterSetIdx] = new AlfCovariance[numClasses];
          for (int k = 0; k < numClasses; k++)
          {
            m_alfCovariance[compIdx][i][j][fixedFilterSetIdx][k].create(m_filterShapes[chType][i].numCoeff);
          }
        }             
      }
    }
#else
    m_alfCovariance[compIdx] = new AlfCovariance**[m_filterShapes[chType].size()];

    for( int i = 0; i != m_filterShapes[chType].size(); i++ )
    {
      m_alfCovariance[compIdx][i] = new AlfCovariance*[m_numCTUsInPic];
      for( int j = 0; j < m_numCTUsInPic; j++ )
      {
        m_alfCovariance[compIdx][i][j] = new AlfCovariance[numClasses];
        for( int k = 0; k < numClasses; k++ )
        {
          m_alfCovariance[compIdx][i][j][k].create( m_filterShapes[chType][i].numCoeff );
        }
      }
    }
#endif
  }

  for( int i = 0; i != m_filterShapes[COMPONENT_Y].size(); i++ )
  {
    for (int j = 0; j <= MAX_NUM_ALF_CLASSES + 1; j++)
    {
      m_alfCovarianceMerged[i][j].create( m_filterShapes[COMPONENT_Y][i].numCoeff );
    }
  }
#if ALF_IMPROVEMENT
  m_filterCoeffSet = new int*[std::max(MAX_NUM_ALF_ALTERNATIVES_LUMA * MAX_NUM_ALF_CLASSES, MAX_NUM_ALF_ALTERNATIVES_CHROMA)];
  m_filterClippSet = new int*[std::max(MAX_NUM_ALF_ALTERNATIVES_LUMA * MAX_NUM_ALF_CLASSES, MAX_NUM_ALF_ALTERNATIVES_CHROMA)];
#else
  m_filterCoeffSet = new int*[std::max(MAX_NUM_ALF_CLASSES, MAX_NUM_ALF_ALTERNATIVES_CHROMA)];
  m_filterClippSet = new int*[std::max(MAX_NUM_ALF_CLASSES, MAX_NUM_ALF_ALTERNATIVES_CHROMA)];
#endif
  m_diffFilterCoeff = new int*[MAX_NUM_ALF_CLASSES];
  for( int i = 0; i < MAX_NUM_ALF_CLASSES; i++ )
  {
#if ALF_IMPROVEMENT
    for (int j = 0; j < MAX_NUM_ALF_ALTERNATIVES_LUMA; j++)
    {
      m_filterCoeffSet[i * MAX_NUM_ALF_ALTERNATIVES_LUMA + j] = new int[MAX_NUM_ALF_LUMA_COEFF];
      m_filterClippSet[i * MAX_NUM_ALF_ALTERNATIVES_LUMA + j] = new int[MAX_NUM_ALF_LUMA_COEFF];
    }
#else
    m_filterCoeffSet[i] = new int[MAX_NUM_ALF_LUMA_COEFF];
    m_filterClippSet[i] = new int[MAX_NUM_ALF_LUMA_COEFF];
#endif
    m_diffFilterCoeff[i] = new int[MAX_NUM_ALF_LUMA_COEFF];
  }
#if !ALF_IMPROVEMENT
  for( int i = 0; i < NUM_FIXED_FILTER_SETS; i++ )
  {
    m_ctbDistortionFixedFilter[i] = new double[m_numCTUsInPic];
  }
#else
  for (int i = 0; i < NUM_FIXED_FILTER_SETS; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      m_ctbDistortionFixedFilter[i][j] = new double[m_numCTUsInPic];
    }
  }
#endif
  for( int i = 0; i < ALF_CTB_MAX_NUM_APS; i++ )
  {
#if ALF_IMPROVEMENT
    for( int j = 0; j< MAX_NUM_ALF_ALTERNATIVES_LUMA; j++ )
    {
      for( int k = 0; k < 2; k++ )
      {
        m_distCtbApsLuma[i][j][k] = new double[m_numCTUsInPic];
      }
    }
#else
    m_distCtbApsLuma[i] = new double[m_numCTUsInPic];
#endif
  }
#if ALF_IMPROVEMENT
  for( int i = 0; i < MAX_NUM_ALF_ALTERNATIVES_LUMA; i++ )
  {
    for (int k = 0; k < 2; k++)
    {
      m_distCtbLumaNewFilt[i][k] = new double[m_numCTUsInPic];
    }
  }
#else
  m_distCtbLumaNewFilt = new double[m_numCTUsInPic];
#endif

  for (int comp = 0; comp < MAX_NUM_COMPONENT; comp++)
  {
    m_ctbDistortionUnfilter[comp] = new double[m_numCTUsInPic];
  }
  m_alfCtbFilterSetIndexTmp.resize(m_numCTUsInPic);
#if !ALF_IMPROVEMENT
  memset(m_clipDefaultEnc, 0, sizeof(m_clipDefaultEnc));
#endif
  m_apsIdCcAlfStart[0] = (int) MAX_NUM_APS;
  m_apsIdCcAlfStart[1] = (int) MAX_NUM_APS;
  for( int compIdx = 1; compIdx < MAX_NUM_COMPONENT; compIdx++ )
  {
    int numFilters = MAX_NUM_CC_ALF_FILTERS;
    m_alfCovarianceCcAlf[compIdx-1] = new AlfCovariance**[m_filterShapesCcAlf[compIdx-1].size()];
    m_alfCovarianceFrameCcAlf[compIdx-1] = new AlfCovariance*[m_filterShapesCcAlf[compIdx-1].size()];
    for( int i = 0; i != m_filterShapesCcAlf[compIdx-1].size(); i++ )
    {
      m_alfCovarianceFrameCcAlf[compIdx - 1][i] = new AlfCovariance[numFilters];
      for (int k = 0; k < numFilters; k++)
      {
        m_alfCovarianceFrameCcAlf[compIdx - 1][i][k].create(m_filterShapesCcAlf[compIdx - 1][i].numCoeff);
      }

      m_alfCovarianceCcAlf[compIdx - 1][i] = new AlfCovariance *[numFilters];
      for (int j = 0; j < numFilters; j++)
      {
        m_alfCovarianceCcAlf[compIdx - 1][i][j] = new AlfCovariance[m_numCTUsInPic];
        for (int k = 0; k < m_numCTUsInPic; k++)
        {
          m_alfCovarianceCcAlf[compIdx - 1][i][j][k].create(m_filterShapesCcAlf[compIdx - 1][i].numCoeff);
        }
      }
    }
  }
  m_trainingCovControl   = new uint8_t[m_numCTUsInPic];
  for ( int i = 0; i < MAX_NUM_CC_ALF_FILTERS; i++ )
  {
    m_trainingDistortion[i] = new uint64_t[m_numCTUsInPic];
  }
  m_filterControl         = new uint8_t[m_numCTUsInPic];
  m_bestFilterControl     = new uint8_t[m_numCTUsInPic];
  uint32_t area           = (picWidth >> getComponentScaleX(COMPONENT_Cb,chromaFormatIDC))*(picHeight >> getComponentScaleY(COMPONENT_Cb,chromaFormatIDC));
  m_bufOrigin             = ( Pel* ) xMalloc( Pel, area );
  m_buf                   = new PelBuf( m_bufOrigin, picWidth >> getComponentScaleX(COMPONENT_Cb,chromaFormatIDC), picWidth >> getComponentScaleX(COMPONENT_Cb,chromaFormatIDC), picHeight >> getComponentScaleY(COMPONENT_Cb,chromaFormatIDC) );
  m_lumaSwingGreaterThanThresholdCount = new uint64_t[m_numCTUsInPic];
  m_chromaSampleCountNearMidPoint = new uint64_t[m_numCTUsInPic];
}

void EncAdaptiveLoopFilter::destroy()
{
  if (!m_created)
  {
    return;
  }
  for( int channelIdx = 0; channelIdx < MAX_NUM_CHANNEL_TYPE; channelIdx++ )
  {
    if( m_alfCovarianceFrame[channelIdx] )
    {
      ChannelType chType = (ChannelType)channelIdx;
      int numClasses = channelIdx ? 1 : MAX_NUM_ALF_CLASSES;

      for (int i = 0; i != m_filterShapes[chType].size(); i++)
      {
#if ALF_IMPROVEMENT
        if (m_alfCovarianceFrame[channelIdx][i] == nullptr)
        {
          continue;
        }
#endif
        for( int k = 0; k < numClasses; k++ )
        {
          m_alfCovarianceFrame[channelIdx][i][k].destroy();
        }
        delete[] m_alfCovarianceFrame[channelIdx][i];
        m_alfCovarianceFrame[channelIdx][i] = nullptr;
      }
      delete[] m_alfCovarianceFrame[channelIdx];
      m_alfCovarianceFrame[channelIdx] = nullptr;
    }

  }

  for( int compIdx = 0; compIdx < MAX_NUM_COMPONENT; compIdx++ )
  {
    if( m_ctuEnableFlagTmp[compIdx] )
    {
      delete[] m_ctuEnableFlagTmp[compIdx];
      m_ctuEnableFlagTmp[compIdx] = nullptr;
    }

    if( m_ctuEnableFlagTmp2[compIdx] )
    {
      delete[] m_ctuEnableFlagTmp2[compIdx];
      m_ctuEnableFlagTmp2[compIdx] = nullptr;
    }

    if( m_ctuAlternativeTmp[compIdx] )
    {
      delete[] m_ctuAlternativeTmp[compIdx];
      m_ctuAlternativeTmp[compIdx] = nullptr;
    }

    if( m_alfCovariance[compIdx] )
    {
      ChannelType chType = toChannelType( ComponentID( compIdx ) );
      int numClasses = compIdx ? 1 : MAX_NUM_ALF_CLASSES;

      for( int i = 0; i != m_filterShapes[chType].size(); i++ )
      {
#if ALF_IMPROVEMENT
        if (m_alfCovariance[compIdx][i] == nullptr)
        {
          continue;
        }
        int numFixedFilterSet = ( m_filterShapes[chType][i].filterType == ALF_FILTER_9_EXT || m_filterShapes[chType][i].filterType == ALF_FILTER_EXT ) ? 2 : 1;
#endif
        for( int j = 0; j < m_numCTUsInPic; j++ )
        {
#if ALF_IMPROVEMENT
          for( int fixedFilterSetIdx = 0; fixedFilterSetIdx < numFixedFilterSet; fixedFilterSetIdx++ )
          {
            for( int k = 0; k < numClasses; k++ )
            {
              m_alfCovariance[compIdx][i][j][fixedFilterSetIdx][k].destroy();
            }
            delete[] m_alfCovariance[compIdx][i][j][fixedFilterSetIdx];
            m_alfCovariance[compIdx][i][j][fixedFilterSetIdx] = nullptr;
          }
#else
          for( int k = 0; k < numClasses; k++ )
          {
            m_alfCovariance[compIdx][i][j][k].destroy();
          }
#endif
          delete[] m_alfCovariance[compIdx][i][j];
          m_alfCovariance[compIdx][i][j] = nullptr;

        }
        delete[] m_alfCovariance[compIdx][i];
        m_alfCovariance[compIdx][i] = nullptr;

      }
      delete[] m_alfCovariance[compIdx];
      m_alfCovariance[compIdx] = nullptr;
    }
  }

  for( int i = 0; i != m_filterShapes[COMPONENT_Y].size(); i++ )
  {
    for (int j = 0; j <= MAX_NUM_ALF_CLASSES + 1; j++)
    {
      m_alfCovarianceMerged[i][j].destroy();
    }
  }

  if( m_filterCoeffSet )
  {
#if ALF_IMPROVEMENT
    for( int i = 0; i < MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_ALTERNATIVES_LUMA; i++ )
#else
    for( int i = 0; i < MAX_NUM_ALF_CLASSES; i++ )
#endif
    {
      delete[] m_filterCoeffSet[i];
      m_filterCoeffSet[i] = nullptr;
    }
    delete[] m_filterCoeffSet;
    m_filterCoeffSet = nullptr;
  }

  if( m_filterClippSet )
  {
#if ALF_IMPROVEMENT
    for ( int i = 0; i < MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_ALTERNATIVES_LUMA; i++ )
#else
    for ( int i = 0; i < MAX_NUM_ALF_CLASSES; i++ )
#endif
    {
      delete[] m_filterClippSet[i];
      m_filterClippSet[i] = nullptr;
    }
    delete[] m_filterClippSet;
    m_filterClippSet = nullptr;
  }

  if( m_diffFilterCoeff )
  {
    for( int i = 0; i < MAX_NUM_ALF_CLASSES; i++ )
    {
      delete[] m_diffFilterCoeff[i];
      m_diffFilterCoeff[i] = nullptr;
    }
    delete[] m_diffFilterCoeff;
    m_diffFilterCoeff = nullptr;
  }
#if ALF_IMPROVEMENT
  for( int i = 0; i < NUM_FIXED_FILTER_SETS; i++ )
  {
    for( int j = 0; j < 2; j++ )
    {
      if( m_ctbDistortionFixedFilter[i][j] )
      {
        delete[] m_ctbDistortionFixedFilter[i][j];
        m_ctbDistortionFixedFilter[i][j] = nullptr;
      }
    }
  }
#else
  for( int i = 0; i < NUM_FIXED_FILTER_SETS; i++ )
  {
    if( m_ctbDistortionFixedFilter[i] )
    {
      delete[] m_ctbDistortionFixedFilter[i];
      m_ctbDistortionFixedFilter[i] = nullptr;
    }
  }
#endif
  for( int i = 0; i < ALF_CTB_MAX_NUM_APS; i++ )
  {
#if ALF_IMPROVEMENT
    for( int j = 0; j< MAX_NUM_ALF_ALTERNATIVES_LUMA; j++ )
    {
      if( m_distCtbApsLuma[i][j] )
      {
        for( int k = 0; k < 2; k++)
        {
          if( m_distCtbApsLuma[i][j][k] )
          {
            delete[] m_distCtbApsLuma[i][j][k];
            m_distCtbApsLuma[i][j][k] = nullptr;
          }
        }
      }
    }
#else
    if( m_distCtbApsLuma[i] )
    {
      delete[] m_distCtbApsLuma[i];
      m_distCtbApsLuma[i] = nullptr;
    }
#endif
  }
#if ALF_IMPROVEMENT
  for( int i = 0; i < MAX_NUM_ALF_ALTERNATIVES_LUMA; i++ )
  {
    if( m_distCtbLumaNewFilt[i] )
    {
      for( int k = 0; k < 2; k++)
      {
        if( m_distCtbLumaNewFilt[i][k] )
        {
          delete[] m_distCtbLumaNewFilt[i][k];
          m_distCtbLumaNewFilt[i][k] = nullptr;
        }
      }
    }
  }
#else
  if( m_distCtbLumaNewFilt )
  {
    delete[] m_distCtbLumaNewFilt;
    m_distCtbLumaNewFilt = nullptr;
  }
#endif

  for (int comp = 0; comp < MAX_NUM_COMPONENT; comp++)
  {
    delete[] m_ctbDistortionUnfilter[comp];
    m_ctbDistortionUnfilter[comp] = nullptr;
  }

  for (int compIdx = 1; compIdx < MAX_NUM_COMPONENT; compIdx++)
  {
    int numFilters = MAX_NUM_CC_ALF_FILTERS;
    if (m_alfCovarianceFrameCcAlf[compIdx - 1])
    {
      for (int i = 0; i != m_filterShapesCcAlf[compIdx - 1].size(); i++)
      {
        for (int k = 0; k < numFilters; k++)
        {
          m_alfCovarianceFrameCcAlf[compIdx - 1][i][k].destroy();
        }
        delete[] m_alfCovarianceFrameCcAlf[compIdx - 1][i];
      }
      delete[] m_alfCovarianceFrameCcAlf[compIdx - 1];
      m_alfCovarianceFrameCcAlf[compIdx - 1] = nullptr;
    }

    if (m_alfCovarianceCcAlf[compIdx - 1])
    {
      for (int i = 0; i != m_filterShapesCcAlf[compIdx - 1].size(); i++)
      {
        for (int j = 0; j < numFilters; j++)
        {
          for (int k = 0; k < m_numCTUsInPic; k++)
          {
            m_alfCovarianceCcAlf[compIdx - 1][i][j][k].destroy();
          }
          delete[] m_alfCovarianceCcAlf[compIdx - 1][i][j];
        }
        delete[] m_alfCovarianceCcAlf[compIdx - 1][i];
      }
      delete[] m_alfCovarianceCcAlf[compIdx - 1];
      m_alfCovarianceCcAlf[compIdx - 1] = nullptr;
    }
  }

  if (m_trainingCovControl)
  {
    delete[] m_trainingCovControl;
    m_trainingCovControl = nullptr;
  }

  for ( int i = 0; i < MAX_NUM_CC_ALF_FILTERS; i++ )
  {
    if (m_trainingDistortion[i])
    {
      delete[] m_trainingDistortion[i];
      m_trainingDistortion[i] = nullptr;
    }
  }

  if (m_filterControl)
  {
    delete[] m_filterControl;
    m_filterControl = nullptr;
  }

  if (m_bestFilterControl)
  {
    delete[] m_bestFilterControl;
    m_bestFilterControl = nullptr;
  }

  if (m_bufOrigin)
  {
    xFree(m_bufOrigin);
    m_bufOrigin = nullptr;
  }

  if (m_buf)
  {
    delete m_buf;
    m_buf = nullptr;
  }

  if (m_lumaSwingGreaterThanThresholdCount)
  {
    delete[] m_lumaSwingGreaterThanThresholdCount;
    m_lumaSwingGreaterThanThresholdCount = nullptr;
  }
  if (m_chromaSampleCountNearMidPoint)
  {
    delete[] m_chromaSampleCountNearMidPoint;
    m_chromaSampleCountNearMidPoint = nullptr;
  }

  AdaptiveLoopFilter::destroy();
}

void EncAdaptiveLoopFilter::initCABACEstimator( CABACEncoder* cabacEncoder, CtxCache* ctxCache, Slice* pcSlice
, ParameterSetMap<APS>* apsMap )
{
  m_apsMap = apsMap;
  m_CABACEstimator = cabacEncoder->getCABACEstimator( pcSlice->getSPS() );
  m_CtxCache = ctxCache;
  m_CABACEstimator->initCtxModels( *pcSlice );
  m_CABACEstimator->resetBits();
}

void EncAdaptiveLoopFilter::xSetupCcAlfAPS( CodingStructure &cs )
{
  if (m_ccAlfFilterParam.ccAlfFilterEnabled[COMPONENT_Cb - 1])
  {
    int  ccAlfCbApsId = cs.slice->getTileGroupCcAlfCbApsId();
    APS* aps = m_apsMap->getPS((cs.slice->getTileGroupCcAlfCbApsId() << NUM_APS_TYPE_LEN) + ALF_APS);
    if (aps == NULL)
    {
      aps = m_apsMap->allocatePS((ccAlfCbApsId << NUM_APS_TYPE_LEN) + ALF_APS);
      aps->setTemporalId(cs.slice->getTLayer());
    }
    aps->getCcAlfAPSParam().ccAlfFilterEnabled[COMPONENT_Cb - 1] = 1;
    aps->getCcAlfAPSParam().ccAlfFilterCount[COMPONENT_Cb - 1] = m_ccAlfFilterParam.ccAlfFilterCount[COMPONENT_Cb - 1];
    for ( int filterIdx = 0; filterIdx < MAX_NUM_CC_ALF_FILTERS; filterIdx++ )
    {
      aps->getCcAlfAPSParam().ccAlfFilterIdxEnabled[COMPONENT_Cb - 1][filterIdx] =
        m_ccAlfFilterParam.ccAlfFilterIdxEnabled[COMPONENT_Cb - 1][filterIdx];
      memcpy(aps->getCcAlfAPSParam().ccAlfCoeff[COMPONENT_Cb - 1][filterIdx],
             m_ccAlfFilterParam.ccAlfCoeff[COMPONENT_Cb - 1][filterIdx], sizeof(short) * MAX_NUM_CC_ALF_CHROMA_COEFF);
    }
    aps->setAPSId(ccAlfCbApsId);
    aps->setAPSType(ALF_APS);
    if (m_reuseApsId[COMPONENT_Cb - 1] < 0)
    {
      aps->getCcAlfAPSParam().newCcAlfFilter[COMPONENT_Cb - 1] = 1;
      m_apsMap->setChangedFlag((ccAlfCbApsId << NUM_APS_TYPE_LEN) + ALF_APS, true);
      aps->setTemporalId(cs.slice->getTLayer());
    }
    cs.slice->setTileGroupCcAlfCbEnabledFlag(true);
  }
  else
  {
    cs.slice->setTileGroupCcAlfCbEnabledFlag(false);
  }
  if (m_ccAlfFilterParam.ccAlfFilterEnabled[COMPONENT_Cr - 1])
  {
    int  ccAlfCrApsId = cs.slice->getTileGroupCcAlfCrApsId();
    APS* aps = m_apsMap->getPS((cs.slice->getTileGroupCcAlfCrApsId() << NUM_APS_TYPE_LEN) + ALF_APS);
    if (aps == NULL)
    {
      aps = m_apsMap->allocatePS((ccAlfCrApsId << NUM_APS_TYPE_LEN) + ALF_APS);
      aps->setTemporalId(cs.slice->getTLayer());
    }
    aps->getCcAlfAPSParam().ccAlfFilterEnabled[COMPONENT_Cr - 1] = 1;
    aps->getCcAlfAPSParam().ccAlfFilterCount[COMPONENT_Cr - 1] = m_ccAlfFilterParam.ccAlfFilterCount[COMPONENT_Cr - 1];
    for ( int filterIdx = 0; filterIdx < MAX_NUM_CC_ALF_FILTERS; filterIdx++ )
    {
      aps->getCcAlfAPSParam().ccAlfFilterIdxEnabled[COMPONENT_Cr - 1][filterIdx] =
        m_ccAlfFilterParam.ccAlfFilterIdxEnabled[COMPONENT_Cr - 1][filterIdx];
      memcpy(aps->getCcAlfAPSParam().ccAlfCoeff[COMPONENT_Cr - 1][filterIdx],
             m_ccAlfFilterParam.ccAlfCoeff[COMPONENT_Cr - 1][filterIdx], sizeof(short) * MAX_NUM_CC_ALF_CHROMA_COEFF);
    }
    aps->setAPSId(ccAlfCrApsId);
    if (m_reuseApsId[COMPONENT_Cr - 1] < 0)
    {
      aps->getCcAlfAPSParam().newCcAlfFilter[COMPONENT_Cr - 1] = 1;
      m_apsMap->setChangedFlag((ccAlfCrApsId << NUM_APS_TYPE_LEN) + ALF_APS, true);
      aps->setTemporalId(cs.slice->getTLayer());
    }
    aps->setAPSType(ALF_APS);
    cs.slice->setTileGroupCcAlfCrEnabledFlag(true);
  }
  else
  {
    cs.slice->setTileGroupCcAlfCrEnabledFlag(false);
  }
}

void EncAdaptiveLoopFilter::ALFProcess( CodingStructure& cs, const double *lambdas
#if ENABLE_QPA
                                        , const double lambdaChromaWeight
#endif
                                        , Picture* pcPic, uint32_t numSliceSegments
)
{
  int layerIdx = cs.vps == nullptr ? 0 : cs.vps->getGeneralLayerIdx( cs.slice->getPic()->layerId );

  // IRAP AU is assumed
  if( !layerIdx && ( cs.slice->getPendingRasInit() || cs.slice->isIDRorBLA() ) )
  {
    memset( cs.slice->getAlfAPSs(), 0, sizeof( *cs.slice->getAlfAPSs() )*ALF_CTB_MAX_NUM_APS );
    m_apsIdStart = ALF_CTB_MAX_NUM_APS;

    m_apsMap->clear();
    for( int i = 0; i < ALF_CTB_MAX_NUM_APS; i++ )
    {
      APS* alfAPS = m_apsMap->getPS( ( i << NUM_APS_TYPE_LEN ) + ALF_APS );
      m_apsMap->clearChangedFlag( ( i << NUM_APS_TYPE_LEN ) + ALF_APS );
      if( alfAPS )
      {
        alfAPS->getAlfAPSParam().reset();
        alfAPS->getCcAlfAPSParam().reset();
        alfAPS = nullptr;
      }
    }
  }

#if ALF_IMPROVEMENT
  cs.slice->setTileGroupAlfFixedFilterSetIdx(-1);
#endif
  AlfParam alfParam;
  alfParam.reset();
  const TempCtx  ctxStart( m_CtxCache, AlfCtx( m_CABACEstimator->getCtx() ) );

  const TempCtx ctxStartCcAlf( m_CtxCache, SubCtx( Ctx::CcAlfFilterControlFlag, m_CABACEstimator->getCtx() ) );

  // set available filter shapes
  alfParam.filterShapes = m_filterShapes;

  // set clipping range
  m_clpRngs = cs.slice->getClpRngs();

  // set CTU ALF enable flags, it was already reset before ALF process
  for( int compIdx = 0; compIdx < MAX_NUM_COMPONENT; compIdx++ )
  {
    m_ctuEnableFlag[compIdx] = cs.picture->getAlfCtuEnableFlag( compIdx );
    m_ctuAlternative[compIdx] = cs.picture->getAlfCtuAlternativeData( compIdx );
  }


  // reset ALF parameters
  alfParam.reset();
  int shiftLuma = 2 * DISTORTION_PRECISION_ADJUSTMENT( m_inputBitDepth[CHANNEL_TYPE_LUMA] );
  int shiftChroma = 2 * DISTORTION_PRECISION_ADJUSTMENT( m_inputBitDepth[CHANNEL_TYPE_CHROMA] );
  m_lambda[COMPONENT_Y] = lambdas[COMPONENT_Y] * double( 1 << shiftLuma );
  m_lambda[COMPONENT_Cb] = lambdas[COMPONENT_Cb] * double( 1 << shiftChroma );
  m_lambda[COMPONENT_Cr] = lambdas[COMPONENT_Cr] * double( 1 << shiftChroma );

#if ALF_SAO_TRUE_ORG
  PelUnitBuf orgYuv = cs.getTrueOrgBuf();
#else
  PelUnitBuf orgYuv = cs.getOrgBuf();
#endif

  m_tempBuf.copyFrom( cs.getRecoBuf() );
  PelUnitBuf recYuv = m_tempBuf.getBuf( cs.area );
#if ALF_IMPROVEMENT
  recYuv.extendBorderPel( MAX_FILTER_LENGTH_FIXED >> 1 );
#else
  recYuv.extendBorderPel( MAX_ALF_FILTER_LENGTH >> 1 );
#endif

  // derive classification
  const CPelBuf& recLuma = recYuv.get( COMPONENT_Y );
  const PreCalcValues& pcv = *cs.pcv;
  bool clipTop = false, clipBottom = false, clipLeft = false, clipRight = false;
  int numHorVirBndry = 0, numVerVirBndry = 0;
  int horVirBndryPos[] = { 0, 0, 0 };
  int verVirBndryPos[] = { 0, 0, 0 };

  for( int yPos = 0; yPos < pcv.lumaHeight; yPos += pcv.maxCUHeight )
  {
    for( int xPos = 0; xPos < pcv.lumaWidth; xPos += pcv.maxCUWidth )
    {
      const int width = ( xPos + pcv.maxCUWidth > pcv.lumaWidth ) ? ( pcv.lumaWidth - xPos ) : pcv.maxCUWidth;
      const int height = ( yPos + pcv.maxCUHeight > pcv.lumaHeight ) ? ( pcv.lumaHeight - yPos ) : pcv.maxCUHeight;
      int rasterSliceAlfPad = 0;
      if( isCrossedByVirtualBoundaries( cs, xPos, yPos, width, height, clipTop, clipBottom, clipLeft, clipRight, numHorVirBndry, numVerVirBndry, horVirBndryPos, verVirBndryPos, rasterSliceAlfPad ) )
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
            buf.copyFrom( recYuv.subBuf( UnitArea( cs.area.chromaFormat, Area( xStart - ( clipL ? 0 : MAX_ALF_PADDING_SIZE ), yStart - ( clipT ? 0 : MAX_ALF_PADDING_SIZE ), wBuf, hBuf ) ) ) );
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

            const Area blkSrc( 0, 0, w, h );
            const Area blkDst( xStart, yStart, w, h );
            deriveClassification( m_classifier, buf.get( COMPONENT_Y ), blkDst, blkSrc 
#if ALF_IMPROVEMENT
              , cs, -1
#endif
            );

            xStart = xEnd;
          }

          yStart = yEnd;
        }
      }
      else
      {
        Area blk( xPos, yPos, width, height );
        deriveClassification( m_classifier, recLuma, blk, blk 
#if ALF_IMPROVEMENT
          , cs, -1
#endif
        );
      }
    }
  }

  // get CTB stats for filtering
  if( m_alfWSSD )
  {
    deriveStatsForFiltering<true>( orgYuv, recYuv, cs );
  }
  else
  {
    deriveStatsForFiltering<false>( orgYuv, recYuv, cs );
  }

  for (int ctbIdx = 0; ctbIdx < m_numCTUsInPic; ctbIdx++)
  {
    cs.slice->getPic()->getAlfCtbFilterIndex()[ctbIdx] = NUM_FIXED_FILTER_SETS;
  }
  // consider using new filter (only)
  alfParam.newFilterFlag[CHANNEL_TYPE_LUMA] = true;
  alfParam.newFilterFlag[CHANNEL_TYPE_CHROMA] = true;
  cs.slice->setTileGroupNumAps(1); // Only new filter for RD cost optimization
  // derive filter (luma)
  alfEncoder( cs, alfParam, orgYuv, recYuv, cs.getRecoBuf(), CHANNEL_TYPE_LUMA
#if ENABLE_QPA
            , lambdaChromaWeight
#endif
            );

  // derive filter (chroma)
  if (isChromaEnabled(cs.pcv->chrFormat))
  {
    alfEncoder( cs, alfParam, orgYuv, recYuv, cs.getRecoBuf(), CHANNEL_TYPE_CHROMA
#if ENABLE_QPA
              , lambdaChromaWeight
#endif
              );
  }

  // let alfEncoderCtb decide now
  alfParam.newFilterFlag[CHANNEL_TYPE_LUMA] = false;
  alfParam.newFilterFlag[CHANNEL_TYPE_CHROMA] = false;
  cs.slice->setTileGroupNumAps(0);
  m_CABACEstimator->getCtx() = AlfCtx(ctxStart);
  alfEncoderCtb(cs, alfParam
#if ENABLE_QPA
    , lambdaChromaWeight
#endif
  );

  for (int s = 0; s < numSliceSegments; s++)
  {
    if (pcPic->slices[s]->isLossless())
    {
      for (uint32_t ctuIdx = 0; ctuIdx < pcPic->slices[s]->getNumCtuInSlice(); ctuIdx++)
      {
        uint32_t ctuRsAddr = pcPic->slices[s]->getCtuAddrInSlice(ctuIdx);
        m_ctuEnableFlag[COMPONENT_Y][ctuRsAddr] = 0;
        m_ctuEnableFlag[COMPONENT_Cb][ctuRsAddr] = 0;
        m_ctuEnableFlag[COMPONENT_Cr][ctuRsAddr] = 0;
      }
    }
  }

  alfReconstructor(cs, recYuv);

  // Do not transmit CC ALF if it is unchanged
  if (cs.slice->getTileGroupAlfEnabledFlag(COMPONENT_Y))
  {
    for (int32_t lumaAlfApsId : cs.slice->getTileGroupApsIdLuma())
    {
      APS* aps = (lumaAlfApsId >= 0) ? m_apsMap->getPS((lumaAlfApsId << NUM_APS_TYPE_LEN) + ALF_APS) : nullptr;
      if (aps && m_apsMap->getChangedFlag((lumaAlfApsId << NUM_APS_TYPE_LEN) + ALF_APS))
      {
        aps->getCcAlfAPSParam().newCcAlfFilter[0] = false;
          aps->getCcAlfAPSParam().newCcAlfFilter[1] = false;
      }
    }
  }
  int chromaAlfApsId = ( cs.slice->getTileGroupAlfEnabledFlag(COMPONENT_Cb) || cs.slice->getTileGroupAlfEnabledFlag(COMPONENT_Cr) ) ? cs.slice->getTileGroupApsIdChroma() : -1;
  APS* aps = (chromaAlfApsId >= 0) ? m_apsMap->getPS((chromaAlfApsId << NUM_APS_TYPE_LEN) + ALF_APS) : nullptr;
  if (aps && m_apsMap->getChangedFlag((chromaAlfApsId << NUM_APS_TYPE_LEN) + ALF_APS))
  {
    aps->getCcAlfAPSParam().newCcAlfFilter[0] = false;
    aps->getCcAlfAPSParam().newCcAlfFilter[1] = false;
  }

  if (!cs.slice->getSPS()->getCCALFEnabledFlag())
  {
    return;
  }

  m_tempBuf.get(COMPONENT_Cb).copyFrom(cs.getRecoBuf().get(COMPONENT_Cb));
  m_tempBuf.get(COMPONENT_Cr).copyFrom(cs.getRecoBuf().get(COMPONENT_Cr));
  recYuv = m_tempBuf.getBuf(cs.area);
  recYuv.extendBorderPel(MAX_ALF_FILTER_LENGTH >> 1);
  
  if( m_alfWSSD )
  {
    deriveStatsForCcAlfFiltering<true>( orgYuv, recYuv, COMPONENT_Cb, m_numCTUsInWidth, ( 0 + 1 ), cs );
    deriveStatsForCcAlfFiltering<true>( orgYuv, recYuv, COMPONENT_Cr, m_numCTUsInWidth, ( 0 + 1 ), cs );
  }
  else
  {
    deriveStatsForCcAlfFiltering<false>( orgYuv, recYuv, COMPONENT_Cb, m_numCTUsInWidth, ( 0 + 1 ), cs );
    deriveStatsForCcAlfFiltering<false>( orgYuv, recYuv, COMPONENT_Cr, m_numCTUsInWidth, ( 0 + 1 ), cs );
  }

  initDistortionCcalf();

  m_CABACEstimator->getCtx() = SubCtx(Ctx::CcAlfFilterControlFlag, ctxStartCcAlf);
  deriveCcAlfFilter(cs, COMPONENT_Cb, orgYuv, recYuv, cs.getRecoBuf());
  m_CABACEstimator->getCtx() = SubCtx(Ctx::CcAlfFilterControlFlag, ctxStartCcAlf);
  deriveCcAlfFilter(cs, COMPONENT_Cr, orgYuv, recYuv, cs.getRecoBuf());

  xSetupCcAlfAPS(cs);

  for (int compIdx = 1; compIdx < getNumberValidComponents(cs.pcv->chrFormat); compIdx++)
  {
    ComponentID compID     = ComponentID(compIdx);
    if (m_ccAlfFilterParam.ccAlfFilterEnabled[compIdx - 1])
    {
      applyCcAlfFilter(cs, compID, cs.getRecoBuf().get(compID), recYuv, m_ccAlfFilterControl[compIdx - 1],
                       m_ccAlfFilterParam.ccAlfCoeff[compIdx - 1], -1);
    }
  }
}

double EncAdaptiveLoopFilter::deriveCtbAlfEnableFlags( CodingStructure& cs, const int iShapeIdx, ChannelType channel,
#if ENABLE_QPA
                                                       const double chromaWeight,
#endif
                                                       const int numClasses, const int numCoeff, double& distUnfilter
#if ALF_IMPROVEMENT
                                                       , int fixedFilterSetIdx
#endif
)
{
  TempCtx        ctxTempStart( m_CtxCache );
  TempCtx        ctxTempBest( m_CtxCache );
  TempCtx        ctxTempAltStart( m_CtxCache );
  TempCtx        ctxTempAltBest( m_CtxCache );
  const ComponentID compIDFirst = isLuma( channel ) ? COMPONENT_Y : COMPONENT_Cb;
  const ComponentID compIDLast = isLuma( channel ) ? COMPONENT_Y : COMPONENT_Cr;
#if ALF_IMPROVEMENT
  const int numAlts = isLuma(channel) ? m_alfParamTemp.numAlternativesLuma : m_alfParamTemp.numAlternativesChroma;
#else
  const int numAlts = isLuma( channel ) ? 1 : m_alfParamTemp.numAlternativesChroma;
#endif
  double cost = 0;
  distUnfilter = 0;

  setEnableFlag(m_alfParamTemp, channel, true);
#if ENABLE_QPA
  CHECK ((chromaWeight > 0.0) && (cs.slice->getFirstCtuRsAddrInSlice() != 0), "incompatible start CTU address, must be 0");
#endif

  reconstructCoeff(m_alfParamTemp, channel, true, isLuma(channel));
#if ALF_IMPROVEMENT
  for( int altIdx = 0; altIdx < (isLuma(channel) ? MAX_NUM_ALF_ALTERNATIVES_LUMA : MAX_NUM_ALF_ALTERNATIVES_CHROMA); altIdx++ )
#else
  for( int altIdx = 0; altIdx < (isLuma(channel) ? 1 : MAX_NUM_ALF_ALTERNATIVES_CHROMA); altIdx++)
#endif
  {
    for (int classIdx = 0; classIdx < (isLuma(channel) ? MAX_NUM_ALF_CLASSES : 1); classIdx++)
    {
      for (int i = 0; i < (isLuma(channel) ? MAX_NUM_ALF_LUMA_COEFF : MAX_NUM_ALF_CHROMA_COEFF); i++)
      {
#if ALF_IMPROVEMENT
        m_filterCoeffSet[altIdx*numClasses + classIdx][i] = isLuma(channel) ? m_coeffFinal[altIdx][classIdx * MAX_NUM_ALF_LUMA_COEFF + i] : m_chromaCoeffFinal[altIdx][i];
        m_filterClippSet[altIdx*numClasses + classIdx][i] = isLuma(channel) ? m_clippFinal[altIdx][classIdx * MAX_NUM_ALF_LUMA_COEFF + i] : m_chromaClippFinal[altIdx][i];
#else
        m_filterCoeffSet[isLuma(channel) ? classIdx : altIdx][i] = isLuma(channel) ? m_coeffFinal[classIdx * MAX_NUM_ALF_LUMA_COEFF + i] : m_chromaCoeffFinal[altIdx][i];
        m_filterClippSet[isLuma(channel) ? classIdx : altIdx][i] = isLuma(channel) ? m_clippFinal[classIdx * MAX_NUM_ALF_LUMA_COEFF + i] : m_chromaClippFinal[altIdx][i];
#endif
      }
    }
  }

  for( int ctuIdx = 0; ctuIdx < m_numCTUsInPic; ctuIdx++ )
  {
    for( int compID = compIDFirst; compID <= compIDLast; compID++ )
    {
#if ENABLE_QPA
      const double ctuLambda = chromaWeight > 0.0 ? (isLuma (channel) ? cs.picture->m_uEnerHpCtu[ctuIdx] : cs.picture->m_uEnerHpCtu[ctuIdx] / chromaWeight) : m_lambda[compID];
#else
      const double ctuLambda = m_lambda[compID];
#endif

      double distUnfilterCtu = m_ctbDistortionUnfilter[compID][ctuIdx];

      ctxTempStart = AlfCtx( m_CABACEstimator->getCtx() );
      m_CABACEstimator->resetBits();
      m_ctuEnableFlag[compID][ctuIdx] = 1;

      m_CABACEstimator->codeAlfCtuEnableFlag( cs, ctuIdx, compID, &m_alfParamTemp );

      if( isLuma( channel ) )
      {
        // Evaluate cost of signaling filter set index for convergence of filters enabled flag / filter derivation
        assert( cs.slice->getPic()->getAlfCtbFilterIndex()[ctuIdx] == NUM_FIXED_FILTER_SETS );
        assert( cs.slice->getTileGroupNumAps() == 1 );
        m_CABACEstimator->codeAlfCtuFilterIndex(cs, ctuIdx, &m_alfParamTemp.enabledFlag[COMPONENT_Y]);
      }
      double costOn = distUnfilterCtu + ctuLambda * FRAC_BITS_SCALE * m_CABACEstimator->getEstFracBits();

      ctxTempBest = AlfCtx( m_CABACEstimator->getCtx() );

      if( isLuma( channel ) )
      {
#if ALF_IMPROVEMENT
        double bestAltCost = MAX_DOUBLE;
        int bestAltIdx = -1;
        ctxTempAltStart = AlfCtx(ctxTempBest);
        for( int altIdx = 0; altIdx < numAlts; ++altIdx )
        {
          if( altIdx )
          {
            m_CABACEstimator->getCtx() = AlfCtx( ctxTempAltStart );
          }

          m_CABACEstimator->resetBits();
          m_ctuAlternative[compID][ctuIdx] = altIdx;
          m_CABACEstimator->codeAlfCtuAlternative(cs, ctuIdx, compID, &m_alfParamTemp, numAlts);
          double r_altCost = ctuLambda * FRAC_BITS_SCALE * m_CABACEstimator->getEstFracBits();
          double altDist = getFilteredDistortion(m_alfCovariance[compID][iShapeIdx][ctuIdx][fixedFilterSetIdx], numClasses, m_alfParamTemp.numLumaFilters[altIdx] - 1, numCoeff, altIdx);
          double altCost = altDist + r_altCost;

          if (altCost < bestAltCost)
          {
            bestAltCost = altCost;
            bestAltIdx = altIdx;
            ctxTempBest = AlfCtx(m_CABACEstimator->getCtx());
          }
        }
        m_ctuAlternative[compID][ctuIdx] = bestAltIdx;
        costOn += bestAltCost;
#else
        costOn += getFilteredDistortion( m_alfCovariance[compID][iShapeIdx][ctuIdx], numClasses, m_alfParamTemp.numLumaFilters - 1, numCoeff );
#endif
      }
      else
      {
        double bestAltCost = MAX_DOUBLE;
        int bestAltIdx = -1;
        ctxTempAltStart = AlfCtx( ctxTempBest );
        for( int altIdx = 0; altIdx < numAlts; ++altIdx )
        {
          if( altIdx )
          {
            m_CABACEstimator->getCtx() = AlfCtx( ctxTempAltStart );
          }

          m_CABACEstimator->resetBits();
          m_ctuAlternative[compID][ctuIdx] = altIdx;
          m_CABACEstimator->codeAlfCtuAlternative( cs, ctuIdx, compID, &m_alfParamTemp );
          double r_altCost = ctuLambda * FRAC_BITS_SCALE * m_CABACEstimator->getEstFracBits();

          double altDist = 0.;
#if ALF_IMPROVEMENT
          altDist += m_alfCovariance[compID][iShapeIdx][ctuIdx][0][0].calcErrorForCoeffs(  m_filterClippSet[altIdx], m_filterCoeffSet[altIdx], numCoeff, m_NUM_BITS );
#else
          altDist += m_alfCovariance[compID][iShapeIdx][ctuIdx][0].calcErrorForCoeffs(  m_filterClippSet[altIdx], m_filterCoeffSet[altIdx], numCoeff, m_NUM_BITS );
#endif

          double altCost = altDist + r_altCost;
          if( altCost < bestAltCost )
          {
            bestAltCost = altCost;
            bestAltIdx = altIdx;
            ctxTempBest = AlfCtx( m_CABACEstimator->getCtx() );
          }
        }
        m_ctuAlternative[compID][ctuIdx] = bestAltIdx;
        costOn += bestAltCost;
      }

      m_CABACEstimator->getCtx() = AlfCtx( ctxTempStart );
      m_CABACEstimator->resetBits();
      m_ctuEnableFlag[compID][ctuIdx] = 0;
      m_CABACEstimator->codeAlfCtuEnableFlag( cs, ctuIdx, compID, &m_alfParamTemp);
      double costOff = distUnfilterCtu + ctuLambda * FRAC_BITS_SCALE * m_CABACEstimator->getEstFracBits();

      if( costOn < costOff )
      {
        cost += costOn;
        m_CABACEstimator->getCtx() = AlfCtx( ctxTempBest );
        m_ctuEnableFlag[compID][ctuIdx] = 1;

      }
      else
      {
        cost += costOff;
        m_ctuEnableFlag[compID][ctuIdx] = 0;
        distUnfilter += distUnfilterCtu;
      }
    }
  }

  if( isChroma( channel ) )
  {
    setEnableFlag(m_alfParamTemp, channel, m_ctuEnableFlag);
  }

  return cost;
}

void EncAdaptiveLoopFilter::alfEncoder( CodingStructure& cs, AlfParam& alfParam, const PelUnitBuf& orgUnitBuf, const PelUnitBuf& recExtBuf, const PelUnitBuf& recBuf, const ChannelType channel
#if ENABLE_QPA
                                      , const double lambdaChromaWeight // = 0.0
#endif
                                      )
{
  const TempCtx  ctxStart( m_CtxCache, AlfCtx( m_CABACEstimator->getCtx() ) );
  TempCtx        ctxBest( m_CtxCache );

  double costMin = MAX_DOUBLE;
  std::vector<AlfFilterShape>& alfFilterShape = alfParam.filterShapes[channel];
  m_bitsNewFilter[channel] = 0;
  const int numClasses = isLuma( channel ) ? MAX_NUM_ALF_CLASSES : 1;
  int uiCoeffBits = 0;

  for( int iShapeIdx = 0; iShapeIdx < alfFilterShape.size(); iShapeIdx++ )
  {
#if ALF_IMPROVEMENT
    if (m_filterTypeTest[channel][alfFilterShape[iShapeIdx].filterType] == false)
    {
      continue;
    }
    int numFixedFilterSet = ( alfFilterShape[iShapeIdx].filterType == ALF_FILTER_EXT || alfFilterShape[iShapeIdx].filterType == ALF_FILTER_9_EXT ) ? 2 : 1;
#endif
      m_alfParamTemp = alfParam;
      //1. get unfiltered distortion
      if( isChroma( channel ) )
      {
        m_alfParamTemp.numAlternativesChroma = 1;
      }
#if ALF_IMPROVEMENT
      else
      {
        m_alfParamTemp.numAlternativesLuma = 1;
      }
#endif

      double cost = isLuma( channel ) ? m_unFiltDistCompnent[COMPONENT_Y] : m_unFiltDistCompnent[COMPONENT_Cb] + m_unFiltDistCompnent[COMPONENT_Cr];
      cost /= 1.001; // slight preference for unfiltered choice

      if( cost < costMin )
      {
        costMin = cost;
        setEnableFlag( alfParam, channel, false );

        // no CABAC signalling
        ctxBest = AlfCtx( ctxStart );
        setCtuEnableFlag( m_ctuEnableFlagTmp, channel, 0 );
        if( isChroma( channel ) )
        {
          setCtuAlternativeChroma( m_ctuAlternativeTmp, 0 );
        }
#if ALF_IMPROVEMENT
        else
        {
          std::fill_n( m_ctuAlternativeTmp[COMPONENT_Y], m_numCTUsInPic, 0 );
        }
#endif
      }
#if ALF_IMPROVEMENT
      const int nonLinearFlagMax = 1;
#else
       const int nonLinearFlagMax =
        ( isLuma( channel ) ? m_encCfg->getUseNonLinearAlfLuma() : m_encCfg->getUseNonLinearAlfChroma()) // For Chroma non linear flag is check for each alternative filter
        ? 2 : 1;
#endif
    for( int nonLinearFlag = 0; nonLinearFlag < nonLinearFlagMax; nonLinearFlag++ )
    {
#if ALF_IMPROVEMENT
      for( int fixedFilterSetIdx = 0; fixedFilterSetIdx < numFixedFilterSet; fixedFilterSetIdx++ )
      {
        for( int numAlternatives = 1; numAlternatives <= (isLuma(channel) ? std::max(1, std::min(m_numCTUsInPic, MAX_NUM_ALF_ALTERNATIVES_LUMA)) : getMaxNumAlternativesChroma()); numAlternatives++ )
#else
        for( int numAlternatives = isLuma( channel ) ? 1 : getMaxNumAlternativesChroma(); numAlternatives > 0; numAlternatives-- )
#endif
        {
          if( isChroma( channel ) )
          {
            m_alfParamTemp.numAlternativesChroma = numAlternatives;
          }
#if ALF_IMPROVEMENT
          else
          {
            m_alfParamTemp.numAlternativesLuma = numAlternatives;
          }
#endif
          //2. all CTUs are on
#if ALF_IMPROVEMENT
          m_alfParamTemp.filterType[channel] = alfFilterShape[iShapeIdx].filterType;
#endif
          setEnableFlag( m_alfParamTemp, channel, true );
#if !ALF_IMPROVEMENT
          m_alfParamTemp.nonLinearFlag[channel] = nonLinearFlag;
#endif
          m_CABACEstimator->getCtx() = AlfCtx( ctxStart );
          setCtuEnableFlag( m_ctuEnableFlag, channel, 1 );

          // all alternatives are on
          if( isChroma( channel ) )
          {
            initCtuAlternativeChroma( m_ctuAlternative );
          }
#if ALF_IMPROVEMENT
          else
          {
            initCtuAlternativeLuma( m_ctuAlternative );
          }
          bool update = false;
#endif
          cost = getFilterCoeffAndCost( cs, 0, channel, true, iShapeIdx, uiCoeffBits
#if ALF_IMPROVEMENT
          , fixedFilterSetIdx
#endif
          );

          if( cost < costMin )
          {
#if ALF_IMPROVEMENT
            update = true;
#endif
            costMin = cost;
            m_bitsNewFilter[channel] = uiCoeffBits;
            copyAlfParam( alfParam, m_alfParamTemp, channel );
            ctxBest = AlfCtx( m_CABACEstimator->getCtx() );
            setCtuEnableFlag( m_ctuEnableFlagTmp, channel, 1 );

            if( isChroma( channel ) )
            {
              copyCtuAlternativeChroma( m_ctuAlternativeTmp, m_ctuAlternative );
            }
#if ALF_IMPROVEMENT
            else
            {
              std::copy_n( m_ctuAlternative[COMPONENT_Y], m_numCTUsInPic, m_ctuAlternativeTmp[COMPONENT_Y] );
            }
#endif
          }

          //3. CTU decision
          double distUnfilter = 0;
#if ALF_IMPROVEMENT
          double prevItCost = cost;
#else
          double prevItCost = MAX_DOUBLE;
#endif
          const int iterNum = isLuma(channel) ? (2 * 4 + 1) : (2 * (2 + m_alfParamTemp.numAlternativesChroma - 1) + 1);

          for( int iter = 0; iter < iterNum; iter++ )
          {
            if ((iter & 0x01) == 0)
            {
              m_CABACEstimator->getCtx() = AlfCtx(ctxStart);
              cost = m_lambda[channel] * uiCoeffBits;
              cost += deriveCtbAlfEnableFlags(cs, iShapeIdx, channel,
#if ENABLE_QPA
                                              lambdaChromaWeight,
#endif
                                              numClasses, alfFilterShape[iShapeIdx].numCoeff, distUnfilter
#if ALF_IMPROVEMENT
                                              , fixedFilterSetIdx
#endif
                                              );
              if (cost < costMin)
              {
#if ALF_IMPROVEMENT
                update = true;
#endif
                costMin = cost;
                ctxBest = AlfCtx(m_CABACEstimator->getCtx());
                copyCtuEnableFlag(m_ctuEnableFlagTmp, m_ctuEnableFlag, channel);
                if( isChroma( channel ) )
                {
                  copyCtuAlternativeChroma( m_ctuAlternativeTmp, m_ctuAlternative );
                }
                else
                {
#if ALF_IMPROVEMENT
                  std::copy_n( m_ctuAlternative[COMPONENT_Y], m_numCTUsInPic, m_ctuAlternativeTmp[COMPONENT_Y] );
#endif
                }
                m_bitsNewFilter[channel] = uiCoeffBits;
                copyAlfParam(alfParam, m_alfParamTemp, channel);
              }
              else if( cost >= prevItCost )
              {
                // High probability that we have converged or we are diverging
                break;
              }
              prevItCost = cost;
            }
            else
            {
              // unfiltered distortion is added due to some CTBs may not use filter
              // no need to reset CABAC here, since uiCoeffBits is not affected
              /*cost = */getFilterCoeffAndCost( cs, distUnfilter, channel, true, iShapeIdx, uiCoeffBits 
#if ALF_IMPROVEMENT
                , fixedFilterSetIdx
#endif
              );
            }
          }//for iter
#if ALF_IMPROVEMENT
          if (update == false && numAlternatives > 1)
          {
            break;
          }
#endif
       }// for number of alternatives and reset ctu params and filters
#if ALF_IMPROVEMENT
      }//for (int fixedFilterSetIdx = 0; fixedFilterSetIdx < numFixedFilterSet; fixedFilterSetIdx++)
#endif
    }// for nonLineaFlag
  }//for shapeIdx
  m_CABACEstimator->getCtx() = AlfCtx( ctxBest );
  if( isChroma( channel ) )
  {
    copyCtuAlternativeChroma( m_ctuAlternative, m_ctuAlternativeTmp );
  }
  else
  {
#if ALF_IMPROVEMENT
    std::copy_n( m_ctuAlternativeTmp[COMPONENT_Y], m_numCTUsInPic, m_ctuAlternative[COMPONENT_Y] );
#endif

  }
  copyCtuEnableFlag( m_ctuEnableFlag, m_ctuEnableFlagTmp, channel );
}

void EncAdaptiveLoopFilter::copyAlfParam( AlfParam& alfParamDst, AlfParam& alfParamSrc, ChannelType channel )
{
  if( isLuma( channel ) )
  {
    memcpy( &alfParamDst, &alfParamSrc, sizeof( AlfParam ) );
  }
  else
  {
    alfParamDst.enabledFlag[COMPONENT_Cb] = alfParamSrc.enabledFlag[COMPONENT_Cb];
    alfParamDst.enabledFlag[COMPONENT_Cr] = alfParamSrc.enabledFlag[COMPONENT_Cr];
    alfParamDst.numAlternativesChroma = alfParamSrc.numAlternativesChroma;
#if ALF_IMPROVEMENT
    alfParamDst.filterType[CHANNEL_TYPE_CHROMA] = alfParamSrc.filterType[CHANNEL_TYPE_CHROMA];
    memcpy( alfParamDst.nonLinearFlag[CHANNEL_TYPE_CHROMA], alfParamSrc.nonLinearFlag[CHANNEL_TYPE_CHROMA], sizeof( alfParamDst.nonLinearFlag[CHANNEL_TYPE_CHROMA] ) );
#else
    alfParamDst.nonLinearFlag[CHANNEL_TYPE_CHROMA] = alfParamSrc.nonLinearFlag[CHANNEL_TYPE_CHROMA];
#endif
    memcpy( alfParamDst.chromaCoeff, alfParamSrc.chromaCoeff, sizeof( alfParamDst.chromaCoeff ) );
    memcpy( alfParamDst.chromaClipp, alfParamSrc.chromaClipp, sizeof( alfParamDst.chromaClipp ) );
  }
}

double EncAdaptiveLoopFilter::getFilterCoeffAndCost( CodingStructure& cs, double distUnfilter, ChannelType channel, bool bReCollectStat, int iShapeIdx, int& uiCoeffBits, 
#if ALF_IMPROVEMENT
  int fixedFilterSetIdx,
#endif
  bool onlyFilterCost )
{
  double dist = distUnfilter;
  uiCoeffBits = 0;
  AlfFilterShape& alfFilterShape = m_alfParamTemp.filterShapes[channel][iShapeIdx];
  //get filter coeff
  if( isLuma( channel ) )
  {
#if ALF_IMPROVEMENT
    for (int altIdx = 0; altIdx < m_alfParamTemp.numAlternativesLuma; ++altIdx)
    {
      //collect stat based on CTU decision
      if ( bReCollectStat )
      {
        getFrameStats( channel, iShapeIdx, altIdx 
#if ALF_IMPROVEMENT
          , fixedFilterSetIdx
#endif
        );
      }
      assert(alfFilterShape.numCoeff == m_alfCovarianceFrame[channel][iShapeIdx][0].numCoeff);

      AlfParam bestSliceParam;
      double bestCost = MAX_DOUBLE;
      double bestDist = MAX_DOUBLE;
      int bestCoeffBits = 0;
      const int nonLinearFlagMax = m_encCfg->getUseNonLinearAlfLuma() ? 2 : 1;
      for (int nonLinearFlag = 0; nonLinearFlag < nonLinearFlagMax; nonLinearFlag++)
      {
        m_alfParamTemp.nonLinearFlag[channel][altIdx] = nonLinearFlag;
        std::fill_n(m_alfClipMerged[iShapeIdx][0][0], MAX_NUM_ALF_LUMA_COEFF*MAX_NUM_ALF_CLASSES*MAX_NUM_ALF_CLASSES, m_alfParamTemp.nonLinearFlag[channel][altIdx] ? AlfNumClippingValues[CHANNEL_TYPE_LUMA] / 2 : 0);

        // Reset Merge Tmp Cov
        m_alfCovarianceMerged[iShapeIdx][MAX_NUM_ALF_CLASSES].reset(AlfNumClippingValues[channel]);
        m_alfCovarianceMerged[iShapeIdx][MAX_NUM_ALF_CLASSES + 1].reset(AlfNumClippingValues[channel]);
        int curCoeffBits;
        double curDist = mergeFiltersAndCost(m_alfParamTemp, alfFilterShape, m_alfCovarianceFrame[channel][iShapeIdx], m_alfCovarianceMerged[iShapeIdx], m_alfClipMerged[iShapeIdx], curCoeffBits, altIdx);
        double cost = curDist + m_lambda[channel] * curCoeffBits;
        if (cost < bestCost)
        {
          bestCost = cost;
          bestDist = curDist;
          bestCoeffBits = curCoeffBits;
          bestSliceParam = m_alfParamTemp;
        }
      }//for (int nonLinearFlag)
      uiCoeffBits += bestCoeffBits;
      dist += bestDist;
      m_alfParamTemp = bestSliceParam;
    } //for (int altIdx
    uiCoeffBits += lengthUvlc(m_alfParamTemp.numAlternativesLuma - 1);
    uiCoeffBits += m_alfParamTemp.numAlternativesLuma; // non-linear flags   
#else
    //collect stat based on CTU decision
    if( bReCollectStat )
    {
      getFrameStats( channel, iShapeIdx, 0 );
    }
    std::fill_n(m_alfClipMerged[iShapeIdx][0][0], MAX_NUM_ALF_LUMA_COEFF*MAX_NUM_ALF_CLASSES*MAX_NUM_ALF_CLASSES, m_alfParamTemp.nonLinearFlag[channel] ? AlfNumClippingValues[CHANNEL_TYPE_LUMA] / 2 : 0);
    // Reset Merge Tmp Cov
    m_alfCovarianceMerged[iShapeIdx][MAX_NUM_ALF_CLASSES].reset(AlfNumClippingValues[channel]);
    m_alfCovarianceMerged[iShapeIdx][MAX_NUM_ALF_CLASSES + 1].reset(AlfNumClippingValues[channel]);
    //distortion
    dist += mergeFiltersAndCost( m_alfParamTemp, alfFilterShape, m_alfCovarianceFrame[channel][iShapeIdx], m_alfCovarianceMerged[iShapeIdx], m_alfClipMerged[iShapeIdx], uiCoeffBits );
#endif
  }
  else
  {
    //distortion
    CHECK( m_alfParamTemp.numAlternativesChroma < 1, "Wrong number of m_alfParamTemp.numAlternativesChroma" );

    for( int altIdx = 0; altIdx < m_alfParamTemp.numAlternativesChroma; ++altIdx )
    {
      //collect stat based on CTU decision
      if ( bReCollectStat )
      {
        getFrameStats( channel, iShapeIdx, altIdx 
#if ALF_IMPROVEMENT
          , fixedFilterSetIdx
#endif
        );
      }
      assert(alfFilterShape.numCoeff == m_alfCovarianceFrame[channel][iShapeIdx][0].numCoeff);

      AlfParam bestSliceParam;
      double bestCost = MAX_DOUBLE;
      double bestDist = MAX_DOUBLE;
      int bestCoeffBits = 0;
      const int nonLinearFlagMax = m_encCfg->getUseNonLinearAlfChroma() ? 2 : 1;

      for( int nonLinearFlag = 0; nonLinearFlag < nonLinearFlagMax; nonLinearFlag++ )
      {
#if ALF_IMPROVEMENT
        m_alfParamTemp.nonLinearFlag[channel][altIdx] = nonLinearFlag;
#else
        int currentNonLinearFlag = m_alfParamTemp.nonLinearFlag[channel] ? 1 : 0;
        if (nonLinearFlag != currentNonLinearFlag)
        {
          continue;
        }
#endif

        std::fill_n( m_filterClippSet[altIdx], MAX_NUM_ALF_CHROMA_COEFF, nonLinearFlag ? AlfNumClippingValues[CHANNEL_TYPE_CHROMA] / 2 : 0 );
        double dist = m_alfCovarianceFrame[channel][iShapeIdx][0].pixAcc + deriveCoeffQuant( m_filterClippSet[altIdx], m_filterCoeffSet[altIdx], m_alfCovarianceFrame[channel][iShapeIdx][0], alfFilterShape, m_NUM_BITS, nonLinearFlag );

        for( int i = 0; i < MAX_NUM_ALF_CHROMA_COEFF; i++ )
        {
          m_alfParamTemp.chromaCoeff[altIdx][i] = m_filterCoeffSet[altIdx][i];
          m_alfParamTemp.chromaClipp[altIdx][i] = m_filterClippSet[altIdx][i];
        }
        int coeffBits = getChromaCoeffRate( m_alfParamTemp, altIdx );
        double cost = dist + m_lambda[channel] * coeffBits;
        if( cost < bestCost )
        {
          bestCost = cost;
          bestDist = dist;
          bestCoeffBits = coeffBits;
          bestSliceParam = m_alfParamTemp;
        }
      }
      uiCoeffBits += bestCoeffBits;
      dist += bestDist;
      m_alfParamTemp = bestSliceParam;
    }
    uiCoeffBits += lengthUvlc( m_alfParamTemp.numAlternativesChroma-1 );
#if ALF_IMPROVEMENT
    uiCoeffBits += m_alfParamTemp.numAlternativesChroma; // non-linear flags
#else
    uiCoeffBits++;
#endif
  }
  if (onlyFilterCost)
  {
    return dist + m_lambda[channel] * uiCoeffBits;
  }
  double rate = uiCoeffBits;
  m_CABACEstimator->resetBits();
  m_CABACEstimator->codeAlfCtuEnableFlags( cs, channel, &m_alfParamTemp);
  for( int ctuIdx = 0; ctuIdx < m_numCTUsInPic; ctuIdx++ )
  {
    if( isLuma( channel ) )
    {
      // Evaluate cost of signaling filter set index for convergence of filters enabled flag / filter derivation
      assert( cs.slice->getPic()->getAlfCtbFilterIndex()[ctuIdx] == NUM_FIXED_FILTER_SETS );
      assert( cs.slice->getTileGroupNumAps() == 1 );
      m_CABACEstimator->codeAlfCtuFilterIndex(cs, ctuIdx, &m_alfParamTemp.enabledFlag[COMPONENT_Y]);
#if ALF_IMPROVEMENT
      m_CABACEstimator->codeAlfCtuAlternative(cs, ctuIdx, COMPONENT_Y, &m_alfParamTemp, m_alfParamTemp.numAlternativesLuma);
#endif
    }
  }
  m_CABACEstimator->codeAlfCtuAlternatives( cs, channel, &m_alfParamTemp );
  rate += FRAC_BITS_SCALE * m_CABACEstimator->getEstFracBits();
  return dist + m_lambda[channel] * rate;
}

int EncAdaptiveLoopFilter::getChromaCoeffRate( AlfParam& alfParam, int altIdx )
{
  int iBits = 0;
#if ALF_IMPROVEMENT
  AlfFilterShape alfShape = m_filterShapes[CHANNEL_TYPE_CHROMA][m_filterTypeToStatIndex[CHANNEL_TYPE_CHROMA][alfParam.filterType[CHANNEL_TYPE_CHROMA]]];
#else
  AlfFilterShape alfShape(5);
#endif
  // Filter coefficients
#if ALF_IMPROVEMENT
  for( int orderIdx = 0; orderIdx < alfShape.numOrder; orderIdx++ )
  {
    int minBits = MAX_INT;
    int startIdx = orderIdx == 0 ? 0 : alfShape.indexSecOrder;
    int endIdx = orderIdx == 0 ? alfShape.indexSecOrder : alfShape.numCoeff - 1;
    for( int k = 0; k < 4; k++ )
    {
      int curBits = 0;
      for( int coeffIdx = startIdx; coeffIdx < endIdx; coeffIdx++ )
      {
        curBits += EncAdaptiveLoopFilter::lengthGolomb( alfParam.chromaCoeff[altIdx][coeffIdx], k + ( orderIdx == 0 ? alfShape.offset0 : ALF_ORDER ) );
      }
      if (curBits < minBits)
      {
        minBits = curBits;
      }
    }
    iBits += minBits + 2;
  }
#else
  for( int i = 0; i < alfShape.numCoeff - 1; i++ )
  {
    iBits += lengthUvlc( abs( alfParam.chromaCoeff[ altIdx ][ i ] ) );  // alf_coeff_chroma[altIdx][i]
    if( ( alfParam.chromaCoeff[ altIdx ][ i ] ) != 0 )
      iBits += 1;
  }
#endif
#if ALF_IMPROVEMENT
  if( m_alfParamTemp.nonLinearFlag[CHANNEL_TYPE_CHROMA][altIdx] )
#else
  if( m_alfParamTemp.nonLinearFlag[CHANNEL_TYPE_CHROMA] )
#endif
  {
    for (int i = 0; i < alfShape.numCoeff - 1; i++)
    {
      if( !abs( alfParam.chromaCoeff[altIdx][i] ) )
      {
        alfParam.chromaClipp[altIdx][i] = 0;
      }
    }
    iBits += ((alfShape.numCoeff - 1) << 1);
  }
  return iBits;
}

double EncAdaptiveLoopFilter::getUnfilteredDistortion( AlfCovariance* cov, ChannelType channel )
{
  double dist = 0;
  if( isLuma( channel ) )
  {
    dist = getUnfilteredDistortion( cov, MAX_NUM_ALF_CLASSES );
  }
  else
  {
    dist = getUnfilteredDistortion( cov, 1 );
  }
  return dist;
}

double EncAdaptiveLoopFilter::getUnfilteredDistortion( AlfCovariance* cov, const int numClasses )
{
  double dist = 0;
  for( int classIdx = 0; classIdx < numClasses; classIdx++ )
  {
    dist += cov[classIdx].pixAcc;
  }
  return dist;
}

double EncAdaptiveLoopFilter::getFilteredDistortion( AlfCovariance* cov, const int numClasses, const int numFiltersMinus1, const int numCoeff 
#if ALF_IMPROVEMENT
  , int altIdx
#endif
)
{
  double dist = 0;

  for( int classIdx = 0; classIdx < numClasses; classIdx++ )
  {
#if ALF_IMPROVEMENT
    dist += cov[classIdx].calcErrorForCoeffs(m_filterClippSet[altIdx * numClasses + classIdx], m_filterCoeffSet[altIdx * numClasses + classIdx], numCoeff, m_NUM_BITS);
#else
    dist += cov[classIdx].calcErrorForCoeffs(m_filterClippSet[classIdx], m_filterCoeffSet[classIdx], numCoeff, m_NUM_BITS);
#endif
  }

  return dist;
}

double EncAdaptiveLoopFilter::mergeFiltersAndCost( AlfParam& alfParam, AlfFilterShape& alfShape, AlfCovariance* covFrame, AlfCovariance* covMerged, int clipMerged[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF], int& uiCoeffBits 
#if ALF_IMPROVEMENT
  , int altIdx
#endif
)
{
  int numFiltersBest = 0;
  int numFilters = MAX_NUM_ALF_CLASSES;
  bool   codedVarBins[MAX_NUM_ALF_CLASSES];
  double errorForce0CoeffTab[MAX_NUM_ALF_CLASSES][2];
  double cost, cost0, dist, distForce0, costMin = MAX_DOUBLE;
  int coeffBits, coeffBitsForce0;

  mergeClasses( alfShape, covFrame, covMerged, clipMerged, MAX_NUM_ALF_CLASSES, m_filterIndices 
#if ALF_IMPROVEMENT
    , altIdx
#endif
  );

  while( numFilters >= 1 )
  {
    dist = deriveFilterCoeffs(covFrame, covMerged, clipMerged, alfShape, m_filterIndices[numFilters - 1], numFilters, errorForce0CoeffTab, alfParam
#if ALF_IMPROVEMENT
      , m_alfParamTemp.nonLinearFlag[CHANNEL_TYPE_LUMA][altIdx]
#endif
    );
    // filter coeffs are stored in m_filterCoeffSet
    distForce0 = getDistForce0( alfShape, numFilters, errorForce0CoeffTab, codedVarBins 
#if ALF_IMPROVEMENT
      , altIdx
#endif
    );
    coeffBits = deriveFilterCoefficientsPredictionMode( alfShape, m_filterCoeffSet, m_diffFilterCoeff, numFilters );
    coeffBitsForce0 = getCostFilterCoeffForce0( alfShape, m_filterCoeffSet, numFilters, codedVarBins 
#if ALF_IMPROVEMENT
      , altIdx
#endif
    );

    cost = dist + m_lambda[COMPONENT_Y] * coeffBits;
    cost0 = distForce0 + m_lambda[COMPONENT_Y] * coeffBitsForce0;

    if( cost0 < cost )
    {
      cost = cost0;
    }

    if( cost <= costMin )
    {
      costMin = cost;
      numFiltersBest = numFilters;
    }
    numFilters--;
  }

  dist = deriveFilterCoeffs( covFrame, covMerged, clipMerged, alfShape, m_filterIndices[numFiltersBest - 1], numFiltersBest, errorForce0CoeffTab, alfParam 
#if ALF_IMPROVEMENT
    , m_alfParamTemp.nonLinearFlag[CHANNEL_TYPE_LUMA][altIdx]
#endif
  );
  coeffBits = deriveFilterCoefficientsPredictionMode( alfShape, m_filterCoeffSet, m_diffFilterCoeff, numFiltersBest );
  distForce0 = getDistForce0( alfShape, numFiltersBest, errorForce0CoeffTab, codedVarBins 
#if ALF_IMPROVEMENT
    , altIdx
#endif
  );
  coeffBitsForce0 = getCostFilterCoeffForce0( alfShape, m_filterCoeffSet, numFiltersBest, codedVarBins 
#if ALF_IMPROVEMENT
    , altIdx
#endif
  );

  cost = dist + m_lambda[COMPONENT_Y] * coeffBits;
  cost0 = distForce0 + m_lambda[COMPONENT_Y] * coeffBitsForce0;
#if ALF_IMPROVEMENT
  alfParam.numLumaFilters[altIdx] = numFiltersBest;
#else
  alfParam.numLumaFilters = numFiltersBest;
#endif
  double distReturn;
  if (cost <= cost0)
  {
    distReturn = dist;
    alfParam.alfLumaCoeffDeltaFlag = 0;
    uiCoeffBits = coeffBits;
  }
  else
  {
    distReturn = distForce0;
    alfParam.alfLumaCoeffDeltaFlag = 1;
    uiCoeffBits = coeffBitsForce0;
    memcpy( alfParam.alfLumaCoeffFlag, codedVarBins, sizeof( codedVarBins ) );

    for( int varInd = 0; varInd < numFiltersBest; varInd++ )
    {
      if( codedVarBins[varInd] == 0 )
      {
        memset( m_filterCoeffSet[varInd], 0, sizeof( int )*MAX_NUM_ALF_LUMA_COEFF );
        memset( m_filterClippSet[varInd], 0, sizeof( int )*MAX_NUM_ALF_LUMA_COEFF );
      }
    }
  }
#if ALF_IMPROVEMENT
  for (int ind = 0; ind < alfParam.numLumaFilters[altIdx]; ++ind)
  {
    for (int i = 0; i < alfShape.numCoeff; i++)
    {
      alfParam.lumaCoeff[altIdx][ind * MAX_NUM_ALF_LUMA_COEFF + i] = m_filterCoeffSet[ind][i];
      alfParam.lumaClipp[altIdx][ind * MAX_NUM_ALF_LUMA_COEFF + i] = m_filterClippSet[ind][i];
    }
  }
  memcpy( alfParam.filterCoeffDeltaIdx[altIdx], m_filterIndices[numFiltersBest - 1], sizeof(short) * MAX_NUM_ALF_CLASSES );
  uiCoeffBits += getNonFilterCoeffRate(alfParam, altIdx);
#else
  for( int ind = 0; ind < alfParam.numLumaFilters; ++ind )
  {
    for( int i = 0; i < alfShape.numCoeff; i++ )
    {
      alfParam.lumaCoeff[ind * MAX_NUM_ALF_LUMA_COEFF + i] = m_filterCoeffSet[ind][i];
      alfParam.lumaClipp[ind * MAX_NUM_ALF_LUMA_COEFF + i] = m_filterClippSet[ind][i];
    }
  }

  memcpy( alfParam.filterCoeffDeltaIdx, m_filterIndices[numFiltersBest - 1], sizeof( short ) * MAX_NUM_ALF_CLASSES );
  uiCoeffBits += getNonFilterCoeffRate( alfParam );
#endif
  return distReturn;
}

int EncAdaptiveLoopFilter::getNonFilterCoeffRate( AlfParam& alfParam 
#if ALF_IMPROVEMENT
  , int altIdx
#endif
)
{
#if ALF_IMPROVEMENT
  CHECK( alfParam.numLumaFilters[altIdx] < 1, "Wrong number of alfParam.numLumaFilters[altIdx]" );
#else
  CHECK( alfParam.numLumaFilters < 1, "Wrong number of alfParam.numLumaFilters[altIdx]" );
#endif

  int len = 0   // alf_coefficients_delta_flag
#if ALF_IMPROVEMENT
    + lengthUvlc(alfParam.numLumaFilters[altIdx] - 1);  // alf_luma_num_filters_signalled_minus1   ue(v)
  if (alfParam.numLumaFilters[altIdx] > 1)
  {
    const int coeffLength = ceilLog2(alfParam.numLumaFilters[altIdx]);
#else
          + 2                                          // slice_alf_chroma_idc                     u(2)
          + lengthUvlc (alfParam.numLumaFilters - 1);  // alf_luma_num_filters_signalled_minus1   ue(v)

  if( alfParam.numLumaFilters > 1 )
  {
    const int coeffLength = ceilLog2(alfParam.numLumaFilters);
#endif
    for( int i = 0; i < MAX_NUM_ALF_CLASSES; i++ )
    {
      len += coeffLength;                              // alf_luma_coeff_delta_idx   u(v)
    }
  }
  return len;
}


int EncAdaptiveLoopFilter::getCostFilterCoeffForce0( AlfFilterShape& alfShape, int **pDiffQFilterCoeffIntPP, const int numFilters, bool* codedVarBins 
#if ALF_IMPROVEMENT
  , int altIdx
#endif
)
{
  int len = 0;
  // Filter coefficients
#if ALF_IMPROVEMENT
  for (int orderIdx = 0; orderIdx < alfShape.numOrder; orderIdx++)
  {
    int minBits = MAX_INT;
    int startIdx = orderIdx == 0 ? 0 : alfShape.indexSecOrder;
    int endIdx = orderIdx == 0 ? alfShape.indexSecOrder : alfShape.numCoeff - 1;
    for (int k = 0; k < 4; k++)
    {
      int curBits = 0;
      for (int filtIdx = 0; filtIdx < numFilters; filtIdx++)
      {
        for (int coeffIdx = startIdx; coeffIdx < endIdx; coeffIdx++)
        {
          if (codedVarBins[filtIdx])
          {
            curBits += lengthGolomb(pDiffQFilterCoeffIntPP[filtIdx][coeffIdx], k + (orderIdx == 0 ? alfShape.offset0 : ALF_ORDER));
          }
          else
          {
            curBits += lengthGolomb(0, k + (orderIdx == 0 ? alfShape.offset0 : ALF_ORDER));
          }
        }
      }
      if (curBits < minBits)
      {
        minBits = curBits;
      }
    }
    len += minBits + 2;
  }
#else
  for( int ind = 0; ind < numFilters; ++ind )
  {
    if( codedVarBins[ind] )
    {
      for( int i = 0; i < alfShape.numCoeff - 1; i++ )
      {
        len += lengthUvlc( abs( pDiffQFilterCoeffIntPP[ ind ][ i ] ) ); // alf_coeff_luma_delta[i][j]
        if( ( abs( pDiffQFilterCoeffIntPP[ ind ][ i ] ) != 0 ) )
          len += 1;
      }
    }
    else
    {
      for (int i = 0; i < alfShape.numCoeff - 1; i++)
      {
        len += lengthUvlc( 0 ); // alf_coeff_luma_delta[i][j]
      }
    }
  }
#endif
#if ALF_IMPROVEMENT
  if( m_alfParamTemp.nonLinearFlag[CHANNEL_TYPE_LUMA][altIdx] )
#else
  if( m_alfParamTemp.nonLinearFlag[CHANNEL_TYPE_LUMA] )
#endif
  {
    for (int ind = 0; ind < numFilters; ++ind)
    {
      for (int i = 0; i < alfShape.numCoeff - 1; i++)
      {
        if (!abs(pDiffQFilterCoeffIntPP[ind][i]))
        {
          m_filterClippSet[ind][i] = 0;
        }
        len += 2;
      }
    }
  }

  return len;
}

int EncAdaptiveLoopFilter::deriveFilterCoefficientsPredictionMode( AlfFilterShape& alfShape, int **filterSet, int** filterCoeffDiff, const int numFilters )
{
  return (m_alfParamTemp.nonLinearFlag[CHANNEL_TYPE_LUMA] ? getCostFilterClipp(alfShape, filterSet, numFilters) : 0) + getCostFilterCoeff(alfShape, filterSet, numFilters);
}

int EncAdaptiveLoopFilter::getCostFilterCoeff( AlfFilterShape& alfShape, int **pDiffQFilterCoeffIntPP, const int numFilters )
{
  return lengthFilterCoeffs( alfShape, numFilters, pDiffQFilterCoeffIntPP );  // alf_coeff_luma_delta[i][j];
}

int EncAdaptiveLoopFilter::getCostFilterClipp( AlfFilterShape& alfShape, int **pDiffQFilterCoeffIntPP, const int numFilters )
{
  for (int filterIdx = 0; filterIdx < numFilters; ++filterIdx)
  {
    for (int i = 0; i < alfShape.numCoeff - 1; i++)
    {
      if (!abs(pDiffQFilterCoeffIntPP[filterIdx][i]))
      {
        m_filterClippSet[filterIdx][i] = 0;
      }
    }
  }
  return (numFilters * (alfShape.numCoeff - 1)) << 1;
}

int EncAdaptiveLoopFilter::lengthFilterCoeffs( AlfFilterShape& alfShape, const int numFilters, int **FilterCoeff )
{
  int bitCnt = 0;

#if ALF_IMPROVEMENT
  for( int orderIdx = 0; orderIdx < alfShape.numOrder; orderIdx++ )
  {
    int minBits = MAX_INT;
    int startIdx = orderIdx == 0 ? 0 : alfShape.indexSecOrder;
    int endIdx = orderIdx == 0 ? alfShape.indexSecOrder : alfShape.numCoeff - 1;
    for( int k = 0; k < 4; k++ )
    {
      int curBits = 0;
      for( int filtIdx = 0; filtIdx < numFilters; filtIdx++ )
      {
        for( int coeffIdx = startIdx; coeffIdx < endIdx; coeffIdx++ )
        {
          curBits += lengthGolomb( FilterCoeff[filtIdx][coeffIdx], k + (orderIdx == 0 ? alfShape.offset0 : ALF_ORDER ) );
        }
      }
      if (curBits < minBits)
      {
        minBits = curBits;
      }
    }
    bitCnt += minBits + 2;
  }
#else
  for( int ind = 0; ind < numFilters; ++ind )
  {
    for( int i = 0; i < alfShape.numCoeff - 1; i++ )
    {
      bitCnt += lengthUvlc( abs( FilterCoeff[ ind ][ i ] ) );
      if( abs( FilterCoeff[ ind ][ i ] ) != 0 )
        bitCnt += 1;
    }
  }
#endif
  return bitCnt;
}


double EncAdaptiveLoopFilter::getDistForce0( AlfFilterShape& alfShape, const int numFilters, double errorTabForce0Coeff[MAX_NUM_ALF_CLASSES][2], bool* codedVarBins 
#if ALF_IMPROVEMENT
  , int altIdx
#endif
)
{
  int bitsVarBin[MAX_NUM_ALF_CLASSES];

#if ALF_IMPROVEMENT
  int bestK[2] = { 0 };
  for( int orderIdx = 0; orderIdx < alfShape.numOrder; orderIdx++)
  {
    int minBits = MAX_INT;
    int startIdx = orderIdx == 0 ? 0 : alfShape.indexSecOrder;
    int endIdx = orderIdx == 0 ? alfShape.indexSecOrder : alfShape.numCoeff - 1;
    for( int k = 0; k < 4; k++ )
    {
      int curBits = 0;
      for( int filtIdx = 0; filtIdx < numFilters; filtIdx++ )
      {
        for( int coeffIdx = startIdx; coeffIdx < endIdx; coeffIdx++ )
        {
          curBits += lengthGolomb( m_filterCoeffSet[filtIdx][coeffIdx], k + ( orderIdx == 0 ? alfShape.offset0 : ALF_ORDER ) );
        }
      }
      if( curBits < minBits )
      {
        bestK[orderIdx] = k;
        minBits = curBits;
      }
    }
  }
  bestK[0] += alfShape.offset0;
  bestK[1] += ALF_ORDER;
#endif
  for( int ind = 0; ind < numFilters; ++ind )
  {
    bitsVarBin[ind] = 0;
    for( int i = 0; i < alfShape.numCoeff - 1; i++ )
    {
#if ALF_IMPROVEMENT
      if (i < alfShape.indexSecOrder)
      {
        bitsVarBin[ind] += lengthGolomb( m_filterCoeffSet[ind][i], bestK[0] );
      }
      else
      {
        bitsVarBin[ind] += lengthGolomb(m_filterCoeffSet[ind][i], bestK[1]);
      }
#else
      bitsVarBin[ ind ] += lengthUvlc( abs( m_filterCoeffSet[ ind ][ i ] ) );
      if( abs( m_filterCoeffSet[ ind ][ i ] ) != 0 )
        bitsVarBin[ ind ] += 1;
#endif
    }
  }

  int zeroBitsVarBin = 0;
  for (int i = 0; i < alfShape.numCoeff - 1; i++)
  {
#if ALF_IMPROVEMENT
    if (i < alfShape.indexSecOrder)
    {
      zeroBitsVarBin += lengthGolomb(0, bestK[0]);
    }
    else
    {
      zeroBitsVarBin += lengthGolomb(0, bestK[1]);
    }
#else
    zeroBitsVarBin += lengthUvlc( 0 );
#endif
  }
#if ALF_IMPROVEMENT
  if( m_alfParamTemp.nonLinearFlag[CHANNEL_TYPE_LUMA][altIdx] )
#else
  if( m_alfParamTemp.nonLinearFlag[CHANNEL_TYPE_LUMA] )
#endif
  {
    for (int ind = 0; ind < numFilters; ++ind)
    {
      for (int i = 0; i < alfShape.numCoeff - 1; i++)
      {
        if (!abs(m_filterCoeffSet[ind][i]))
        {
          m_filterClippSet[ind][i] = 0;
        }
      }
    }
  }

  double distForce0 = getDistCoeffForce0( codedVarBins, errorTabForce0Coeff, bitsVarBin, zeroBitsVarBin, numFilters);

  return distForce0;
}
double EncAdaptiveLoopFilter::getDistCoeffForce0( bool* codedVarBins, double errorForce0CoeffTab[MAX_NUM_ALF_CLASSES][2], int* bitsVarBin, int zeroBitsVarBin, const int numFilters)
{
  double distForce0 = 0;
  std::memset( codedVarBins, 0, sizeof( *codedVarBins ) * MAX_NUM_ALF_CLASSES );

  for( int filtIdx = 0; filtIdx < numFilters; filtIdx++ )
  {
    double costDiff = (errorForce0CoeffTab[filtIdx][0] + m_lambda[COMPONENT_Y] * zeroBitsVarBin) - (errorForce0CoeffTab[filtIdx][1] + m_lambda[COMPONENT_Y] * bitsVarBin[filtIdx]);
    codedVarBins[filtIdx] = costDiff > 0 ? true : false;
    distForce0 += errorForce0CoeffTab[filtIdx][codedVarBins[filtIdx] ? 1 : 0];
  }

  return distForce0;
}

int EncAdaptiveLoopFilter::lengthUvlc( int uiCode )
{
  int uiLength = 1;
  int uiTemp = ++uiCode;

  CHECK( !uiTemp, "Integer overflow" );

  while( 1 != uiTemp )
  {
    uiTemp >>= 1;
    uiLength += 2;
  }
  // Take care of cases where uiLength > 32
  return ( uiLength >> 1 ) + ( ( uiLength + 1 ) >> 1 );
}

#if ALF_IMPROVEMENT
int EncAdaptiveLoopFilter::lengthGolomb(int coeffVal, int k, bool signed_coeff)
{
  int numBins = 0;
  unsigned int symbol = abs(coeffVal);
  while (symbol >= (unsigned int)(1 << k))
  {
    numBins++;
    symbol -= 1 << k;
    k++;
  }
  numBins += (k + 1);
  if (signed_coeff && coeffVal != 0)
  {
    numBins++;
  }
  return numBins;
}
#endif

double EncAdaptiveLoopFilter::deriveFilterCoeffs( AlfCovariance* cov, AlfCovariance* covMerged, int clipMerged[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF], AlfFilterShape& alfShape, short* filterIndices, int numFilters, double errorTabForce0Coeff[MAX_NUM_ALF_CLASSES][2], AlfParam& alfParam 
#if ALF_IMPROVEMENT
  , bool nonLinear
#endif
)
{
  double error = 0.0;
  AlfCovariance& tmpCov = covMerged[MAX_NUM_ALF_CLASSES];

  for( int filtIdx = 0; filtIdx < numFilters; filtIdx++ )
  {
    tmpCov.reset();
    bool found_clip = false;
    for( int classIdx = 0; classIdx < MAX_NUM_ALF_CLASSES; classIdx++ )
    {
      if( filterIndices[classIdx] == filtIdx )
      {
        tmpCov += cov[classIdx];
        if( !found_clip )
        {
          found_clip = true; // clip should be at the adress of shortest one
          memcpy(m_filterClippSet[filtIdx], clipMerged[numFilters-1][classIdx], sizeof(int[MAX_NUM_ALF_LUMA_COEFF]));
        }
      }
    }

    // Find coeffcients
    assert(alfShape.numCoeff == tmpCov.numCoeff);
#if ALF_IMPROVEMENT
    errorTabForce0Coeff[filtIdx][1] = tmpCov.pixAcc + deriveCoeffQuant( m_filterClippSet[filtIdx], m_filterCoeffSet[filtIdx], tmpCov, alfShape, m_NUM_BITS, nonLinear );
#else
    errorTabForce0Coeff[filtIdx][1] = tmpCov.pixAcc + deriveCoeffQuant( m_filterClippSet[filtIdx], m_filterCoeffSet[filtIdx], tmpCov, alfShape, m_NUM_BITS, false );
#endif
    errorTabForce0Coeff[filtIdx][0] = tmpCov.pixAcc;
    error += errorTabForce0Coeff[filtIdx][1];
  }
  return error;
}

double EncAdaptiveLoopFilter::deriveCoeffQuant( int *filterClipp, int *filterCoeffQuant, const AlfCovariance& cov, const AlfFilterShape& shape, const int bitDepth, const bool optimizeClip )
{
  const int factor = 1 << ( bitDepth - 1 );
  const int max_value = factor - 1;
  const int min_value = -factor + 1;

  const int numCoeff = shape.numCoeff;
  double    filterCoeff[MAX_NUM_ALF_LUMA_COEFF];

  cov.optimizeFilter( shape, filterClipp, filterCoeff, optimizeClip );
  roundFiltCoeff( filterCoeffQuant, filterCoeff, numCoeff, factor );

  for ( int i = 0; i < numCoeff - 1; i++ )
  {
    filterCoeffQuant[i] = std::min( max_value, std::max( min_value, filterCoeffQuant[i] ) );
  }
  filterCoeffQuant[numCoeff - 1] = 0;

  int modified=1;

  double errRef=cov.calcErrorForCoeffs( filterClipp, filterCoeffQuant, numCoeff, bitDepth );
  while( modified )
  {
    modified=0;
    for( int sign: {1, -1} )
    {
      double errMin = MAX_DOUBLE;
      int minInd = -1;

      for( int k = 0; k < numCoeff-1; k++ )
      {
        if( filterCoeffQuant[k] - sign > max_value || filterCoeffQuant[k] - sign < min_value )
        {
          continue;
        }

        filterCoeffQuant[k] -= sign;

        double error = cov.calcErrorForCoeffs( filterClipp, filterCoeffQuant, numCoeff, bitDepth );
        if( error < errMin )
        {
          errMin = error;
          minInd = k;
        }
        filterCoeffQuant[k] += sign;
      }
      if( errMin < errRef )
      {
        filterCoeffQuant[minInd] -= sign;
        modified++;
        errRef = errMin;
      }
    }
  }

  return errRef;
}

void EncAdaptiveLoopFilter::roundFiltCoeff( int *filterCoeffQuant, double *filterCoeff, const int numCoeff, const int factor )
{
  for( int i = 0; i < numCoeff; i++ )
  {
    int sign = filterCoeff[i] > 0 ? 1 : -1;
    filterCoeffQuant[i] = int( filterCoeff[i] * sign * factor + 0.5 ) * sign;
  }
}

void EncAdaptiveLoopFilter::roundFiltCoeffCCALF(int16_t *filterCoeffQuant, double *filterCoeff, const int numCoeff,
                                                const int factor)
{
  for( int i = 0; i < numCoeff; i++ )
  {
    int sign = filterCoeff[i] > 0 ? 1 : -1;
    double best_err = 128.0*128.0;
    int best_index = 0;
    for(int k = 0; k < CCALF_CANDS_COEFF_NR; k++)
    {
      double err = (filterCoeff[i] * sign * factor - CCALF_SMALL_TAB[k]);
      err = err*err;
      if(err < best_err)
      {
        best_err = err;
        best_index = k;
      }
    }
    filterCoeffQuant[i] = CCALF_SMALL_TAB[best_index] * sign;
  }
}

void EncAdaptiveLoopFilter::mergeClasses( const AlfFilterShape& alfShape, AlfCovariance* cov, AlfCovariance* covMerged, int clipMerged[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF], const int numClasses, short filterIndices[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES] 
#if ALF_IMPROVEMENT
  , const int altIdx
#endif
)
{
  int     tmpClip[MAX_NUM_ALF_LUMA_COEFF];
  int     bestMergeClip[MAX_NUM_ALF_LUMA_COEFF];
  double  err[MAX_NUM_ALF_CLASSES];
  double  bestMergeErr = std::numeric_limits<double>::max();
  bool    availableClass[MAX_NUM_ALF_CLASSES];
  uint8_t indexList[MAX_NUM_ALF_CLASSES];
  uint8_t indexListTemp[MAX_NUM_ALF_CLASSES];
  int numRemaining = numClasses;

  memset( filterIndices, 0, sizeof( short ) * MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_CLASSES );

  for( int i = 0; i < numClasses; i++ )
  {
    filterIndices[numRemaining - 1][i] = i;
    indexList[i] = i;
    availableClass[i] = true;
    covMerged[i] = cov[i];
#if ALF_IMPROVEMENT
    covMerged[i].numBins = m_alfParamTemp.nonLinearFlag[CHANNEL_TYPE_LUMA][altIdx] ? AlfNumClippingValues[COMPONENT_Y] : 1;
#else
    covMerged[i].numBins = m_alfParamTemp.nonLinearFlag[CHANNEL_TYPE_LUMA] ? AlfNumClippingValues[COMPONENT_Y] : 1;
#endif
  }

  // Try merging different covariance matrices

  // temporal AlfCovariance structure is allocated as the last element in covMerged array, the size of covMerged is MAX_NUM_ALF_CLASSES + 1
  AlfCovariance& tmpCov = covMerged[MAX_NUM_ALF_CLASSES];
#if ALF_IMPROVEMENT
  tmpCov.numBins = m_alfParamTemp.nonLinearFlag[CHANNEL_TYPE_LUMA][altIdx] ? AlfNumClippingValues[COMPONENT_Y] : 1;
#else
  tmpCov.numBins = m_alfParamTemp.nonLinearFlag[CHANNEL_TYPE_LUMA] ? AlfNumClippingValues[COMPONENT_Y] : 1;
#endif

#if ALF_IMPROVEMENT
  for (int i = 0; i < MAX_NUM_ALF_CLASSES; i++)
  {
    for (int j = 0; j < MAX_NUM_ALF_CLASSES; j++)
    {
      classChanged[i][j] = true;
      memset(clipHistory[i][j], 0, sizeof(int)*MAX_NUM_ALF_LUMA_COEFF);
    }
  }
  memset(errorHistory, 0, sizeof(double)*MAX_NUM_ALF_CLASSES*MAX_NUM_ALF_CLASSES);
#endif

  // init Clip
  for( int i = 0; i < numClasses; i++ )
  {
#if ALF_IMPROVEMENT
    std::fill_n(clipMerged[numRemaining - 1][i], MAX_NUM_ALF_LUMA_COEFF, m_alfParamTemp.nonLinearFlag[CHANNEL_TYPE_LUMA][altIdx] ? AlfNumClippingValues[CHANNEL_TYPE_LUMA] / 2 : 0);
    if ( m_alfParamTemp.nonLinearFlag[CHANNEL_TYPE_LUMA][altIdx] )
#else
    std::fill_n(clipMerged[numRemaining-1][i], MAX_NUM_ALF_LUMA_COEFF, m_alfParamTemp.nonLinearFlag[CHANNEL_TYPE_LUMA] ? AlfNumClippingValues[CHANNEL_TYPE_LUMA] / 2 : 0);
    if ( m_alfParamTemp.nonLinearFlag[CHANNEL_TYPE_LUMA] )
#endif
    {
      err[i] = covMerged[i].optimizeFilterClip( alfShape, clipMerged[numRemaining-1][i] );
    }
    else
    {
      err[i] = covMerged[i].calculateError( clipMerged[numRemaining-1][i] );
    }
  }

  while( numRemaining >= 2 )
  {
    double errorMin = std::numeric_limits<double>::max();
    int bestToMergeIdx1 = 0, bestToMergeIdx2 = 1;

    for( int i = 0; i < numClasses - 1; i++ )
    {
      if( availableClass[i] )
      {
        for( int j = i + 1; j < numClasses; j++ )
        {
          if( availableClass[j] )
          {
            double error1 = err[i];
            double error2 = err[j];

            tmpCov.add( covMerged[i], covMerged[j] );
            for( int l = 0; l < MAX_NUM_ALF_LUMA_COEFF; ++l )
            {
              tmpClip[l] = (clipMerged[numRemaining-1][i][l] + clipMerged[numRemaining-1][j][l] + 1 ) >> 1;
            }
#if ALF_IMPROVEMENT
            double errorMerged = 0;
            if (classChanged[i][j])
            {
              errorMerged = tmpCov.calculateError(tmpClip);
              classChanged[i][j] = false;
              errorHistory[i][j] = errorMerged;
              memcpy(clipHistory[i][j], tmpClip, sizeof(tmpClip));
            }
            else
            {
              errorMerged = errorHistory[i][j];
              memcpy(tmpClip, clipHistory[i][j], sizeof(tmpClip));
            }
#else
            double errorMerged = m_alfParamTemp.nonLinearFlag[CHANNEL_TYPE_LUMA] ? tmpCov.optimizeFilterClip(alfShape, tmpClip) : tmpCov.calculateError(tmpClip);
#endif

            double error = errorMerged - error1 - error2;

            if( error < errorMin )
            {
              bestMergeErr = errorMerged;
              memcpy(bestMergeClip, tmpClip, sizeof(bestMergeClip));
              errorMin = error;
              bestToMergeIdx1 = i;
              bestToMergeIdx2 = j;
            }
          }
        }
      }
    }

    covMerged[bestToMergeIdx1] += covMerged[bestToMergeIdx2];
    memcpy(clipMerged[numRemaining-2], clipMerged[numRemaining-1], sizeof(int[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF]));
    memcpy(clipMerged[numRemaining-2][bestToMergeIdx1], bestMergeClip, sizeof(bestMergeClip));
    err[bestToMergeIdx1] = bestMergeErr;
    availableClass[bestToMergeIdx2] = false;
#if ALF_IMPROVEMENT
    for (int i = 0; i < numClasses; i++)
    {
      classChanged[bestToMergeIdx1][i] = classChanged[i][bestToMergeIdx1] = true;
    }
#endif

    for( int i = 0; i < numClasses; i++ )
    {
      if( indexList[i] == bestToMergeIdx2 )
      {
        indexList[i] = bestToMergeIdx1;
      }
    }

    numRemaining--;
    if( numRemaining <= numClasses )
    {
      std::memcpy( indexListTemp, indexList, sizeof( uint8_t ) * numClasses );

      bool exist = false;
      int ind = 0;

      for( int j = 0; j < numClasses; j++ )
      {
        exist = false;
        for( int i = 0; i < numClasses; i++ )
        {
          if( indexListTemp[i] == j )
          {
            exist = true;
            break;
          }
        }

        if( exist )
        {
          for( int i = 0; i < numClasses; i++ )
          {
            if( indexListTemp[i] == j )
            {
              filterIndices[numRemaining - 1][i] = ind;
              indexListTemp[i] = -1;
            }
          }
          ind++;
        }
      }
    }
  }
}

void EncAdaptiveLoopFilter::getFrameStats( ChannelType channel, int iShapeIdx, int altIdx 
#if ALF_IMPROVEMENT
  , int fixedFilterSetIdx
#endif
)
{
  int numClasses = isLuma( channel ) ? MAX_NUM_ALF_CLASSES : 1;

  for (int i = 0; i < numClasses; i++)
  {
    m_alfCovarianceFrame[channel][iShapeIdx][i].reset(AlfNumClippingValues[channel]);
  }
  if (isLuma(channel))
  {
#if ALF_IMPROVEMENT
    getFrameStat( m_alfCovarianceFrame[CHANNEL_TYPE_LUMA][iShapeIdx], m_alfCovariance[COMPONENT_Y][iShapeIdx], m_ctuEnableFlag[COMPONENT_Y], m_ctuAlternative[COMPONENT_Y], numClasses, altIdx, fixedFilterSetIdx );
#else
    getFrameStat( m_alfCovarianceFrame[CHANNEL_TYPE_LUMA][iShapeIdx], m_alfCovariance[COMPONENT_Y][iShapeIdx], m_ctuEnableFlag[COMPONENT_Y], nullptr, numClasses, altIdx );
#endif
  }
  else
  {
    getFrameStat( m_alfCovarianceFrame[CHANNEL_TYPE_CHROMA][iShapeIdx], m_alfCovariance[COMPONENT_Cb][iShapeIdx], m_ctuEnableFlag[COMPONENT_Cb], m_ctuAlternative[COMPONENT_Cb], numClasses, altIdx 
#if ALF_IMPROVEMENT
      , fixedFilterSetIdx
#endif
    );
    getFrameStat( m_alfCovarianceFrame[CHANNEL_TYPE_CHROMA][iShapeIdx], m_alfCovariance[COMPONENT_Cr][iShapeIdx], m_ctuEnableFlag[COMPONENT_Cr], m_ctuAlternative[COMPONENT_Cr], numClasses, altIdx 
#if ALF_IMPROVEMENT
      , fixedFilterSetIdx
#endif
    );
  }
}

#if ALF_IMPROVEMENT
void EncAdaptiveLoopFilter::getFrameStat( AlfCovariance* frameCov, AlfCovariance*** ctbCov, uint8_t* ctbEnableFlags, uint8_t* ctbAltIdx, const int numClasses, int altIdx, int fixedFilterSetIdx)
#else
void EncAdaptiveLoopFilter::getFrameStat( AlfCovariance* frameCov, AlfCovariance** ctbCov, uint8_t* ctbEnableFlags, uint8_t* ctbAltIdx, const int numClasses, int altIdx )
#endif
{
#if !ALF_IMPROVEMENT
  const ChannelType channel = (!ctbAltIdx ? CHANNEL_TYPE_LUMA : CHANNEL_TYPE_CHROMA);
#endif
  for( int ctuIdx = 0; ctuIdx < m_numCTUsInPic; ctuIdx++ )
  {
    if( ctbEnableFlags[ctuIdx]  )
    {
      for( int classIdx = 0; classIdx < numClasses; classIdx++ )
      {
#if ALF_IMPROVEMENT
        if (altIdx == ctbAltIdx[ctuIdx])
        {
          frameCov[classIdx] += ctbCov[ctuIdx][fixedFilterSetIdx][classIdx];
        }
#else
        if( isLuma( channel ) || altIdx == ctbAltIdx[ctuIdx] )
        {
          frameCov[classIdx] += ctbCov[ctuIdx][classIdx];
        }
#endif
      }
    }
  }
}

template<bool alfWSSD>
void EncAdaptiveLoopFilter::deriveStatsForFiltering( PelUnitBuf& orgYuv, PelUnitBuf& recYuv, CodingStructure& cs )
{
  int ctuRsAddr = 0;
  const int numberOfComponents = getNumberValidComponents( m_chromaFormat );

  // init CTU stats buffers
  for( int compIdx = 0; compIdx < numberOfComponents; compIdx++ )
  {
    const ComponentID compID = ComponentID( compIdx );
    const int numClasses = isLuma( compID ) ? MAX_NUM_ALF_CLASSES : 1;

    for( int shape = 0; shape != m_filterShapes[toChannelType( compID )].size(); shape++ )
    {
#if ALF_IMPROVEMENT
      if (m_alfCovariance[compIdx][shape] == nullptr)
      {
        continue;
      }
      int numFixedFilterSet = ( m_filterShapes[toChannelType(compID)][shape].filterType == ALF_FILTER_EXT || m_filterShapes[toChannelType(compID)][shape].filterType == ALF_FILTER_9_EXT ) ? 2 : 1;
#endif
      for( int classIdx = 0; classIdx < numClasses; classIdx++ )
      {
        for( int ctuIdx = 0; ctuIdx < m_numCTUsInPic; ctuIdx++ )
        {
#if ALF_IMPROVEMENT
          for( int fixedFilterSetIdx = 0; fixedFilterSetIdx < numFixedFilterSet; fixedFilterSetIdx++ )
          {
            m_alfCovariance[compIdx][shape][ctuIdx][fixedFilterSetIdx][classIdx].reset(AlfNumClippingValues[toChannelType(compID)]);
          }
#else
          m_alfCovariance[compIdx][shape][ctuIdx][classIdx].reset(AlfNumClippingValues[toChannelType( compID )]);
#endif
        }
      }
    }
  }

  // init Frame stats buffers
  const int numberOfChannels = getNumberValidChannels( m_chromaFormat );
  for( int channelIdx = 0; channelIdx < numberOfChannels; channelIdx++ )
  {
    const ChannelType channelID = ChannelType( channelIdx );
    const int numClasses = isLuma( channelID ) ? MAX_NUM_ALF_CLASSES : 1;

    for( int shape = 0; shape != m_filterShapes[channelIdx].size(); shape++ )
    {
#if ALF_IMPROVEMENT
      if (m_alfCovarianceFrame[channelIdx][shape] == nullptr)
      {
        continue;
      }
#endif
      for( int classIdx = 0; classIdx < numClasses; classIdx++ )
      {
        m_alfCovarianceFrame[channelIdx][shape][classIdx].reset(AlfNumClippingValues[channelID]);
      }
    }
  }

  const PreCalcValues& pcv = *cs.pcv;
  bool clipTop = false, clipBottom = false, clipLeft = false, clipRight = false;
  int numHorVirBndry = 0, numVerVirBndry = 0;
  int horVirBndryPos[] = { 0, 0, 0 };
  int verVirBndryPos[] = { 0, 0, 0 };

  for( int yPos = 0; yPos < m_picHeight; yPos += m_maxCUHeight )
  {
    for( int xPos = 0; xPos < m_picWidth; xPos += m_maxCUWidth )
    {
      const int width = ( xPos + m_maxCUWidth > m_picWidth ) ? ( m_picWidth - xPos ) : m_maxCUWidth;
      const int height = ( yPos + m_maxCUHeight > m_picHeight ) ? ( m_picHeight - yPos ) : m_maxCUHeight;
      int rasterSliceAlfPad = 0;

      if( isCrossedByVirtualBoundaries( cs, xPos, yPos, width, height, clipTop, clipBottom, clipLeft, clipRight, numHorVirBndry, numVerVirBndry, horVirBndryPos, verVirBndryPos, rasterSliceAlfPad ) )
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
            const int wBuf = w + (clipL ? 0 : MAX_ALF_PADDING_SIZE) + (clipR ? 0 : MAX_ALF_PADDING_SIZE);
            const int hBuf = h + (clipT ? 0 : MAX_ALF_PADDING_SIZE) + (clipB ? 0 : MAX_ALF_PADDING_SIZE);
            PelUnitBuf recBuf = m_tempBuf2.subBuf( UnitArea( cs.area.chromaFormat, Area( 0, 0, wBuf, hBuf ) ) );
            recBuf.copyFrom( recYuv.subBuf( UnitArea( cs.area.chromaFormat, Area( xStart - (clipL ? 0 : MAX_ALF_PADDING_SIZE), yStart - (clipT ? 0 : MAX_ALF_PADDING_SIZE), wBuf, hBuf ) ) ) );
            // pad top-left unavailable samples for raster slice
            if ( xStart == xPos && yStart == yPos && ( rasterSliceAlfPad & 1 ) )
            {
              recBuf.padBorderPel( MAX_ALF_PADDING_SIZE, 1 );
            }

            // pad bottom-right unavailable samples for raster slice
            if ( xEnd == xPos + width && yEnd == yPos + height && ( rasterSliceAlfPad & 2 ) )
            {
              recBuf.padBorderPel( MAX_ALF_PADDING_SIZE, 2 );
            }
            recBuf.extendBorderPel( MAX_ALF_PADDING_SIZE );
            recBuf = recBuf.subBuf( UnitArea ( cs.area.chromaFormat, Area( clipL ? 0 : MAX_ALF_PADDING_SIZE, clipT ? 0 : MAX_ALF_PADDING_SIZE, w, h ) ) );

            const UnitArea area( m_chromaFormat, Area( 0, 0, w, h ) );
            const UnitArea areaDst( m_chromaFormat, Area( xStart, yStart, w, h ) );
            for( int compIdx = 0; compIdx < numberOfComponents; compIdx++ )
            {
              const ComponentID compID = ComponentID( compIdx );
              const CompArea& compArea = area.block( compID );

              int  recStride = recBuf.get( compID ).stride;
              Pel* rec = recBuf.get( compID ).bufAt( compArea );

              int  orgStride = orgYuv.get(compID).stride;
              Pel* org = orgYuv.get(compID).bufAt(xStart >> ::getComponentScaleX(compID, m_chromaFormat), yStart >> ::getComponentScaleY(compID, m_chromaFormat));
              ChannelType chType = toChannelType( compID );

              for( int shape = 0; shape != m_filterShapes[chType].size(); shape++ )
              {
              const CompArea& compAreaDst = areaDst.block( compID );
#if ALF_IMPROVEMENT
              for (int fixedFilterSetIdx = 0; fixedFilterSetIdx < ((m_filterShapes[chType][shape].filterType == ALF_FILTER_EXT || m_filterShapes[chType][shape].filterType == ALF_FILTER_9_EXT) ? 2 : 1); fixedFilterSetIdx++)
              {
                getBlkStats<alfWSSD>( m_alfCovariance[compIdx][shape][ctuRsAddr][fixedFilterSetIdx], m_filterShapes[chType][shape], compIdx ? nullptr : m_classifier, org, orgStride, rec, recStride, compAreaDst, compArea, chType, fixedFilterSetIdx);
              }
#else
              getBlkStats( m_alfCovariance[compIdx][shape][ctuRsAddr], m_filterShapes[chType][shape], compIdx ? nullptr : m_classifier, org, orgStride, rec, recStride, compAreaDst, compArea, chType, ((compIdx == 0) ? m_alfVBLumaCTUHeight : m_alfVBChmaCTUHeight), (compIdx == 0) ? m_alfVBLumaPos : m_alfVBChmaPos );
#endif
              }
            }

            xStart = xEnd;
          }

          yStart = yEnd;
        }

        for( int compIdx = 0; compIdx < numberOfComponents; compIdx++ )
        {
          const ComponentID compID = ComponentID( compIdx );

          ChannelType chType = toChannelType( compID );

          for( int shape = 0; shape != m_filterShapes[chType].size(); shape++ )
          {
            const int numClasses = isLuma( compID ) ? MAX_NUM_ALF_CLASSES : 1;

            for( int classIdx = 0; classIdx < numClasses; classIdx++ )
            {
#if ALF_IMPROVEMENT
              m_alfCovarianceFrame[chType][shape][isLuma( compID ) ? classIdx : 0] += m_alfCovariance[compIdx][shape][ctuRsAddr][0][classIdx];
#else
              m_alfCovarianceFrame[chType][shape][isLuma( compID ) ? classIdx : 0] += m_alfCovariance[compIdx][shape][ctuRsAddr][classIdx];
#endif
            }
          }
        }
      }
      else
      {
        const UnitArea area( m_chromaFormat, Area( xPos, yPos, width, height ) );

        for( int compIdx = 0; compIdx < numberOfComponents; compIdx++ )
        {
          const ComponentID compID = ComponentID( compIdx );
          const CompArea &  compArea = area.block( compID );

          int  recStride = recYuv.get( compID ).stride;
          Pel *rec = recYuv.get( compID ).bufAt( compArea );

          int  orgStride = orgYuv.get( compID ).stride;
          Pel *org = orgYuv.get( compID ).bufAt( compArea );

          ChannelType chType = toChannelType( compID );
          const int numClasses = isLuma( compID ) ? MAX_NUM_ALF_CLASSES : 1;

          for( int shape = 0; shape != m_filterShapes[chType].size(); shape++ )
          {
#if ALF_IMPROVEMENT
            if( m_filterTypeTest[chType][m_filterShapes[chType][shape].filterType] == false )
            {
              continue;
            }
            for( int fixedFilterSetIdx = 0; fixedFilterSetIdx < ((m_filterShapes[chType][shape].filterType == ALF_FILTER_EXT || m_filterShapes[chType][shape].filterType == ALF_FILTER_9_EXT) ? 2 : 1); fixedFilterSetIdx++ )
            {
              getBlkStats<alfWSSD>( m_alfCovariance[compIdx][shape][ctuRsAddr][fixedFilterSetIdx], m_filterShapes[chType][shape], compIdx ? nullptr : m_classifier, org, orgStride, rec, recStride, compArea, compArea, chType, fixedFilterSetIdx);
            }
#else
            getBlkStats( m_alfCovariance[compIdx][shape][ctuRsAddr], m_filterShapes[chType][shape], compIdx ? nullptr : m_classifier, org, orgStride, rec, recStride, compArea, compArea, chType, ( ( compIdx == 0 ) ? m_alfVBLumaCTUHeight : m_alfVBChmaCTUHeight ), ( compIdx == 0 ) ? m_alfVBLumaPos : m_alfVBChmaPos );
#endif        

          for( int classIdx = 0; classIdx < numClasses; classIdx++ )
          {
#if ALF_IMPROVEMENT
            m_alfCovarianceFrame[chType][shape][isLuma( compID ) ? classIdx : 0] += m_alfCovariance[compIdx][shape][ctuRsAddr][0][classIdx];
#else
            m_alfCovarianceFrame[chType][shape][isLuma( compID ) ? classIdx : 0] += m_alfCovariance[compIdx][shape][ctuRsAddr][classIdx];
#endif
          }
          }
        }
      }
      ctuRsAddr++;
    }
  }

  initDistortion(
#if ALF_IMPROVEMENT
    cs
#endif
  );
}

#if ALF_IMPROVEMENT
template<bool m_alfWSSD>
void EncAdaptiveLoopFilter::getBlkStats( AlfCovariance* alfCovariance, const AlfFilterShape& shape, AlfClassifier** classifier, const Pel* org, const int orgStride, const Pel* rec, const int recStride, const CompArea& areaDst, const CompArea& area, const ChannelType channel, int fixedFilterSetIdx )
#else
void EncAdaptiveLoopFilter::getBlkStats( AlfCovariance* alfCovariance, const AlfFilterShape& shape, AlfClassifier** classifier, Pel* org, const int orgStride, Pel* rec, const int recStride, const CompArea& areaDst, const CompArea& area, const ChannelType channel, int vbCTUHeight, int vbPos )
#endif
{
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  Pel ELocal[MAX_NUM_ALF_LUMA_COEFF][MaxAlfNumClippingValues];
#else
  int ELocal[MAX_NUM_ALF_LUMA_COEFF][MaxAlfNumClippingValues];
#endif

  const int numBins = AlfNumClippingValues[channel];
  int transposeIdx = 0;
  int classIdx = 0;

  for( int i = 0; i < area.height; i++ )
  {
#if !ALF_IMPROVEMENT
    int vbDistance = ((areaDst.y + i) % vbCTUHeight) - vbPos;
#endif
    for( int j = 0; j < area.width; j++ )
    {
#if ALF_IMPROVEMENT
      if( classifier && (classifier[areaDst.y + i][areaDst.x + j]>>2) == m_ALF_UNUSED_CLASSIDX && (classifier[areaDst.y + i][areaDst.x + j] & 0x3) == m_ALF_UNUSED_TRANSPOSIDX )
#else
      if( classifier && classifier[areaDst.y + i][areaDst.x + j].classIdx == m_ALF_UNUSED_CLASSIDX && classifier[areaDst.y + i][areaDst.x + j].transposeIdx == m_ALF_UNUSED_TRANSPOSIDX )
#endif
      {
        continue;
      }
      std::memset( ELocal, 0, sizeof( ELocal ) );
      if( classifier )
      {
        AlfClassifier& cl = classifier[areaDst.y + i][areaDst.x + j];
#if ALF_IMPROVEMENT
        transposeIdx = cl & 0x3;
        classIdx = cl >> 2;
#else
        transposeIdx = cl.transposeIdx;
        classIdx = cl.classIdx;
#endif
      }

      double weight = 1.0;
      if( m_alfWSSD )
      {
        weight = m_lumaLevelToWeightPLUT[org[j]];
      }
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
      Intermediate_Int yLocal = org[j] - rec[j];
#else
      int yLocal = org[j] - rec[j];
#endif

#if ALF_IMPROVEMENT
      calcCovariance( ELocal, rec + j, recStride, shape, transposeIdx, channel, ( shape.filterType == ALF_FILTER_9_EXT || shape.filterType == ALF_FILTER_EXT ) ? m_fixFilterResult : nullptr, Position(area.x + j ,area.y + i), fixedFilterSetIdx );
#else
      calcCovariance( ELocal, rec + j, recStride, shape, transposeIdx, channel, vbDistance );
#endif

      for( int k = 0; k < shape.numCoeff; k++ )
      {
        for( int l = k; l < shape.numCoeff; l++ )
        {
          for( int b0 = 0; b0 < numBins; b0++ )
          {
            for( int b1 = 0; b1 < numBins; b1++ )
            {
              if (m_alfWSSD)
              {
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
                alfCovariance[classIdx].E[b0][b1][k][l] += weight * (ELocal[k][b0] * (double)ELocal[l][b1]);
#else
                alfCovariance[classIdx].E[b0][b1][k][l] += weight * (double)(ELocal[k][b0] * ELocal[l][b1]);
#endif
              }
              else
              {
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
                alfCovariance[classIdx].E[b0][b1][k][l] += ELocal[k][b0] * (double)ELocal[l][b1];
#else
                alfCovariance[classIdx].E[b0][b1][k][l] += ELocal[k][b0] * ELocal[l][b1];
#endif
              }
            }
          }
        }
        for( int b = 0; b < numBins; b++ )
        {
          if (m_alfWSSD)
          {
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
            alfCovariance[classIdx].y[b][k] += weight * (ELocal[k][b] * (double)yLocal);
#else
            alfCovariance[classIdx].y[b][k] += weight * (double)(ELocal[k][b] * yLocal);
#endif
          }
          else
          {
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
            alfCovariance[classIdx].y[b][k] += ELocal[k][b] * (double)yLocal;
#else
            alfCovariance[classIdx].y[b][k] += ELocal[k][b] * yLocal;
#endif
          }
        }
      }
      if (m_alfWSSD)
      {
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
        alfCovariance[classIdx].pixAcc += weight * (yLocal * (double)yLocal);
#else
        alfCovariance[classIdx].pixAcc += weight * (double)(yLocal * yLocal);
#endif
      }
      else
      {
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
        alfCovariance[classIdx].pixAcc += yLocal * (double)yLocal;
#else
        alfCovariance[classIdx].pixAcc += yLocal * yLocal;
#endif
      }
    }
    org += orgStride;
    rec += recStride;
  }

  int numClasses = classifier ? MAX_NUM_ALF_CLASSES : 1;
  for( classIdx = 0; classIdx < numClasses; classIdx++ )
  {
    for( int k = 1; k < shape.numCoeff; k++ )
    {
      for( int l = 0; l < k; l++ )
      {
        for( int b0 = 0; b0 < numBins; b0++ )
        {
          for( int b1 = 0; b1 < numBins; b1++ )
          {
            alfCovariance[classIdx].E[b0][b1][k][l] = alfCovariance[classIdx].E[b1][b0][l][k];
          }
        }
      }
    }
  }
}

#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
#if ALF_IMPROVEMENT
void EncAdaptiveLoopFilter::calcCovariance( Pel ELocal[MAX_NUM_ALF_LUMA_COEFF][MaxAlfNumClippingValues], const Pel *rec, const int stride, const AlfFilterShape& shape, const int transposeIdx, const ChannelType channel, Pel ***fixedFitlerResults, Position pos, int fixedFilterSetIdx )
#else
void EncAdaptiveLoopFilter::calcCovariance( Pel ELocal[MAX_NUM_ALF_LUMA_COEFF][MaxAlfNumClippingValues], const Pel *rec, const int stride, const AlfFilterShape& shape, const int transposeIdx, const ChannelType channel, int vbDistance )
#endif
#else
#if ALF_IMPROVEMENT
void EncAdaptiveLoopFilter::calcCovariance( int ELocal[MAX_NUM_ALF_LUMA_COEFF][MaxAlfNumClippingValues], const Pel *rec, const int stride, const AlfFilterShape& shape, const int transposeIdx, const ChannelType channel, Pel *fixedFitlerResults, int fixedFilterSetIdx )
#else
void EncAdaptiveLoopFilter::calcCovariance( int ELocal[MAX_NUM_ALF_LUMA_COEFF][MaxAlfNumClippingValues], const Pel *rec, const int stride, const AlfFilterShape& shape, const int transposeIdx, const ChannelType channel, int vbDistance )
#endif
#endif
{
#if !ALF_IMPROVEMENT
  int clipTopRow = -4;
  int clipBotRow = 4;
  if (vbDistance >= -3 && vbDistance < 0)
  {
    clipBotRow = -vbDistance - 1;
    clipTopRow = -clipBotRow; // symmetric
  }
  else if (vbDistance >= 0 && vbDistance < 3)
  {
    clipTopRow = -vbDistance;
    clipBotRow = -clipTopRow; // symmetric
  }
#endif
  const int *filterPattern = shape.pattern.data();
  const int halfFilterLength = shape.filterLength >> 1;
  const Pel* clip = m_alfClippingValues[channel];
  const int numBins = AlfNumClippingValues[channel];

  int k = 0;

#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  const Pel curr = rec[0];
#else
  const short curr = rec[0];
#endif

  if( transposeIdx == 0 )
  {
    for( int i = -halfFilterLength; i < 0; i++ )
    {
#if ALF_IMPROVEMENT
      const Pel* rec0 = rec + i * stride;
      const Pel* rec1 = rec - i * stride;
#else
      const Pel* rec0 = rec + std::max(i, clipTopRow) * stride;
      const Pel* rec1 = rec - std::max(i, -clipBotRow) * stride;
#endif
      for( int j = -halfFilterLength - i; j <= halfFilterLength + i; j++, k++ )
      {
        for( int b = 0; b < numBins; b++ )
        {
          ELocal[filterPattern[k]][b] += clipALF(clip[b], curr, rec0[j], rec1[-j]);
        }
      }
    }
    for( int j = -halfFilterLength; j < 0; j++, k++ )
    {
      for( int b = 0; b < numBins; b++ )
      {
        ELocal[filterPattern[k]][b] += clipALF(clip[b], curr, rec[j], rec[-j]);
      }
    }
  }
  else if( transposeIdx == 1 )
  {
    for( int j = -halfFilterLength; j < 0; j++ )
    {
      const Pel* rec0 = rec + j;
      const Pel* rec1 = rec - j;
      for (int i = -halfFilterLength - j; i <= halfFilterLength + j; i++, k++)
      {
        for (int b = 0; b < numBins; b++)
        {
#if ALF_IMPROVEMENT
          ELocal[filterPattern[k]][b] += clipALF(clip[b], curr, rec0[i * stride], rec1[-i * stride]);
#else
          ELocal[filterPattern[k]][b] += clipALF(clip[b], curr, rec0[std::max(i, clipTopRow) * stride], rec1[-std::max(i, -clipBotRow) * stride]);
#endif
        }
      }
    }
    for (int i = -halfFilterLength; i < 0; i++, k++)
    {
      for (int b = 0; b < numBins; b++)
      {
#if ALF_IMPROVEMENT
        ELocal[filterPattern[k]][b] += clipALF(clip[b], curr, rec[i * stride], rec[-i * stride]);
#else
        ELocal[filterPattern[k]][b] += clipALF(clip[b], curr, rec[std::max(i, clipTopRow) * stride], rec[-std::max(i, -clipBotRow) * stride]);
#endif
      }
    }
  }
  else if( transposeIdx == 2 )
  {
    for( int i = -halfFilterLength; i < 0; i++ )
    {
#if ALF_IMPROVEMENT
      const Pel* rec0 = rec + i * stride;
      const Pel* rec1 = rec - i * stride;
#else
      const Pel* rec0 = rec + std::max(i, clipTopRow) * stride;
      const Pel* rec1 = rec - std::max(i, -clipBotRow) * stride;
#endif

      for( int j = halfFilterLength + i; j >= -halfFilterLength - i; j--, k++ )
      {
        for( int b = 0; b < numBins; b++ )
        {
          ELocal[filterPattern[k]][b] += clipALF(clip[b], curr, rec0[j], rec1[-j]);
        }
      }
    }
    for( int j = -halfFilterLength; j < 0; j++, k++ )
    {
      for( int b = 0; b < numBins; b++ )
      {
        ELocal[filterPattern[k]][b] += clipALF(clip[b], curr, rec[j], rec[-j]);
      }
    }
  }
  else
  {
    for( int j = -halfFilterLength; j < 0; j++ )
    {
      const Pel* rec0 = rec + j;
      const Pel* rec1 = rec - j;
      for (int i = halfFilterLength + j; i >= -halfFilterLength - j; i--, k++)
      {
        for (int b = 0; b < numBins; b++)
        {
#if ALF_IMPROVEMENT
          ELocal[filterPattern[k]][b] += clipALF(clip[b], curr, rec0[i * stride], rec1[-i * stride]);
#else
          ELocal[filterPattern[k]][b] += clipALF(clip[b], curr, rec0[std::max(i, clipTopRow) * stride], rec1[-std::max(i, -clipBotRow) * stride]);
#endif
        }
      }
    }
    for (int i = -halfFilterLength; i < 0; i++, k++)
    {
      for (int b = 0; b < numBins; b++)
      {
#if ALF_IMPROVEMENT
        ELocal[filterPattern[k]][b] += clipALF(clip[b], curr, rec[i * stride], rec[-i * stride]);
#else
        ELocal[filterPattern[k]][b] += clipALF(clip[b], curr, rec[std::max(i, clipTopRow) * stride], rec[-std::max(i, -clipBotRow) * stride]);
#endif
      }
    }
  }
#if ALF_IMPROVEMENT
  if( shape.filterType == ALF_FILTER_9_EXT )
  {
    for( int i = filterPattern[k]; i < filterPattern[k] + EXT_LENGTH; i++ )
    {
      for ( int b = 0; b < numBins; b++ )
      {
        ELocal[i][b] += clipALF(clip[b], curr, fixedFitlerResults[(i - filterPattern[k]) * EXT_LENGTH + fixedFilterSetIdx][pos.y][pos.x]);
      }
    }

    for( int b = 0; b < numBins; b++ )
    {
      ELocal[shape.numCoeff - 1][b] += curr;
    }
  }
  else if( shape.filterType == ALF_FILTER_EXT )
  {
    for( int i = 0; i < EXT_LENGTH; i++ )
    {
      for( int b = 0; b < numBins; b++ )
      {
        ELocal[i][b] += clipALF(clip[b], curr, fixedFitlerResults[i * EXT_LENGTH + fixedFilterSetIdx][pos.y][pos.x]);
      }
    }

    for( int b = 0; b < numBins; b++ )
    {
      ELocal[shape.numCoeff - 1][b] += curr;
    }
  }
  else
#endif
  for( int b = 0; b < numBins; b++ )
  {
    ELocal[filterPattern[k]][b] += curr;
  }
}



void EncAdaptiveLoopFilter::setEnableFlag( AlfParam& alfSlicePara, ChannelType channel, bool val )
{
  if( channel == CHANNEL_TYPE_LUMA )
  {
    alfSlicePara.enabledFlag[COMPONENT_Y] = val;
  }
  else
  {
    alfSlicePara.enabledFlag[COMPONENT_Cb] = alfSlicePara.enabledFlag[COMPONENT_Cr] = val;
  }
}

void EncAdaptiveLoopFilter::setEnableFlag( AlfParam& alfSlicePara, ChannelType channel, uint8_t** ctuFlags )
{
  const ComponentID compIDFirst = isLuma( channel ) ? COMPONENT_Y : COMPONENT_Cb;
  const ComponentID compIDLast = isLuma( channel ) ? COMPONENT_Y : COMPONENT_Cr;
  for( int compId = compIDFirst; compId <= compIDLast; compId++ )
  {
    alfSlicePara.enabledFlag[compId] = false;
    for( int i = 0; i < m_numCTUsInPic; i++ )
    {
      if( ctuFlags[compId][i] )
      {
        alfSlicePara.enabledFlag[compId] = true;
        break;
      }
    }
  }
}

void EncAdaptiveLoopFilter::copyCtuEnableFlag( uint8_t** ctuFlagsDst, uint8_t** ctuFlagsSrc, ChannelType channel )
{
  if( isLuma( channel ) )
  {
    memcpy( ctuFlagsDst[COMPONENT_Y], ctuFlagsSrc[COMPONENT_Y], sizeof( uint8_t ) * m_numCTUsInPic );
  }
  else
  {
    memcpy( ctuFlagsDst[COMPONENT_Cb], ctuFlagsSrc[COMPONENT_Cb], sizeof( uint8_t ) * m_numCTUsInPic );
    memcpy( ctuFlagsDst[COMPONENT_Cr], ctuFlagsSrc[COMPONENT_Cr], sizeof( uint8_t ) * m_numCTUsInPic );
  }
}

void EncAdaptiveLoopFilter::setCtuEnableFlag( uint8_t** ctuFlags, ChannelType channel, uint8_t val )
{
  if( isLuma( channel ) )
  {
    memset( ctuFlags[COMPONENT_Y], val, sizeof( uint8_t ) * m_numCTUsInPic );
  }
  else
  {
    memset( ctuFlags[COMPONENT_Cb], val, sizeof( uint8_t ) * m_numCTUsInPic );
    memset( ctuFlags[COMPONENT_Cr], val, sizeof( uint8_t ) * m_numCTUsInPic );
  }
}

std::vector<int> EncAdaptiveLoopFilter::getAvaiApsIdsLuma(CodingStructure& cs, int &newApsId)
{
  APS** apss = cs.slice->getAlfAPSs();
  for (int i = 0; i < ALF_CTB_MAX_NUM_APS; i++)
  {
    apss[i] = m_apsMap->getPS((i << NUM_APS_TYPE_LEN) + ALF_APS);
  }

  std::vector<int> result;
  int apsIdChecked = 0, curApsId = m_apsIdStart;
  if (curApsId < ALF_CTB_MAX_NUM_APS)
  {
    while (apsIdChecked < ALF_CTB_MAX_NUM_APS && !cs.slice->isIntra() && result.size() < ALF_CTB_MAX_NUM_APS && !cs.slice->getPendingRasInit() && !cs.slice->isIDRorBLA())
    {
      APS* curAPS = cs.slice->getAlfAPSs()[curApsId];

      if( curAPS && curAPS->getLayerId() == cs.slice->getPic()->layerId && curAPS->getTemporalId() <= cs.slice->getTLayer() && curAPS->getAlfAPSParam().newFilterFlag[CHANNEL_TYPE_LUMA] )
      {
        result.push_back(curApsId);
      }
      apsIdChecked++;
      curApsId = (curApsId + 1) % ALF_CTB_MAX_NUM_APS;
    }
  }
  cs.slice->setTileGroupNumAps((int)result.size());
  cs.slice->setAlfAPSs(result);
  newApsId = m_apsIdStart - 1;
  if (newApsId < 0)
  {
    newApsId = ALF_CTB_MAX_NUM_APS - 1;
  }
  CHECK(newApsId >= ALF_CTB_MAX_NUM_APS, "Wrong APS index assignment in getAvaiApsIdsLuma");
  return result;
}
void  EncAdaptiveLoopFilter::initDistortion(
#if ALF_IMPROVEMENT
  CodingStructure& cs
#endif
)
{
  for (int comp = 0; comp < MAX_NUM_COMPONENT; comp++)
  {
    m_unFiltDistCompnent[comp] = 0.0;

    for (int ctbIdx = 0; ctbIdx < m_numCTUsInPic; ctbIdx++)
    {
#if ALF_IMPROVEMENT
      for (int shapeIdx = 0; shapeIdx < m_filterShapes[toChannelType((ComponentID)comp)].size(); shapeIdx++)
      {
        if (m_filterTypeTest[toChannelType((ComponentID)comp)][m_filterShapes[toChannelType((ComponentID)comp)][shapeIdx].filterType])
        {
          m_ctbDistortionUnfilter[comp][ctbIdx] = getUnfilteredDistortion(m_alfCovariance[comp][shapeIdx][ctbIdx][0], comp == 0 ? MAX_NUM_ALF_CLASSES : 1);
          m_unFiltDistCompnent[comp] += m_ctbDistortionUnfilter[comp][ctbIdx];
          break;
        }
      }
#else
      m_ctbDistortionUnfilter[comp][ctbIdx] = getUnfilteredDistortion(m_alfCovariance[comp][0][ctbIdx], comp == 0 ? MAX_NUM_ALF_CLASSES : 1);
      m_unFiltDistCompnent[comp] += m_ctbDistortionUnfilter[comp][ctbIdx];
#endif
    }
  }

#if ALF_IMPROVEMENT
#if ALF_SAO_TRUE_ORG
  PelUnitBuf orgYuv = cs.getTrueOrgBuf();
#else
  PelUnitBuf orgYuv = cs.getOrgBuf();
#endif
  int orgStride = orgYuv.get(COMPONENT_Y).stride;
  int ctbIdx = 0;
  for (int yPos = 0; yPos < m_picHeight; yPos += m_maxCUHeight)
  {
    for (int xPos = 0; xPos < m_picWidth; xPos += m_maxCUWidth)
    {
      const int width = (xPos + m_maxCUWidth > m_picWidth) ? (m_picWidth - xPos) : m_maxCUWidth;
      const int height = (yPos + m_maxCUHeight > m_picHeight) ? (m_picHeight - yPos) : m_maxCUHeight;
      const Area blk(xPos, yPos, width, height);
      for (int classifierIdx = 0; classifierIdx < NUM_FIXED_FILTER_SETS; classifierIdx++)
      {
        for (int filterSetIdx = 0; filterSetIdx < 2; filterSetIdx++)
        {
          m_ctbDistortionFixedFilter[classifierIdx][filterSetIdx][ctbIdx] = 0;
          Pel *org = orgYuv.get(COMPONENT_Y).bufAt(blk);
          int fixedFilter = classifierIdx * 2 + filterSetIdx;
          for (int y = 0; y < height; y++)
          {
            for (int x = 0; x < width; x++)
            {
              m_ctbDistortionFixedFilter[classifierIdx][filterSetIdx][ctbIdx] += (org[x] - m_fixFilterResult[fixedFilter][y + yPos][x + xPos]) *  (org[x] - m_fixFilterResult[fixedFilter][y + yPos][x + xPos]);
            }
            org += orgStride;
          }
}
      }
      ctbIdx++;
    }
  }
#else
  for( int ctbIdx = 0; ctbIdx < m_numCTUsInPic; ctbIdx++ )
  {
    for( int filterSetIdx = 0; filterSetIdx < NUM_FIXED_FILTER_SETS; filterSetIdx++ )    
    {
      m_ctbDistortionFixedFilter[filterSetIdx][ctbIdx] = m_ctbDistortionUnfilter[COMPONENT_Y][ctbIdx];
      for( int classIdx = 0; classIdx < MAX_NUM_ALF_CLASSES; classIdx++ )
      {
        int filterIdx = m_classToFilterMapping[filterSetIdx][classIdx];
        m_ctbDistortionFixedFilter[filterSetIdx][ctbIdx] += m_alfCovariance[COMPONENT_Y][0][ctbIdx][classIdx].calcErrorForCoeffs( m_clipDefaultEnc, m_fixedFilterSetCoeff[filterIdx], MAX_NUM_ALF_LUMA_COEFF, m_NUM_BITS );
      }
    }
  }  
#endif
}

void  EncAdaptiveLoopFilter::initDistortionCcalf()
{
  for (int comp = 1; comp < MAX_NUM_COMPONENT; comp++)
  {
    for (int ctbIdx = 0; ctbIdx < m_numCTUsInPic; ctbIdx++)
    {
      m_ctbDistortionUnfilter[comp][ctbIdx] = m_alfCovarianceCcAlf[comp - 1][0][0][ctbIdx].pixAcc;
    }
  }
}

void  EncAdaptiveLoopFilter::getDistNewFilter( AlfParam& alfParam )
{
  reconstructCoeff( alfParam, CHANNEL_TYPE_LUMA, true, true );
#if ALF_IMPROVEMENT
  AlfFilterType filterTypeCtb = alfParam.filterType[CHANNEL_TYPE_LUMA];
  int numFixedFilterSet = (filterTypeCtb == ALF_FILTER_EXT || filterTypeCtb == ALF_FILTER_9_EXT) ? 2 : 1;
  for( int altIdx = 0; altIdx < alfParam.numAlternativesLuma; altIdx++ )
  {
    for( int classIdx = 0; classIdx < MAX_NUM_ALF_CLASSES; classIdx++ )
    {
      for( int coeff = 0; coeff < MAX_NUM_ALF_LUMA_COEFF; coeff++ )
      {
        m_filterTmp[coeff] = m_coeffFinal[altIdx][classIdx * MAX_NUM_ALF_LUMA_COEFF + coeff];
        m_clipTmp[coeff]   = m_clippFinal[altIdx][classIdx * MAX_NUM_ALF_LUMA_COEFF + coeff];
      }
      for( int ctbIdx = 0; ctbIdx < m_numCTUsInPic; ctbIdx++ )
      {
        for( int fixedFilterSetIdx = 0; fixedFilterSetIdx < numFixedFilterSet; fixedFilterSetIdx++ )
        {
          if (classIdx == 0)
          {
            m_distCtbLumaNewFilt[altIdx][fixedFilterSetIdx][ctbIdx] = m_ctbDistortionUnfilter[COMPONENT_Y][ctbIdx];
          }
          m_distCtbLumaNewFilt[altIdx][fixedFilterSetIdx][ctbIdx] += m_alfCovariance[COMPONENT_Y][m_filterTypeToStatIndex[CHANNEL_TYPE_LUMA][filterTypeCtb]][ctbIdx][fixedFilterSetIdx][classIdx].calcErrorForCoeffs(m_clipTmp, m_filterTmp, m_filterShapes[CHANNEL_TYPE_LUMA][m_filterTypeToStatIndex[CHANNEL_TYPE_LUMA][filterTypeCtb]].numCoeff, m_NUM_BITS);
        }
      }
    }
  }
#else
  for( int classIdx = 0; classIdx < MAX_NUM_ALF_CLASSES; classIdx++ )
  {
    for( int coeff = 0; coeff < MAX_NUM_ALF_LUMA_COEFF; coeff++ )
    {
      m_filterTmp[coeff] = m_coeffFinal[classIdx * MAX_NUM_ALF_LUMA_COEFF + coeff];
      m_clipTmp[coeff] = m_clippFinal[classIdx * MAX_NUM_ALF_LUMA_COEFF + coeff];
    }
    for( int ctbIdx = 0; ctbIdx < m_numCTUsInPic; ctbIdx++ )
    {
      if( classIdx == 0 )
      {
        m_distCtbLumaNewFilt[ctbIdx] = m_ctbDistortionUnfilter[COMPONENT_Y][ctbIdx];
      }
      m_distCtbLumaNewFilt[ctbIdx] += m_alfCovariance[COMPONENT_Y][0][ctbIdx][classIdx].calcErrorForCoeffs(m_clipTmp, m_filterTmp, MAX_NUM_ALF_LUMA_COEFF, m_NUM_BITS);
    }
  }
#endif
}
void  EncAdaptiveLoopFilter::getDistApsFilter( CodingStructure& cs, std::vector<int> apsIds )
{
  APS** aps = cs.slice->getAlfAPSs();
  AlfParam alfParamTmp;
  APS* curAPS;
  for( int apsIdx = 0; apsIdx < apsIds.size(); apsIdx++ )
  {
    int apsId = apsIds[apsIdx];
    curAPS = aps[apsId];
    alfParamTmp = curAPS->getAlfAPSParam();
    reconstructCoeff( alfParamTmp, CHANNEL_TYPE_LUMA, true, true );
#if ALF_IMPROVEMENT
    AlfFilterType filterTypeCtb = alfParamTmp.filterType[CHANNEL_TYPE_LUMA];
    m_filterTypeApsLuma[apsIdx] = filterTypeCtb;
    m_numLumaAltAps[apsIdx] = alfParamTmp.numAlternativesLuma;
    int numFixedFilterSet = ( filterTypeCtb == ALF_FILTER_EXT || filterTypeCtb == ALF_FILTER_9_EXT ) ? 2 : 1;
    for( int altIdx = 0; altIdx < alfParamTmp.numAlternativesLuma; altIdx++ )
    {
      for( int classIdx = 0; classIdx < MAX_NUM_ALF_CLASSES; classIdx++ )
      {
        for ( int coeff = 0; coeff < MAX_NUM_ALF_LUMA_COEFF; coeff++ )
        {
          m_filterTmp[coeff] = m_coeffFinal[altIdx][classIdx * MAX_NUM_ALF_LUMA_COEFF + coeff];
          m_clipTmp[coeff]   = m_clippFinal[altIdx][classIdx * MAX_NUM_ALF_LUMA_COEFF + coeff];
        }
        for( int ctbIdx = 0; ctbIdx < m_numCTUsInPic; ctbIdx++ )
        {
          for( int fixedFilterSetIdx = 0; fixedFilterSetIdx < numFixedFilterSet; fixedFilterSetIdx++ )
          {
            if( classIdx == 0 )
            {
              m_distCtbApsLuma[apsIdx][altIdx][fixedFilterSetIdx][ctbIdx] = m_ctbDistortionUnfilter[COMPONENT_Y][ctbIdx];
            }
            m_distCtbApsLuma[apsIdx][altIdx][fixedFilterSetIdx][ctbIdx] += m_alfCovariance[COMPONENT_Y][m_filterTypeToStatIndex[CHANNEL_TYPE_LUMA][filterTypeCtb]][ctbIdx][fixedFilterSetIdx][classIdx].calcErrorForCoeffs(m_clipTmp, m_filterTmp, m_filterShapes[CHANNEL_TYPE_LUMA][m_filterTypeToStatIndex[CHANNEL_TYPE_LUMA][filterTypeCtb]].numCoeff, m_NUM_BITS);
          }          
        }
      }
    }
#else
    for( int classIdx = 0; classIdx < MAX_NUM_ALF_CLASSES; classIdx++ )
    {
      for( int coeff = 0; coeff < MAX_NUM_ALF_LUMA_COEFF; coeff++ )
      {
        m_filterTmp[coeff] = m_coeffFinal[classIdx * MAX_NUM_ALF_LUMA_COEFF + coeff];
        m_clipTmp[coeff] = m_clippFinal[classIdx * MAX_NUM_ALF_LUMA_COEFF + coeff];
      }
      for( int ctbIdx = 0; ctbIdx < m_numCTUsInPic; ctbIdx++ )
      {
        if( classIdx == 0 )
        {
          m_distCtbApsLuma[apsIdx][ctbIdx] = m_ctbDistortionUnfilter[COMPONENT_Y][ctbIdx];
        }
        m_distCtbApsLuma[apsIdx][ctbIdx] += m_alfCovariance[COMPONENT_Y][0][ctbIdx][classIdx].calcErrorForCoeffs(m_clipTmp, m_filterTmp, MAX_NUM_ALF_LUMA_COEFF, m_NUM_BITS);
      }
    }
#endif
  }
}

void  EncAdaptiveLoopFilter::alfEncoderCtb(CodingStructure& cs, AlfParam& alfParamNewFilters
#if ENABLE_QPA
  , const double lambdaChromaWeight
#endif
)
{
  TempCtx        ctxStart(m_CtxCache, AlfCtx(m_CABACEstimator->getCtx()));
  TempCtx        ctxBest(m_CtxCache);
  TempCtx        ctxTempStart(m_CtxCache);
  TempCtx        ctxTempBest(m_CtxCache);
  TempCtx        ctxTempAltStart(m_CtxCache);
  TempCtx        ctxTempAltBest(m_CtxCache);
  AlfParam  alfParamNewFiltersBest = alfParamNewFilters;
  bool     hasNewFilters[2] = { alfParamNewFilters.enabledFlag[COMPONENT_Y] , alfParamNewFilters.enabledFlag[COMPONENT_Cb] || alfParamNewFilters.enabledFlag[COMPONENT_Cr] };
  m_alfParamTemp = alfParamNewFilters;
  APS**          apss = cs.slice->getAlfAPSs();
  short*     alfCtbFilterSetIndex = cs.picture->getAlfCtbFilterIndex();

  //luma
  setCtuEnableFlag(m_ctuEnableFlag, CHANNEL_TYPE_LUMA, 1);
#if ALF_IMPROVEMENT
  std::fill_n(m_ctuAlternative[COMPONENT_Y], m_numCTUsInPic, 0);
  if (m_alfParamTemp.numAlternativesLuma < 1)
  {
    m_alfParamTemp.numAlternativesLuma = 1;
  }
#endif
  double costOff = m_unFiltDistCompnent[COMPONENT_Y];
  setCtuEnableFlag(m_ctuEnableFlag, CHANNEL_TYPE_LUMA, 0);
  int newApsId;
  std::vector<int> apsIds = getAvaiApsIdsLuma(cs, newApsId);
  std::vector<int> bestApsIds;
  double costMin = MAX_DOUBLE;

  getDistApsFilter( cs, apsIds );
  int numLoops = hasNewFilters[CHANNEL_TYPE_LUMA] ? 2 : 1;
#if ALF_IMPROVEMENT
  AlfFilterType filterTypeNewFilter = ALF_NUM_OF_FILTER_TYPES;
  int numAlfLumaNew = alfParamNewFilters.numAlternativesLuma;
  int bestFixedFilterSetIdx = -1;
  for (int fixedFilterSetIdx = 0; fixedFilterSetIdx < 2; fixedFilterSetIdx++)
  {
#endif
    for (int useNewFilter = 0; useNewFilter < numLoops; useNewFilter++)
    {
      int bitsNewFilter = 0;
      if (useNewFilter == 1)
      {
        if (!hasNewFilters[CHANNEL_TYPE_LUMA])
        {
          continue;
        }
        else
        {
          bitsNewFilter = m_bitsNewFilter[CHANNEL_TYPE_LUMA];
          getDistNewFilter(alfParamNewFilters);
#if ALF_IMPROVEMENT
          filterTypeNewFilter = alfParamNewFilters.filterType[CHANNEL_TYPE_LUMA];
          numAlfLumaNew = alfParamNewFilters.numAlternativesLuma;
#endif
        }
      }
      int numIter = useNewFilter ? 2 : 1;
      for (int numTemporalAps = 0; numTemporalAps <= apsIds.size(); numTemporalAps++)
      {
        if (numTemporalAps + useNewFilter >= ALF_CTB_MAX_NUM_APS)
        {
          continue;
        }
        cs.slice->setTileGroupNumAps(numTemporalAps + useNewFilter);
        int numFilterSet = NUM_FIXED_FILTER_SETS + numTemporalAps + useNewFilter;
        if (numTemporalAps == apsIds.size() && numTemporalAps > 0 && useNewFilter && newApsId == apsIds.back()) //last temporalAPS is occupied by new filter set and this temporal APS becomes unavailable
        {
          continue;
        }
        for (int iter = 0; iter < numIter; iter++)
        {
          m_alfParamTemp = alfParamNewFilters;
          m_alfParamTemp.enabledFlag[CHANNEL_TYPE_LUMA] = true;
          double curCost = 3 * m_lambda[CHANNEL_TYPE_LUMA];
          if (iter > 0)  //re-derive new filter-set
          {
            double dDistOrgNewFilter = 0;
            int blocksUsingNewFilter = 0;
            for (int ctbIdx = 0; ctbIdx < m_numCTUsInPic; ctbIdx++)
            {
              if (m_ctuEnableFlag[COMPONENT_Y][ctbIdx] && alfCtbFilterSetIndex[ctbIdx] != NUM_FIXED_FILTER_SETS)
              {
                m_ctuEnableFlag[COMPONENT_Y][ctbIdx] = 0;
              }
              else if (m_ctuEnableFlag[COMPONENT_Y][ctbIdx] && alfCtbFilterSetIndex[ctbIdx] == NUM_FIXED_FILTER_SETS)
              {
                blocksUsingNewFilter++;
#if ALF_IMPROVEMENT
                int altIdx = m_ctuAlternative[COMPONENT_Y][ctbIdx];
#endif
                dDistOrgNewFilter += m_ctbDistortionUnfilter[COMPONENT_Y][ctbIdx];
                for (int classIdx = 0; classIdx < MAX_NUM_ALF_CLASSES; classIdx++)
                {
#if ALF_IMPROVEMENT
                  short* pCoeff = m_coeffFinal[altIdx];
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
                  Pel* pClipp = m_clippFinal[altIdx];
#else
                  short* pClipp = m_clippFinal[altIdx];
#endif
#else
                  short* pCoeff = m_coeffFinal;
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
                  Pel* pClipp = m_clippFinal;
#else
                  short* pClipp = m_clippFinal;
#endif
#endif
                  for (int i = 0; i < MAX_NUM_ALF_LUMA_COEFF; i++)
                  {
                    m_filterTmp[i] = pCoeff[classIdx * MAX_NUM_ALF_LUMA_COEFF + i];
                    m_clipTmp[i] = pClipp[classIdx * MAX_NUM_ALF_LUMA_COEFF + i];
                  }
#if ALF_IMPROVEMENT
                  dDistOrgNewFilter += m_alfCovariance[COMPONENT_Y][m_filterTypeToStatIndex[CHANNEL_TYPE_LUMA][filterTypeNewFilter]][ctbIdx][fixedFilterSetIdx][classIdx].calcErrorForCoeffs(m_clipTmp, m_filterTmp, m_filterShapes[CHANNEL_TYPE_LUMA][m_filterTypeToStatIndex[CHANNEL_TYPE_LUMA][filterTypeNewFilter]].numCoeff, m_NUM_BITS);
#else
                  dDistOrgNewFilter += m_alfCovariance[COMPONENT_Y][0][ctbIdx][classIdx].calcErrorForCoeffs(m_clipTmp, m_filterTmp, MAX_NUM_ALF_LUMA_COEFF, m_NUM_BITS);
#endif
                }
              }
            }
            if (blocksUsingNewFilter > 0 && blocksUsingNewFilter < m_numCTUsInPic)
            {
#if ALF_IMPROVEMENT
              int bitsNewFilterTempLuma = 0, bitsTemp = 0;
              double err = 0.0;
              double costNewFilter = MAX_DOUBLE;
              for (int shapeIdx = 0; shapeIdx < m_filterShapes[CHANNEL_TYPE_LUMA].size(); shapeIdx++)
              {
                if (m_filterTypeTest[CHANNEL_TYPE_LUMA][m_filterShapes[CHANNEL_TYPE_LUMA][shapeIdx].filterType] == false)
                {
                  continue;
                }
                m_alfParamTemp.filterType[CHANNEL_TYPE_LUMA] = m_filterShapes[CHANNEL_TYPE_LUMA][shapeIdx].filterType;
                err = getFilterCoeffAndCost(cs, 0, CHANNEL_TYPE_LUMA, true, shapeIdx, bitsTemp, fixedFilterSetIdx, true);
                if (err < costNewFilter)
                {
                  costNewFilter = err;
                  bitsNewFilterTempLuma = bitsTemp;
                  m_alfParamTempNL = m_alfParamTemp;
                }
              }
              err = costNewFilter;
              m_alfParamTemp = m_alfParamTempNL;
#else
              int bitNL[2] = { 0, 0 };
              double errNL[2] = { 0.0, 0.0 };
              m_alfParamTemp.nonLinearFlag[CHANNEL_TYPE_LUMA] = 1;
              if (m_encCfg->getUseNonLinearAlfLuma())
              {
                errNL[1] = getFilterCoeffAndCost(cs, 0, CHANNEL_TYPE_LUMA, true, 0, bitNL[1], true);
                m_alfParamTempNL = m_alfParamTemp;
              }
              else
              {
                errNL[1] = MAX_DOUBLE;
              }
              m_alfParamTemp.nonLinearFlag[CHANNEL_TYPE_LUMA] = 0;
              errNL[0] = getFilterCoeffAndCost(cs, 0, CHANNEL_TYPE_LUMA, true, 0, bitNL[0], true);

              int bitsNewFilterTempLuma = bitNL[0];
              double err = errNL[0];
              if (errNL[1] < errNL[0])
              {
                err = errNL[1];
                bitsNewFilterTempLuma = bitNL[1];
                m_alfParamTemp = m_alfParamTempNL;
              }
#endif
              if (dDistOrgNewFilter + m_lambda[CHANNEL_TYPE_LUMA] * m_bitsNewFilter[CHANNEL_TYPE_LUMA] < err) //re-derived filter is not good, skip
              {
                continue;
              }

              getDistNewFilter(m_alfParamTemp);
              bitsNewFilter = bitsNewFilterTempLuma;
#if ALF_IMPROVEMENT
              filterTypeNewFilter = m_alfParamTemp.filterType[CHANNEL_TYPE_LUMA];
              numAlfLumaNew = m_alfParamTemp.numAlternativesLuma;
#endif
            }
            else //no blocks using new filter, skip
            {
              continue;
            }
          }

          m_CABACEstimator->getCtx() = ctxStart;

          for (int ctbIdx = 0; ctbIdx < m_numCTUsInPic; ctbIdx++)
          {
            double distUnfilterCtb = m_ctbDistortionUnfilter[COMPONENT_Y][ctbIdx];

            //ctb on
            m_ctuEnableFlag[COMPONENT_Y][ctbIdx] = 1;
            double         costOn = MAX_DOUBLE;
            ctxTempStart = AlfCtx(m_CABACEstimator->getCtx());
            int iBestFilterSetIdx = 0;
#if ALF_IMPROVEMENT
            int iBestAltIdx = 0;
#endif

            for (int filterSetIdx = 0; filterSetIdx < numFilterSet; filterSetIdx++)  // to select best APS index / fixed filter index
            {
#if ALF_IMPROVEMENT
              int iterAltLuma = 1;
              if (filterSetIdx >= NUM_FIXED_FILTER_SETS)
              {
                if (useNewFilter && filterSetIdx == NUM_FIXED_FILTER_SETS)
                {
                  iterAltLuma = numAlfLumaNew;
                }
                else if (useNewFilter)
                {
                  iterAltLuma = m_numLumaAltAps[filterSetIdx - 1 - NUM_FIXED_FILTER_SETS];
                }
                else
                {
                  iterAltLuma = m_numLumaAltAps[filterSetIdx - NUM_FIXED_FILTER_SETS];
                }
              }
              for (int altIdx = 0; altIdx < iterAltLuma; altIdx++)     // to select best filter set index inside one APS
              {
#endif
                //rate
                m_CABACEstimator->getCtx() = AlfCtx(ctxTempStart);
                m_CABACEstimator->resetBits();
                m_CABACEstimator->codeAlfCtuEnableFlag(cs, ctbIdx, COMPONENT_Y, &m_alfParamTemp);
                alfCtbFilterSetIndex[ctbIdx] = filterSetIdx;
                m_CABACEstimator->codeAlfCtuFilterIndex(cs, ctbIdx, &m_alfParamTemp.enabledFlag[COMPONENT_Y]);
#if ALF_IMPROVEMENT
                m_ctuAlternative[COMPONENT_Y][ctbIdx] = altIdx;
                if (filterSetIdx >= NUM_FIXED_FILTER_SETS)
                {
                  m_CABACEstimator->codeAlfCtuAlternative(cs, ctbIdx, COMPONENT_Y, &m_alfParamTemp, iterAltLuma);
                }
#endif
                double rateOn = FRAC_BITS_SCALE * m_CABACEstimator->getEstFracBits();
                //distortion
                double dist;
#if ALF_IMPROVEMENT
                if (filterSetIdx < NUM_FIXED_FILTER_SETS)
                {
                  dist = m_ctbDistortionFixedFilter[filterSetIdx][fixedFilterSetIdx][ctbIdx];
                }
                else if (useNewFilter && filterSetIdx == NUM_FIXED_FILTER_SETS)
                {
                  dist = m_distCtbLumaNewFilt[altIdx][fixedFilterSetIdx][ctbIdx];
                }
                else
                {
                  dist = m_distCtbApsLuma[filterSetIdx - NUM_FIXED_FILTER_SETS - useNewFilter][altIdx][fixedFilterSetIdx][ctbIdx];
                }
#else
                if (filterSetIdx < NUM_FIXED_FILTER_SETS)
                {
                  dist = m_ctbDistortionFixedFilter[filterSetIdx][ctbIdx];
                }
                else if (useNewFilter && filterSetIdx == NUM_FIXED_FILTER_SETS)
                {
                  dist = m_distCtbLumaNewFilt[ctbIdx];
                }
                else
                {
                  dist = m_distCtbApsLuma[filterSetIdx - NUM_FIXED_FILTER_SETS - useNewFilter][ctbIdx];
                }
#endif

                //cost
                double costOnTmp = dist + m_lambda[COMPONENT_Y] * rateOn;
                if (costOnTmp < costOn)
                {
                  ctxTempBest = AlfCtx(m_CABACEstimator->getCtx());
                  costOn = costOnTmp;
                  iBestFilterSetIdx = filterSetIdx;
#if ALF_IMPROVEMENT
                  iBestAltIdx = altIdx;
#endif
                }
#if ALF_IMPROVEMENT
              }
#endif
            }

            //ctb off
            m_ctuEnableFlag[COMPONENT_Y][ctbIdx] = 0;
            //rate
            m_CABACEstimator->getCtx() = AlfCtx(ctxTempStart);
            m_CABACEstimator->resetBits();
            m_CABACEstimator->codeAlfCtuEnableFlag(cs, ctbIdx, COMPONENT_Y, &m_alfParamTemp);
            //cost
            double costOff = distUnfilterCtb + m_lambda[COMPONENT_Y] * FRAC_BITS_SCALE * m_CABACEstimator->getEstFracBits();

            if (costOn < costOff)
            {
              m_CABACEstimator->getCtx() = AlfCtx(ctxTempBest);
              m_ctuEnableFlag[COMPONENT_Y][ctbIdx] = 1;
              alfCtbFilterSetIndex[ctbIdx] = iBestFilterSetIdx;
#if ALF_IMPROVEMENT
              m_ctuAlternative[COMPONENT_Y][ctbIdx] = iBestAltIdx;
#endif
              curCost += costOn;
            }
            else
            {
              m_ctuEnableFlag[COMPONENT_Y][ctbIdx] = 0;
              curCost += costOff;
            }
          } //for(ctbIdx)
#if ALF_IMPROVEMENT
          int tmpBits = bitsNewFilter + 3 * (numFilterSet - NUM_FIXED_FILTER_SETS) + 1;
#else
          int tmpBits = bitsNewFilter + 3 * (numFilterSet - NUM_FIXED_FILTER_SETS);
#endif
          curCost += tmpBits * m_lambda[COMPONENT_Y];
          if (curCost < costMin)
          {
            costMin = curCost;
            bestApsIds.resize(numFilterSet - NUM_FIXED_FILTER_SETS);
            for (int i = 0; i < bestApsIds.size(); i++)
            {
              if (i == 0 && useNewFilter)
              {
                bestApsIds[i] = newApsId;
              }
              else
              {
                bestApsIds[i] = apsIds[i - useNewFilter];
              }
            }
#if ALF_IMPROVEMENT
            bestFixedFilterSetIdx = fixedFilterSetIdx;
#endif
            alfParamNewFiltersBest = m_alfParamTemp;
            ctxBest = AlfCtx(m_CABACEstimator->getCtx());
            copyCtuEnableFlag(m_ctuEnableFlagTmp, m_ctuEnableFlag, CHANNEL_TYPE_LUMA);
            for (int ctuIdx = 0; ctuIdx < m_numCTUsInPic; ctuIdx++)
            {
              m_alfCtbFilterSetIndexTmp[ctuIdx] = alfCtbFilterSetIndex[ctuIdx];
#if ALF_IMPROVEMENT
              m_ctuAlternativeTmp[COMPONENT_Y][ctuIdx] = m_ctuAlternative[COMPONENT_Y][ctuIdx];
#endif
            }
            alfParamNewFiltersBest.newFilterFlag[CHANNEL_TYPE_LUMA] = useNewFilter;
          }
        }//for (int iter = 0; iter < numIter; iter++)
      }// for (int numTemporalAps = 0; numTemporalAps < apsIds.size(); numTemporalAps++)
    }//for (int useNewFilter = 0; useNewFilter <= 1; useNewFilter++)
#if ALF_IMPROVEMENT
  }//for(fixedFilterSetIdx)
#endif

  cs.slice->setTileGroupCcAlfCbApsId(newApsId);
  cs.slice->setTileGroupCcAlfCrApsId(newApsId);

  if (costOff <= costMin)
  {
    cs.slice->resetTileGroupAlfEnabledFlag();
    cs.slice->setTileGroupNumAps(0);
    setCtuEnableFlag(m_ctuEnableFlag, CHANNEL_TYPE_LUMA, 0);
    setCtuEnableFlag(m_ctuEnableFlag, CHANNEL_TYPE_CHROMA, 0);
    return;
  }
  else
  {
#if ALF_IMPROVEMENT
    cs.slice->setTileGroupAlfFixedFilterSetIdx( bestFixedFilterSetIdx );
#endif
    cs.slice->setTileGroupAlfEnabledFlag(COMPONENT_Y, true);
    cs.slice->setTileGroupNumAps((int)bestApsIds.size());
    cs.slice->setAlfAPSs(bestApsIds);
    copyCtuEnableFlag(m_ctuEnableFlag, m_ctuEnableFlagTmp, CHANNEL_TYPE_LUMA);
    for (int ctuIdx = 0; ctuIdx < m_numCTUsInPic; ctuIdx++)
    {
      alfCtbFilterSetIndex[ctuIdx] = m_alfCtbFilterSetIndexTmp[ctuIdx];
#if ALF_IMPROVEMENT
      m_ctuAlternative[COMPONENT_Y][ctuIdx] = m_ctuAlternativeTmp[COMPONENT_Y][ctuIdx];
#endif
    }
    if (alfParamNewFiltersBest.newFilterFlag[CHANNEL_TYPE_LUMA])
    {
      APS* newAPS = m_apsMap->getPS((newApsId << NUM_APS_TYPE_LEN) + ALF_APS);
      if (newAPS == NULL)
      {
        newAPS = m_apsMap->allocatePS((newApsId << NUM_APS_TYPE_LEN) + ALF_APS);
        newAPS->setAPSId(newApsId);
        newAPS->setAPSType(ALF_APS);
      }
      newAPS->setAlfAPSParam(alfParamNewFiltersBest);
      newAPS->setTemporalId( cs.slice->getTLayer() );
      newAPS->getAlfAPSParam().newFilterFlag[CHANNEL_TYPE_CHROMA] = false;
      m_apsMap->setChangedFlag((newApsId << NUM_APS_TYPE_LEN) + ALF_APS);
      m_apsIdStart = newApsId;
    }

    std::vector<int> apsIds = cs.slice->getTileGroupApsIdLuma();
    for (int i = 0; i < (int)cs.slice->getTileGroupNumAps(); i++)
    {
      apss[apsIds[i]] = m_apsMap->getPS((apsIds[i] << NUM_APS_TYPE_LEN) + ALF_APS);
    }
  }

  //chroma
  if( isChromaEnabled( cs.pcv->chrFormat ) )
  {
    m_alfParamTemp = alfParamNewFiltersBest;
    if( m_alfParamTemp.numAlternativesChroma < 1 )
    {
      m_alfParamTemp.numAlternativesChroma = 1;
    }
    setCtuAlternativeChroma( m_ctuAlternative, 0 );
    setCtuEnableFlag( m_ctuEnableFlag, CHANNEL_TYPE_CHROMA, 1 );
    costOff = m_unFiltDistCompnent[COMPONENT_Cb] + m_unFiltDistCompnent[COMPONENT_Cr];
    costMin = MAX_DOUBLE;
    m_CABACEstimator->getCtx() = AlfCtx( ctxBest );
    ctxStart = AlfCtx( m_CABACEstimator->getCtx() );
    int newApsIdChroma = -1;
    if( alfParamNewFiltersBest.newFilterFlag[CHANNEL_TYPE_LUMA] && ( alfParamNewFiltersBest.enabledFlag[COMPONENT_Cb] || alfParamNewFiltersBest.enabledFlag[COMPONENT_Cr] ) )
    {
      newApsIdChroma = newApsId;
    }
    else if( alfParamNewFiltersBest.enabledFlag[COMPONENT_Cb] || alfParamNewFiltersBest.enabledFlag[COMPONENT_Cr] )
    {
      int curId = m_apsIdStart;
      while( newApsIdChroma < 0 )
      {
        curId--;
        if( curId < 0 )
        {
          curId = ALF_CTB_MAX_NUM_APS - 1;
        }
        if( std::find( bestApsIds.begin(), bestApsIds.end(), curId ) == bestApsIds.end() )
        {
          newApsIdChroma = curId;
        }
      }
    }

    for( int curApsId = 0; curApsId < ALF_CTB_MAX_NUM_APS; curApsId++ )
    {
      if( ( cs.slice->getPendingRasInit() || cs.slice->isIDRorBLA() || cs.slice->isIntra() )
          && curApsId != newApsIdChroma )
      {
        continue;
      }
      APS *curAPS = m_apsMap->getPS( ( curApsId << NUM_APS_TYPE_LEN ) + ALF_APS );

      if( curAPS && curAPS->getLayerId() != cs.slice->getPic()->layerId )
      {
        continue;
      }

      double curCost = m_lambda[CHANNEL_TYPE_CHROMA] * 3;
      if( curApsId == newApsIdChroma )
      {
        m_alfParamTemp = alfParamNewFilters;
        curCost += m_lambda[CHANNEL_TYPE_CHROMA] * m_bitsNewFilter[CHANNEL_TYPE_CHROMA];
      }
      else if( curAPS && curAPS->getTemporalId() <= cs.slice->getTLayer() && curAPS->getAlfAPSParam().newFilterFlag[CHANNEL_TYPE_CHROMA] )
      {
        m_alfParamTemp = curAPS->getAlfAPSParam();
      }
      else
      {
        continue;
      }
      reconstructCoeff( m_alfParamTemp, CHANNEL_TYPE_CHROMA, true, true );
#if ALF_IMPROVEMENT
      AlfFilterType filterTypeChroma = m_alfParamTemp.filterType[CHANNEL_TYPE_CHROMA];
#endif
      m_CABACEstimator->getCtx() = AlfCtx( ctxStart );
      for( int compId = 1; compId < MAX_NUM_COMPONENT; compId++ )
      {
        m_alfParamTemp.enabledFlag[compId] = true;
        for( int ctbIdx = 0; ctbIdx < m_numCTUsInPic; ctbIdx++ )
        {
          double distUnfilterCtu = m_ctbDistortionUnfilter[compId][ctbIdx];
          //cost on
          m_ctuEnableFlag[compId][ctbIdx] = 1;
          ctxTempStart = AlfCtx( m_CABACEstimator->getCtx() );
          //rate
          m_CABACEstimator->getCtx() = AlfCtx( ctxTempStart );
          m_CABACEstimator->resetBits();
          //ctb flag
          m_CABACEstimator->codeAlfCtuEnableFlag( cs, ctbIdx, compId, &m_alfParamTemp );
          double rateOn = FRAC_BITS_SCALE * m_CABACEstimator->getEstFracBits();
#if ENABLE_QPA
          const double ctuLambda =
            lambdaChromaWeight > 0.0 ? cs.picture->m_uEnerHpCtu[ctbIdx] / lambdaChromaWeight : m_lambda[compId];
#else
          const double ctuLambda = m_lambda[compId];
#endif
          double dist = MAX_DOUBLE;
          int    numAlts = m_alfParamTemp.numAlternativesChroma;
          ctxTempBest = AlfCtx( m_CABACEstimator->getCtx() );
          double bestAltRate = 0;
          double bestAltCost = MAX_DOUBLE;
          int    bestAltIdx = -1;
          ctxTempAltStart = AlfCtx( ctxTempBest );
          for( int altIdx = 0; altIdx < numAlts; ++altIdx )
          {
            if( altIdx )
            {
              m_CABACEstimator->getCtx() = AlfCtx( ctxTempAltStart );
            }
            m_CABACEstimator->resetBits();
            m_ctuAlternative[compId][ctbIdx] = altIdx;
            m_CABACEstimator->codeAlfCtuAlternative( cs, ctbIdx, compId, &m_alfParamTemp );
            double altRate = FRAC_BITS_SCALE * m_CABACEstimator->getEstFracBits();
            double r_altCost = ctuLambda * altRate;

            // distortion
            for( int i = 0; i < MAX_NUM_ALF_CHROMA_COEFF; i++ )
            {
              m_filterTmp[i] = m_chromaCoeffFinal[altIdx][i];
              m_clipTmp[i] = m_chromaClippFinal[altIdx][i];
            }
#if ALF_IMPROVEMENT
            double altDist = m_alfCovariance[compId][m_filterTypeToStatIndex[CHANNEL_TYPE_CHROMA][filterTypeChroma]][ctbIdx][0][0].calcErrorForCoeffs( m_clipTmp, m_filterTmp, m_filterShapes[CHANNEL_TYPE_CHROMA][m_filterTypeToStatIndex[CHANNEL_TYPE_CHROMA][filterTypeChroma]].numCoeff, m_NUM_BITS );
#else            
            double altDist = m_alfCovariance[compId][0][ctbIdx][0].calcErrorForCoeffs( m_clipTmp, m_filterTmp, MAX_NUM_ALF_CHROMA_COEFF, m_NUM_BITS );
#endif
            double altCost = altDist + r_altCost;
            if( altCost < bestAltCost )
            {
              bestAltCost = altCost;
              bestAltIdx = altIdx;
              bestAltRate = altRate;
              ctxTempBest = AlfCtx( m_CABACEstimator->getCtx() );
              dist = altDist;
            }
          }
          m_ctuAlternative[compId][ctbIdx] = bestAltIdx;
          rateOn += bestAltRate;
          dist += distUnfilterCtu;
          // cost
          double costOn = dist + ctuLambda * rateOn;
          // cost off
          m_ctuEnableFlag[compId][ctbIdx] = 0;
          // rate
          m_CABACEstimator->getCtx() = AlfCtx( ctxTempStart );
          m_CABACEstimator->resetBits();
          m_CABACEstimator->codeAlfCtuEnableFlag( cs, ctbIdx, compId, &m_alfParamTemp );
          // cost
          double costOff = distUnfilterCtu + m_lambda[compId] * FRAC_BITS_SCALE * m_CABACEstimator->getEstFracBits();
          if( costOn < costOff )
          {
            m_CABACEstimator->getCtx() = AlfCtx( ctxTempBest );
            m_ctuEnableFlag[compId][ctbIdx] = 1;
            curCost += costOn;
          }
          else
          {
            m_ctuEnableFlag[compId][ctbIdx] = 0;
            curCost += costOff;
          }
        }
      }
      // chroma idc
      setEnableFlag( m_alfParamTemp, CHANNEL_TYPE_CHROMA, m_ctuEnableFlag );

      if( curCost < costMin )
      {
        costMin = curCost;
        cs.slice->setTileGroupApsIdChroma( curApsId );
        cs.slice->setTileGroupAlfEnabledFlag( COMPONENT_Cb, m_alfParamTemp.enabledFlag[COMPONENT_Cb] );
        cs.slice->setTileGroupAlfEnabledFlag( COMPONENT_Cr, m_alfParamTemp.enabledFlag[COMPONENT_Cr] );
        copyCtuEnableFlag( m_ctuEnableFlagTmp, m_ctuEnableFlag, CHANNEL_TYPE_CHROMA );
        copyCtuAlternativeChroma( m_ctuAlternativeTmp, m_ctuAlternative );
      }
    }

    if( newApsIdChroma >= 0 )
    {
      cs.slice->setTileGroupCcAlfCbApsId( newApsIdChroma );
      cs.slice->setTileGroupCcAlfCrApsId( newApsIdChroma );
    }
    if( costOff < costMin )
    {
      cs.slice->setTileGroupAlfEnabledFlag( COMPONENT_Cb, false );
      cs.slice->setTileGroupAlfEnabledFlag( COMPONENT_Cr, false );
      setCtuEnableFlag( m_ctuEnableFlag, CHANNEL_TYPE_CHROMA, 0 );
    }
    else
    {
      copyCtuEnableFlag( m_ctuEnableFlag, m_ctuEnableFlagTmp, CHANNEL_TYPE_CHROMA );
      copyCtuAlternativeChroma( m_ctuAlternative, m_ctuAlternativeTmp );
      if( cs.slice->getTileGroupApsIdChroma() == newApsIdChroma )  //new filter
      {
        APS* newAPS = m_apsMap->getPS( ( newApsIdChroma << NUM_APS_TYPE_LEN ) + ALF_APS );
        if( newAPS == NULL )
        {
          newAPS = m_apsMap->allocatePS( ( newApsIdChroma << NUM_APS_TYPE_LEN ) + ALF_APS );
          newAPS->setAPSType( ALF_APS );
          newAPS->setAPSId( newApsIdChroma );
          newAPS->getAlfAPSParam().reset();
        }
        newAPS->getAlfAPSParam().newFilterFlag[CHANNEL_TYPE_CHROMA] = true;
        if( !alfParamNewFiltersBest.newFilterFlag[CHANNEL_TYPE_LUMA] )
        {
          newAPS->getAlfAPSParam().newFilterFlag[CHANNEL_TYPE_LUMA] = false;
        }
        newAPS->getAlfAPSParam().numAlternativesChroma = alfParamNewFilters.numAlternativesChroma;
#if ALF_IMPROVEMENT
        for( int altIdx = 0; altIdx < MAX_NUM_ALF_ALTERNATIVES_CHROMA; ++altIdx )
        {
          newAPS->getAlfAPSParam().nonLinearFlag[CHANNEL_TYPE_CHROMA][altIdx] = alfParamNewFilters.nonLinearFlag[CHANNEL_TYPE_CHROMA][altIdx];
        }
        newAPS->getAlfAPSParam().filterType[CHANNEL_TYPE_CHROMA] = alfParamNewFilters.filterType[CHANNEL_TYPE_CHROMA];
#else
        newAPS->getAlfAPSParam().nonLinearFlag[CHANNEL_TYPE_CHROMA] = alfParamNewFilters.nonLinearFlag[CHANNEL_TYPE_CHROMA];
#endif
        newAPS->setTemporalId( cs.slice->getTLayer() );
        for( int altIdx = 0; altIdx < MAX_NUM_ALF_ALTERNATIVES_CHROMA; ++altIdx )
        {
          for( int i = 0; i < MAX_NUM_ALF_CHROMA_COEFF; i++ )
          {
            newAPS->getAlfAPSParam().chromaCoeff[altIdx][i] = alfParamNewFilters.chromaCoeff[altIdx][i];
            newAPS->getAlfAPSParam().chromaClipp[altIdx][i] = alfParamNewFilters.chromaClipp[altIdx][i];
          }
        }
        m_apsMap->setChangedFlag( ( newApsIdChroma << NUM_APS_TYPE_LEN ) + ALF_APS );
        m_apsIdStart = newApsIdChroma;
      }
      apss[cs.slice->getTileGroupApsIdChroma()] = m_apsMap->getPS( ( cs.slice->getTileGroupApsIdChroma() << NUM_APS_TYPE_LEN ) + ALF_APS );
    }
  }
}

void EncAdaptiveLoopFilter::alfReconstructor(CodingStructure& cs, const PelUnitBuf& recExtBuf)
{
  if (!cs.slice->getTileGroupAlfEnabledFlag(COMPONENT_Y))
  {
    return;
  }
  reconstructCoeffAPSs(cs, true, cs.slice->getTileGroupAlfEnabledFlag(COMPONENT_Cb) || cs.slice->getTileGroupAlfEnabledFlag(COMPONENT_Cr), false);
  short* alfCtuFilterIndex = cs.slice->getPic()->getAlfCtbFilterIndex();
  PelUnitBuf& recBuf = cs.getRecoBufRef();
  const PreCalcValues& pcv = *cs.pcv;

#if ALF_IMPROVEMENT
  int fixedFilterSetIdx = cs.slice->getTileGroupAlfFixedFilterSetIdx();
#endif

  int ctuIdx = 0;
  bool clipTop = false, clipBottom = false, clipLeft = false, clipRight = false;
  int numHorVirBndry = 0, numVerVirBndry = 0;
  int horVirBndryPos[] = { 0, 0, 0 };
  int verVirBndryPos[] = { 0, 0, 0 };
  for (int yPos = 0; yPos < pcv.lumaHeight; yPos += pcv.maxCUHeight)
  {
    for (int xPos = 0; xPos < pcv.lumaWidth; xPos += pcv.maxCUWidth)
    {
      const int width = (xPos + pcv.maxCUWidth > pcv.lumaWidth) ? (pcv.lumaWidth - xPos) : pcv.maxCUWidth;
      const int height = (yPos + pcv.maxCUHeight > pcv.lumaHeight) ? (pcv.lumaHeight - yPos) : pcv.maxCUHeight;

      bool ctuEnableFlag = m_ctuEnableFlag[COMPONENT_Y][ctuIdx];
      for (int compIdx = 1; compIdx < MAX_NUM_COMPONENT; compIdx++)
      {
        ctuEnableFlag |= m_ctuEnableFlag[compIdx][ctuIdx] > 0;
      }
      int rasterSliceAlfPad = 0;
      if ( ctuEnableFlag && isCrossedByVirtualBoundaries( cs, xPos, yPos, width, height, clipTop, clipBottom, clipLeft, clipRight, numHorVirBndry, numVerVirBndry, horVirBndryPos, verVirBndryPos, rasterSliceAlfPad ) )
      {
        int yStart = yPos;
        for (int i = 0; i <= numHorVirBndry; i++)
        {
          const int yEnd = i == numHorVirBndry ? yPos + height : horVirBndryPos[i];
          const int h = yEnd - yStart;
          const bool clipT = (i == 0 && clipTop) || (i > 0) || (yStart == 0);
          const bool clipB = (i == numHorVirBndry && clipBottom) || (i < numHorVirBndry ) || (yEnd == pcv.lumaHeight);
          int xStart = xPos;
          for (int j = 0; j <= numVerVirBndry; j++)
          {
            const int xEnd = j == numVerVirBndry ? xPos + width : verVirBndryPos[j];
            const int w = xEnd - xStart;
            const bool clipL = (j == 0 && clipLeft) || (j > 0) || (xStart == 0);
            const bool clipR = (j == numVerVirBndry && clipRight) || (j < numVerVirBndry ) || (xEnd == pcv.lumaWidth);
            const int wBuf = w + (clipL ? 0 : MAX_ALF_PADDING_SIZE) + (clipR ? 0 : MAX_ALF_PADDING_SIZE);
            const int hBuf = h + (clipT ? 0 : MAX_ALF_PADDING_SIZE) + (clipB ? 0 : MAX_ALF_PADDING_SIZE);
            PelUnitBuf buf = m_tempBuf2.subBuf(UnitArea(cs.area.chromaFormat, Area(0, 0, wBuf, hBuf)));
            buf.copyFrom(recExtBuf.subBuf(UnitArea(cs.area.chromaFormat, Area(xStart - (clipL ? 0 : MAX_ALF_PADDING_SIZE), yStart - (clipT ? 0 : MAX_ALF_PADDING_SIZE), wBuf, hBuf))));
            // pad top-left unavailable samples for raster slice
            if ( xStart == xPos && yStart == yPos && ( rasterSliceAlfPad & 1 ) )
            {
              buf.padBorderPel( MAX_ALF_PADDING_SIZE, 1 );
            }

            // pad bottom-right unavailable samples for raster slice
            if ( xEnd == xPos + width && yEnd == yPos + height && ( rasterSliceAlfPad & 2 ) )
            {
              buf.padBorderPel( MAX_ALF_PADDING_SIZE, 2 );
            }
            buf.extendBorderPel(MAX_ALF_PADDING_SIZE);
            buf = buf.subBuf(UnitArea(cs.area.chromaFormat, Area(clipL ? 0 : MAX_ALF_PADDING_SIZE, clipT ? 0 : MAX_ALF_PADDING_SIZE, w, h)));


            if (m_ctuEnableFlag[COMPONENT_Y][ctuIdx])
            {
              const Area blkSrc(0, 0, w, h);
              const Area blkDst(xStart, yStart, w, h);
              short filterSetIndex = alfCtuFilterIndex[ctuIdx];
              short *coeff;
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
              Pel *clip;
#else
              short *clip;
#endif
#if ALF_IMPROVEMENT
              if (filterSetIndex < NUM_FIXED_FILTER_SETS)
              {
                copyFixedFilterResults(recBuf, blkDst, COMPONENT_Y, m_fixFilterResult, fixedFilterSetIdx, filterSetIndex);
              }
              else
              {
                const int alt_num = m_ctuAlternative[COMPONENT_Y][ctuIdx];
                coeff = m_coeffApsLuma[filterSetIndex - NUM_FIXED_FILTER_SETS][alt_num];
                clip = m_clippApsLuma[filterSetIndex - NUM_FIXED_FILTER_SETS][alt_num];
                AlfFilterType filterTypeCtb = m_filterTypeApsLuma[filterSetIndex - NUM_FIXED_FILTER_SETS];
                alfFiltering( m_classifier, recBuf, buf, blkDst, blkSrc, COMPONENT_Y, coeff, clip, m_clpRngs.comp[COMPONENT_Y], cs, filterTypeCtb, m_fixFilterResult, fixedFilterSetIdx );
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
              m_filter7x7Blk( m_classifier, recBuf, buf, blkDst, blkSrc, COMPONENT_Y, coeff, clip, m_clpRngs.comp[COMPONENT_Y], cs, m_alfVBLumaCTUHeight, m_alfVBLumaPos );
#endif
            }

            for (int compIdx = 1; compIdx < MAX_NUM_COMPONENT; compIdx++)
            {
              ComponentID compID = ComponentID(compIdx);
              const int chromaScaleX = getComponentScaleX(compID, recBuf.chromaFormat);
              const int chromaScaleY = getComponentScaleY(compID, recBuf.chromaFormat);
              if (m_ctuEnableFlag[compIdx][ctuIdx])
              {
                const Area blkSrc(0, 0, w >> chromaScaleX, h >> chromaScaleY);
                const Area blkDst(xStart >> chromaScaleX, yStart >> chromaScaleY, w >> chromaScaleX, h >> chromaScaleY);
                const int alt_num = m_ctuAlternative[compID][ctuIdx];
#if ALF_IMPROVEMENT
                alfFiltering( m_classifier, recBuf, buf, blkDst, blkSrc, compID, m_chromaCoeffFinal[alt_num], m_chromaClippFinal[alt_num], m_clpRngs.comp[compIdx], cs, m_filterTypeApsChroma, nullptr, -1 );
#else
                m_filter5x5Blk( m_classifier, recBuf, buf, blkDst, blkSrc, compID, m_chromaCoeffFinal[alt_num], m_chromaClippFinal[alt_num], m_clpRngs.comp[compIdx], cs, m_alfVBChmaCTUHeight, m_alfVBChmaPos );
#endif
              }
            }

            xStart = xEnd;
          }

          yStart = yEnd;
        }
      }
      else
      {

      const UnitArea area(cs.area.chromaFormat, Area(xPos, yPos, width, height));
      if (m_ctuEnableFlag[COMPONENT_Y][ctuIdx])
      {
        Area blk(xPos, yPos, width, height);
        short filterSetIndex = alfCtuFilterIndex[ctuIdx];
        short *coeff;
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
          Pel *clip;
#else        
        short *clip;
#endif
#if ALF_IMPROVEMENT
        if( filterSetIndex < NUM_FIXED_FILTER_SETS )
        {
          copyFixedFilterResults(recBuf, blk, COMPONENT_Y, m_fixFilterResult, fixedFilterSetIdx, filterSetIndex);
        }
        else
        {
          const int alt_num = m_ctuAlternative[COMPONENT_Y][ctuIdx];
          coeff = m_coeffApsLuma[filterSetIndex - NUM_FIXED_FILTER_SETS][alt_num];
          clip = m_clippApsLuma[filterSetIndex - NUM_FIXED_FILTER_SETS][alt_num];
          AlfFilterType filterTypeCtb = m_filterTypeApsLuma[filterSetIndex - NUM_FIXED_FILTER_SETS];
          alfFiltering(m_classifier, recBuf, recExtBuf, blk, blk, COMPONENT_Y, coeff, clip, m_clpRngs.comp[COMPONENT_Y], cs, filterTypeCtb, m_fixFilterResult, fixedFilterSetIdx);
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
        m_filter7x7Blk( m_classifier, recBuf, recExtBuf, blk, blk, COMPONENT_Y, coeff, clip, m_clpRngs.comp[COMPONENT_Y], cs, m_alfVBLumaCTUHeight, m_alfVBLumaPos );
#endif
      }

      for (int compIdx = 1; compIdx < MAX_NUM_COMPONENT; compIdx++)
      {
        ComponentID compID = ComponentID(compIdx);
        const int chromaScaleX = getComponentScaleX(compID, recBuf.chromaFormat);
        const int chromaScaleY = getComponentScaleY(compID, recBuf.chromaFormat);
        if (m_ctuEnableFlag[compIdx][ctuIdx])
        {
          Area blk(xPos >> chromaScaleX, yPos >> chromaScaleY, width >> chromaScaleX, height >> chromaScaleY);
          const int alt_num = m_ctuAlternative[compID][ctuIdx];
#if ALF_IMPROVEMENT
          alfFiltering( m_classifier, recBuf, recExtBuf, blk, blk, compID, m_chromaCoeffFinal[alt_num], m_chromaClippFinal[alt_num], m_clpRngs.comp[compIdx], cs, m_filterTypeApsChroma, nullptr, -1 );
#else
          m_filter5x5Blk( m_classifier, recBuf, recExtBuf, blk, blk, compID, m_chromaCoeffFinal[alt_num], m_chromaClippFinal[alt_num], m_clpRngs.comp[compIdx], cs, m_alfVBChmaCTUHeight, m_alfVBChmaPos );
#endif
        }
      }
      }
      ctuIdx++;
    }
  }
}

void EncAdaptiveLoopFilter::copyCtuAlternativeChroma( uint8_t* ctuAltsDst[MAX_NUM_COMPONENT], uint8_t* ctuAltsSrc[MAX_NUM_COMPONENT] )
{
  std::copy_n( ctuAltsSrc[COMPONENT_Cb], m_numCTUsInPic, ctuAltsDst[COMPONENT_Cb] );
  std::copy_n( ctuAltsSrc[COMPONENT_Cr], m_numCTUsInPic, ctuAltsDst[COMPONENT_Cr] );
}

void EncAdaptiveLoopFilter::setCtuAlternativeChroma( uint8_t* ctuAlts[MAX_NUM_COMPONENT], uint8_t val )
{
  std::fill_n( ctuAlts[COMPONENT_Cb], m_numCTUsInPic, val );
  std::fill_n( ctuAlts[COMPONENT_Cr], m_numCTUsInPic, val );
}

#if ALF_IMPROVEMENT
void EncAdaptiveLoopFilter::initCtuAlternativeLuma(uint8_t* ctuAlts[MAX_NUM_COMPONENT])
{
  uint8_t altIdx = 0;
  for (int ctuIdx = 0; ctuIdx < m_numCTUsInPic; ++ctuIdx)
  {
    ctuAlts[COMPONENT_Y][ctuIdx] = altIdx;
    if ((ctuIdx + 1) * m_alfParamTemp.numAlternativesLuma >= (altIdx + 1)*m_numCTUsInPic)
      ++altIdx;
  }
}
#endif

void EncAdaptiveLoopFilter::initCtuAlternativeChroma( uint8_t* ctuAlts[MAX_NUM_COMPONENT] )
{
  uint8_t altIdx = 0;
  for( int ctuIdx = 0; ctuIdx < m_numCTUsInPic; ++ctuIdx )
  {
    ctuAlts[COMPONENT_Cb][ctuIdx] = altIdx;
    if( (ctuIdx+1) * m_alfParamTemp.numAlternativesChroma >= (altIdx+1)*m_numCTUsInPic )
      ++altIdx;
  }
  altIdx = 0;
  for( int ctuIdx = 0; ctuIdx < m_numCTUsInPic; ++ctuIdx )
  {
    ctuAlts[COMPONENT_Cr][ctuIdx] = altIdx;
    if( (ctuIdx+1) * m_alfParamTemp.numAlternativesChroma >= (altIdx+1)*m_numCTUsInPic )
      ++altIdx;
  }
}

int EncAdaptiveLoopFilter::getMaxNumAlternativesChroma( )
{
  return std::min<int>( m_numCTUsInPic * 2, m_encCfg->getMaxNumAlfAlternativesChroma() );
}

int EncAdaptiveLoopFilter::getCoeffRateCcAlf(short chromaCoeff[MAX_NUM_CC_ALF_FILTERS][MAX_NUM_CC_ALF_CHROMA_COEFF], bool filterEnabled[MAX_NUM_CC_ALF_FILTERS], uint8_t filterCount, ComponentID compID)
{
  int bits = 0;

  if ( filterCount > 0 )
  {
    bits += lengthUvlc(filterCount - 1);
    int signaledFilterCount = 0;
    for ( int filterIdx=0; filterIdx<MAX_NUM_CC_ALF_FILTERS; filterIdx++ )
    {
      if (filterEnabled[filterIdx])
      {
        AlfFilterShape alfShape(size_CC_ALF);
        // Filter coefficients
        for (int i = 0; i < alfShape.numCoeff - 1; i++)
        {
          bits += CCALF_BITS_PER_COEFF_LEVEL + (chromaCoeff[filterIdx][i] == 0 ? 0 : 1);
        }

        signaledFilterCount++;
      }
    }
    CHECK(signaledFilterCount != filterCount, "Number of filter signaled not same as indicated");
  }

  return bits;
}

void EncAdaptiveLoopFilter::deriveCcAlfFilterCoeff( ComponentID compID, const PelUnitBuf& recYuv, const PelUnitBuf& recYuvExt, short filterCoeff[MAX_NUM_CC_ALF_FILTERS][MAX_NUM_CC_ALF_CHROMA_COEFF], const uint8_t filterIdx )
{
  int forward_tab[CCALF_CANDS_COEFF_NR * 2 - 1] = {0};
  for (int i = 0; i < CCALF_CANDS_COEFF_NR; i++)
  {
    forward_tab[CCALF_CANDS_COEFF_NR - 1 + i] = CCALF_SMALL_TAB[i];
    forward_tab[CCALF_CANDS_COEFF_NR - 1 - i] = (-1) * CCALF_SMALL_TAB[i];
  }
  using TE = double[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF];
  using Ty = double[MAX_NUM_ALF_LUMA_COEFF];

  double filterCoeffDbl[MAX_NUM_CC_ALF_CHROMA_COEFF];
  int16_t filterCoeffInt[MAX_NUM_CC_ALF_CHROMA_COEFF];

  std::fill_n(filterCoeffInt, MAX_NUM_CC_ALF_CHROMA_COEFF, 0);

  TE        kE;
  Ty        ky;
  const int size = m_filterShapesCcAlf[compID - 1][0].numCoeff - 1;

  for (int k = 0; k < size; k++)
  {
    ky[k] = m_alfCovarianceFrameCcAlf[compID - 1][0][filterIdx].y[0][k];
    for (int l = 0; l < size; l++)
    {
      kE[k][l] = m_alfCovarianceFrameCcAlf[compID - 1][0][filterIdx].E[0][0][k][l];
    }
  }

  m_alfCovarianceFrameCcAlf[compID - 1][0][filterIdx].gnsSolveByChol(kE, ky, filterCoeffDbl, size);
  roundFiltCoeffCCALF(filterCoeffInt, filterCoeffDbl, size, (1 << m_scaleBits));

  for (int k = 0; k < size; k++)
  {
    CHECK( filterCoeffInt[k] < -(1 << CCALF_DYNAMIC_RANGE), "this is not possible: filterCoeffInt[k] <  -(1 << CCALF_DYNAMIC_RANGE)");
    CHECK( filterCoeffInt[k] > (1 << CCALF_DYNAMIC_RANGE), "this is not possible: filterCoeffInt[k] >  (1 << CCALF_DYNAMIC_RANGE)");
  }

  // Refine quanitzation
  int modified       = 1;
  double errRef      = m_alfCovarianceFrameCcAlf[compID - 1][0][filterIdx].calcErrorForCcAlfCoeffs(filterCoeffInt, size, (m_scaleBits+1));
  while (modified)
  {
    modified = 0;
    for (int delta : { 1, -1 })
    {
      double errMin = MAX_DOUBLE;
      int    idxMin = -1;
      int minIndex = -1;

      for (int k = 0; k < size; k++)
      {
        int org_idx = -1;
        for (int i = 0; i < CCALF_CANDS_COEFF_NR * 2 - 1; i++)
        {
          if (forward_tab[i] == filterCoeffInt[k])
          {
            org_idx = i;
            break;
          }
        }
        CHECK( org_idx < 0, "this is wrong, does not find coeff from forward_tab");
        if ( (org_idx - delta < 0) || (org_idx - delta >= CCALF_CANDS_COEFF_NR * 2 - 1) )
          continue;

        filterCoeffInt[k] = forward_tab[org_idx - delta];
        double error = m_alfCovarianceFrameCcAlf[compID - 1][0][filterIdx].calcErrorForCcAlfCoeffs(filterCoeffInt, size, (m_scaleBits+1));
        if( error < errMin )
        {
          errMin = error;
          idxMin = k;
          minIndex = org_idx;
        }
        filterCoeffInt[k] = forward_tab[org_idx];
      }
      if (errMin < errRef)
      {
        minIndex -= delta;
        CHECK( minIndex < 0, "this is wrong, index - delta < 0");
        CHECK( minIndex >= CCALF_CANDS_COEFF_NR * 2 - 1, "this is wrong, index - delta >= CCALF_CANDS_COEFF_NR * 2 - 1");
        filterCoeffInt[idxMin] = forward_tab[minIndex];
        modified++;
        errRef = errMin;
      }
    }
  }

  for (int k = 0; k < (size + 1); k++)
  {
    CHECK((filterCoeffInt[k] < -(1 << CCALF_DYNAMIC_RANGE)) || (filterCoeffInt[k] > (1 << CCALF_DYNAMIC_RANGE)), "Exceeded valid range for CC ALF coefficient");
    filterCoeff[filterIdx][k] = filterCoeffInt[k];
  }
}

void EncAdaptiveLoopFilter::determineControlIdcValues(CodingStructure &cs, const ComponentID compID, const PelBuf *buf,
                                                      const int ctuWidthC, const int ctuHeightC, const int picWidthC,
                                                      const int picHeightC, double **unfilteredDistortion,
                                                      uint64_t *trainingDistortion[MAX_NUM_CC_ALF_FILTERS],
                                                      uint64_t *lumaSwingGreaterThanThresholdCount,
                                                      uint64_t *chromaSampleCountNearMidPoint,
                                                      bool reuseTemporalFilterCoeff, uint8_t *trainingCovControl,
                                                      uint8_t *filterControl, uint64_t &curTotalDistortion,
                                                      double &curTotalRate, bool filterEnabled[MAX_NUM_CC_ALF_FILTERS],
                                                      uint8_t  mapFilterIdxToFilterIdc[MAX_NUM_CC_ALF_FILTERS + 1],
                                                      uint8_t &ccAlfFilterCount)
{
  bool curFilterEnabled[MAX_NUM_CC_ALF_FILTERS];
  std::fill_n(curFilterEnabled, MAX_NUM_CC_ALF_FILTERS, false);

#if MAX_NUM_CC_ALF_FILTERS>1
  FilterIdxCount filterIdxCount[MAX_NUM_CC_ALF_FILTERS];
  for (int i = 0; i < MAX_NUM_CC_ALF_FILTERS; i++)
  {
    filterIdxCount[i].count     = 0;
    filterIdxCount[i].filterIdx = i;
  }

  double prevRate = curTotalRate;
#endif

  TempCtx ctxInitial(m_CtxCache);
  TempCtx ctxBest(m_CtxCache);
  TempCtx ctxStart(m_CtxCache);
  ctxInitial = SubCtx(Ctx::CcAlfFilterControlFlag, m_CABACEstimator->getCtx());
  ctxBest    = SubCtx(Ctx::CcAlfFilterControlFlag, m_CABACEstimator->getCtx());

  int ctuIdx = 0;
  for (int yCtu = 0; yCtu < buf->height; yCtu += ctuHeightC)
  {
    for (int xCtu = 0; xCtu < buf->width; xCtu += ctuWidthC)
    {
      uint64_t ssd;
      double   rate;
      double   cost;

      uint64_t bestSSD       = MAX_UINT64;
      double   bestRate      = MAX_DOUBLE;
      double   bestCost      = MAX_DOUBLE;
      uint8_t  bestFilterIdc = 0;
      uint8_t  bestFilterIdx = 0;
      const uint32_t thresholdS = std::min<int>(buf->height - yCtu, ctuHeightC) << getComponentScaleY(COMPONENT_Cb, m_chromaFormat);
      const uint32_t numberOfChromaSamples = std::min<int>(buf->height - yCtu, ctuHeightC) * std::min<int>(buf->width - xCtu, ctuWidthC);
      const uint32_t thresholdC = (numberOfChromaSamples >> 2);

      m_CABACEstimator->getCtx() = ctxBest;
      ctxStart                   = SubCtx(Ctx::CcAlfFilterControlFlag, m_CABACEstimator->getCtx());

      for (int filterIdx = 0; filterIdx <= MAX_NUM_CC_ALF_FILTERS; filterIdx++)
      {
        uint8_t filterIdc = mapFilterIdxToFilterIdc[filterIdx];
        if (filterIdx < MAX_NUM_CC_ALF_FILTERS && !filterEnabled[filterIdx])
        {
          continue;
        }

        if (filterIdx == MAX_NUM_CC_ALF_FILTERS)
        {
          ssd = (uint64_t)unfilteredDistortion[compID][ctuIdx];   // restore saved distortion computation
        }
        else
        {
          ssd = trainingDistortion[filterIdx][ctuIdx];
        }
        m_CABACEstimator->getCtx() = ctxStart;
        m_CABACEstimator->resetBits();
        const Position lumaPos = Position({ xCtu << getComponentScaleX(compID, cs.pcv->chrFormat),
          yCtu << getComponentScaleY(compID, cs.pcv->chrFormat) });
        m_CABACEstimator->codeCcAlfFilterControlIdc(filterIdc, cs, compID, ctuIdx, filterControl, lumaPos,
                                                    ccAlfFilterCount);
        rate = FRAC_BITS_SCALE * m_CABACEstimator->getEstFracBits();
        cost = rate * m_lambda[compID] + ssd;

        bool limitationExceeded = false;
        if (m_limitCcAlf && filterIdx < MAX_NUM_CC_ALF_FILTERS)
        {
          limitationExceeded = limitationExceeded || (lumaSwingGreaterThanThresholdCount[ctuIdx] >= thresholdS);
          limitationExceeded = limitationExceeded || (chromaSampleCountNearMidPoint[ctuIdx] >= thresholdC);
        }
        if (cost < bestCost && !limitationExceeded)
        {
          bestCost      = cost;
          bestRate      = rate;
          bestSSD       = ssd;
          bestFilterIdc = filterIdc;
          bestFilterIdx = filterIdx;

          ctxBest = SubCtx(Ctx::CcAlfFilterControlFlag, m_CABACEstimator->getCtx());

          trainingCovControl[ctuIdx] = (filterIdx == MAX_NUM_CC_ALF_FILTERS) ? 0 : (filterIdx + 1);
          filterControl[ctuIdx]      = (filterIdx == MAX_NUM_CC_ALF_FILTERS) ? 0 : (filterIdx + 1);
        }
      }
      if (bestFilterIdc != 0)
      {
        curFilterEnabled[bestFilterIdx] = true;
#if MAX_NUM_CC_ALF_FILTERS>1
        filterIdxCount[bestFilterIdx].count++;
#endif
      }
      curTotalRate += bestRate;
      curTotalDistortion += bestSSD;
      ctuIdx++;
    }
  }

#if MAX_NUM_CC_ALF_FILTERS>1
  if (!reuseTemporalFilterCoeff)
  {
    std::copy_n(curFilterEnabled, MAX_NUM_CC_ALF_FILTERS, filterEnabled);

    std::stable_sort(filterIdxCount, filterIdxCount + MAX_NUM_CC_ALF_FILTERS, compareCounts);

    int filterIdc = 1;
    ccAlfFilterCount = 0;
    for ( FilterIdxCount &s : filterIdxCount )
    {
      const int filterIdx = s.filterIdx;
      if (filterEnabled[filterIdx])
      {
        mapFilterIdxToFilterIdc[filterIdx] = filterIdc;
        filterIdc++;
        ccAlfFilterCount++;
      }
    }

    curTotalRate = prevRate;
    m_CABACEstimator->getCtx() = ctxInitial;
    m_CABACEstimator->resetBits();
    int ctuIdx = 0;
    for (int y = 0; y < buf->height; y += ctuHeightC)
    {
      for (int x = 0; x < buf->width; x += ctuWidthC)
      {
        const int filterIdxPlus1 = filterControl[ctuIdx];

        const Position lumaPos = Position(
                                          { x << getComponentScaleX(compID, cs.pcv->chrFormat), y << getComponentScaleY(compID, cs.pcv->chrFormat) });

        m_CABACEstimator->codeCcAlfFilterControlIdc(filterIdxPlus1 == 0 ? 0
                                                    : mapFilterIdxToFilterIdc[filterIdxPlus1 - 1],
                                                    cs, compID, ctuIdx, filterControl, lumaPos, ccAlfFilterCount);

        ctuIdx++;
      }
    }
    curTotalRate += FRAC_BITS_SCALE*m_CABACEstimator->getEstFracBits();
  }
#endif

  // restore for next iteration
  m_CABACEstimator->getCtx() = ctxInitial;
}

std::vector<int> EncAdaptiveLoopFilter::getAvailableCcAlfApsIds(CodingStructure& cs, ComponentID compID)
{
  APS** apss = cs.slice->getAlfAPSs();
  for (int i = 0; i < ALF_CTB_MAX_NUM_APS; i++)
  {
    apss[i] = m_apsMap->getPS((i << NUM_APS_TYPE_LEN) + ALF_APS);
  }

  std::vector<int> result;
  int apsIdChecked = 0, curApsId = m_apsIdStart;
  if (curApsId < ALF_CTB_MAX_NUM_APS)
  {
    while (apsIdChecked < ALF_CTB_MAX_NUM_APS && !cs.slice->isIntra() && result.size() < ALF_CTB_MAX_NUM_APS && !cs.slice->getPendingRasInit() && !cs.slice->isIDRorBLA())
    {
      APS* curAPS = cs.slice->getAlfAPSs()[curApsId];
      if (curAPS && curAPS->getLayerId() == cs.slice->getPic()->layerId
          && curAPS->getTemporalId() <= cs.slice->getTLayer() && curAPS->getCcAlfAPSParam().newCcAlfFilter[compID - 1])
      {
        result.push_back(curApsId);
      }
      apsIdChecked++;
      curApsId = (curApsId + 1) % ALF_CTB_MAX_NUM_APS;
    }
  }
  return result;
}

void EncAdaptiveLoopFilter::getFrameStatsCcalf(ComponentID compIdx, int filterIdc)
{
        int ctuRsAddr = 0;
  const int filterIdx = filterIdc - 1;

  // init Frame stats buffers
  for (int shape = 0; shape != m_filterShapesCcAlf[compIdx - 1].size(); shape++)
  {
    m_alfCovarianceFrameCcAlf[compIdx - 1][shape][filterIdx].reset();
  }

  for (int yPos = 0; yPos < m_picHeight; yPos += m_maxCUHeight)
  {
    for (int xPos = 0; xPos < m_picWidth; xPos += m_maxCUWidth)
    {
      if (m_trainingCovControl[ctuRsAddr] == filterIdc)
      {
        for (int shape = 0; shape != m_filterShapesCcAlf[compIdx - 1].size(); shape++)
        {
          m_alfCovarianceFrameCcAlf[compIdx - 1][shape][filterIdx] +=
            m_alfCovarianceCcAlf[compIdx - 1][shape][0][ctuRsAddr];
        }
      }
      ctuRsAddr++;
    }
  }
}

void EncAdaptiveLoopFilter::deriveCcAlfFilter( CodingStructure& cs, ComponentID compID, const PelUnitBuf& orgYuv, const PelUnitBuf& tempDecYuvBuf, const PelUnitBuf& dstYuv )
{
  if (!cs.slice->getTileGroupAlfEnabledFlag(COMPONENT_Y))
  {
    m_ccAlfFilterParam.ccAlfFilterEnabled[compID - 1] = false;
    return;
  }

  m_limitCcAlf = m_encCfg->getBaseQP() >= m_encCfg->getCCALFQpThreshold();
  if (m_limitCcAlf && cs.slice->getSliceQp() <= m_encCfg->getBaseQP() + 1)
  {
    m_ccAlfFilterParam.ccAlfFilterEnabled[compID - 1] = false;
    return;
  }

  uint8_t bestMapFilterIdxToFilterIdc[MAX_NUM_CC_ALF_FILTERS+1];
  const int scaleX               = getComponentScaleX(compID, cs.pcv->chrFormat);
  const int scaleY               = getComponentScaleY(compID, cs.pcv->chrFormat);
  const int ctuWidthC            = cs.pcv->maxCUWidth >> scaleX;
  const int ctuHeightC           = cs.pcv->maxCUHeight >> scaleY;
  const int picWidthC            = cs.pcv->lumaWidth >> scaleX;
  const int picHeightC           = cs.pcv->lumaHeight >> scaleY;
  const int maxTrainingIterCount = 15;

  if (m_limitCcAlf)
  {
    countLumaSwingGreaterThanThreshold(dstYuv.get(COMPONENT_Y).bufAt(0, 0), dstYuv.get(COMPONENT_Y).stride, dstYuv.get(COMPONENT_Y).height, dstYuv.get(COMPONENT_Y).width, cs.pcv->maxCUWidthLog2, cs.pcv->maxCUHeightLog2, m_lumaSwingGreaterThanThresholdCount, m_numCTUsInWidth);
  }
  if (m_limitCcAlf)
  {
    countChromaSampleValueNearMidPoint(dstYuv.get(compID).bufAt(0, 0), dstYuv.get(compID).stride, dstYuv.get(compID).height, dstYuv.get(compID).width, cs.pcv->maxCUWidthLog2 - scaleX, cs.pcv->maxCUHeightLog2 - scaleY, m_chromaSampleCountNearMidPoint, m_numCTUsInWidth);
  }

  for ( int filterIdx = 0; filterIdx <= MAX_NUM_CC_ALF_FILTERS; filterIdx++ )
  {
    if ( filterIdx < MAX_NUM_CC_ALF_FILTERS)
    {
      memset( m_bestFilterCoeffSet[filterIdx], 0, sizeof(m_bestFilterCoeffSet[filterIdx]) );
      bestMapFilterIdxToFilterIdc[filterIdx] = filterIdx + 1;
    }
    else
    {
      bestMapFilterIdxToFilterIdc[filterIdx] = 0;
    }
  }
  memset(m_bestFilterControl, 0, sizeof(uint8_t) * m_numCTUsInPic);
  int ccalfReuseApsId      = -1;
  m_reuseApsId[compID - 1] = -1;

  const TempCtx ctxStartCcAlfFilterControlFlag  ( m_CtxCache, SubCtx( Ctx::CcAlfFilterControlFlag, m_CABACEstimator->getCtx() ) );

  // compute cost of not filtering
  uint64_t unfilteredDistortion = 0;
  for (int ctbIdx = 0; ctbIdx < m_numCTUsInPic; ctbIdx++)
  {
    unfilteredDistortion += (uint64_t)m_alfCovarianceCcAlf[compID - 1][0][0][ctbIdx].pixAcc;
  }

  double bestUnfilteredTotalCost = 1 * m_lambda[compID] + unfilteredDistortion;   // 1 bit is for gating flag

  bool             ccAlfFilterIdxEnabled[MAX_NUM_CC_ALF_FILTERS];
  short            ccAlfFilterCoeff[MAX_NUM_CC_ALF_FILTERS][MAX_NUM_CC_ALF_CHROMA_COEFF];
  uint8_t          ccAlfFilterCount             = MAX_NUM_CC_ALF_FILTERS;
  double bestFilteredTotalCost        = MAX_DOUBLE;
  bool   bestreuseTemporalFilterCoeff = false;
  std::vector<int> apsIds             = getAvailableCcAlfApsIds(cs, compID);

  for (int testFilterIdx = 0; testFilterIdx < ( apsIds.size() + 1 ); testFilterIdx++ )
  {
    bool referencingExistingAps   = (testFilterIdx < apsIds.size()) ? true : false;
    int maxNumberOfFiltersBeingTested = MAX_NUM_CC_ALF_FILTERS - (testFilterIdx - static_cast<int>(apsIds.size()));

    if (maxNumberOfFiltersBeingTested < 0)
    {
      maxNumberOfFiltersBeingTested = 1;
    }

    {
      // Instead of rewriting the control buffer for every training iteration just keep a mapping from filterIdx to filterIdc
      uint8_t mapFilterIdxToFilterIdc[MAX_NUM_CC_ALF_FILTERS + 1];
      for (int filterIdx = 0; filterIdx <= MAX_NUM_CC_ALF_FILTERS; filterIdx++)
      {
        if (filterIdx == MAX_NUM_CC_ALF_FILTERS)
        {
          mapFilterIdxToFilterIdc[filterIdx] = 0;
        }
        else
        {
          mapFilterIdxToFilterIdc[filterIdx] = filterIdx + 1;
        }
      }

      // initialize filters
      for ( int filterIdx = 0; filterIdx < MAX_NUM_CC_ALF_FILTERS; filterIdx++ )
      {
        ccAlfFilterIdxEnabled[filterIdx] = false;
        memset(ccAlfFilterCoeff[filterIdx], 0, sizeof(ccAlfFilterCoeff[filterIdx]));
      }
      if ( referencingExistingAps )
      {
        maxNumberOfFiltersBeingTested = m_apsMap->getPS((apsIds[testFilterIdx] << NUM_APS_TYPE_LEN) + ALF_APS)->getCcAlfAPSParam().ccAlfFilterCount[compID - 1];
        ccAlfFilterCount = maxNumberOfFiltersBeingTested;
        for (int filterIdx = 0; filterIdx < maxNumberOfFiltersBeingTested; filterIdx++)
        {
          ccAlfFilterIdxEnabled[filterIdx] = true;
          memcpy(ccAlfFilterCoeff[filterIdx], m_ccAlfFilterParam.ccAlfCoeff[compID - 1][filterIdx],
                 sizeof(ccAlfFilterCoeff[filterIdx]));
        }
        memcpy( ccAlfFilterCoeff, m_apsMap->getPS((apsIds[testFilterIdx] << NUM_APS_TYPE_LEN) + ALF_APS)->getCcAlfAPSParam().ccAlfCoeff[compID - 1], sizeof(ccAlfFilterCoeff) );
      }
      else
      {
        for (int i = 0; i < maxNumberOfFiltersBeingTested; i++)
        {
          ccAlfFilterIdxEnabled[i] = true;
        }
        ccAlfFilterCount = maxNumberOfFiltersBeingTested;
      }

      // initialize
      int controlIdx = 0;
      const int columnSize = ( m_buf->width / maxNumberOfFiltersBeingTested);
      for (int y = 0; y < m_buf->height; y += ctuHeightC)
      {
        for (int x = 0; x < m_buf->width; x += ctuWidthC)
        {
          m_trainingCovControl[controlIdx] = ( x / columnSize ) + 1;
          controlIdx++;
        }
      }

      // compute cost of filtering
      int    trainingIterCount = 0;
      bool   keepTraining      = true;
      bool   improvement       = false;
      double prevTotalCost     = MAX_DOUBLE;
      while (keepTraining)
      {
        improvement = false;
        for (int filterIdx = 0; filterIdx < maxNumberOfFiltersBeingTested; filterIdx++)
        {
          if (ccAlfFilterIdxEnabled[filterIdx])
          {
            if (!referencingExistingAps)
            {
              getFrameStatsCcalf(compID, (filterIdx + 1));
              deriveCcAlfFilterCoeff(compID, dstYuv, tempDecYuvBuf, ccAlfFilterCoeff, filterIdx);
            }
            const int numCoeff  = m_filterShapesCcAlf[compID - 1][0].numCoeff - 1;
            int log2BlockWidth  = cs.pcv->maxCUWidthLog2 - scaleX;
            int log2BlockHeight = cs.pcv->maxCUHeightLog2 - scaleY;
            for (int y = 0; y < m_buf->height; y += (1 << log2BlockHeight))
            {
              for (int x = 0; x < m_buf->width; x += (1 << log2BlockWidth))
              {
                int ctuIdx = (y >> log2BlockHeight) * m_numCTUsInWidth + (x >> log2BlockWidth);
                m_trainingDistortion[filterIdx][ctuIdx] =
                  int(m_ctbDistortionUnfilter[compID][ctuIdx]
                      + m_alfCovarianceCcAlf[compID - 1][0][0][ctuIdx].calcErrorForCcAlfCoeffs(
                        ccAlfFilterCoeff[filterIdx], numCoeff, m_scaleBits + 1));
              }
            }
          }
        }

        m_CABACEstimator->getCtx() = ctxStartCcAlfFilterControlFlag;

        uint64_t curTotalDistortion = 0;
        double curTotalRate = 0;
        determineControlIdcValues(cs, compID, m_buf, ctuWidthC, ctuHeightC, picWidthC, picHeightC,
                                  m_ctbDistortionUnfilter, m_trainingDistortion,
                                  m_lumaSwingGreaterThanThresholdCount,
                                  m_chromaSampleCountNearMidPoint,
                                  (referencingExistingAps == true),
                                  m_trainingCovControl, m_filterControl, curTotalDistortion, curTotalRate,
                                  ccAlfFilterIdxEnabled, mapFilterIdxToFilterIdc, ccAlfFilterCount);

        // compute coefficient coding bit cost
        if (ccAlfFilterCount > 0)
        {
          if (referencingExistingAps)
          {
            curTotalRate += 1 + 3; // +1 for enable flag, +3 APS ID in slice header
          }
          else
          {
            curTotalRate += getCoeffRateCcAlf(ccAlfFilterCoeff, ccAlfFilterIdxEnabled, ccAlfFilterCount, compID) + 1
            + 9;   // +1 for the enable flag, +9 3-bit for APS ID in slice header, 5-bit for APS ID in APS, a 1-bit
            // new filter flags (ignore shared cost such as other new-filter flags/NALU header/RBSP
            // terminating bit/byte alignment bits)
          }

          double curTotalCost = curTotalRate * m_lambda[compID] + curTotalDistortion;

          if (curTotalCost < prevTotalCost)
          {
            prevTotalCost = curTotalCost;
            improvement = true;
          }

          if (curTotalCost < bestFilteredTotalCost)
          {
            bestFilteredTotalCost = curTotalCost;
            memcpy(m_bestFilterIdxEnabled, ccAlfFilterIdxEnabled, sizeof(ccAlfFilterIdxEnabled));
            memcpy(m_bestFilterCoeffSet, ccAlfFilterCoeff, sizeof(ccAlfFilterCoeff));
            memcpy(m_bestFilterControl, m_filterControl, sizeof(uint8_t) * m_numCTUsInPic);
            m_bestFilterCount = ccAlfFilterCount;
            ccalfReuseApsId = referencingExistingAps ? apsIds[testFilterIdx] : -1;
            memcpy(bestMapFilterIdxToFilterIdc, mapFilterIdxToFilterIdc, sizeof(mapFilterIdxToFilterIdc));
          }
        }

        trainingIterCount++;
        if (!improvement || trainingIterCount > maxTrainingIterCount || referencingExistingAps)
        {
          keepTraining = false;
        }
      }
    }
  }

  if (bestUnfilteredTotalCost < bestFilteredTotalCost)
  {
    memset(m_bestFilterControl, 0, sizeof(uint8_t) * m_numCTUsInPic);
  }

  // save best coeff and control
  bool atleastOneBlockUndergoesFitlering = false;
  for (int controlIdx = 0; m_bestFilterCount > 0 && controlIdx < m_numCTUsInPic; controlIdx++)
  {
    if (m_bestFilterControl[controlIdx])
    {
      atleastOneBlockUndergoesFitlering = true;
      break;
    }
  }
  m_ccAlfFilterParam.numberValidComponents          = getNumberValidComponents(m_chromaFormat);
  m_ccAlfFilterParam.ccAlfFilterEnabled[compID - 1] = atleastOneBlockUndergoesFitlering;
  if (atleastOneBlockUndergoesFitlering)
  {
    // update the filter control indicators
    if (bestreuseTemporalFilterCoeff!=1)
    {
      short storedBestFilterCoeffSet[MAX_NUM_CC_ALF_FILTERS][MAX_NUM_CC_ALF_CHROMA_COEFF];
      for (int filterIdx=0; filterIdx<MAX_NUM_CC_ALF_FILTERS; filterIdx++)
      {
        memcpy(storedBestFilterCoeffSet[filterIdx], m_bestFilterCoeffSet[filterIdx], sizeof(m_bestFilterCoeffSet[filterIdx]));
      }
      memcpy(m_filterControl, m_bestFilterControl, sizeof(uint8_t) * m_numCTUsInPic);

      int filterCount = 0;
      for ( int filterIdx = 0; filterIdx < MAX_NUM_CC_ALF_FILTERS; filterIdx++ )
      {
        uint8_t curFilterIdc = bestMapFilterIdxToFilterIdc[filterIdx];
        if (m_bestFilterIdxEnabled[filterIdx])
        {
          for (int controlIdx = 0; controlIdx < m_numCTUsInPic; controlIdx++)
          {
            if (m_filterControl[controlIdx] == (filterIdx+1) )
            {
              m_bestFilterControl[controlIdx] = curFilterIdc;
            }
          }
          memcpy( m_bestFilterCoeffSet[curFilterIdc-1], storedBestFilterCoeffSet[filterIdx], sizeof(storedBestFilterCoeffSet[filterIdx]) );
          filterCount++;
        }
        m_bestFilterIdxEnabled[filterIdx] = ( filterIdx < m_bestFilterCount ) ? true : false;
      }
      CHECK( filterCount != m_bestFilterCount, "Number of filters enabled did not match the filter count");
    }

    m_ccAlfFilterParam.ccAlfFilterCount[compID - 1] = m_bestFilterCount;
    // cleanup before copying
    memset(m_ccAlfFilterControl[compID - 1], 0, sizeof(uint8_t) * m_numCTUsInPic);
    for ( int filterIdx = 0; filterIdx < MAX_NUM_CC_ALF_FILTERS; filterIdx++ )
    {
      memset(m_ccAlfFilterParam.ccAlfCoeff[compID - 1][filterIdx], 0,
             sizeof(m_ccAlfFilterParam.ccAlfCoeff[compID - 1][filterIdx]));
    }
    memset(m_ccAlfFilterParam.ccAlfFilterIdxEnabled[compID - 1], false,
           sizeof(m_ccAlfFilterParam.ccAlfFilterIdxEnabled[compID - 1]));
    for ( int filterIdx = 0; filterIdx < m_bestFilterCount; filterIdx++ )
    {
      m_ccAlfFilterParam.ccAlfFilterIdxEnabled[compID - 1][filterIdx] = m_bestFilterIdxEnabled[filterIdx];
      memcpy(m_ccAlfFilterParam.ccAlfCoeff[compID - 1][filterIdx], m_bestFilterCoeffSet[filterIdx],
             sizeof(m_bestFilterCoeffSet[filterIdx]));
    }
    memcpy(m_ccAlfFilterControl[compID - 1], m_bestFilterControl, sizeof(uint8_t) * m_numCTUsInPic);
    if ( ccalfReuseApsId >= 0 )
    {
      m_reuseApsId[compID - 1] = ccalfReuseApsId;
      if (compID == COMPONENT_Cb)
      {
        cs.slice->setTileGroupCcAlfCbApsId(ccalfReuseApsId);
      }
      else
      {
        cs.slice->setTileGroupCcAlfCrApsId(ccalfReuseApsId);
      }
    }
  }
}

template<bool alfWSSD>
void EncAdaptiveLoopFilter::deriveStatsForCcAlfFiltering( const PelUnitBuf &orgYuv, const PelUnitBuf &recYuv,
                                                          const int compIdx, const int maskStride,
                                                          const uint8_t filterIdc, CodingStructure &cs )
{
  const int filterIdx = filterIdc - 1;

  // init CTU stats buffers
  for( int shape = 0; shape != m_filterShapesCcAlf[compIdx - 1].size(); shape++ )
  {
    for( int ctuIdx = 0; ctuIdx < m_numCTUsInPic; ctuIdx++ )
    {
      m_alfCovarianceCcAlf[compIdx - 1][shape][filterIdx][ctuIdx].reset();
    }
  }

  // init Frame stats buffers
  for( int shape = 0; shape != m_filterShapesCcAlf[compIdx - 1].size(); shape++ )
  {
    m_alfCovarianceFrameCcAlf[compIdx - 1][shape][filterIdx].reset();
  }

  int                  ctuRsAddr = 0;
  const PreCalcValues &pcv = *cs.pcv;
  bool                 clipTop = false, clipBottom = false, clipLeft = false, clipRight = false;
  int                  numHorVirBndry = 0, numVerVirBndry = 0;
  int                  horVirBndryPos[] = { 0, 0, 0 };
  int                  verVirBndryPos[] = { 0, 0, 0 };

  for( int yPos = 0; yPos < m_picHeight; yPos += m_maxCUHeight )
  {
    for( int xPos = 0; xPos < m_picWidth; xPos += m_maxCUWidth )
    {
      const int width = ( xPos + m_maxCUWidth > m_picWidth ) ? ( m_picWidth - xPos ) : m_maxCUWidth;
      const int height = ( yPos + m_maxCUHeight > m_picHeight ) ? ( m_picHeight - yPos ) : m_maxCUHeight;
      int       rasterSliceAlfPad = 0;
      if( isCrossedByVirtualBoundaries( cs, xPos, yPos, width, height, clipTop, clipBottom, clipLeft, clipRight,
                                        numHorVirBndry, numVerVirBndry, horVirBndryPos, verVirBndryPos,
                                        rasterSliceAlfPad ) )
      {
        int yStart = yPos;
        for( int i = 0; i <= numHorVirBndry; i++ )
        {
          const int  yEnd = i == numHorVirBndry ? yPos + height : horVirBndryPos[i];
          const int  h = yEnd - yStart;
          const bool clipT = ( i == 0 && clipTop ) || ( i > 0 ) || ( yStart == 0 );
          const bool clipB = ( i == numHorVirBndry && clipBottom ) || ( i < numHorVirBndry ) || ( yEnd == pcv.lumaHeight );
          int        xStart = xPos;
          for( int j = 0; j <= numVerVirBndry; j++ )
          {
            const int  xEnd = j == numVerVirBndry ? xPos + width : verVirBndryPos[j];
            const int  w = xEnd - xStart;
            const bool clipL = ( j == 0 && clipLeft ) || ( j > 0 ) || ( xStart == 0 );
            const bool clipR = ( j == numVerVirBndry && clipRight ) || ( j < numVerVirBndry ) || ( xEnd == pcv.lumaWidth );
            const int  wBuf = w + ( clipL ? 0 : MAX_ALF_PADDING_SIZE ) + ( clipR ? 0 : MAX_ALF_PADDING_SIZE );
            const int  hBuf = h + ( clipT ? 0 : MAX_ALF_PADDING_SIZE ) + ( clipB ? 0 : MAX_ALF_PADDING_SIZE );
            PelUnitBuf recBuf = m_tempBuf2.subBuf( UnitArea( cs.area.chromaFormat, Area( 0, 0, wBuf, hBuf ) ) );
            recBuf.copyFrom( recYuv.subBuf(
              UnitArea( cs.area.chromaFormat, Area( xStart - ( clipL ? 0 : MAX_ALF_PADDING_SIZE ),
                                                    yStart - ( clipT ? 0 : MAX_ALF_PADDING_SIZE ), wBuf, hBuf ) ) ) );
            // pad top-left unavailable samples for raster slice
            if( xStart == xPos && yStart == yPos && ( rasterSliceAlfPad & 1 ) )
            {
              recBuf.padBorderPel( MAX_ALF_PADDING_SIZE, 1 );
            }

            // pad bottom-right unavailable samples for raster slice
            if( xEnd == xPos + width && yEnd == yPos + height && ( rasterSliceAlfPad & 2 ) )
            {
              recBuf.padBorderPel( MAX_ALF_PADDING_SIZE, 2 );
            }
            recBuf.extendBorderPel( MAX_ALF_PADDING_SIZE );
            recBuf = recBuf.subBuf( UnitArea(
              cs.area.chromaFormat, Area( clipL ? 0 : MAX_ALF_PADDING_SIZE, clipT ? 0 : MAX_ALF_PADDING_SIZE, w, h ) ) );

            const UnitArea area( m_chromaFormat, Area( 0, 0, w, h ) );
            const UnitArea areaDst( m_chromaFormat, Area( xStart, yStart, w, h ) );

            const ComponentID compID = ComponentID( compIdx );

            for( int shape = 0; shape != m_filterShapesCcAlf[compIdx - 1].size(); shape++ )
            {
              getBlkStatsCcAlf<alfWSSD>( m_alfCovarianceCcAlf[compIdx - 1][0][filterIdx][ctuRsAddr],
                                         m_filterShapesCcAlf[compIdx - 1][shape], orgYuv, recBuf, areaDst, area, compID, yPos );
              m_alfCovarianceFrameCcAlf[compIdx - 1][shape][filterIdx] +=
                m_alfCovarianceCcAlf[compIdx - 1][shape][filterIdx][ctuRsAddr];
            }

            xStart = xEnd;
          }

          yStart = yEnd;
        }
      }
      else
      {
        const UnitArea area( m_chromaFormat, Area( xPos, yPos, width, height ) );

        const ComponentID compID = ComponentID( compIdx );

        for( int shape = 0; shape != m_filterShapesCcAlf[compIdx - 1].size(); shape++ )
        {
          getBlkStatsCcAlf<alfWSSD>( m_alfCovarianceCcAlf[compIdx - 1][0][filterIdx][ctuRsAddr],
                                     m_filterShapesCcAlf[compIdx - 1][shape], orgYuv, recYuv, area, area, compID, yPos );
                                     m_alfCovarianceFrameCcAlf[compIdx - 1][shape][filterIdx] +=
                                     m_alfCovarianceCcAlf[compIdx - 1][shape][filterIdx][ctuRsAddr];
        }
      }
      ctuRsAddr++;
    }
  }
}

template<bool m_alfWSSD>
void EncAdaptiveLoopFilter::getBlkStatsCcAlf(AlfCovariance &alfCovariance, const AlfFilterShape &shape,
                                             const PelUnitBuf &orgYuv, const PelUnitBuf &recYuv,
                                             const UnitArea &areaDst, const UnitArea &area, const ComponentID compID,
                                             const int yPos)
{
  const int numberOfComponents = getNumberValidComponents( m_chromaFormat );
  const CompArea &compArea           = areaDst.block(compID);
  int  recStride[MAX_NUM_COMPONENT];
  const Pel* rec[MAX_NUM_COMPONENT];
  for ( int cIdx = 0; cIdx < numberOfComponents; cIdx++ )
  {
    recStride[cIdx] = recYuv.get(ComponentID(cIdx)).stride;
    rec[cIdx] = recYuv.get(ComponentID(cIdx)).bufAt(isLuma(ComponentID(cIdx)) ? area.lumaPos() : area.chromaPos());
  }

  int        orgStride = orgYuv.get(compID).stride;
  const Pel *org       = orgYuv.get(compID).bufAt(compArea);
  const int  numBins   = 1;

#if !ALF_IMPROVEMENT
  int vbCTUHeight = m_alfVBLumaCTUHeight;
  int vbPos       = m_alfVBLumaPos;
  if ((yPos + m_maxCUHeight) >= m_picHeight)
  {
    vbPos = m_picHeight;
  }
#endif

#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  Pel ELocal[MAX_NUM_CC_ALF_CHROMA_COEFF][1];
#else
  int ELocal[MAX_NUM_CC_ALF_CHROMA_COEFF][1];
#endif

  for (int i = 0; i < compArea.height; i++)
  {
#if !ALF_IMPROVEMENT
    int vbDistance = ((i << getComponentScaleY(compID, m_chromaFormat)) % vbCTUHeight) - vbPos;
    const bool skipThisRow = getComponentScaleY(compID, m_chromaFormat) == 0 && (vbDistance == 0 || vbDistance == 1);
    for (int j = 0; j < compArea.width && (!skipThisRow); j++)
#else
    for (int j = 0; j < compArea.width; j++)
#endif
    {
      std::memset(ELocal, 0, sizeof(ELocal));

      double weight = 1.0;
      if (m_alfWSSD)
      {
        weight = m_lumaLevelToWeightPLUT[org[j]];
      }

#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
      Intermediate_Int yLocal = org[j] - rec[compID][j];
#else
      int yLocal = org[j] - rec[compID][j];
#endif

#if ALF_IMPROVEMENT
      calcCovarianceCcAlf( ELocal, rec[COMPONENT_Y] + ( j << getComponentScaleX( compID, m_chromaFormat ) ), recStride[COMPONENT_Y], shape );
#else
      calcCovarianceCcAlf( ELocal, rec[COMPONENT_Y] + ( j << getComponentScaleX(compID, m_chromaFormat)), recStride[COMPONENT_Y], shape, vbDistance );
#endif

      for( int k = 0; k < (shape.numCoeff - 1); k++ )
      {
        for( int l = k; l < (shape.numCoeff - 1); l++ )
        {
          for( int b0 = 0; b0 < numBins; b0++ )
          {
            for (int b1 = 0; b1 < numBins; b1++)
            {
              if (m_alfWSSD)
              {
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
                alfCovariance.E[b0][b1][k][l] += weight * (ELocal[k][b0] * (double)ELocal[l][b1]);
#else
                alfCovariance.E[b0][b1][k][l] += weight * (double) (ELocal[k][b0] * ELocal[l][b1]);
#endif
              }
              else
              {
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
                alfCovariance.E[b0][b1][k][l] += ELocal[k][b0] * (double)ELocal[l][b1];
#else
                alfCovariance.E[b0][b1][k][l] += ELocal[k][b0] * ELocal[l][b1];
#endif
              }
            }
          }
        }
        for (int b = 0; b < numBins; b++)
        {
          if (m_alfWSSD)
          {
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
            alfCovariance.y[b][k] += weight * (ELocal[k][b] * (double)yLocal);
#else
            alfCovariance.y[b][k] += weight * (double) (ELocal[k][b] * yLocal);
#endif
          }
          else
          {
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
            alfCovariance.y[b][k] += ELocal[k][b] * (double)yLocal;
#else
            alfCovariance.y[b][k] += ELocal[k][b] * yLocal;
#endif
          }
        }
      }
      if (m_alfWSSD)
      {
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
        alfCovariance.pixAcc += weight * (yLocal * (double)yLocal);
#else
        alfCovariance.pixAcc += weight * (double) (yLocal * yLocal);
#endif
      }
      else
      {
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
        alfCovariance.pixAcc += yLocal * (double)yLocal;
#else
        alfCovariance.pixAcc += yLocal * yLocal;
#endif
      }
    }
    org += orgStride;
    for (int srcCIdx = 0; srcCIdx < numberOfComponents; srcCIdx++)
    {
      ComponentID srcCompID = ComponentID(srcCIdx);
      if (toChannelType(srcCompID) == toChannelType(compID))
      {
        rec[srcCIdx] += recStride[srcCIdx];
      }
      else
      {
        if (isLuma(compID))
        {
          rec[srcCIdx] += (recStride[srcCIdx] >> getComponentScaleY(srcCompID, m_chromaFormat));
        }
        else
        {
          rec[srcCIdx] += (recStride[srcCIdx] << getComponentScaleY(compID, m_chromaFormat));
        }
      }
    }
  }

  for (int k = 1; k < (MAX_NUM_CC_ALF_CHROMA_COEFF - 1); k++)
  {
    for (int l = 0; l < k; l++)
    {
      for (int b0 = 0; b0 < numBins; b0++)
      {
        for (int b1 = 0; b1 < numBins; b1++)
        {
          alfCovariance.E[b0][b1][k][l] = alfCovariance.E[b1][b0][l][k];
        }
      }
    }
  }
}
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
#if ALF_IMPROVEMENT
void EncAdaptiveLoopFilter::calcCovarianceCcAlf( Pel ELocal[MAX_NUM_CC_ALF_CHROMA_COEFF][1], const Pel *rec, const int stride, const AlfFilterShape& shape )
#else
void EncAdaptiveLoopFilter::calcCovarianceCcAlf(Pel ELocal[MAX_NUM_CC_ALF_CHROMA_COEFF][1], const Pel *rec, const int stride, const AlfFilterShape& shape, int vbDistance)
#endif
#else
#if ALF_IMPROVEMENT
void EncAdaptiveLoopFilter::calcCovarianceCcAlf( int ELocal[MAX_NUM_CC_ALF_CHROMA_COEFF][1], const Pel *rec, const int stride, const AlfFilterShape& shape )
#else
void EncAdaptiveLoopFilter::calcCovarianceCcAlf(int ELocal[MAX_NUM_CC_ALF_CHROMA_COEFF][1], const Pel *rec, const int stride, const AlfFilterShape& shape, int vbDistance)
#endif
#endif
{
  CHECK(shape.filterType != CC_ALF, "Bad CC ALF shape");

  const Pel *recYM1 = rec - 1 * stride;
  const Pel *recY0  = rec;
  const Pel *recYP1 = rec + 1 * stride;
  const Pel *recYP2 = rec + 2 * stride;

#if !ALF_IMPROVEMENT
  if (vbDistance == -2 || vbDistance == +1)
  {
    recYP2 = recYP1;
  }
  else if (vbDistance == -1 || vbDistance == 0)
  {
    recYM1 = recY0;
    recYP2 = recYP1 = recY0;
  }
#endif

  for (int b = 0; b < 1; b++)
  {
    const Pel centerValue = recY0[+0];
    ELocal[0][b] += recYM1[+0] - centerValue;
    ELocal[1][b] += recY0[-1] - centerValue;
    ELocal[2][b] += recY0[+1] - centerValue;
    ELocal[3][b] += recYP1[-1] - centerValue;
    ELocal[4][b] += recYP1[+0] - centerValue;
    ELocal[5][b] += recYP1[+1] - centerValue;
    ELocal[6][b] += recYP2[+0] - centerValue;
  }
}

void EncAdaptiveLoopFilter::countLumaSwingGreaterThanThreshold(const Pel* luma, int lumaStride, int height, int width, int log2BlockWidth, int log2BlockHeight, uint64_t* lumaSwingGreaterThanThresholdCount, int lumaCountStride)
{
  const int lumaBitDepth = m_inputBitDepth[CH_L];
  const int threshold = (1 << ( m_inputBitDepth[CH_L] - 2 )) - 1;

  // 3x4 Diamond
  int xSupport[] = {  0, -1, 0, 1, -1, 0, 1, 0 };
  int ySupport[] = { -1,  0, 0, 0,  1, 1, 1, 2 };

  for (int y = 0; y < height; y += (1 << log2BlockHeight))
  {
    for (int x = 0; x < width; x += (1 << log2BlockWidth))
    {
      lumaSwingGreaterThanThresholdCount[(y >> log2BlockHeight) * lumaCountStride + (x >> log2BlockWidth)] = 0;

      for (int yOff = 0; yOff < (1 << log2BlockHeight); yOff++)
      {
        for (int xOff = 0; xOff < (1 << log2BlockWidth); xOff++)
        {
          if ((y + yOff) >= (height - 2) || (x + xOff) >= (width - 1) || (y + yOff) < 1 || (x + xOff) < 1) // only consider samples that are fully supported by picture
          {
            continue;
          }

          int minVal = ((1 << lumaBitDepth) - 1);
          int maxVal = 0;
          for (int i = 0; i < 8; i++)
          {
            Pel p = luma[(yOff + ySupport[i]) * lumaStride + x + xOff + xSupport[i]];

            if ( p < minVal )
            {
              minVal = p;
            }
            if ( p > maxVal )
            {
              maxVal = p;
            }
          }

          if ((maxVal - minVal) > threshold)
          {
            lumaSwingGreaterThanThresholdCount[(y >> log2BlockHeight) * lumaCountStride + (x >> log2BlockWidth)]++;
          }
        }
      }
    }
    luma += (lumaStride << log2BlockHeight);
  }
}

void EncAdaptiveLoopFilter::countChromaSampleValueNearMidPoint(const Pel* chroma, int chromaStride, int height, int width, int log2BlockWidth, int log2BlockHeight, uint64_t* chromaSampleCountNearMidPoint, int chromaSampleCountNearMidPointStride)
{
  const int midPoint  = (1 << m_inputBitDepth[CH_C]) >> 1;
  const int threshold = 16;

  for (int y = 0; y < height; y += (1 << log2BlockHeight))
  {
    for (int x = 0; x < width; x += (1 << log2BlockWidth))
    {
      chromaSampleCountNearMidPoint[(y >> log2BlockHeight)* chromaSampleCountNearMidPointStride + (x >> log2BlockWidth)] = 0;

      for (int yOff = 0; yOff < (1 << log2BlockHeight); yOff++)
      {
        for (int xOff = 0; xOff < (1 << log2BlockWidth); xOff++)
        {
          if ((y + yOff) >= height || (x + xOff) >= width)
          {
            continue;
          }

          int distanceToMidPoint = abs(chroma[yOff * chromaStride + x + xOff] - midPoint);
          if (distanceToMidPoint < threshold)
          {
            chromaSampleCountNearMidPoint[(y >> log2BlockHeight)* chromaSampleCountNearMidPointStride + (x >> log2BlockWidth)]++;
          }
        }
      }
    }
    chroma += (chromaStride << log2BlockHeight);
  }
}

