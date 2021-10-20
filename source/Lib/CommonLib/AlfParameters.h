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

/** \file     AlfParameters.h
    \brief    Define types for storing ALF parameters
*/

#ifndef __ALFPARAMETERS__
#define __ALFPARAMETERS__

#include <vector>
#include "CommonDef.h"

//! \ingroup AlfParameters
//! \{

enum AlfFilterType
{
  ALF_FILTER_5,
  ALF_FILTER_7,
  CC_ALF,

#if ALF_IMPROVEMENT
  ALF_FILTER_9,
  ALF_FILTER_9_EXT,
  ALF_FILTER_EXT,
#endif
  ALF_NUM_OF_FILTER_TYPES
};


static const int size_CC_ALF = -1;
#if ALF_IMPROVEMENT
static const int size_ALF_FILTER_9_EXT = -2;
static const int size_ALF_FILTER_EXT = -3;
#endif

struct AlfFilterShape
{
  AlfFilterShape( int size )
    : filterLength( size ),
    numCoeff( size * size / 4 + 1 ),
    filterSize( size * size / 2 + 1 )
#if ALF_IMPROVEMENT
    , numOrder( 0 )
    , indexSecOrder( 0 )
    , offset0( 0 )
#endif
  {
    if( size == 5 )
    {
      pattern = {
                 0,
             1,  2,  3,
         4,  5,  6,  5,  4,
             3,  2,  1,
                 0
      };

      weights = {
                 2,
              2, 2, 2,
           2, 2, 1, 1
      };

      filterType = ALF_FILTER_5;
#if ALF_IMPROVEMENT
      numOrder = 1;
      indexSecOrder = numCoeff - 1;
      offset0 = 0;
#endif
    }
    else if( size == 7 )
    {
      pattern = {
                     0,
                 1,  2,  3,
             4,  5,  6,  7,  8,
         9, 10, 11, 12, 11, 10, 9,
             8,  7,  6,  5,  4,
                 3,  2,  1,
                     0
      };

      weights = {
                    2,
                2,  2,  2,
            2,  2,  2,  2,  2,
        2,  2,  2,  1,  1
      };

      filterType = ALF_FILTER_7;
#if ALF_IMPROVEMENT
      numOrder = 1;
      indexSecOrder = numCoeff - 1;
      offset0 = 0;
#endif
    }
#if ALF_IMPROVEMENT
    else if( size == 9 )
    {
      pattern = {
                     0,
                 1,  2,  3,
             4,  5,  6,  7,  8,
         9, 10, 11, 12, 13, 14, 15,
    16, 17, 18, 19, 20, 19, 18, 17, 16,
        15, 14, 13, 12, 11, 10,  9, 
             8,  7,  6,  5,  4,
                 3,  2,  1,
                     0
      };

      weights = {
                    2,
                2,  2,  2,
            2,  2,  2,  2,  2,
        2,  2,  2,  2,  2,  2,  2,
    2,  2,  2,  2,  1,  1
      };

      filterType = ALF_FILTER_9;
#if ALF_IMPROVEMENT
      numOrder = 1;
      indexSecOrder = numCoeff - 1;
      offset0 = 0;
#endif
    }
    else if( size == size_ALF_FILTER_EXT )
    {
      size = filterLength = 0;
      numCoeff = filterSize = EXT_LENGTH;
      filterType = ALF_FILTER_EXT;
#if ALF_IMPROVEMENT
      numOrder = 1;
      indexSecOrder = numCoeff - 1;
      offset0 = ALF_ORDER;
#endif
    }
    else if( size == size_ALF_FILTER_9_EXT )
    {
      size = 9;
      numCoeff = MAX_NUM_ALF_LUMA_COEFF;
      filterSize = MAX_NUM_ALF_LUMA_COEFF;
      filterLength = 9;
      filterType = ALF_FILTER_9_EXT;
      pattern = {
                     0,
                 1,  2,  3,
             4,  5,  6,  7,  8,
         9, 10, 11, 12, 13, 14, 15,
    16, 17, 18, 19, 20, 19, 18, 17, 16,
        15, 14, 13, 12, 11, 10,  9,
             8,  7,  6,  5,  4,
                 3,  2,  1,
                     0
      };
#if ALF_IMPROVEMENT
      numOrder = 2;
      indexSecOrder = 20;
      offset0 = 0;
#endif
    }
#endif
    else if( size == size_CC_ALF )
    {
      size = 4;
      filterLength = 8;
      numCoeff = 8;
      filterSize = 8;
      filterType   = CC_ALF;
    }
    else
    {
      filterType = ALF_NUM_OF_FILTER_TYPES;
      CHECK( 0, "Wrong ALF filter shape" );
    }
  }

  AlfFilterType filterType;
  int filterLength;
  int numCoeff;      //TO DO: check whether we need both numCoeff and filterSize
  int filterSize;
  std::vector<int> pattern;
  std::vector<int> weights;
#if ALF_IMPROVEMENT
  int numOrder;
  int indexSecOrder;
  int offset0;
#endif
};

struct AlfParam
{
  bool                         enabledFlag[MAX_NUM_COMPONENT];                          // alf_slice_enable_flag, alf_chroma_idc
#if ALF_IMPROVEMENT
  AlfFilterType                filterType[MAX_NUM_CHANNEL_TYPE];
  bool                         nonLinearFlag[MAX_NUM_CHANNEL_TYPE][32]; // alf_[luma/chroma]_clip_flag
  int                          numAlternativesLuma;
  short                        lumaCoeff[MAX_NUM_ALF_ALTERNATIVES_LUMA][MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF]; // alf_coeff_luma_delta[i][j]
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  Pel                          lumaClipp[MAX_NUM_ALF_ALTERNATIVES_LUMA][MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF]; // alf_clipp_luma_[i][j]
#else  
  short                        lumaClipp[MAX_NUM_ALF_ALTERNATIVES_LUMA][MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF]; // alf_clipp_luma_[i][j]
#endif  
  int                          numLumaFilters[MAX_NUM_ALF_ALTERNATIVES_LUMA];                                          // number_of_filters_minus1 + 1
  short                        filterCoeffDeltaIdx[MAX_NUM_ALF_ALTERNATIVES_LUMA][MAX_NUM_ALF_CLASSES];                // filter_coeff_delta[i]// number_of_filters_minus1 + 1
#else
  bool                         nonLinearFlag[MAX_NUM_CHANNEL_TYPE];                     // alf_[luma/chroma]_clip_flag
  short                        lumaCoeff[MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF]; // alf_coeff_luma_delta[i][j]
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  Pel                          lumaClipp[MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF]; // alf_clipp_luma_[i][j]
#else
  short                        lumaClipp[MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF]; // alf_clipp_luma_[i][j]
#endif
  int                          numLumaFilters;                                          // number_of_filters_minus1 + 1
  short                        filterCoeffDeltaIdx[MAX_NUM_ALF_CLASSES];                // filter_coeff_delta[i]
#endif
  int                          numAlternativesChroma;                                                  // alf_chroma_num_alts_minus_one + 1
  short                        chromaCoeff[MAX_NUM_ALF_ALTERNATIVES_CHROMA][MAX_NUM_ALF_CHROMA_COEFF]; // alf_coeff_chroma[i]
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
  Pel                          chromaClipp[MAX_NUM_ALF_ALTERNATIVES_CHROMA][MAX_NUM_ALF_CHROMA_COEFF]; // alf_clipp_chroma[i]
#else
  short                        chromaClipp[MAX_NUM_ALF_ALTERNATIVES_CHROMA][MAX_NUM_ALF_CHROMA_COEFF]; // alf_clipp_chroma[i]
#endif
  bool                         alfLumaCoeffFlag[MAX_NUM_ALF_CLASSES];                   // alf_luma_coeff_flag[i]
  bool                         alfLumaCoeffDeltaFlag;                                   // alf_luma_coeff_delta_flag
  std::vector<AlfFilterShape>* filterShapes;
  bool                         newFilterFlag[MAX_NUM_CHANNEL_TYPE];

  AlfParam()
  {
    reset();
  }

  void reset()
  {
    std::memset( enabledFlag, false, sizeof( enabledFlag ) );
    std::memset( nonLinearFlag, false, sizeof( nonLinearFlag ) );
    std::memset( lumaCoeff, 0, sizeof( lumaCoeff ) );
    std::memset( lumaClipp, 0, sizeof( lumaClipp ) );
    numAlternativesChroma = 1;
    std::memset( chromaCoeff, 0, sizeof( chromaCoeff ) );
    std::memset( chromaClipp, 0, sizeof( chromaClipp ) );
    std::memset( filterCoeffDeltaIdx, 0, sizeof( filterCoeffDeltaIdx ) );
    std::memset( alfLumaCoeffFlag, true, sizeof( alfLumaCoeffFlag ) );
#if ALF_IMPROVEMENT
    std::memset(filterType, false, sizeof(filterType));
    numAlternativesLuma = 1;
    for (int i = 0; i < MAX_NUM_ALF_ALTERNATIVES_LUMA; i++)
    {
      numLumaFilters[i] = 1;
    }
#else
    numLumaFilters = 1;
#endif
    alfLumaCoeffDeltaFlag = false;
    memset(newFilterFlag, 0, sizeof(newFilterFlag));
  }

  const AlfParam& operator = ( const AlfParam& src )
  {
    std::memcpy( enabledFlag, src.enabledFlag, sizeof( enabledFlag ) );
    std::memcpy( nonLinearFlag, src.nonLinearFlag, sizeof( nonLinearFlag ) );
    std::memcpy( lumaCoeff, src.lumaCoeff, sizeof( lumaCoeff ) );
    std::memcpy( lumaClipp, src.lumaClipp, sizeof( lumaClipp ) );
    numAlternativesChroma = src.numAlternativesChroma;
    std::memcpy( chromaCoeff, src.chromaCoeff, sizeof( chromaCoeff ) );
    std::memcpy( chromaClipp, src.chromaClipp, sizeof( chromaClipp ) );
    std::memcpy( filterCoeffDeltaIdx, src.filterCoeffDeltaIdx, sizeof( filterCoeffDeltaIdx ) );
    std::memcpy( alfLumaCoeffFlag, src.alfLumaCoeffFlag, sizeof( alfLumaCoeffFlag ) );
#if ALF_IMPROVEMENT
    std::memcpy(filterType, src.filterType, sizeof(filterType));
    numAlternativesLuma = src.numAlternativesLuma;
    std::memcpy(numLumaFilters, src.numLumaFilters, sizeof(numLumaFilters));
#else
    numLumaFilters = src.numLumaFilters;
#endif
    alfLumaCoeffDeltaFlag = src.alfLumaCoeffDeltaFlag;
    filterShapes = src.filterShapes;
    std::memcpy(newFilterFlag, src.newFilterFlag, sizeof(newFilterFlag));
    return *this;
  }

  bool operator==( const AlfParam& other )
  {
    if( memcmp( enabledFlag, other.enabledFlag, sizeof( enabledFlag ) ) )
    {
      return false;
    }
    if( memcmp( nonLinearFlag, other.nonLinearFlag, sizeof( nonLinearFlag ) ) )
    {
      return false;
    }
#if ALF_IMPROVEMENT
    if( memcmp( filterType, other.filterType, sizeof( filterType ) ) )
    {
      return false;
    }
#endif
    if( memcmp( lumaCoeff, other.lumaCoeff, sizeof( lumaCoeff ) ) )
    {
      return false;
    }
    if( memcmp( lumaClipp, other.lumaClipp, sizeof( lumaClipp ) ) )
    {
      return false;
    }
    if( memcmp( chromaCoeff, other.chromaCoeff, sizeof( chromaCoeff ) ) )
    {
      return false;
    }
    if( memcmp( chromaClipp, other.chromaClipp, sizeof( chromaClipp ) ) )
    {
      return false;
    }
    if( memcmp( filterCoeffDeltaIdx, other.filterCoeffDeltaIdx, sizeof( filterCoeffDeltaIdx ) ) )
    {
      return false;
    }
    if( memcmp( alfLumaCoeffFlag, other.alfLumaCoeffFlag, sizeof( alfLumaCoeffFlag ) ) )
    {
      return false;
    }
    if( memcmp( newFilterFlag, other.newFilterFlag, sizeof( newFilterFlag ) ) )
    {
      return false;
    }
    if( numAlternativesChroma != other.numAlternativesChroma )
    {
      return false;
    }
#if ALF_IMPROVEMENT
    if (numAlternativesLuma != other.numAlternativesLuma)
    {
      return false;
    }
    if (memcmp(numLumaFilters, other.numLumaFilters, sizeof(numLumaFilters)))
    {
      return false;
    }
#else
    if( numLumaFilters != other.numLumaFilters )
    {
      return false;
    }
#endif
    if( alfLumaCoeffDeltaFlag != other.alfLumaCoeffDeltaFlag )
    {
      return false;
    }

    return true;
  }

  bool operator!=( const AlfParam& other )
  {
    return !( *this == other );
  }
};

struct CcAlfFilterParam
{
  bool    ccAlfFilterEnabled[2];
  bool    ccAlfFilterIdxEnabled[2][MAX_NUM_CC_ALF_FILTERS];
  uint8_t ccAlfFilterCount[2];
  short   ccAlfCoeff[2][MAX_NUM_CC_ALF_FILTERS][MAX_NUM_CC_ALF_CHROMA_COEFF];
  int     newCcAlfFilter[2];
  int     numberValidComponents;
  CcAlfFilterParam()
  {
    reset();
  }
  void reset()
  {
    std::memset( ccAlfFilterEnabled, false, sizeof( ccAlfFilterEnabled ) );
    std::memset( ccAlfFilterIdxEnabled, false, sizeof( ccAlfFilterIdxEnabled ) );
    std::memset( ccAlfCoeff, 0, sizeof( ccAlfCoeff ) );
    ccAlfFilterCount[0] = ccAlfFilterCount[1] = MAX_NUM_CC_ALF_FILTERS;
    numberValidComponents = 3;
    newCcAlfFilter[0] = newCcAlfFilter[1] = 0;
  }
  const CcAlfFilterParam& operator = ( const CcAlfFilterParam& src )
  {
    std::memcpy( ccAlfFilterEnabled, src.ccAlfFilterEnabled, sizeof( ccAlfFilterEnabled ) );
    std::memcpy( ccAlfFilterIdxEnabled, src.ccAlfFilterIdxEnabled, sizeof( ccAlfFilterIdxEnabled ) );
    std::memcpy( ccAlfCoeff, src.ccAlfCoeff, sizeof( ccAlfCoeff ) );
    ccAlfFilterCount[0] = src.ccAlfFilterCount[0];
    ccAlfFilterCount[1] = src.ccAlfFilterCount[1];
    numberValidComponents = src.numberValidComponents;
    newCcAlfFilter[0] = src.newCcAlfFilter[0];
    newCcAlfFilter[1] = src.newCcAlfFilter[1];

    return *this;
  }
};

#if JVET_W0066_CCSAO
struct CcSaoComParam
{
  bool     enabled   [MAX_NUM_COMPONENT];
  uint8_t  setNum    [MAX_NUM_COMPONENT];
  bool     setEnabled[MAX_NUM_COMPONENT][MAX_CCSAO_SET_NUM];
  uint16_t candPos   [MAX_NUM_COMPONENT][MAX_CCSAO_SET_NUM][MAX_NUM_LUMA_COMP];
  uint16_t bandNum   [MAX_NUM_COMPONENT][MAX_CCSAO_SET_NUM][MAX_NUM_COMPONENT];
  short    offset    [MAX_NUM_COMPONENT][MAX_CCSAO_SET_NUM][MAX_CCSAO_CLASS_NUM];
  CcSaoComParam()
  {
    reset();
  }
  void reset()
  {
    std::memset( enabled,    false, sizeof( enabled    ) );
    std::memset( setNum,         0, sizeof( setNum     ) );
    std::memset( setEnabled, false, sizeof( setEnabled ) );
    std::memset( candPos,        0, sizeof( candPos    ) );
    std::memset( bandNum,        0, sizeof( bandNum    ) );
    std::memset( offset,         0, sizeof( offset     ) );
  }
  void reset(ComponentID compID)
  {
    enabled[compID] = false;
    setNum [compID] = 0;
    std::memset( setEnabled[compID], false, sizeof( setEnabled[compID]) );
    std::memset( candPos   [compID],     0, sizeof( candPos   [compID]) );
    std::memset( bandNum   [compID],     0, sizeof( bandNum   [compID]) );
    std::memset( offset    [compID],     0, sizeof( offset    [compID]) );
  }
  const CcSaoComParam& operator = ( const CcSaoComParam& src )
  {
    std::memcpy( enabled,    src.enabled,    sizeof( enabled    ) );
    std::memcpy( setNum,     src.setNum,     sizeof( setNum     ) );
    std::memcpy( setEnabled, src.setEnabled, sizeof( setEnabled ) );
    std::memcpy( candPos,    src.candPos,    sizeof( candPos    ) );
    std::memcpy( bandNum,    src.bandNum,    sizeof( bandNum    ) );
    std::memcpy( offset,     src.offset,     sizeof( offset     ) );
    return *this;
  }
};
#endif
//! \}

#endif  // end of #ifndef  __ALFPARAMETERS__
