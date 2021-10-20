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

/** \file     SampleAdaptiveOffset.cpp
    \brief    sample adaptive offset class
*/

#include "SampleAdaptiveOffset.h"

#include "UnitTools.h"
#include "UnitPartitioner.h"
#include "CodingStructure.h"
#include "CommonLib/dtrace_codingstruct.h"
#include "CommonLib/dtrace_buffer.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//! \ingroup CommonLib
//! \{

SAOOffset::SAOOffset()
{
  reset();
}

SAOOffset::~SAOOffset()
{

}

void SAOOffset::reset()
{
  modeIdc = SAO_MODE_OFF;
  typeIdc = -1;
  typeAuxInfo = -1;
  ::memset(offset, 0, sizeof(int)* MAX_NUM_SAO_CLASSES);
}

const SAOOffset& SAOOffset::operator= (const SAOOffset& src)
{
  modeIdc = src.modeIdc;
  typeIdc = src.typeIdc;
  typeAuxInfo = src.typeAuxInfo;
  ::memcpy(offset, src.offset, sizeof(int)* MAX_NUM_SAO_CLASSES);

  return *this;
}


SAOBlkParam::SAOBlkParam()
{
  reset();
}

SAOBlkParam::~SAOBlkParam()
{

}

void SAOBlkParam::reset()
{
  for(int compIdx = 0; compIdx < MAX_NUM_COMPONENT; compIdx++)
  {
    offsetParam[compIdx].reset();
  }
}

const SAOBlkParam& SAOBlkParam::operator= (const SAOBlkParam& src)
{
  for(int compIdx = 0; compIdx < MAX_NUM_COMPONENT; compIdx++)
  {
    offsetParam[compIdx] = src.offsetParam[compIdx];
  }
  return *this;
}




SampleAdaptiveOffset::SampleAdaptiveOffset()
{
  m_numberOfComponents = 0;
}


SampleAdaptiveOffset::~SampleAdaptiveOffset()
{
  destroy();

  m_signLineBuf1.clear();
  m_signLineBuf2.clear();
}

void SampleAdaptiveOffset::create( int picWidth, int picHeight, ChromaFormat format, uint32_t maxCUWidth, uint32_t maxCUHeight, uint32_t maxCUDepth, uint32_t lumaBitShift, uint32_t chromaBitShift )
{
  //temporary picture buffer
  UnitArea picArea(format, Area(0, 0, picWidth, picHeight));

  m_tempBuf.destroy();
  m_tempBuf.create( picArea );

  //bit-depth related
  for(int compIdx = 0; compIdx < MAX_NUM_COMPONENT; compIdx++)
  {
    m_offsetStepLog2  [compIdx] = isLuma(ComponentID(compIdx))? lumaBitShift : chromaBitShift;
  }
  m_numberOfComponents = getNumberValidComponents(format);

#if JVET_W0066_CCSAO
  if (m_created)
  {
    return;
  }
  m_created = true;

  m_ccSaoBuf.destroy();
  m_ccSaoBuf.create(format, Area(0, 0, picWidth, picHeight), maxCUWidth, MAX_CCSAO_FILTER_LENGTH >> 1, 0, false);

  m_picWidth = picWidth;
  m_picHeight = picHeight;
  m_maxCUWidth = maxCUWidth;
  m_maxCUHeight = maxCUHeight;

  m_numCTUsInWidth = ( m_picWidth / m_maxCUWidth ) + ( ( m_picWidth % m_maxCUWidth ) ? 1 : 0 );
  m_numCTUsInHeight = ( m_picHeight / m_maxCUHeight ) + ( ( m_picHeight % m_maxCUHeight ) ? 1 : 0 );
  m_numCTUsInPic = m_numCTUsInHeight * m_numCTUsInWidth;

  for (int compIdx = 0; compIdx < MAX_NUM_COMPONENT; compIdx++)
  {
    m_ccSaoControl[compIdx] = new uint8_t[m_numCTUsInPic];
    ::memset(m_ccSaoControl[compIdx], 0, sizeof(uint8_t) * m_numCTUsInPic);
  }
#endif
}

void SampleAdaptiveOffset::destroy()
{
  m_tempBuf.destroy();

#if JVET_W0066_CCSAO
  if (!m_created)
  {
    return;
  }
  m_created = false;

  m_ccSaoBuf.destroy();

  for (int compIdx = 0; compIdx < MAX_NUM_COMPONENT; compIdx++)
  {
    if (m_ccSaoControl[compIdx]) 
    { 
      delete [] m_ccSaoControl[compIdx];
      m_ccSaoControl[compIdx] = nullptr;
    }
  }
#endif
}

void SampleAdaptiveOffset::invertQuantOffsets(ComponentID compIdx, int typeIdc, int typeAuxInfo, int* dstOffsets, int* srcOffsets)
{
  int codedOffset[MAX_NUM_SAO_CLASSES];

  ::memcpy(codedOffset, srcOffsets, sizeof(int)*MAX_NUM_SAO_CLASSES);
  ::memset(dstOffsets, 0, sizeof(int)*MAX_NUM_SAO_CLASSES);

  if(typeIdc == SAO_TYPE_START_BO)
  {
    for(int i=0; i< 4; i++)
    {
      dstOffsets[(typeAuxInfo+ i)%NUM_SAO_BO_CLASSES] = codedOffset[(typeAuxInfo+ i)%NUM_SAO_BO_CLASSES]*(1<<m_offsetStepLog2[compIdx]);
    }
  }
  else //EO
  {
    for(int i=0; i< NUM_SAO_EO_CLASSES; i++)
    {
      dstOffsets[i] = codedOffset[i] *(1<<m_offsetStepLog2[compIdx]);
    }
    CHECK(dstOffsets[SAO_CLASS_EO_PLAIN] != 0, "EO offset is not '0'"); //keep EO plain offset as zero
  }

}

int SampleAdaptiveOffset::getMergeList(CodingStructure& cs, int ctuRsAddr, SAOBlkParam* blkParams, SAOBlkParam* mergeList[NUM_SAO_MERGE_TYPES])
{
  const PreCalcValues& pcv = *cs.pcv;

  int ctuX = ctuRsAddr % pcv.widthInCtus;
  int ctuY = ctuRsAddr / pcv.widthInCtus;
  const CodingUnit& cu = *cs.getCU(Position(ctuX*pcv.maxCUWidth, ctuY*pcv.maxCUHeight), CH_L);
  int mergedCTUPos;
  int numValidMergeCandidates = 0;

  for(int mergeType=0; mergeType< NUM_SAO_MERGE_TYPES; mergeType++)
  {
    SAOBlkParam* mergeCandidate = NULL;

    switch(mergeType)
    {
    case SAO_MERGE_ABOVE:
      {
        if(ctuY > 0)
        {
          mergedCTUPos = ctuRsAddr- pcv.widthInCtus;
          if(cs.getCURestricted(Position(ctuX*pcv.maxCUWidth, (ctuY-1)*pcv.maxCUHeight), cu, cu.chType))
          {
            mergeCandidate = &(blkParams[mergedCTUPos]);
          }
        }
      }
      break;
    case SAO_MERGE_LEFT:
      {
        if(ctuX > 0)
        {
          mergedCTUPos = ctuRsAddr- 1;
          if(cs.getCURestricted(Position((ctuX-1)*pcv.maxCUWidth, ctuY*pcv.maxCUHeight), cu, cu.chType))
          {
            mergeCandidate = &(blkParams[mergedCTUPos]);
          }
        }
      }
      break;
    default:
      {
        THROW("not a supported merge type");
      }
    }

    mergeList[mergeType]=mergeCandidate;
    if (mergeCandidate != NULL)
    {
      numValidMergeCandidates++;
    }
  }

  return numValidMergeCandidates;
}


void SampleAdaptiveOffset::reconstructBlkSAOParam(SAOBlkParam& recParam, SAOBlkParam* mergeList[NUM_SAO_MERGE_TYPES])
{
  const int numberOfComponents = m_numberOfComponents;
  for(int compIdx = 0; compIdx < numberOfComponents; compIdx++)
  {
    const ComponentID component = ComponentID(compIdx);
    SAOOffset& offsetParam = recParam[component];

    if(offsetParam.modeIdc == SAO_MODE_OFF)
    {
      continue;
    }

    switch(offsetParam.modeIdc)
    {
    case SAO_MODE_NEW:
      {
        invertQuantOffsets(component, offsetParam.typeIdc, offsetParam.typeAuxInfo, offsetParam.offset, offsetParam.offset);
      }
      break;
    case SAO_MODE_MERGE:
      {
        SAOBlkParam* mergeTarget = mergeList[offsetParam.typeIdc];
        CHECK(mergeTarget == NULL, "Merge target does not exist");

        offsetParam = (*mergeTarget)[component];
      }
      break;
    default:
      {
        THROW("Not a supported mode");
      }
    }
  }
}

void SampleAdaptiveOffset::xReconstructBlkSAOParams(CodingStructure& cs, SAOBlkParam* saoBlkParams)
{
  for(uint32_t compIdx = 0; compIdx < MAX_NUM_COMPONENT; compIdx++)
  {
    m_picSAOEnabled[compIdx] = false;
  }

  const uint32_t numberOfComponents = getNumberValidComponents(cs.pcv->chrFormat);

  for(int ctuRsAddr=0; ctuRsAddr< cs.pcv->sizeInCtus; ctuRsAddr++)
  {
    SAOBlkParam* mergeList[NUM_SAO_MERGE_TYPES] = { NULL };
    getMergeList(cs, ctuRsAddr, saoBlkParams, mergeList);

    reconstructBlkSAOParam(saoBlkParams[ctuRsAddr], mergeList);

    for(uint32_t compIdx = 0; compIdx < numberOfComponents; compIdx++)
    {
      if(saoBlkParams[ctuRsAddr][compIdx].modeIdc != SAO_MODE_OFF)
      {
        m_picSAOEnabled[compIdx] = true;
      }
    }
  }
}


void SampleAdaptiveOffset::offsetBlock(const int channelBitDepth, const ClpRng& clpRng, int typeIdx, int* offset
                                          , const Pel* srcBlk, Pel* resBlk, int srcStride, int resStride,  int width, int height
                                          , bool isLeftAvail,  bool isRightAvail, bool isAboveAvail, bool isBelowAvail, bool isAboveLeftAvail, bool isAboveRightAvail, bool isBelowLeftAvail, bool isBelowRightAvail
                                          , bool isCtuCrossedByVirtualBoundaries, int horVirBndryPos[], int verVirBndryPos[], int numHorVirBndry, int numVerVirBndry
  )
{
  int x,y, startX, startY, endX, endY, edgeType;
  int firstLineStartX, firstLineEndX, lastLineStartX, lastLineEndX;
  int8_t signLeft, signRight, signDown;

  const Pel* srcLine = srcBlk;
        Pel* resLine = resBlk;

  switch(typeIdx)
  {
  case SAO_TYPE_EO_0:
    {
      offset += 2;
      startX = isLeftAvail ? 0 : 1;
      endX   = isRightAvail ? width : (width -1);
      for (y=0; y< height; y++)
      {
        signLeft = (int8_t)sgn(srcLine[startX] - srcLine[startX-1]);
        for (x=startX; x< endX; x++)
        {
          signRight = (int8_t)sgn(srcLine[x] - srcLine[x+1]);
          if (isCtuCrossedByVirtualBoundaries && isProcessDisabled(x, y, numVerVirBndry, 0, verVirBndryPos, horVirBndryPos))
          {
            signLeft = -signRight;
            continue;
          }
          edgeType =  signRight + signLeft;
          signLeft  = -signRight;

          resLine[x] = ClipPel<int>( srcLine[x] + offset[edgeType], clpRng);
        }
        srcLine  += srcStride;
        resLine += resStride;
      }

    }
    break;
  case SAO_TYPE_EO_90:
    {
      offset += 2;
      int8_t *signUpLine = &m_signLineBuf1[0];

      startY = isAboveAvail ? 0 : 1;
      endY   = isBelowAvail ? height : height-1;
      if (!isAboveAvail)
      {
        srcLine += srcStride;
        resLine += resStride;
      }

      const Pel* srcLineAbove= srcLine- srcStride;
      for (x=0; x< width; x++)
      {
        signUpLine[x] = (int8_t)sgn(srcLine[x] - srcLineAbove[x]);
      }

      const Pel* srcLineBelow;
      for (y=startY; y<endY; y++)
      {
        srcLineBelow= srcLine+ srcStride;

        for (x=0; x< width; x++)
        {
          signDown  = (int8_t)sgn(srcLine[x] - srcLineBelow[x]);
          if (isCtuCrossedByVirtualBoundaries && isProcessDisabled(x, y, 0, numHorVirBndry, verVirBndryPos, horVirBndryPos))
          {
            signUpLine[x] = -signDown;
            continue;
          }
          edgeType = signDown + signUpLine[x];
          signUpLine[x]= -signDown;

          resLine[x] = ClipPel<int>(srcLine[x] + offset[edgeType], clpRng);
        }
        srcLine += srcStride;
        resLine += resStride;
      }

    }
    break;
  case SAO_TYPE_EO_135:
    {
      offset += 2;
      int8_t *signUpLine, *signDownLine, *signTmpLine;

      signUpLine  = &m_signLineBuf1[0];
      signDownLine= &m_signLineBuf2[0];

      startX = isLeftAvail ? 0 : 1 ;
      endX   = isRightAvail ? width : (width-1);

      //prepare 2nd line's upper sign
      const Pel* srcLineBelow= srcLine+ srcStride;
      for (x=startX; x< endX+1; x++)
      {
        signUpLine[x] = (int8_t)sgn(srcLineBelow[x] - srcLine[x- 1]);
      }

      //1st line
      const Pel* srcLineAbove= srcLine- srcStride;
      firstLineStartX = isAboveLeftAvail ? 0 : 1;
      firstLineEndX   = isAboveAvail? endX: 1;
      for(x= firstLineStartX; x< firstLineEndX; x++)
      {
        if (isCtuCrossedByVirtualBoundaries && isProcessDisabled(x, 0, numVerVirBndry, numHorVirBndry, verVirBndryPos, horVirBndryPos))
        {
          continue;
        }
        edgeType  =  sgn(srcLine[x] - srcLineAbove[x- 1]) - signUpLine[x+1];

        resLine[x] = ClipPel<int>( srcLine[x] + offset[edgeType], clpRng);
      }
      srcLine  += srcStride;
      resLine  += resStride;


      //middle lines
      for (y= 1; y< height-1; y++)
      {
        srcLineBelow= srcLine+ srcStride;

        for (x=startX; x<endX; x++)
        {
          signDown =  (int8_t)sgn(srcLine[x] - srcLineBelow[x+ 1]);
          if (isCtuCrossedByVirtualBoundaries && isProcessDisabled(x, y, numVerVirBndry, numHorVirBndry, verVirBndryPos, horVirBndryPos))
          {
            signDownLine[x + 1] = -signDown;
            continue;
          }
          edgeType =  signDown + signUpLine[x];
          resLine[x] = ClipPel<int>( srcLine[x] + offset[edgeType], clpRng);

          signDownLine[x+1] = -signDown;
        }
        signDownLine[startX] = (int8_t)sgn(srcLineBelow[startX] - srcLine[startX-1]);

        signTmpLine  = signUpLine;
        signUpLine   = signDownLine;
        signDownLine = signTmpLine;

        srcLine += srcStride;
        resLine += resStride;
      }

      //last line
      srcLineBelow= srcLine+ srcStride;
      lastLineStartX = isBelowAvail ? startX : (width -1);
      lastLineEndX   = isBelowRightAvail ? width : (width -1);
      for(x= lastLineStartX; x< lastLineEndX; x++)
      {
        if (isCtuCrossedByVirtualBoundaries && isProcessDisabled(x, height - 1, numVerVirBndry, numHorVirBndry, verVirBndryPos, horVirBndryPos))
        {
          continue;
        }
        edgeType =  sgn(srcLine[x] - srcLineBelow[x+ 1]) + signUpLine[x];
        resLine[x] = ClipPel<int>( srcLine[x] + offset[edgeType], clpRng);

      }
    }
    break;
  case SAO_TYPE_EO_45:
    {
      offset += 2;
      int8_t *signUpLine = &m_signLineBuf1[1];

      startX = isLeftAvail ? 0 : 1;
      endX   = isRightAvail ? width : (width -1);

      //prepare 2nd line upper sign
      const Pel* srcLineBelow= srcLine+ srcStride;
      for (x=startX-1; x< endX; x++)
      {
        signUpLine[x] = (int8_t)sgn(srcLineBelow[x] - srcLine[x+1]);
      }


      //first line
      const Pel* srcLineAbove= srcLine- srcStride;
      firstLineStartX = isAboveAvail ? startX : (width -1 );
      firstLineEndX   = isAboveRightAvail ? width : (width-1);
      for(x= firstLineStartX; x< firstLineEndX; x++)
      {
        if (isCtuCrossedByVirtualBoundaries && isProcessDisabled(x, 0, numVerVirBndry, numHorVirBndry, verVirBndryPos, horVirBndryPos))
        {
          continue;
        }
        edgeType = sgn(srcLine[x] - srcLineAbove[x+1]) -signUpLine[x-1];
        resLine[x] = ClipPel<int>(srcLine[x] + offset[edgeType], clpRng);
      }
      srcLine += srcStride;
      resLine += resStride;

      //middle lines
      for (y= 1; y< height-1; y++)
      {
        srcLineBelow= srcLine+ srcStride;

        for(x= startX; x< endX; x++)
        {
          signDown =  (int8_t)sgn(srcLine[x] - srcLineBelow[x-1]);
          if (isCtuCrossedByVirtualBoundaries && isProcessDisabled(x, y, numVerVirBndry, numHorVirBndry, verVirBndryPos, horVirBndryPos))
          {
            signUpLine[x - 1] = -signDown;
            continue;
          }
          edgeType =  signDown + signUpLine[x];
          resLine[x] = ClipPel<int>(srcLine[x] + offset[edgeType], clpRng);
          signUpLine[x-1] = -signDown;
        }
        signUpLine[endX-1] = (int8_t)sgn(srcLineBelow[endX-1] - srcLine[endX]);
        srcLine  += srcStride;
        resLine += resStride;
      }

      //last line
      srcLineBelow= srcLine+ srcStride;
      lastLineStartX = isBelowLeftAvail ? 0 : 1;
      lastLineEndX   = isBelowAvail ? endX : 1;
      for(x= lastLineStartX; x< lastLineEndX; x++)
      {
        if (isCtuCrossedByVirtualBoundaries && isProcessDisabled(x, height - 1, numVerVirBndry, numHorVirBndry, verVirBndryPos, horVirBndryPos))
        {
          continue;
        }
        edgeType = sgn(srcLine[x] - srcLineBelow[x-1]) + signUpLine[x];
        resLine[x] = ClipPel<int>(srcLine[x] + offset[edgeType], clpRng);

      }
    }
    break;
  case SAO_TYPE_BO:
    {
      const int shiftBits = channelBitDepth - NUM_SAO_BO_CLASSES_LOG2;
      for (y=0; y< height; y++)
      {
        for (x=0; x< width; x++)
        {
          resLine[x] = ClipPel<int>(srcLine[x] + offset[srcLine[x] >> shiftBits], clpRng );
        }
        srcLine += srcStride;
        resLine += resStride;
      }
    }
    break;
  default:
    {
      THROW("Not a supported SAO types\n");
    }
  }
}

#if JVET_V0094_BILATERAL_FILTER
void SampleAdaptiveOffset::offsetBlockNoClip(const int channelBitDepth, const ClpRng& clpRng, int typeIdx, int* offset
                                             , const Pel* srcBlk, Pel* resBlk, int srcStride, int resStride,  int width, int height
                                             , bool isLeftAvail,  bool isRightAvail, bool isAboveAvail, bool isBelowAvail, bool isAboveLeftAvail, bool isAboveRightAvail, bool isBelowLeftAvail, bool isBelowRightAvail)
{
  int x,y, startX, startY, endX, endY, edgeType;
  int firstLineStartX, firstLineEndX, lastLineStartX, lastLineEndX;
  int8_t signLeft, signRight, signDown;
  
  const Pel* srcLine = srcBlk;
  Pel* resLine = resBlk;
  
  switch(typeIdx)
  {
    case SAO_TYPE_EO_0:
    {
      offset += 2;
      startX = isLeftAvail ? 0 : 1;
      endX   = isRightAvail ? width : (width -1);
      for (y=0; y< height; y++)
      {
        signLeft = (int8_t)sgn(srcLine[startX] - srcLine[startX-1]);
        for (x=startX; x< endX; x++)
        {
          signRight = (int8_t)sgn(srcLine[x] - srcLine[x+1]);
          edgeType =  signRight + signLeft;
          signLeft  = -signRight;
          
          resLine[x] = srcLine[x] + offset[edgeType];
        }
        srcLine  += srcStride;
        resLine += resStride;
      }
      
    }
      break;
    case SAO_TYPE_EO_90:
    {
      offset += 2;
      int8_t *signUpLine = &m_signLineBuf1[0];
      
      startY = isAboveAvail ? 0 : 1;
      endY   = isBelowAvail ? height : height-1;
      if (!isAboveAvail)
      {
        srcLine += srcStride;
        resLine += resStride;
      }
      
      const Pel* srcLineAbove= srcLine- srcStride;
      for (x=0; x< width; x++)
      {
        signUpLine[x] = (int8_t)sgn(srcLine[x] - srcLineAbove[x]);
      }
      
      const Pel* srcLineBelow;
      for (y=startY; y<endY; y++)
      {
        srcLineBelow= srcLine+ srcStride;
        
        for (x=0; x< width; x++)
        {
          signDown  = (int8_t)sgn(srcLine[x] - srcLineBelow[x]);
          edgeType = signDown + signUpLine[x];
          signUpLine[x]= -signDown;
          
          resLine[x] = srcLine[x] + offset[edgeType];
        }
        srcLine += srcStride;
        resLine += resStride;
      }
      
    }
      break;
    case SAO_TYPE_EO_135:
    {
      offset += 2;
      int8_t *signUpLine, *signDownLine, *signTmpLine;
      
      signUpLine  = &m_signLineBuf1[0];
      signDownLine= &m_signLineBuf2[0];
      
      startX = isLeftAvail ? 0 : 1 ;
      endX   = isRightAvail ? width : (width-1);
      
      //prepare 2nd line's upper sign
      const Pel* srcLineBelow= srcLine+ srcStride;
      for (x=startX; x< endX+1; x++)
      {
        signUpLine[x] = (int8_t)sgn(srcLineBelow[x] - srcLine[x- 1]);
      }
      
      //1st line
      const Pel* srcLineAbove= srcLine- srcStride;
      firstLineStartX = isAboveLeftAvail ? 0 : 1;
      firstLineEndX   = isAboveAvail? endX: 1;
      for(x= firstLineStartX; x< firstLineEndX; x++)
      {
        edgeType  =  sgn(srcLine[x] - srcLineAbove[x- 1]) - signUpLine[x+1];
        
        resLine[x] = srcLine[x] + offset[edgeType];
      }
      srcLine  += srcStride;
      resLine  += resStride;
      
      
      //middle lines
      for (y= 1; y< height-1; y++)
      {
        srcLineBelow= srcLine+ srcStride;
        
        for (x=startX; x<endX; x++)
        {
          signDown =  (int8_t)sgn(srcLine[x] - srcLineBelow[x+ 1]);
          edgeType =  signDown + signUpLine[x];
          resLine[x] = srcLine[x] + offset[edgeType];
          
          signDownLine[x+1] = -signDown;
        }
        signDownLine[startX] = (int8_t)sgn(srcLineBelow[startX] - srcLine[startX-1]);
        
        signTmpLine  = signUpLine;
        signUpLine   = signDownLine;
        signDownLine = signTmpLine;
        
        srcLine += srcStride;
        resLine += resStride;
      }
      
      //last line
      srcLineBelow= srcLine+ srcStride;
      lastLineStartX = isBelowAvail ? startX : (width -1);
      lastLineEndX   = isBelowRightAvail ? width : (width -1);
      for(x= lastLineStartX; x< lastLineEndX; x++)
      {
        edgeType =  sgn(srcLine[x] - srcLineBelow[x+ 1]) + signUpLine[x];
        resLine[x] = srcLine[x] + offset[edgeType];
        
      }
    }
      break;
    case SAO_TYPE_EO_45:
    {
      offset += 2;
      int8_t *signUpLine = &m_signLineBuf1[1];
      
      startX = isLeftAvail ? 0 : 1;
      endX   = isRightAvail ? width : (width -1);
      
      //prepare 2nd line upper sign
      const Pel* srcLineBelow= srcLine+ srcStride;
      for (x=startX-1; x< endX; x++)
      {
        signUpLine[x] = (int8_t)sgn(srcLineBelow[x] - srcLine[x+1]);
      }
      
      
      //first line
      const Pel* srcLineAbove= srcLine- srcStride;
      firstLineStartX = isAboveAvail ? startX : (width -1 );
      firstLineEndX   = isAboveRightAvail ? width : (width-1);
      for(x= firstLineStartX; x< firstLineEndX; x++)
      {
        edgeType = sgn(srcLine[x] - srcLineAbove[x+1]) -signUpLine[x-1];
        resLine[x] = srcLine[x] + offset[edgeType];
      }
      srcLine += srcStride;
      resLine += resStride;
      
      //middle lines
      for (y= 1; y< height-1; y++)
      {
        srcLineBelow= srcLine+ srcStride;
        
        for(x= startX; x< endX; x++)
        {
          signDown =  (int8_t)sgn(srcLine[x] - srcLineBelow[x-1]);
          edgeType =  signDown + signUpLine[x];
          resLine[x] = srcLine[x] + offset[edgeType];
          signUpLine[x-1] = -signDown;
        }
        signUpLine[endX-1] = (int8_t)sgn(srcLineBelow[endX-1] - srcLine[endX]);
        srcLine  += srcStride;
        resLine += resStride;
      }
      
      //last line
      srcLineBelow= srcLine+ srcStride;
      lastLineStartX = isBelowLeftAvail ? 0 : 1;
      lastLineEndX   = isBelowAvail ? endX : 1;
      for(x= lastLineStartX; x< lastLineEndX; x++)
      {
        edgeType = sgn(srcLine[x] - srcLineBelow[x-1]) + signUpLine[x];
        resLine[x] = srcLine[x] + offset[edgeType];
        
      }
    }
      break;
    case SAO_TYPE_BO:
    {
      const int shiftBits = channelBitDepth - NUM_SAO_BO_CLASSES_LOG2;
      for (y=0; y< height; y++)
      {
        for (x=0; x< width; x++)
        {
          resLine[x] = srcLine[x] + offset[srcLine[x] >> shiftBits];
        }
        srcLine += srcStride;
        resLine += resStride;
      }
    }
      break;
    default:
    {
      THROW("Not a supported SAO types\n");
    }
  }
}
#endif

void SampleAdaptiveOffset::offsetCTU( const UnitArea& area, const CPelUnitBuf& src, PelUnitBuf& res, SAOBlkParam& saoblkParam, CodingStructure& cs)
{
  const uint32_t numberOfComponents = getNumberValidComponents( area.chromaFormat );
  bool bAllOff=true;
  for( uint32_t compIdx = 0; compIdx < numberOfComponents; compIdx++)
  {
    if (saoblkParam[compIdx].modeIdc != SAO_MODE_OFF)
    {
      bAllOff=false;
    }
  }
  if (bAllOff)
  {
    return;
  }

  bool isLeftAvail, isRightAvail, isAboveAvail, isBelowAvail, isAboveLeftAvail, isAboveRightAvail, isBelowLeftAvail, isBelowRightAvail;

  //block boundary availability
  deriveLoopFilterBoundaryAvailibility(cs, area.Y(), isLeftAvail,isRightAvail,isAboveAvail,isBelowAvail,isAboveLeftAvail,isAboveRightAvail,isBelowLeftAvail,isBelowRightAvail);

  const size_t lineBufferSize = area.Y().width + 1;
  if (m_signLineBuf1.size() < lineBufferSize)
  {
    m_signLineBuf1.resize(lineBufferSize);
    m_signLineBuf2.resize(lineBufferSize);
  }

  int numHorVirBndry = 0, numVerVirBndry = 0;
  int horVirBndryPos[] = { -1,-1,-1 };
  int verVirBndryPos[] = { -1,-1,-1 };
  int horVirBndryPosComp[] = { -1,-1,-1 };
  int verVirBndryPosComp[] = { -1,-1,-1 };
  bool isCtuCrossedByVirtualBoundaries = isCrossedByVirtualBoundaries(area.Y().x, area.Y().y, area.Y().width, area.Y().height, numHorVirBndry, numVerVirBndry, horVirBndryPos, verVirBndryPos, cs.picHeader );
  for(int compIdx = 0; compIdx < numberOfComponents; compIdx++)
  {
    const ComponentID compID = ComponentID(compIdx);
    const CompArea& compArea = area.block(compID);
    SAOOffset& ctbOffset     = saoblkParam[compIdx];

    if(ctbOffset.modeIdc != SAO_MODE_OFF)
    {
      int  srcStride    = src.get(compID).stride;
      const Pel* srcBlk = src.get(compID).bufAt(compArea);
      int  resStride    = res.get(compID).stride;
      Pel* resBlk       = res.get(compID).bufAt(compArea);
      for (int i = 0; i < numHorVirBndry; i++)
      {
        horVirBndryPosComp[i] = (horVirBndryPos[i] >> ::getComponentScaleY(compID, area.chromaFormat)) - compArea.y;
      }
      for (int i = 0; i < numVerVirBndry; i++)
      {
        verVirBndryPosComp[i] = (verVirBndryPos[i] >> ::getComponentScaleX(compID, area.chromaFormat)) - compArea.x;
      }

      offsetBlock( cs.sps->getBitDepth(toChannelType(compID)),
                   cs.slice->clpRng(compID),
                   ctbOffset.typeIdc, ctbOffset.offset
                  , srcBlk, resBlk, srcStride, resStride, compArea.width, compArea.height
                  , isLeftAvail, isRightAvail
                  , isAboveAvail, isBelowAvail
                  , isAboveLeftAvail, isAboveRightAvail
                  , isBelowLeftAvail, isBelowRightAvail
                  , isCtuCrossedByVirtualBoundaries, horVirBndryPosComp, verVirBndryPosComp, numHorVirBndry, numVerVirBndry
                  );
    }
  } //compIdx
}

#if JVET_V0094_BILATERAL_FILTER
void SampleAdaptiveOffset::offsetCTUnoClip( const UnitArea& area, const CPelUnitBuf& src, PelUnitBuf& res, SAOBlkParam& saoblkParam, CodingStructure& cs)
{
  const uint32_t numberOfComponents = getNumberValidComponents( area.chromaFormat );
  bool bAllOff=true;
  for( uint32_t compIdx = 0; compIdx < numberOfComponents; compIdx++)
  {
    if (saoblkParam[compIdx].modeIdc != SAO_MODE_OFF)
    {
      bAllOff=false;
    }
  }
  if (bAllOff)
  {
    return;
  }
  
  bool isLeftAvail, isRightAvail, isAboveAvail, isBelowAvail, isAboveLeftAvail, isAboveRightAvail, isBelowLeftAvail, isBelowRightAvail;
  
  //block boundary availability
  deriveLoopFilterBoundaryAvailibility(cs, area.Y(), isLeftAvail,isRightAvail,isAboveAvail,isBelowAvail,isAboveLeftAvail,isAboveRightAvail,isBelowLeftAvail,isBelowRightAvail);
  
  const size_t lineBufferSize = area.Y().width + 1;
  if (m_signLineBuf1.size() < lineBufferSize)
  {
    m_signLineBuf1.resize(lineBufferSize);
    m_signLineBuf2.resize(lineBufferSize);
  }
  
#if !JVET_W0066_CCSAO
  int numHorVirBndry = 0, numVerVirBndry = 0;
  int horVirBndryPos[] = { -1,-1,-1 };
  int verVirBndryPos[] = { -1,-1,-1 };
  int horVirBndryPosComp[] = { -1,-1,-1 };
  int verVirBndryPosComp[] = { -1,-1,-1 };
  bool isCtuCrossedByVirtualBoundaries = isCrossedByVirtualBoundaries(area.Y().x, area.Y().y, area.Y().width, area.Y().height, numHorVirBndry, numVerVirBndry, horVirBndryPos, verVirBndryPos, cs.picHeader);
#endif
  for(int compIdx = 0; compIdx < numberOfComponents; compIdx++)
  {
    const ComponentID compID = ComponentID(compIdx);
    const CompArea& compArea = area.block(compID);
    SAOOffset& ctbOffset     = saoblkParam[compIdx];
    
    if(ctbOffset.modeIdc != SAO_MODE_OFF)
    {
      int  srcStride    = src.get(compID).stride;
      const Pel* srcBlk = src.get(compID).bufAt(compArea);
      int  resStride    = res.get(compID).stride;
      Pel* resBlk       = res.get(compID).bufAt(compArea);
#if !JVET_W0066_CCSAO
      for (int i = 0; i < numHorVirBndry; i++)
      {
        horVirBndryPosComp[i] = (horVirBndryPos[i] >> ::getComponentScaleY(compID, area.chromaFormat)) - compArea.y;
      }
      for (int i = 0; i < numVerVirBndry; i++)
      {
        verVirBndryPosComp[i] = (verVirBndryPos[i] >> ::getComponentScaleX(compID, area.chromaFormat)) - compArea.x;
      }
#endif
#if JVET_W0066_CCSAO    
      // Do not clip the final output for both luma and chroma. Clipping is done jontly for SAO/BIF/CCSAO.
      
	    offsetBlockNoClip(cs.sps->getBitDepth(toChannelType(compID)),
		    cs.slice->clpRng(compID),
		    ctbOffset.typeIdc, ctbOffset.offset
		    , srcBlk, resBlk, srcStride, resStride, compArea.width, compArea.height
		    , isLeftAvail, isRightAvail
		    , isAboveAvail, isBelowAvail
		    , isAboveLeftAvail, isAboveRightAvail
		    , isBelowLeftAvail, isBelowRightAvail
		  //                 , isCtuCrossedByVirtualBoundaries, horVirBndryPosComp, verVirBndryPosComp, numHorVirBndry, numVerVirBndry
	    );    
#else
      if(compID == COMPONENT_Y)
      {
        // If it is luma we should not clip, since we will clip
        // after BIF has been added.

        offsetBlockNoClip( cs.sps->getBitDepth(toChannelType(compID)),
                    cs.slice->clpRng(compID),
                    ctbOffset.typeIdc, ctbOffset.offset
                    , srcBlk, resBlk, srcStride, resStride, compArea.width, compArea.height
                    , isLeftAvail, isRightAvail
                    , isAboveAvail, isBelowAvail
                    , isAboveLeftAvail, isAboveRightAvail
                    , isBelowLeftAvail, isBelowRightAvail
   //                 , isCtuCrossedByVirtualBoundaries, horVirBndryPosComp, verVirBndryPosComp, numHorVirBndry, numVerVirBndry
                    );
        }
      else
      {
        // If it is chroma we should clip as normal, since
        // chroma is not bilaterally filtered.
        offsetBlock( cs.sps->getBitDepth(toChannelType(compID)),
                    cs.slice->clpRng(compID),
                    ctbOffset.typeIdc, ctbOffset.offset
                    , srcBlk, resBlk, srcStride, resStride, compArea.width, compArea.height
                    , isLeftAvail, isRightAvail
                    , isAboveAvail, isBelowAvail
                    , isAboveLeftAvail, isAboveRightAvail
                    , isBelowLeftAvail, isBelowRightAvail
                    , isCtuCrossedByVirtualBoundaries, horVirBndryPosComp, verVirBndryPosComp, numHorVirBndry, numVerVirBndry
                    );

      }
      
#endif
    }
  } //compIdx
}
#endif

void SampleAdaptiveOffset::SAOProcess( CodingStructure& cs, SAOBlkParam* saoBlkParams
                                      )
{
  CHECK(!saoBlkParams, "No parameters present");

#if JVET_V0094_BILATERAL_FILTER
  // In code without BIF, SAOProcess would not be run if 'SAO=0'.
  // However, in the BIF-enabled code, we still might go here if 'SAO=0' and 'BIF=1'.
  // Hence we must check getSAOEnabledFlag() for some of the function calls.
  if( cs.sps->getSAOEnabledFlag() )
  {
    xReconstructBlkSAOParams(cs, saoBlkParams);
  }
#else
  xReconstructBlkSAOParams(cs, saoBlkParams);
#endif

  const uint32_t numberOfComponents = getNumberValidComponents(cs.area.chromaFormat);
  bool bAllDisabled = true;
#if JVET_V0094_BILATERAL_FILTER
  // If 'SAO=0' we would not normally get here. However, now we might get
  // here if 'SAO=0' and 'BIF=1'. Hence we should only run this if
  // getSAOEnabledFlag() is true. Note that if getSAOEnabledFlag() is false,
  // we will not run the code and bAllDisabled will stay true, which will
  // give the correct behavior.
  if( cs.sps->getSAOEnabledFlag() )
  {
#endif
  for (uint32_t compIdx = 0; compIdx < numberOfComponents; compIdx++)
  {
    if (m_picSAOEnabled[compIdx])
    {
      bAllDisabled = false;
    }
  }
#if JVET_V0094_BILATERAL_FILTER
  }
#endif
  if (bAllDisabled)
  {
#if JVET_V0094_BILATERAL_FILTER
    // Even if we are not doing SAO we might still need to do BIF
    // so we cannot return even if SAO is never used.
    if(!cs.pps->getUseBIF())
    {
      // However, if we are not doing BIF it is safe to return.
      return;
    }
#else
    return;
#endif
  }

  const PreCalcValues& pcv = *cs.pcv;
  PelUnitBuf rec = cs.getRecoBuf();
  m_tempBuf.copyFrom( rec );

  int ctuRsAddr = 0;
  for( uint32_t yPos = 0; yPos < pcv.lumaHeight; yPos += pcv.maxCUHeight )
  {
    for( uint32_t xPos = 0; xPos < pcv.lumaWidth; xPos += pcv.maxCUWidth )
    {
      const uint32_t width  = (xPos + pcv.maxCUWidth  > pcv.lumaWidth)  ? (pcv.lumaWidth - xPos)  : pcv.maxCUWidth;
      const uint32_t height = (yPos + pcv.maxCUHeight > pcv.lumaHeight) ? (pcv.lumaHeight - yPos) : pcv.maxCUHeight;
      const UnitArea area( cs.area.chromaFormat, Area(xPos , yPos, width, height) );

#if JVET_W0066_CCSAO
      // Always do non-clipped version for SAO/BIF, the clipping is done jointly after CCSAO is also applied
      if (!bAllDisabled)
        offsetCTUnoClip(area, m_tempBuf, rec, cs.picture->getSAO()[ctuRsAddr], cs);

      BifParams& bifParams = cs.picture->getBifParam();

      // And now we traverse the CTU to do BIF
      for (auto& currCU : cs.traverseCUs(CS::getArea(cs, area, CH_L), CH_L))
      {
        for (auto& currTU : CU::traverseTUs(currCU))
        {
          bool isInter = (currCU.predMode == MODE_INTER) ? true : false;
          if (bifParams.ctuOn[ctuRsAddr] && ((TU::getCbf(currTU, COMPONENT_Y) || isInter == false) && (currTU.cu->qp > 17)) && (128 > std::max(currTU.lumaSize().width, currTU.lumaSize().height)) && ((isInter == false) || (32 > std::min(currTU.lumaSize().width, currTU.lumaSize().height))))
          {
            m_bilateralFilter.bilateralFilterDiamond5x5NoClip(m_tempBuf, rec, currTU.cu->qp, cs.slice->clpRng(COMPONENT_Y), currTU);
          }
        }
      }
#else
#if JVET_V0094_BILATERAL_FILTER
      if(cs.pps->getUseBIF())
      {
        // We are using BIF, so we run SAO without clipping
        // However, if 'SAO=0', bAllDisabled=true and we should not run offsetCTUnoClip.
        if( !bAllDisabled )
          offsetCTUnoClip( area, m_tempBuf, rec, cs.picture->getSAO()[ctuRsAddr], cs);
        
        // We don't need to clip if SAO was not performed on luma.
        SAOBlkParam mySAOblkParam = cs.picture->getSAO()[ctuRsAddr];
        SAOOffset& myCtbOffset     = mySAOblkParam[0];
        BifParams& bifParams = cs.picture->getBifParam();
        
        bool clipLumaIfNoBilat = false;
        if(!bAllDisabled && myCtbOffset.modeIdc != SAO_MODE_OFF)
          clipLumaIfNoBilat = true;
        
        // And now we traverse the CTU to do BIF
        for (auto &currCU : cs.traverseCUs(CS::getArea(cs, area, CH_L), CH_L))
        {
          for (auto &currTU : CU::traverseTUs(currCU))
          {
            
            bool isInter = (currCU.predMode == MODE_INTER) ? true : false;
            if ( bifParams.ctuOn[ctuRsAddr] && ((TU::getCbf(currTU, COMPONENT_Y) || isInter == false) && (currTU.cu->qp > 17)) && (128 > std::max(currTU.lumaSize().width, currTU.lumaSize().height)) && ((isInter == false) || (32 > std::min(currTU.lumaSize().width, currTU.lumaSize().height))))
            {
              m_bilateralFilter.bilateralFilterDiamond5x5(m_tempBuf, rec, currTU.cu->qp, cs.slice->clpRng(COMPONENT_Y), currTU);
            }
            else
            {
              // We don't need to clip if SAO was not performed on luma.
              if(clipLumaIfNoBilat)
                m_bilateralFilter.clipNotBilaterallyFilteredBlocks(m_tempBuf, rec, cs.slice->clpRng(COMPONENT_Y), currTU);
            }
          }
        }
      }
      else
      {
        // BIF is not used, use old SAO code
        offsetCTU( area, m_tempBuf, rec, cs.picture->getSAO()[ctuRsAddr], cs);
      }
#else
      offsetCTU( area, m_tempBuf, rec, cs.picture->getSAO()[ctuRsAddr], cs);
#endif
#endif
      ctuRsAddr++;
    }
  }

  DTRACE_UPDATE(g_trace_ctx, (std::make_pair("poc", cs.slice->getPOC())));
  DTRACE_PIC_COMP(D_REC_CB_LUMA_SAO, cs, cs.getRecoBuf(), COMPONENT_Y);
  DTRACE_PIC_COMP(D_REC_CB_CHROMA_SAO, cs, cs.getRecoBuf(), COMPONENT_Cb);
  DTRACE_PIC_COMP(D_REC_CB_CHROMA_SAO, cs, cs.getRecoBuf(), COMPONENT_Cr);

  DTRACE    ( g_trace_ctx, D_CRC, "SAO" );
  DTRACE_CRC( g_trace_ctx, D_CRC, cs, cs.getRecoBuf() );

}

#if JVET_W0066_CCSAO
void SampleAdaptiveOffset::CCSAOProcess(CodingStructure& cs)
{
  const uint32_t numberOfComponents = getNumberValidComponents(cs.area.chromaFormat);
  bool bAllDisabled = true;
  for (uint32_t compIdx = 0; compIdx < numberOfComponents; compIdx++)
  {
    if (m_ccSaoComParam.enabled[compIdx])
    {
      bAllDisabled = false;
    }
  }
  if (bAllDisabled)
  {
    return;
  }
  
  const PreCalcValues& pcv = *cs.pcv;
  PelUnitBuf dstYuv = cs.getRecoBuf();
  PelUnitBuf srcYuv = m_ccSaoBuf.getBuf( cs.area );
  srcYuv.extendBorderPel( MAX_CCSAO_FILTER_LENGTH >> 1 );

  applyCcSao(cs, pcv, srcYuv, dstYuv);
}

void SampleAdaptiveOffset::applyCcSao(CodingStructure &cs, const PreCalcValues& pcv, const CPelUnitBuf& srcYuv, PelUnitBuf& dstYuv)
{
  int ctuRsAddr = 0;
  for (uint32_t yPos = 0; yPos < pcv.lumaHeight; yPos += pcv.maxCUHeight)
  {
    for (uint32_t xPos = 0; xPos < pcv.lumaWidth; xPos += pcv.maxCUWidth)
    {
      const uint32_t width  = (xPos + pcv.maxCUWidth  > pcv.lumaWidth ) ? (pcv.lumaWidth  - xPos) : pcv.maxCUWidth;
      const uint32_t height = (yPos + pcv.maxCUHeight > pcv.lumaHeight) ? (pcv.lumaHeight - yPos) : pcv.maxCUHeight;
      const UnitArea area(cs.area.chromaFormat, Area(xPos, yPos, width, height));

      offsetCTUCcSaoNoClip(cs, area, srcYuv, dstYuv, ctuRsAddr);
      ctuRsAddr++;
    }
  }
}

void SampleAdaptiveOffset::jointClipSaoBifCcSao(CodingStructure& cs)
{
  if (!cs.sps->getSAOEnabledFlag() && !cs.pps->getUseBIF() && !cs.sps->getCCSAOEnabledFlag())
    return;

  const PreCalcValues& pcv = *cs.pcv;
  PelUnitBuf dstYuv = cs.getRecoBuf();

  // Iterate all CTUs and check if any of the filters is on for a given component
  int ctuRsAddr = 0;
  for (uint32_t yPos = 0; yPos < pcv.lumaHeight; yPos += pcv.maxCUHeight)
  {
    for (uint32_t xPos = 0; xPos < pcv.lumaWidth; xPos += pcv.maxCUWidth)
    {
      const uint32_t width = (xPos + pcv.maxCUWidth > pcv.lumaWidth) ? (pcv.lumaWidth - xPos) : pcv.maxCUWidth;
      const uint32_t height = (yPos + pcv.maxCUHeight > pcv.lumaHeight) ? (pcv.lumaHeight - yPos) : pcv.maxCUHeight;
      const UnitArea area(cs.area.chromaFormat, Area(xPos, yPos, width, height));
      const uint32_t numberOfComponents = getNumberValidComponents(area.chromaFormat);
      SAOBlkParam mySAOblkParam = cs.picture->getSAO()[ctuRsAddr];

      for (int compIdx = 0; compIdx < numberOfComponents; compIdx++)
      {
        const int setIdc = m_ccSaoControl[compIdx][ctuRsAddr];
        SAOOffset& myCtbOffset = mySAOblkParam[compIdx];
        if ((m_ccSaoComParam.enabled[compIdx] && setIdc != 0) || (myCtbOffset.modeIdc != SAO_MODE_OFF))
        {
          // We definitely need to clip if either SAO or CCSAO is on for the given component of the CTU                  
          clipCTU(cs, dstYuv, area, ComponentID(compIdx));
        }
        else
        {
          // When BIF is on, the luma component might need to be clipped
          if (cs.pps->getUseBIF())
          {
            if (compIdx == COMPONENT_Y)
            {
              BifParams& bifParams = cs.picture->getBifParam();

              // And now we traverse the CTU to do clipping
              for (auto& currCU : cs.traverseCUs(CS::getArea(cs, area, CH_L), CH_L))
              {
                for (auto& currTU : CU::traverseTUs(currCU))
                {
                  bool isInter = (currCU.predMode == MODE_INTER) ? true : false;
                  if (bifParams.ctuOn[ctuRsAddr] && ((TU::getCbf(currTU, COMPONENT_Y) || isInter == false) && (currTU.cu->qp > 17)) && (128 > std::max(currTU.lumaSize().width, currTU.lumaSize().height)) && ((isInter == false) || (32 > std::min(currTU.lumaSize().width, currTU.lumaSize().height))))
                  {
                    m_bilateralFilter.clipNotBilaterallyFilteredBlocks(m_tempBuf, dstYuv, cs.slice->clpRng(COMPONENT_Y), currTU);
                  }
                }
              }
            }
          }
        }
      }
      ctuRsAddr++;
    }
  }
}

void SampleAdaptiveOffset::clipCTU(CodingStructure& cs, PelUnitBuf& dstYuv, const UnitArea& area, const ComponentID compID)
{
  const CompArea &compArea = area.block(compID);
  const uint32_t height = compArea.height;
  const uint32_t width = compArea.width;
  Pel *dst = dstYuv.get(compID).bufAt(area.block(compID));
  int dstStride = dstYuv.get(compID).stride;
  
  for (uint32_t y = 0; y < height; y++)
  {
    for (uint32_t x = 0; x < width; x++)
    {
      // new result = old result (which is SAO-treated already) + clipping
      dst[x] = ClipPel<int>(dst[x], cs.slice->clpRng(compID));
    }
    dst += dstStride;
  }
}

void SampleAdaptiveOffset::offsetCTUCcSaoNoClip(CodingStructure& cs, const UnitArea& area, const CPelUnitBuf& srcYuv, PelUnitBuf& dstYuv, const int ctuRsAddr)
{
  const uint32_t numberOfComponents = getNumberValidComponents(area.chromaFormat);
  bool bAllOff = true;
  for (uint32_t compIdx = 0; compIdx < numberOfComponents; compIdx++)
  {
    if (m_ccSaoComParam.enabled[compIdx])
    {
      bAllOff = false;
    }
  }
  if (bAllOff)
  {
    return;
  }

  bool isLeftAvail, isRightAvail, isAboveAvail, isBelowAvail, isAboveLeftAvail, isAboveRightAvail, isBelowLeftAvail, isBelowRightAvail;
  deriveLoopFilterBoundaryAvailibility(cs, area.Y(), isLeftAvail, isRightAvail, isAboveAvail, isBelowAvail, isAboveLeftAvail, isAboveRightAvail, isBelowLeftAvail, isBelowRightAvail);

  for (int compIdx = 0; compIdx < numberOfComponents; compIdx++)
  {
    if (m_ccSaoComParam.enabled[compIdx])
    {
      const int setIdc = m_ccSaoControl[compIdx][ctuRsAddr];

      if (setIdc != 0)
      {
        const ComponentID compID     = ComponentID(compIdx);
        const CompArea   &compArea   = area.block(compID);
        const int         srcStrideY = srcYuv.get(COMPONENT_Y ).stride;
        const int         srcStrideU = srcYuv.get(COMPONENT_Cb).stride;
        const int         srcStrideV = srcYuv.get(COMPONENT_Cr).stride;
        const int         dstStride  = dstYuv.get(compID      ).stride;
        const Pel        *srcBlkY    = srcYuv.get(COMPONENT_Y ).bufAt(area.block(COMPONENT_Y ));
        const Pel        *srcBlkU    = srcYuv.get(COMPONENT_Cb).bufAt(area.block(COMPONENT_Cb));
        const Pel        *srcBlkV    = srcYuv.get(COMPONENT_Cr).bufAt(area.block(COMPONENT_Cr));
              Pel        *dstBlk     = dstYuv.get(compID      ).bufAt(compArea);
              
        const uint16_t    candPosY   = m_ccSaoComParam.candPos[compIdx][setIdc - 1][COMPONENT_Y ];
        const uint16_t    bandNumY   = m_ccSaoComParam.bandNum[compIdx][setIdc - 1][COMPONENT_Y ];
        const uint16_t    bandNumU   = m_ccSaoComParam.bandNum[compIdx][setIdc - 1][COMPONENT_Cb];
        const uint16_t    bandNumV   = m_ccSaoComParam.bandNum[compIdx][setIdc - 1][COMPONENT_Cr];
        const short      *offset     = m_ccSaoComParam.offset [compIdx][setIdc - 1];

        offsetBlockCcSaoNoClip(compID, cs.sps->getBitDepth(toChannelType(compID)), cs.slice->clpRng(compID)
                             , candPosY, bandNumY, bandNumU, bandNumV
                             , offset
                             , srcBlkY, srcBlkU, srcBlkV, dstBlk
                             , srcStrideY, srcStrideU, srcStrideV, dstStride
                             , compArea.width, compArea.height
                             , isLeftAvail, isRightAvail
                             , isAboveAvail, isBelowAvail
                             , isAboveLeftAvail, isAboveRightAvail
                             , isBelowLeftAvail, isBelowRightAvail
                              );
      }
    }
  }
}

void SampleAdaptiveOffset::offsetCTUCcSao(CodingStructure& cs, const UnitArea& area, const CPelUnitBuf& srcYuv, PelUnitBuf& dstYuv, const int ctuRsAddr)
{
  const uint32_t numberOfComponents = getNumberValidComponents( area.chromaFormat );
  bool bAllOff = true;
  for( uint32_t compIdx = 0; compIdx < numberOfComponents; compIdx++)
  {
    if (m_ccSaoComParam.enabled[compIdx])
    {
      bAllOff = false;
    }
  }
  if (bAllOff)
  {
    return;
  }

  bool isLeftAvail, isRightAvail, isAboveAvail, isBelowAvail, isAboveLeftAvail, isAboveRightAvail, isBelowLeftAvail, isBelowRightAvail;
  deriveLoopFilterBoundaryAvailibility(cs, area.Y(), isLeftAvail,isRightAvail,isAboveAvail,isBelowAvail,isAboveLeftAvail,isAboveRightAvail,isBelowLeftAvail,isBelowRightAvail);

  for(int compIdx = 0; compIdx < numberOfComponents; compIdx++)
  {
    if(m_ccSaoComParam.enabled[compIdx])
    {
      const int setIdc = m_ccSaoControl[compIdx][ctuRsAddr];

      if (setIdc != 0)
      {
        const ComponentID compID     = ComponentID(compIdx);
        const CompArea   &compArea   = area.block(compID);
        const int         srcStrideY = srcYuv.get(COMPONENT_Y ).stride;
        const int         srcStrideU = srcYuv.get(COMPONENT_Cb).stride;
        const int         srcStrideV = srcYuv.get(COMPONENT_Cr).stride;
        const int         dstStride  = dstYuv.get(compID      ).stride;
        const Pel        *srcBlkY    = srcYuv.get(COMPONENT_Y ).bufAt(area.block(COMPONENT_Y ));
        const Pel        *srcBlkU    = srcYuv.get(COMPONENT_Cb).bufAt(area.block(COMPONENT_Cb));
        const Pel        *srcBlkV    = srcYuv.get(COMPONENT_Cr).bufAt(area.block(COMPONENT_Cr));
              Pel        *dstBlk     = dstYuv.get(compID      ).bufAt(compArea);
              
        const uint16_t    candPosY   = m_ccSaoComParam.candPos[compIdx][setIdc - 1][COMPONENT_Y ];
        const uint16_t    bandNumY   = m_ccSaoComParam.bandNum[compIdx][setIdc - 1][COMPONENT_Y ];
        const uint16_t    bandNumU   = m_ccSaoComParam.bandNum[compIdx][setIdc - 1][COMPONENT_Cb];
        const uint16_t    bandNumV   = m_ccSaoComParam.bandNum[compIdx][setIdc - 1][COMPONENT_Cr];
        const short      *offset     = m_ccSaoComParam.offset [compIdx][setIdc - 1];

        offsetBlockCcSao( compID, cs.sps->getBitDepth(toChannelType(compID)), cs.slice->clpRng(compID)
                        , candPosY, bandNumY, bandNumU, bandNumV
                        , offset
                        , srcBlkY, srcBlkU, srcBlkV, dstBlk
                        , srcStrideY, srcStrideU, srcStrideV, dstStride
                        , compArea.width, compArea.height
                        , isLeftAvail, isRightAvail
                        , isAboveAvail, isBelowAvail
                        , isAboveLeftAvail, isAboveRightAvail
                        , isBelowLeftAvail, isBelowRightAvail
                        );
      }
    }
  }
}

void SampleAdaptiveOffset::offsetBlockCcSaoNoClip(const ComponentID compID, const int bitDepth, const ClpRng& clpRng
                                                , const uint16_t candPosY
                                                , const uint16_t bandNumY, const uint16_t bandNumU, const uint16_t bandNumV
                                                , const short* offset
                                                , const Pel* srcY, const Pel* srcU, const Pel* srcV, Pel* dst
                                                , const int srcStrideY, const int srcStrideU, const int srcStrideV, const int dstStride
                                                , const int width, const int height
                                                , bool isLeftAvail, bool isRightAvail, bool isAboveAvail, bool isBelowAvail, bool isAboveLeftAvail, bool isAboveRightAvail, bool isBelowLeftAvail, bool isBelowRightAvail
                                                 )
{
  const int candPosYX = g_ccSaoCandPosX[COMPONENT_Y][candPosY];
  const int candPosYY = g_ccSaoCandPosY[COMPONENT_Y][candPosY];

  switch (compID)
  {
  case COMPONENT_Y:
    {
      for (int y = 0; y < height; y++)
      {
        for (int x = 0; x < width; x++)
        {
          const Pel* colY = srcY +  x + srcStrideY * candPosYY + candPosYX;
          const Pel* colU = srcU + (x >> 1);
          const Pel* colV = srcV + (x >> 1);

          const int bandY    = (*colY * bandNumY) >> bitDepth;
          const int bandU    = (*colU * bandNumU) >> bitDepth;
          const int bandV    = (*colV * bandNumV) >> bitDepth;
          const int bandIdx  = bandY * bandNumU * bandNumV
                             + bandU * bandNumV
                             + bandV;
          const int classIdx = bandIdx;

          //dst[x] = ClipPel<int>(dst[x] + offset[classIdx], clpRng);
          dst[x] = dst[x] + offset[classIdx];
        }

        srcY += srcStrideY;
        srcU += srcStrideU * (y & 0x1);
        srcV += srcStrideV * (y & 0x1);
        dst  += dstStride;
      }
    }
  break;
  case COMPONENT_Cb:
  case COMPONENT_Cr:
    {
      for (int y = 0; y < height; y++)
      {
        for (int x = 0; x < width; x++)
        {
          const Pel* colY = srcY + (x << 1) + srcStrideY * candPosYY + candPosYX;
          const Pel* colU = srcU +  x;
          const Pel* colV = srcV +  x;

          const int bandY    = (*colY * bandNumY) >> bitDepth;
          const int bandU    = (*colU * bandNumU) >> bitDepth;
          const int bandV    = (*colV * bandNumV) >> bitDepth;
          const int bandIdx  = bandY * bandNumU * bandNumV
                             + bandU * bandNumV
                             + bandV;
          const int classIdx = bandIdx;

          //dst[x] = ClipPel<int>(dst[x] + offset[classIdx], clpRng);
          dst[x] = dst[x] + offset[classIdx];
        }
        srcY += srcStrideY << 1;
        srcU += srcStrideU;
        srcV += srcStrideV;
        dst  += dstStride;
      }
    }
  break;
  default:
    {
      THROW("Not a supported CCSAO compID\n");
    }
  }
}

void SampleAdaptiveOffset::offsetBlockCcSao(const ComponentID compID, const int bitDepth, const ClpRng& clpRng
                                          , const uint16_t candPosY
                                          , const uint16_t bandNumY, const uint16_t bandNumU, const uint16_t bandNumV
                                          , const short* offset
                                          , const Pel* srcY, const Pel* srcU, const Pel* srcV, Pel* dst
                                          , const int srcStrideY, const int srcStrideU, const int srcStrideV, const int dstStride
                                          , const int width, const int height
                                          , bool isLeftAvail, bool isRightAvail, bool isAboveAvail, bool isBelowAvail, bool isAboveLeftAvail, bool isAboveRightAvail, bool isBelowLeftAvail, bool isBelowRightAvail
                                           )
{
  const int candPosYX = g_ccSaoCandPosX[COMPONENT_Y][candPosY];
  const int candPosYY = g_ccSaoCandPosY[COMPONENT_Y][candPosY];

  switch(compID)
  {
  case COMPONENT_Y:
    {
      for (int y = 0; y < height; y++)
      {
        for (int x = 0; x < width; x++)
        {
          const Pel *colY = srcY +  x  + srcStrideY * candPosYY + candPosYX;
          const Pel *colU = srcU + (x >> 1);
          const Pel *colV = srcV + (x >> 1);

          const int bandY    = (*colY * bandNumY) >> bitDepth;
          const int bandU    = (*colU * bandNumU) >> bitDepth;
          const int bandV    = (*colV * bandNumV) >> bitDepth;
          const int bandIdx  = bandY * bandNumU * bandNumV
                             + bandU * bandNumV
                             + bandV;
          const int classIdx = bandIdx;

          dst[x] = ClipPel<int>(dst[x] + offset[classIdx], clpRng);
        }

        srcY += srcStrideY;
        srcU += srcStrideU * (y & 0x1);
        srcV += srcStrideV * (y & 0x1);
        dst  += dstStride;
      }
    }
    break;
  case COMPONENT_Cb:
  case COMPONENT_Cr:
    {
      for (int y = 0; y < height; y++)
      {
        for (int x = 0; x < width; x++)
        {
          const Pel *colY = srcY + (x << 1) + srcStrideY * candPosYY + candPosYX;
          const Pel *colU = srcU + x;
          const Pel *colV = srcV + x;

          const int bandY    = (*colY * bandNumY) >> bitDepth;
          const int bandU    = (*colU * bandNumU) >> bitDepth;
          const int bandV    = (*colV * bandNumV) >> bitDepth;
          const int bandIdx  = bandY * bandNumU * bandNumV
                             + bandU * bandNumV
                             + bandV;
          const int classIdx = bandIdx;

          dst[x] = ClipPel<int>(dst[x] + offset[classIdx], clpRng);
        }

        srcY += srcStrideY << 1;
        srcU += srcStrideU;
        srcV += srcStrideV;
        dst  += dstStride;
      }
    }
    break;
  default:
    {
      THROW("Not a supported CCSAO compID\n");
    }
  }
}
#endif

void SampleAdaptiveOffset::deriveLoopFilterBoundaryAvailibility(CodingStructure& cs, const Position &pos,
  bool& isLeftAvail,
  bool& isRightAvail,
  bool& isAboveAvail,
  bool& isBelowAvail,
  bool& isAboveLeftAvail,
  bool& isAboveRightAvail,
  bool& isBelowLeftAvail,
  bool& isBelowRightAvail
  ) const
{
  const int width = cs.pcv->maxCUWidth;
  const int height = cs.pcv->maxCUHeight;
  const CodingUnit* cuCurr = cs.getCU(pos, CH_L);
  const CodingUnit* cuLeft = cs.getCU(pos.offset(-width, 0), CH_L);
  const CodingUnit* cuRight = cs.getCU(pos.offset(width, 0), CH_L);
  const CodingUnit* cuAbove = cs.getCU(pos.offset(0, -height), CH_L);
  const CodingUnit* cuBelow = cs.getCU(pos.offset(0, height), CH_L);
  const CodingUnit* cuAboveLeft = cs.getCU(pos.offset(-width, -height), CH_L);
  const CodingUnit* cuAboveRight = cs.getCU(pos.offset(width, -height), CH_L);
  const CodingUnit* cuBelowLeft = cs.getCU(pos.offset(-width, height), CH_L);
  const CodingUnit* cuBelowRight = cs.getCU(pos.offset(width, height), CH_L);

  // check cross slice flags
  const bool isLoopFilterAcrossSlicePPS = cs.pps->getLoopFilterAcrossSlicesEnabledFlag();
  if (!isLoopFilterAcrossSlicePPS)
  {
    isLeftAvail       = (cuLeft == NULL)       ? false : CU::isSameSlice(*cuCurr, *cuLeft);
    isAboveAvail      = (cuAbove == NULL)      ? false : CU::isSameSlice(*cuCurr, *cuAbove);
    isRightAvail      = (cuRight == NULL)      ? false : CU::isSameSlice(*cuCurr, *cuRight);
    isBelowAvail      = (cuBelow == NULL)      ? false : CU::isSameSlice(*cuCurr, *cuBelow);
    isAboveLeftAvail  = (cuAboveLeft == NULL)  ? false : CU::isSameSlice(*cuCurr, *cuAboveLeft);
    isAboveRightAvail = (cuAboveRight == NULL) ? false : CU::isSameSlice(*cuCurr, *cuAboveRight);
    isBelowLeftAvail  = (cuBelowLeft == NULL)  ? false : CU::isSameSlice(*cuCurr, *cuBelowLeft);
    isBelowRightAvail = (cuBelowRight == NULL) ? false : CU::isSameSlice(*cuCurr, *cuBelowRight);
  }
  else
  {
    isLeftAvail       = (cuLeft != NULL);
    isAboveAvail      = (cuAbove != NULL);
    isRightAvail      = (cuRight != NULL);
    isBelowAvail      = (cuBelow != NULL);
    isAboveLeftAvail  = (cuAboveLeft != NULL);
    isAboveRightAvail = (cuAboveRight != NULL);
    isBelowLeftAvail  = (cuBelowLeft != NULL);
    isBelowRightAvail = (cuBelowRight != NULL);
  }

  // check cross tile flags
  const bool isLoopFilterAcrossTilePPS = cs.pps->getLoopFilterAcrossTilesEnabledFlag();
  if (!isLoopFilterAcrossTilePPS)
  {
    isLeftAvail       = (!isLeftAvail)       ? false : CU::isSameTile(*cuCurr, *cuLeft);
    isAboveAvail      = (!isAboveAvail)      ? false : CU::isSameTile(*cuCurr, *cuAbove);
    isRightAvail      = (!isRightAvail)      ? false : CU::isSameTile(*cuCurr, *cuRight);
    isBelowAvail      = (!isBelowAvail)      ? false : CU::isSameTile(*cuCurr, *cuBelow);
    isAboveLeftAvail  = (!isAboveLeftAvail)  ? false : CU::isSameTile(*cuCurr, *cuAboveLeft);
    isAboveRightAvail = (!isAboveRightAvail) ? false : CU::isSameTile(*cuCurr, *cuAboveRight);
    isBelowLeftAvail  = (!isBelowLeftAvail)  ? false : CU::isSameTile(*cuCurr, *cuBelowLeft);
    isBelowRightAvail = (!isBelowRightAvail) ? false : CU::isSameTile(*cuCurr, *cuBelowRight);
  }

  // check cross subpic flags
  const SubPic& curSubPic = cs.pps->getSubPicFromCU(*cuCurr);
  if (!curSubPic.getloopFilterAcrossEnabledFlag())
  {
    isLeftAvail       = (!isLeftAvail)       ? false : CU::isSameSubPic(*cuCurr, *cuLeft);
    isAboveAvail      = (!isAboveAvail)      ? false : CU::isSameSubPic(*cuCurr, *cuAbove);
    isRightAvail      = (!isRightAvail)      ? false : CU::isSameSubPic(*cuCurr, *cuRight);
    isBelowAvail      = (!isBelowAvail)      ? false : CU::isSameSubPic(*cuCurr, *cuBelow);
    isAboveLeftAvail  = (!isAboveLeftAvail)  ? false : CU::isSameSubPic(*cuCurr, *cuAboveLeft);
    isAboveRightAvail = (!isAboveRightAvail) ? false : CU::isSameSubPic(*cuCurr, *cuAboveRight);
    isBelowLeftAvail  = (!isBelowLeftAvail)  ? false : CU::isSameSubPic(*cuCurr, *cuBelowLeft);
    isBelowRightAvail = (!isBelowRightAvail) ? false : CU::isSameSubPic(*cuCurr, *cuBelowRight);
  }
}

bool SampleAdaptiveOffset::isCrossedByVirtualBoundaries(const int xPos, const int yPos, const int width, const int height, int& numHorVirBndry, int& numVerVirBndry, int horVirBndryPos[], int verVirBndryPos[], const PicHeader* picHeader )
{
  numHorVirBndry = 0; numVerVirBndry = 0;
  if( picHeader->getVirtualBoundariesPresentFlag() )
  {
    for (int i = 0; i < picHeader->getNumHorVirtualBoundaries(); i++)
    {
     if (yPos <= picHeader->getVirtualBoundariesPosY(i) && picHeader->getVirtualBoundariesPosY(i) <= yPos + height)
      {
        horVirBndryPos[numHorVirBndry++] = picHeader->getVirtualBoundariesPosY(i);
      }
    }
    for (int i = 0; i < picHeader->getNumVerVirtualBoundaries(); i++)
    {
      if (xPos <= picHeader->getVirtualBoundariesPosX(i) && picHeader->getVirtualBoundariesPosX(i) <= xPos + width)
      {
        verVirBndryPos[numVerVirBndry++] = picHeader->getVirtualBoundariesPosX(i);
      }
    }
  }
  return numHorVirBndry > 0 || numVerVirBndry > 0 ;
}
//! \}
