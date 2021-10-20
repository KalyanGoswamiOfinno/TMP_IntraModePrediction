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

/** \file     EncCu.cpp
    \brief    Coding Unit (CU) encoder class
*/

#include "EncCu.h"

#include "EncLib.h"
#include "Analyze.h"
#include "AQp.h"

#include "CommonLib/dtrace_codingstruct.h"
#include "CommonLib/Picture.h"
#include "CommonLib/UnitTools.h"
#include "MCTS.h"


#include "CommonLib/dtrace_buffer.h"

#include <stdio.h>
#include <cmath>
#include <algorithm>



//! \ingroup EncoderLib
//! \{

// ====================================================================================================================
EncCu::EncCu() : m_GeoModeTest
{
  GeoMotionInfo(0, 1), GeoMotionInfo(1, 0),GeoMotionInfo(0, 2), GeoMotionInfo(1, 2), GeoMotionInfo(2, 0),
  GeoMotionInfo(2, 1), GeoMotionInfo(0, 3),GeoMotionInfo(1, 3), GeoMotionInfo(2, 3), GeoMotionInfo(3, 0),
  GeoMotionInfo(3, 1), GeoMotionInfo(3, 2),GeoMotionInfo(0, 4), GeoMotionInfo(1, 4), GeoMotionInfo(2, 4),
  GeoMotionInfo(3, 4), GeoMotionInfo(4, 0),GeoMotionInfo(4, 1), GeoMotionInfo(4, 2), GeoMotionInfo(4, 3),
  GeoMotionInfo(0, 5), GeoMotionInfo(1, 5),GeoMotionInfo(2, 5), GeoMotionInfo(3, 5), GeoMotionInfo(4, 5),
  GeoMotionInfo(5, 0), GeoMotionInfo(5, 1),GeoMotionInfo(5, 2), GeoMotionInfo(5, 3), GeoMotionInfo(5, 4)
}
{
#if NON_ADJACENT_MRG_CAND
  int numGeoModeTest = 0;
  for (int i = 1; i < GEO_MAX_NUM_UNI_CANDS; i++)
  {
    for (int j = 0; j < i; j++)
    {
      m_GeoModeTest[numGeoModeTest] = GeoMotionInfo(j, i);
      numGeoModeTest++;
    }
    for (int j = 0; j < i; j++)
    {
      m_GeoModeTest[numGeoModeTest] = GeoMotionInfo(i, j);
      numGeoModeTest++;
    }
  }
#endif
#if JVET_W0097_GPM_MMVD_TM
  fastGpmMmvdSearch = false;
  fastGpmMmvdRelatedCU = false;
  includeMoreMMVDCandFirstPass = false;
  maxNumGPMDirFirstPass = 64;
  numCandPerPar = 5;
#endif
}

void EncCu::create( EncCfg* encCfg )
{
  unsigned      uiMaxWidth    = encCfg->getMaxCUWidth();
  unsigned      uiMaxHeight   = encCfg->getMaxCUHeight();
  ChromaFormat  chromaFormat  = encCfg->getChromaFormatIdc();

  unsigned      numWidths     = gp_sizeIdxInfo->numWidths();
  unsigned      numHeights    = gp_sizeIdxInfo->numHeights();
  m_pTempCS = new CodingStructure**  [numWidths];
  m_pBestCS = new CodingStructure**  [numWidths];
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  m_pTempCS2 = new CodingStructure** [numWidths];
  m_pBestCS2 = new CodingStructure** [numWidths];
#endif
#if ENABLE_OBMC
  m_tempWoOBMCBuffer.create(UnitArea(chromaFormat, Area(0, 0, MAX_CU_SIZE, MAX_CU_SIZE)));
#endif
  for( unsigned w = 0; w < numWidths; w++ )
  {
    m_pTempCS[w] = new CodingStructure*  [numHeights];
    m_pBestCS[w] = new CodingStructure*  [numHeights];
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
    m_pTempCS2[w] = new CodingStructure* [numHeights];
    m_pBestCS2[w] = new CodingStructure* [numHeights];
#endif
    for( unsigned h = 0; h < numHeights; h++ )
    {
      unsigned width  = gp_sizeIdxInfo->sizeFrom( w );
      unsigned height = gp_sizeIdxInfo->sizeFrom( h );

      if( gp_sizeIdxInfo->isCuSize( width ) && gp_sizeIdxInfo->isCuSize( height ) )
      {
        m_pTempCS[w][h] = new CodingStructure( m_unitCache.cuCache, m_unitCache.puCache, m_unitCache.tuCache );
        m_pBestCS[w][h] = new CodingStructure( m_unitCache.cuCache, m_unitCache.puCache, m_unitCache.tuCache );

        m_pTempCS[w][h]->create(chromaFormat, Area(0, 0, width, height), false, (bool)encCfg->getPLTMode());
        m_pBestCS[w][h]->create(chromaFormat, Area(0, 0, width, height), false, (bool)encCfg->getPLTMode());
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
        m_pTempCS2[w][h] = new CodingStructure( m_unitCache.cuCache, m_unitCache.puCache, m_unitCache.tuCache );
        m_pBestCS2[w][h] = new CodingStructure( m_unitCache.cuCache, m_unitCache.puCache, m_unitCache.tuCache );

        m_pTempCS2[w][h]->create(chromaFormat, Area(0, 0, width, height), false, (bool)encCfg->getPLTMode());
        m_pBestCS2[w][h]->create(chromaFormat, Area(0, 0, width, height), false, (bool)encCfg->getPLTMode());
#endif
      }
      else
      {
        m_pTempCS[w][h] = nullptr;
        m_pBestCS[w][h] = nullptr;
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
        m_pTempCS2[w][h] = nullptr;
        m_pBestCS2[w][h] = nullptr;
#endif
      }
    }
  }
#if ENABLE_OBMC
  m_pTempCUWoOBMC = nullptr;
  m_pPredBufWoOBMC = nullptr;

  if (encCfg->getUseOBMC())
  {
    m_pTempCUWoOBMC = new CodingStructure**[numWidths];
    m_pPredBufWoOBMC = new PelStorage*[numWidths];

    for (unsigned w = 0; w < numWidths; w++)
    {
      m_pTempCUWoOBMC[w] = new CodingStructure*[numHeights];
      m_pPredBufWoOBMC[w] = new PelStorage[numHeights];

      for (unsigned h = 0; h < numHeights; h++)
      {
        uint32_t width = gp_sizeIdxInfo->sizeFrom(w);
        uint32_t height = gp_sizeIdxInfo->sizeFrom(h);

        if (gp_sizeIdxInfo->isCuSize(width) && gp_sizeIdxInfo->isCuSize(height))
        {
          m_pTempCUWoOBMC[w][h] = new CodingStructure(m_unitCache.cuCache, m_unitCache.puCache, m_unitCache.tuCache);
          m_pTempCUWoOBMC[w][h]->create(chromaFormat, Area(0, 0, width, height), false, false);
          m_pPredBufWoOBMC[w][h].create(UnitArea(chromaFormat, Area(0, 0, width, height)));
        }
      }
    }
  }
#endif
  m_cuChromaQpOffsetIdxPlus1 = 0;

  unsigned maxDepth = numWidths + numHeights;

  m_modeCtrl = new EncModeCtrlMTnoRQT();

  m_modeCtrl->create( *encCfg );

  for (unsigned ui = 0; ui < MMVD_MRG_MAX_RD_BUF_NUM; ui++)
  {
    m_acMergeBuffer[ui].create( chromaFormat, Area( 0, 0, uiMaxWidth, uiMaxHeight ) );
  }
  for (unsigned ui = 0; ui < MRG_MAX_NUM_CANDS; ui++)
  {
    m_acRealMergeBuffer[ui].create(chromaFormat, Area(0, 0, uiMaxWidth, uiMaxHeight));
#if INTER_LIC || MULTI_HYP_PRED
    m_acRealMergeBuffer[ui+MRG_MAX_NUM_CANDS].create(chromaFormat, Area(0, 0, uiMaxWidth, uiMaxHeight));
#endif
    m_acMergeTmpBuffer[ui].create(chromaFormat, Area(0, 0, uiMaxWidth, uiMaxHeight));
  }

  for( unsigned ui = 0; ui < GEO_MAX_TRY_WEIGHTED_SAD; ui++ )
  {
    m_acGeoWeightedBuffer[ui].create( chromaFormat, Area( 0, 0, uiMaxWidth, uiMaxHeight ) );
  }
#if JVET_W0097_GPM_MMVD_TM
  for (unsigned ui = 0; ui < MRG_MAX_NUM_CANDS; ui++)
  {
    for (unsigned vi = 0; vi < GPM_EXT_MMVD_MAX_REFINE_NUM; vi++)
    {
      m_acGeoMMVDBuffer[ui][vi].create(chromaFormat, Area(0, 0, uiMaxWidth, uiMaxHeight));
      m_acGeoMMVDTmpBuffer[ui][vi].create(chromaFormat, Area(0, 0, uiMaxWidth, uiMaxHeight));
    }
  }
  int sourceWidth = encCfg->getSourceWidth();
  int sourceHeight = encCfg->getSourceHeight();
  fastGpmMmvdSearch = (((encCfg->getIntraPeriod() > 0) && ((sourceWidth * sourceHeight) <= (1920 * 1080))) || ((encCfg->getIntraPeriod() < 0) && ((sourceWidth * sourceHeight) >= (1280 * 720)))) && !encCfg->getIBCMode();
#if TM_MRG
  fastGpmMmvdRelatedCU = ((encCfg->getIntraPeriod() < 0) && ((sourceWidth * sourceHeight) >= (1920 * 1080))) && !encCfg->getIBCMode();
#else
  fastGpmMmvdRelatedCU = ((encCfg->getIntraPeriod() < 0) && ((sourceWidth * sourceHeight) >= (1280 * 720))) && !encCfg->getIBCMode();
#endif

  includeMoreMMVDCandFirstPass = ((encCfg->getIntraPeriod() > 0) || ((encCfg->getIntraPeriod() < 0) && fastGpmMmvdSearch));
  maxNumGPMDirFirstPass = ((encCfg->getIntraPeriod() < 0) ? 50 : (fastGpmMmvdSearch ? 36 : 64));
  numCandPerPar = (fastGpmMmvdSearch ? 4 : 5);
#if TM_MRG
  for (uint16_t ui = 0; ui < GEO_TM_MAX_NUM_CANDS; ui++)
  {
    m_acGeoMergeTmpBuffer[ui].create(chromaFormat, Area(0, 0, uiMaxWidth, uiMaxHeight));
    m_acGeoSADTmpBuffer[ui].create(chromaFormat, Area(0, 0, uiMaxWidth, uiMaxHeight));
  }
#endif
#endif
  m_ciipBuffer[0].create(chromaFormat, Area(0, 0, uiMaxWidth, uiMaxHeight));
  m_ciipBuffer[1].create(chromaFormat, Area(0, 0, uiMaxWidth, uiMaxHeight));

  m_CtxBuffer.resize( maxDepth );
  m_CurrCtx = 0;
}


void EncCu::destroy()
{
  unsigned numWidths  = gp_sizeIdxInfo->numWidths();
  unsigned numHeights = gp_sizeIdxInfo->numHeights();

#if JVET_V0094_BILATERAL_FILTER
  delete m_bilateralFilter;
#endif

#if ENABLE_OBMC
  m_tempWoOBMCBuffer.destroy();
#endif
  for( unsigned w = 0; w < numWidths; w++ )
  {
    for( unsigned h = 0; h < numHeights; h++ )
    {
      if( m_pBestCS[w][h] ) m_pBestCS[w][h]->destroy();
      if( m_pTempCS[w][h] ) m_pTempCS[w][h]->destroy();

      delete m_pBestCS[w][h];
      delete m_pTempCS[w][h];
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
      if( m_pBestCS2[w][h] ) m_pBestCS2[w][h]->destroy();
      if( m_pTempCS2[w][h] ) m_pTempCS2[w][h]->destroy();

      delete m_pBestCS2[w][h];
      delete m_pTempCS2[w][h];
#endif
    }

    delete[] m_pTempCS[w];
    delete[] m_pBestCS[w];
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
    delete[] m_pTempCS2[w];
    delete[] m_pBestCS2[w];
#endif
  }

  delete[] m_pBestCS; m_pBestCS = nullptr;
  delete[] m_pTempCS; m_pTempCS = nullptr;
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  delete[] m_pBestCS2; m_pBestCS2 = nullptr;
  delete[] m_pTempCS2; m_pTempCS2 = nullptr;
#endif
#if REUSE_CU_RESULTS
  if (m_tmpStorageLCU)
  {
    m_tmpStorageLCU->destroy();
    delete m_tmpStorageLCU;  m_tmpStorageLCU = nullptr;
  }
#endif

#if REUSE_CU_RESULTS
  m_modeCtrl->destroy();

#endif
  delete m_modeCtrl;
  m_modeCtrl = nullptr;
#if ENABLE_OBMC
  if (m_pTempCUWoOBMC)
  {
    for (unsigned w = 0; w < numWidths; w++)
    {
      for (unsigned h = 0; h < numHeights; h++)
      {
        if ((w == h) && gp_sizeIdxInfo->isCuSize(gp_sizeIdxInfo->sizeFrom(w)) && gp_sizeIdxInfo->isCuSize(gp_sizeIdxInfo->sizeFrom(h)))
        {
          m_pTempCUWoOBMC[w][h]->destroy();
          delete m_pTempCUWoOBMC[w][h];

          m_pPredBufWoOBMC[w][h].destroy();
        }
      }
      delete[] m_pTempCUWoOBMC[w];
      delete[] m_pPredBufWoOBMC[w];
    }

    delete[] m_pTempCUWoOBMC;
    m_pTempCUWoOBMC = nullptr;

    delete[] m_pPredBufWoOBMC;
    m_pPredBufWoOBMC = nullptr;
  }

#endif
  for (unsigned ui = 0; ui < MMVD_MRG_MAX_RD_BUF_NUM; ui++)
  {
    m_acMergeBuffer[ui].destroy();
  }
  for (unsigned ui = 0; ui < MRG_MAX_NUM_CANDS; ui++)
  {
    m_acRealMergeBuffer[ui].destroy();
#if INTER_LIC || MULTI_HYP_PRED
    m_acRealMergeBuffer[ui+MRG_MAX_NUM_CANDS].destroy();
#endif
    m_acMergeTmpBuffer[ui].destroy();
  }

  for (unsigned ui = 0; ui < GEO_MAX_TRY_WEIGHTED_SAD; ui++)
  {
    m_acGeoWeightedBuffer[ui].destroy();
  }
#if JVET_W0097_GPM_MMVD_TM
  for (unsigned ui = 0; ui < MRG_MAX_NUM_CANDS; ui++)
  {
    for (unsigned vi = 0; vi < GPM_EXT_MMVD_MAX_REFINE_NUM; vi++)
    {
      m_acGeoMMVDBuffer[ui][vi].destroy();
      m_acGeoMMVDTmpBuffer[ui][vi].destroy();
    }
  }
#if TM_MRG
  for (uint16_t ui = 0; ui < GEO_TM_MAX_NUM_CANDS; ui++)
  {
    m_acGeoMergeTmpBuffer[ui].destroy();
    m_acGeoSADTmpBuffer[ui].destroy();
  }
#endif
#endif
  m_ciipBuffer[0].destroy();
  m_ciipBuffer[1].destroy();
}



EncCu::~EncCu()
{
}



/** \param    pcEncLib      pointer of encoder class
 */
void EncCu::init( EncLib* pcEncLib, const SPS& sps PARL_PARAM( const int tId ) )
{
  m_pcEncCfg           = pcEncLib;
  m_pcIntraSearch      = pcEncLib->getIntraSearch( PARL_PARAM0( tId ) );
  m_pcInterSearch      = pcEncLib->getInterSearch( PARL_PARAM0( tId ) );
  m_pcTrQuant          = pcEncLib->getTrQuant( PARL_PARAM0( tId ) );
  m_pcRdCost           = pcEncLib->getRdCost ( PARL_PARAM0( tId ) );
  m_CABACEstimator     = pcEncLib->getCABACEncoder( PARL_PARAM0( tId ) )->getCABACEstimator( &sps );
  m_CABACEstimator->setEncCu(this);
  m_CtxCache           = pcEncLib->getCtxCache( PARL_PARAM0( tId ) );
  m_pcRateCtrl         = pcEncLib->getRateCtrl();
  m_pcSliceEncoder     = pcEncLib->getSliceEncoder();
#if ENABLE_SPLIT_PARALLELISM
  m_pcEncLib           = pcEncLib;
  m_dataId             = tId;
#endif
  m_pcLoopFilter       = pcEncLib->getLoopFilter();

  m_GeoCostList.init(GEO_NUM_PARTITION_MODE, m_pcEncCfg->getMaxNumGeoCand());
  m_AFFBestSATDCost = MAX_DOUBLE;

  DecCu::init( m_pcTrQuant, m_pcIntraSearch, m_pcInterSearch );

  m_modeCtrl->init( m_pcEncCfg, m_pcRateCtrl, m_pcRdCost );

  m_pcInterSearch->setModeCtrl( m_modeCtrl );
  m_modeCtrl->setInterSearch(m_pcInterSearch);
  m_pcIntraSearch->setModeCtrl( m_modeCtrl );

}

// ====================================================================================================================
// Public member functions
// ====================================================================================================================

void EncCu::compressCtu( CodingStructure& cs, const UnitArea& area, const unsigned ctuRsAddr, const int prevQP[], const int currQP[] )
{
  m_modeCtrl->initCTUEncoding( *cs.slice );
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  cs.treeType = TREE_D;
#endif
  cs.slice->m_mapPltCost[0].clear();
  cs.slice->m_mapPltCost[1].clear();
#if ENABLE_SPLIT_PARALLELISM
  if( m_pcEncCfg->getNumSplitThreads() > 1 )
  {
    for( int jId = 1; jId < NUM_RESERVERD_SPLIT_JOBS; jId++ )
    {
      EncCu*            jobEncCu  = m_pcEncLib->getCuEncoder( cs.picture->scheduler.getSplitDataId( jId ) );
      CacheBlkInfoCtrl* cacheCtrl = dynamic_cast< CacheBlkInfoCtrl* >( jobEncCu->m_modeCtrl );
#if REUSE_CU_RESULTS
      BestEncInfoCache* bestCache = dynamic_cast< BestEncInfoCache* >( jobEncCu->m_modeCtrl );
#endif
      SaveLoadEncInfoSbt *sbtCache = dynamic_cast< SaveLoadEncInfoSbt* >( jobEncCu->m_modeCtrl );
      if( cacheCtrl )
      {
        cacheCtrl->init( *cs.slice );
      }
#if REUSE_CU_RESULTS
      if (bestCache)
      {
        bestCache->init(*cs.slice);
      }
#endif
      if (sbtCache)
      {
        sbtCache->init(*cs.slice);
      }
    }
  }

#if REUSE_CU_RESULTS
  if( auto* cacheCtrl = dynamic_cast<BestEncInfoCache*>( m_modeCtrl ) ) { cacheCtrl->tick(); }
#endif
  if( auto* cacheCtrl = dynamic_cast<CacheBlkInfoCtrl*>( m_modeCtrl ) ) { cacheCtrl->tick(); }
#endif
  // init the partitioning manager
  QTBTPartitioner partitioner;
  partitioner.initCtu(area, CH_L, *cs.slice);
  if (m_pcEncCfg->getIBCMode())
  {
    if (area.lx() == 0 && area.ly() == 0)
    {
      m_pcInterSearch->resetIbcSearch();
    }
    m_pcInterSearch->resetCtuRecord();
    m_ctuIbcSearchRangeX = m_pcEncCfg->getIBCLocalSearchRangeX();
    m_ctuIbcSearchRangeY = m_pcEncCfg->getIBCLocalSearchRangeY();
  }
  if (m_pcEncCfg->getIBCMode() && m_pcEncCfg->getIBCHashSearch() && (m_pcEncCfg->getIBCFastMethod() & IBC_FAST_METHOD_ADAPTIVE_SEARCHRANGE))
  {
    const int hashHitRatio = m_ibcHashMap.getHashHitRatio(area.Y()); // in percent
    if (hashHitRatio < 5) // 5%
    {
      m_ctuIbcSearchRangeX >>= 1;
      m_ctuIbcSearchRangeY >>= 1;
    }
    if (cs.slice->getNumRefIdx(REF_PIC_LIST_0) > 0)
    {
      m_ctuIbcSearchRangeX >>= 1;
      m_ctuIbcSearchRangeY >>= 1;
    }
  }
  // init current context pointer
  m_CurrCtx = m_CtxBuffer.data();

  CodingStructure *tempCS = m_pTempCS[gp_sizeIdxInfo->idxFrom( area.lumaSize().width )][gp_sizeIdxInfo->idxFrom( area.lumaSize().height )];
  CodingStructure *bestCS = m_pBestCS[gp_sizeIdxInfo->idxFrom( area.lumaSize().width )][gp_sizeIdxInfo->idxFrom( area.lumaSize().height )];

  cs.initSubStructure(*tempCS, partitioner.chType, partitioner.currArea(), false);
  cs.initSubStructure(*bestCS, partitioner.chType, partitioner.currArea(), false);
  tempCS->currQP[CH_L] = bestCS->currQP[CH_L] =
  tempCS->baseQP       = bestCS->baseQP       = currQP[CH_L];
  tempCS->prevQP[CH_L] = bestCS->prevQP[CH_L] = prevQP[CH_L];

  xCompressCU(tempCS, bestCS, partitioner);

  
  //Kalyan start
   
  for (int i = 0; i < (int) bestCS->cus.size(); i++)
  {

    const bool tpmAllowed = (*cs.sps).getUseIntraTMP() && isLuma(partitioner.chType)
                            && ((bestCS->cus[i]->lfnstIdx == 0) || allowLfnstWithTmp());
    const bool testTpm = tpmAllowed && (bestCS->cus[i]->lwidth() <= (*cs.sps).getIntraTMPMaxSize()
                             && bestCS->cus[i]->lheight() <= (*cs.sps).getIntraTMPMaxSize());
    if (testTpm) 
    {
      CodingUnit      cu_cpy       = *bestCS->cus[i];
      RefTemplateType TemplateType = m_pcIntraSearch->GetRefTemplateType(cu_cpy, cu_cpy.blocks[COMPONENT_Y]);
      if (TemplateType != NO_TEMPLATE)
      {
          //auto pu = cu_cpy.firstPU;
          PredictionUnit pu = *bestCS->pus[i];
          m_pcTrQuant->getTargetTemplate(&cu_cpy, pu.lwidth(), pu.lheight(), TemplateType);
          m_pcTrQuant->candidateSearchIntra(&cu_cpy, pu.lwidth(), pu.lheight(), TemplateType);
          int  foundCandiNum = m_pcTrQuant->getVaildCandiNumTNP();
          if (foundCandiNum > 0)
          {
            printf("\n Kalyan Template is matched for i = %d \t size = %d", i, cu_cpy.lwidth());            

            int pX   = m_pcTrQuant->getTempLibFast().getX();
            int pY   = m_pcTrQuant->getTempLibFast().getY();
            PredictionUnit  pu_temp = *(bestCS->pus[i]);

            int intra_dir = bestCS->pus[i]->intraDir[0];
            //printf("\t intra_dir =%d ", intra_dir);

            Position        xyPos    = pu_temp.lumaPos();
            int             pX_new  = pX + xyPos.x;
            int             pY_new  = pY + xyPos.y;

            try
            {
              const CompArea &area = cu_cpy.blocks[COMPONENT_Y];
              CompArea        tmpArea(COMPONENT_Y, area.chromaFormat, Position(pX_new, pY_new), area.size());
              //PelBuf          tmpLuma = m_tmpStorageLCU->getBuf(tmpArea);
              //CompArea        areaPredReg(COMPONENT_Y, cu_cpy.chromaFormat, area);
              //PelBuf          tmpLuma1 = cs.getPredBuf(areaPredReg);
              //pu.setLumaPos(pX_new, pY_new);
              //pu_temp.setLumaPos(pX_new, pY_new);

              int timd_mode_tmp = m_pcIntraSearch->deriveTimdMode(bestCS->picture->getRecoBuf(tmpArea), area, cu_cpy);
              //printf("\t timd_mode_tmp =%d ", timd_mode_tmp);

              
              //m_pcIntraSearch->predIntraAng(COMPONENT_Y, tmpLuma1, pu);
              //int intra_dir_temp = pu_temp.intraDir[0];
              //int intra_dir_temp = pu.intraDir[0];
              //printf("\t intra_dir_temp =%d ", intra_dir_temp);

              if (intra_dir == timd_mode_tmp)
                printf("\t ------ Matching Happen -----");
              if (bestCS->cus[i]->tmpFlag)
                printf("\t TMP is choosen for i = %d", i);
            }
            catch(exception) 
            {
              continue;
            }
            
            //const CPelBuf &srcBuf = CPelBuf(getPredictorPtr(COMPONENT_Y), srcStride, srcHStride);
            //IntraPrediction::xPredIntraPlanar(srcBuf, piPred);

          }
      }      
    }
    //if (bestCS->cus[i]->tmpFlag)
      //printf("\t Kalyan TMP is choosen for i = %d", i);
    //else
    //{      
      //int intra_dir = bestCS->pus[i]->intraDir[0];
      //printf("intra_dir =%d ", intra_dir);
      
    //}
  }// Kalyan end
  
  cs.slice->m_mapPltCost[0].clear();
  cs.slice->m_mapPltCost[1].clear();
  // all signals were already copied during compression if the CTU was split - at this point only the structures are copied to the top level CS
  const bool copyUnsplitCTUSignals = bestCS->cus.size() == 1;
  cs.useSubStructure(*bestCS, partitioner.chType, CS::getArea(*bestCS, area, partitioner.chType), copyUnsplitCTUSignals,
                     false, false, copyUnsplitCTUSignals, true);

  if (CS::isDualITree (cs) && isChromaEnabled (cs.pcv->chrFormat))
  {
    m_CABACEstimator->getCtx() = m_CurrCtx->start;

    partitioner.initCtu(area, CH_C, *cs.slice);

    cs.initSubStructure(*tempCS, partitioner.chType, partitioner.currArea(), false);
    cs.initSubStructure(*bestCS, partitioner.chType, partitioner.currArea(), false);
    tempCS->currQP[CH_C] = bestCS->currQP[CH_C] =
    tempCS->baseQP       = bestCS->baseQP       = currQP[CH_C];
    tempCS->prevQP[CH_C] = bestCS->prevQP[CH_C] = prevQP[CH_C];

    xCompressCU(tempCS, bestCS, partitioner);

    const bool copyUnsplitCTUSignals = bestCS->cus.size() == 1;
    cs.useSubStructure(*bestCS, partitioner.chType, CS::getArea(*bestCS, area, partitioner.chType),
                       copyUnsplitCTUSignals, false, false, copyUnsplitCTUSignals, true);
  }

  if (m_pcEncCfg->getUseRateCtrl())
  {
    (m_pcRateCtrl->getRCPic()->getLCU(ctuRsAddr)).m_actualMSE = (double)bestCS->dist / (double)m_pcRateCtrl->getRCPic()->getLCU(ctuRsAddr).m_numberOfPixel;
  }
  // reset context states and uninit context pointer
  m_CABACEstimator->getCtx() = m_CurrCtx->start;
  m_CurrCtx                  = 0;


  // Ensure that a coding was found
  // Selected mode's RD-cost must be not MAX_DOUBLE.
  CHECK( bestCS->cus.empty()                                   , "No possible encoding found" );
  CHECK( bestCS->cus[0]->predMode == NUMBER_OF_PREDICTION_MODES, "No possible encoding found" );
  CHECK( bestCS->cost             == MAX_DOUBLE                , "No possible encoding found" );
}

// ====================================================================================================================
// Protected member functions
// ====================================================================================================================

static int xCalcHADs8x8_ISlice(const Pel *piOrg, const int iStrideOrg)
{
  int k, i, j, jj;
  int diff[64], m1[8][8], m2[8][8], m3[8][8], iSumHad = 0;

  for (k = 0; k < 64; k += 8)
  {
    diff[k + 0] = piOrg[0];
    diff[k + 1] = piOrg[1];
    diff[k + 2] = piOrg[2];
    diff[k + 3] = piOrg[3];
    diff[k + 4] = piOrg[4];
    diff[k + 5] = piOrg[5];
    diff[k + 6] = piOrg[6];
    diff[k + 7] = piOrg[7];

    piOrg += iStrideOrg;
  }

  //horizontal
  for (j = 0; j < 8; j++)
  {
    jj = j << 3;
    m2[j][0] = diff[jj    ] + diff[jj + 4];
    m2[j][1] = diff[jj + 1] + diff[jj + 5];
    m2[j][2] = diff[jj + 2] + diff[jj + 6];
    m2[j][3] = diff[jj + 3] + diff[jj + 7];
    m2[j][4] = diff[jj    ] - diff[jj + 4];
    m2[j][5] = diff[jj + 1] - diff[jj + 5];
    m2[j][6] = diff[jj + 2] - diff[jj + 6];
    m2[j][7] = diff[jj + 3] - diff[jj + 7];

    m1[j][0] = m2[j][0] + m2[j][2];
    m1[j][1] = m2[j][1] + m2[j][3];
    m1[j][2] = m2[j][0] - m2[j][2];
    m1[j][3] = m2[j][1] - m2[j][3];
    m1[j][4] = m2[j][4] + m2[j][6];
    m1[j][5] = m2[j][5] + m2[j][7];
    m1[j][6] = m2[j][4] - m2[j][6];
    m1[j][7] = m2[j][5] - m2[j][7];

    m2[j][0] = m1[j][0] + m1[j][1];
    m2[j][1] = m1[j][0] - m1[j][1];
    m2[j][2] = m1[j][2] + m1[j][3];
    m2[j][3] = m1[j][2] - m1[j][3];
    m2[j][4] = m1[j][4] + m1[j][5];
    m2[j][5] = m1[j][4] - m1[j][5];
    m2[j][6] = m1[j][6] + m1[j][7];
    m2[j][7] = m1[j][6] - m1[j][7];
  }

  //vertical
  for (i = 0; i < 8; i++)
  {
    m3[0][i] = m2[0][i] + m2[4][i];
    m3[1][i] = m2[1][i] + m2[5][i];
    m3[2][i] = m2[2][i] + m2[6][i];
    m3[3][i] = m2[3][i] + m2[7][i];
    m3[4][i] = m2[0][i] - m2[4][i];
    m3[5][i] = m2[1][i] - m2[5][i];
    m3[6][i] = m2[2][i] - m2[6][i];
    m3[7][i] = m2[3][i] - m2[7][i];

    m1[0][i] = m3[0][i] + m3[2][i];
    m1[1][i] = m3[1][i] + m3[3][i];
    m1[2][i] = m3[0][i] - m3[2][i];
    m1[3][i] = m3[1][i] - m3[3][i];
    m1[4][i] = m3[4][i] + m3[6][i];
    m1[5][i] = m3[5][i] + m3[7][i];
    m1[6][i] = m3[4][i] - m3[6][i];
    m1[7][i] = m3[5][i] - m3[7][i];

    m2[0][i] = m1[0][i] + m1[1][i];
    m2[1][i] = m1[0][i] - m1[1][i];
    m2[2][i] = m1[2][i] + m1[3][i];
    m2[3][i] = m1[2][i] - m1[3][i];
    m2[4][i] = m1[4][i] + m1[5][i];
    m2[5][i] = m1[4][i] - m1[5][i];
    m2[6][i] = m1[6][i] + m1[7][i];
    m2[7][i] = m1[6][i] - m1[7][i];
  }

  for (i = 0; i < 8; i++)
  {
    for (j = 0; j < 8; j++)
    {
      iSumHad += abs(m2[i][j]);
    }
  }
  iSumHad -= abs(m2[0][0]);
  iSumHad = (iSumHad + 2) >> 2;
  return(iSumHad);
}

int  EncCu::updateCtuDataISlice(const CPelBuf buf)
{
  int  xBl, yBl;
  const int iBlkSize = 8;
  const Pel* pOrgInit = buf.buf;
  int  iStrideOrig = buf.stride;

  int iSumHad = 0;
  for( yBl = 0; ( yBl + iBlkSize ) <= buf.height; yBl += iBlkSize )
  {
    for( xBl = 0; ( xBl + iBlkSize ) <= buf.width; xBl += iBlkSize )
    {
      const Pel* pOrg = pOrgInit + iStrideOrig*yBl + xBl;
      iSumHad += xCalcHADs8x8_ISlice( pOrg, iStrideOrig );
    }
  }
  return( iSumHad );
}

bool EncCu::xCheckBestMode( CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &partitioner, const EncTestMode& encTestMode )
{
  bool bestCSUpdated = false;

  if( !tempCS->cus.empty() )
  {
    if( tempCS->cus.size() == 1 )
    {
      const CodingUnit& cu = *tempCS->cus.front();
      CHECK( cu.skip && !cu.firstPU->mergeFlag, "Skip flag without a merge flag is not allowed!" );
    }

#if WCG_EXT
    DTRACE_BEST_MODE( tempCS, bestCS, m_pcRdCost->getLambda( true ) );
#else
    DTRACE_BEST_MODE( tempCS, bestCS, m_pcRdCost->getLambda() );
#endif

#if MULTI_HYP_PRED
    if (tempCS->sps->getUseInterMultiHyp() && tempCS->slice->isInterB())
    {
      m_baseResultsForMH.insert(m_baseResultsForMH.end(), tempCS->m_MEResults.begin(), tempCS->m_MEResults.end());
#if MULTI_HYP_PRED
      tempCS->m_MEResults.clear(); // avoid duplicate insert
#endif
    }
#endif
    if( m_modeCtrl->useModeResult( encTestMode, tempCS, partitioner ) )
    {

      std::swap( tempCS, bestCS );
      // store temp best CI for next CU coding
      m_CurrCtx->best = m_CABACEstimator->getCtx();
      m_bestModeUpdated = true;
      bestCSUpdated = true;
    }
  }

  // reset context states
  m_CABACEstimator->getCtx() = m_CurrCtx->start;
  return bestCSUpdated;

}

void EncCu::xCompressCU( CodingStructure*& tempCS, CodingStructure*& bestCS, Partitioner& partitioner, double maxCostAllowed )
{
  CHECK(maxCostAllowed < 0, "Wrong value of maxCostAllowed!");
#if ENABLE_SPLIT_PARALLELISM
  CHECK( m_dataId != tempCS->picture->scheduler.getDataId(), "Working in the wrong dataId!" );

  if( m_pcEncCfg->getNumSplitThreads() != 1 && tempCS->picture->scheduler.getSplitJobId() == 0 )
  {
    if( m_modeCtrl->isParallelSplit( *tempCS, partitioner ) )
    {
      m_modeCtrl->setParallelSplit( true );
      xCompressCUParallel( tempCS, bestCS, partitioner );
      return;
    }
  }

#endif
  uint32_t compBegin;
  uint32_t numComp;
  bool jointPLT = false;
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  if (CS::isDualITree(*bestCS))
#else
  if (partitioner.isSepTree( *tempCS ))
#endif
  {
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
    if( !CS::isDualITree(*tempCS) && partitioner.treeType != TREE_D )
    {
      compBegin = COMPONENT_Y;
      numComp = (tempCS->area.chromaFormat != CHROMA_400)?3: 1;
      jointPLT = true;
    }
    else
#endif
    {
    if (isLuma(partitioner.chType))
    {
      compBegin = COMPONENT_Y;
      numComp = 1;
    }
    else
    {
      compBegin = COMPONENT_Cb;
      numComp = 2;
    }
    }
  }
  else
  {
    compBegin = COMPONENT_Y;
    numComp = (tempCS->area.chromaFormat != CHROMA_400) ? 3 : 1;
    jointPLT = true;
  }
  SplitSeries splitmode = -1;
  uint8_t   bestLastPLTSize[MAX_NUM_CHANNEL_TYPE];
  Pel       bestLastPLT[MAX_NUM_COMPONENT][MAXPLTPREDSIZE]; // store LastPLT for
  uint8_t   curLastPLTSize[MAX_NUM_CHANNEL_TYPE];
  Pel       curLastPLT[MAX_NUM_COMPONENT][MAXPLTPREDSIZE]; // store LastPLT if no partition
  for (int i = compBegin; i < (compBegin + numComp); i++)
  {
    ComponentID comID = jointPLT ? (ComponentID)compBegin : ((i > 0) ? COMPONENT_Cb : COMPONENT_Y);
    bestLastPLTSize[comID] = 0;
    curLastPLTSize[comID] = tempCS->prevPLT.curPLTSize[comID];
    memcpy(curLastPLT[i], tempCS->prevPLT.curPLT[i], tempCS->prevPLT.curPLTSize[comID] * sizeof(Pel));
  }

  Slice&   slice      = *tempCS->slice;
  const PPS &pps      = *tempCS->pps;
  const SPS &sps      = *tempCS->sps;
  const uint32_t uiLPelX  = tempCS->area.Y().lumaPos().x;
  const uint32_t uiTPelY  = tempCS->area.Y().lumaPos().y;
#if ENABLE_OBMC
  const unsigned wIdx = gp_sizeIdxInfo->idxFrom(partitioner.currArea().lwidth());
  const unsigned hIdx = gp_sizeIdxInfo->idxFrom(partitioner.currArea().lheight());
#endif
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  const ModeType modeTypeParent  = partitioner.modeType;
  const TreeType treeTypeParent  = partitioner.treeType;
  const ChannelType chTypeParent = partitioner.chType;
#endif
  const UnitArea currCsArea = clipArea( CS::getArea( *bestCS, bestCS->area, partitioner.chType ), *tempCS->picture );
#if MULTI_HYP_PRED
  m_baseResultsForMH.clear();
#endif
#if ENABLE_OBMC
  if (m_pTempCUWoOBMC && !slice.isIntra())
  {
    tempCS->initSubStructure(*m_pTempCUWoOBMC[wIdx][hIdx], partitioner.chType, partitioner.currArea(), false);
  }
#endif
  m_modeCtrl->initCULevel( partitioner, *tempCS );
  if( partitioner.currQtDepth == 0 && partitioner.currMtDepth == 0 && !tempCS->slice->isIntra() && ( sps.getUseSBT() || sps.getUseInterMTS() ) )
  {
    auto slsSbt = dynamic_cast<SaveLoadEncInfoSbt*>( m_modeCtrl );
    int maxSLSize = sps.getUseSBT() ? tempCS->slice->getSPS()->getMaxTbSize() : MTS_INTER_MAX_CU_SIZE;
    slsSbt->resetSaveloadSbt( maxSLSize );
#if ENABLE_SPLIT_PARALLELISM
    CHECK( tempCS->picture->scheduler.getSplitJobId() != 0, "The SBT search reset need to happen in sequential region." );
    if (m_pcEncCfg->getNumSplitThreads() > 1)
    {
      for (int jId = 1; jId < NUM_RESERVERD_SPLIT_JOBS; jId++)
      {
        auto slsSbt = dynamic_cast<SaveLoadEncInfoSbt *>(m_pcEncLib->getCuEncoder(jId)->m_modeCtrl);
        slsSbt->resetSaveloadSbt(maxSLSize);
      }
    }
#endif
  }
  m_sbtCostSave[0] = m_sbtCostSave[1] = MAX_DOUBLE;

  m_CurrCtx->start = m_CABACEstimator->getCtx();

  m_cuChromaQpOffsetIdxPlus1 = 0;

  if( slice.getUseChromaQpAdj() )
  {
    // TODO M0133 : double check encoder decisions with respect to chroma QG detection and actual encode
    int lgMinCuSize = sps.getLog2MinCodingBlockSize() +
      std::max<int>(0, floorLog2(sps.getCTUSize()) - sps.getLog2MinCodingBlockSize() - int(slice.getCuChromaQpOffsetSubdiv() / 2));
    if( partitioner.currQgChromaEnable() )
    {
      m_cuChromaQpOffsetIdxPlus1 = ( ( uiLPelX >> lgMinCuSize ) + ( uiTPelY >> lgMinCuSize ) ) % ( pps.getChromaQpOffsetListLen() + 1 );
    }
  }

  if( !m_modeCtrl->anyMode() )
  {
    m_modeCtrl->finishCULevel( partitioner );
    return;
  }

  DTRACE_UPDATE( g_trace_ctx, std::make_pair( "cux", uiLPelX ) );
  DTRACE_UPDATE( g_trace_ctx, std::make_pair( "cuy", uiTPelY ) );
  DTRACE_UPDATE( g_trace_ctx, std::make_pair( "cuw", tempCS->area.lwidth() ) );
  DTRACE_UPDATE( g_trace_ctx, std::make_pair( "cuh", tempCS->area.lheight() ) );
  DTRACE( g_trace_ctx, D_COMMON, "@(%4d,%4d) [%2dx%2d]\n", tempCS->area.lx(), tempCS->area.ly(), tempCS->area.lwidth(), tempCS->area.lheight() );


  m_pcInterSearch->resetSavedAffineMotion();
#if TM_AMVP
  if (!slice.isIntra())
  {
    m_pcInterSearch->clearTplAmvpBuffer();
  }
#endif
#if INTER_LIC
  m_pcInterSearch->m_fastLicCtrl.init();
#endif

#if MULTI_HYP_PRED
  const UnitArea localUnitArea(tempCS->area.chromaFormat, Area(0, 0, tempCS->area.Y().width, tempCS->area.Y().height));
  m_pcInterSearch->initMHPTmpBuffer(m_acRealMergeBuffer + MRG_MAX_NUM_CANDS, GEO_MAX_NUM_UNI_CANDS,
    m_acGeoWeightedBuffer, GEO_MAX_TRY_WEIGHTED_SAD,
    localUnitArea);
#endif
#if JVET_W0097_GPM_MMVD_TM
  m_mergeCandAvail = false;
#endif
  double bestIntPelCost = MAX_DOUBLE;

  if (tempCS->slice->getSPS()->getUseColorTrans())
  {
    tempCS->tmpColorSpaceCost = MAX_DOUBLE;
    bestCS->tmpColorSpaceCost = MAX_DOUBLE;
    tempCS->firstColorSpaceSelected = true;
    bestCS->firstColorSpaceSelected = true;
  }

  if (tempCS->slice->getSPS()->getUseColorTrans() && !CS::isDualITree(*tempCS))
  {
    tempCS->firstColorSpaceTestOnly = false;
    bestCS->firstColorSpaceTestOnly = false;
    tempCS->tmpColorSpaceIntraCost[0] = MAX_DOUBLE;
    tempCS->tmpColorSpaceIntraCost[1] = MAX_DOUBLE;
    bestCS->tmpColorSpaceIntraCost[0] = MAX_DOUBLE;
    bestCS->tmpColorSpaceIntraCost[1] = MAX_DOUBLE;

    if (tempCS->bestParent && tempCS->bestParent->firstColorSpaceTestOnly)
    {
      tempCS->firstColorSpaceTestOnly = bestCS->firstColorSpaceTestOnly = true;
    }
  }

  do
  {
    for (int i = compBegin; i < (compBegin + numComp); i++)
    {
      ComponentID comID = jointPLT ? (ComponentID)compBegin : ((i > 0) ? COMPONENT_Cb : COMPONENT_Y);
      tempCS->prevPLT.curPLTSize[comID] = curLastPLTSize[comID];
      memcpy(tempCS->prevPLT.curPLT[i], curLastPLT[i], curLastPLTSize[comID] * sizeof(Pel));
    }
    EncTestMode currTestMode = m_modeCtrl->currTestMode();
    currTestMode.maxCostAllowed = maxCostAllowed;
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
    if (pps.getUseDQP() && CS::isDualITree(*tempCS) && isChroma(partitioner.chType))
#else
    if (pps.getUseDQP() && partitioner.isSepTree(*tempCS) && isChroma( partitioner.chType ))
#endif
    {
      const Position chromaCentral(tempCS->area.Cb().chromaPos().offset(tempCS->area.Cb().chromaSize().width >> 1, tempCS->area.Cb().chromaSize().height >> 1));
      const Position lumaRefPos(chromaCentral.x << getComponentScaleX(COMPONENT_Cb, tempCS->area.chromaFormat), chromaCentral.y << getComponentScaleY(COMPONENT_Cb, tempCS->area.chromaFormat));
      const CodingStructure* baseCS = bestCS->picture->cs;
      const CodingUnit* colLumaCu = baseCS->getCU(lumaRefPos, CHANNEL_TYPE_LUMA);

      if (colLumaCu)
      {
        currTestMode.qp = colLumaCu->qp;
      }
    }

#if SHARP_LUMA_DELTA_QP || ENABLE_QPA_SUB_CTU
    if (partitioner.currQgEnable() && (
#if SHARP_LUMA_DELTA_QP
        (m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled()) ||
#endif
#if ENABLE_QPA_SUB_CTU
        (m_pcEncCfg->getUsePerceptQPA() && !m_pcEncCfg->getUseRateCtrl() && pps.getUseDQP())
#else
        false
#endif
      ))
    {
#if ENABLE_SPLIT_PARALLELISM
      CHECK( tempCS->picture->scheduler.getSplitJobId() > 0, "Changing lambda is only allowed in the master thread!" );
#endif
      if (currTestMode.qp >= 0)
      {
        updateLambda (&slice, currTestMode.qp,
 #if WCG_EXT && ER_CHROMA_QP_WCG_PPS
                      m_pcEncCfg->getWCGChromaQPControl().isEnabled(),
 #endif
                      CS::isDualITree (*tempCS) || (partitioner.currDepth == 0));
      }
    }
#endif

    if( currTestMode.type == ETM_INTER_ME )
    {
#if ENABLE_OBMC
      bool tryObmc = true;
#endif

      if( ( currTestMode.opts & ETO_IMV ) != 0 )
      {
        const bool skipAltHpelIF = ( int( ( currTestMode.opts & ETO_IMV ) >> ETO_IMV_SHIFT ) == 4 ) && ( bestIntPelCost > 1.25 * bestCS->cost );
        if (!skipAltHpelIF)
        {
          tempCS->bestCS = bestCS;
#if ENABLE_OBMC
          tryObmc = xCheckRDCostInterIMV(tempCS, bestCS, partitioner, currTestMode, bestIntPelCost);
#else
          xCheckRDCostInterIMV(tempCS, bestCS, partitioner, currTestMode, bestIntPelCost);
#endif
          tempCS->bestCS = nullptr;
        }
      }
      else
      {
        tempCS->bestCS = bestCS;
        xCheckRDCostInter( tempCS, bestCS, partitioner, currTestMode );
        tempCS->bestCS = nullptr;
      }
#if ENABLE_OBMC
      if (tryObmc && tempCS->cus.size() != 0)//todo
      {
        xCheckRDCostInterWoOBMC(tempCS, bestCS, partitioner, currTestMode);
      }
#endif
    }
    else if (currTestMode.type == ETM_HASH_INTER)
    {
      xCheckRDCostHashInter( tempCS, bestCS, partitioner, currTestMode );
    }
#if !MERGE_ENC_OPT
    else if( currTestMode.type == ETM_AFFINE )
    {
      xCheckRDCostAffineMerge2Nx2N( tempCS, bestCS, partitioner, currTestMode );
    }
#endif
#if AFFINE_MMVD && !MERGE_ENC_OPT
    else if (currTestMode.type == ETM_AF_MMVD)
    {
      xCheckRDCostAffineMmvd2Nx2N(tempCS, bestCS, partitioner, currTestMode);
    }
#endif
#if TM_MRG && !MERGE_ENC_OPT
    else if (currTestMode.type == ETM_MERGE_TM)
    {
      xCheckRDCostTMMerge2Nx2N(tempCS, bestCS, partitioner, currTestMode);
    }
#endif
#if REUSE_CU_RESULTS
    else if( currTestMode.type == ETM_RECO_CACHED )
    {
      xReuseCachedResult( tempCS, bestCS, partitioner );
    }
#endif
    else if( currTestMode.type == ETM_MERGE_SKIP )
    {
      xCheckRDCostMerge2Nx2N( tempCS, bestCS, partitioner, currTestMode );
      CodingUnit* cu = bestCS->getCU(partitioner.chType);
      if (cu)
      cu->mmvdSkip = cu->skip == false ? false : cu->mmvdSkip;
    }
    else if( currTestMode.type == ETM_MERGE_GEO )
    {
#if JVET_W0097_GPM_MMVD_TM
      CodedCUInfo    &relatedCU = ((EncModeCtrlMTnoRQT *)m_modeCtrl)->getBlkInfo(partitioner.currArea());
      if (!relatedCU.isGPMTested)
      {
        xCheckRDCostMergeGeoComb2Nx2N(tempCS, bestCS, partitioner, currTestMode);
        relatedCU.isGPMTested = 1;
      }
      else
      {
        xCheckRDCostMergeGeoComb2Nx2N(tempCS, bestCS, partitioner, currTestMode, true);
      }
#else
      xCheckRDCostMergeGeo2Nx2N( tempCS, bestCS, partitioner, currTestMode );
#endif
    }
#if MULTI_HYP_PRED
    else if (currTestMode.type == ETM_INTER_MULTIHYP)
    {
      xCheckRDCostInterMultiHyp2Nx2N(tempCS, bestCS, partitioner, currTestMode);
    }
#endif
    else if( currTestMode.type == ETM_INTRA )
    {
      if (slice.getSPS()->getUseColorTrans() && !CS::isDualITree(*tempCS))
      {
        bool skipSecColorSpace = false;
        skipSecColorSpace = xCheckRDCostIntra(tempCS, bestCS, partitioner, currTestMode, (m_pcEncCfg->getRGBFormatFlag() ? true : false));
        if ((m_pcEncCfg->getCostMode() == COST_LOSSLESS_CODING && slice.isLossless()) && !m_pcEncCfg->getRGBFormatFlag())
        {
          skipSecColorSpace = true;
        }
        if (!skipSecColorSpace && !tempCS->firstColorSpaceTestOnly)
        {
          xCheckRDCostIntra(tempCS, bestCS, partitioner, currTestMode, (m_pcEncCfg->getRGBFormatFlag() ? false : true));
        }

        if (!tempCS->firstColorSpaceTestOnly)
        {
          if (tempCS->tmpColorSpaceIntraCost[0] != MAX_DOUBLE && tempCS->tmpColorSpaceIntraCost[1] != MAX_DOUBLE)
          {
            double skipCostRatio = m_pcEncCfg->getRGBFormatFlag() ? 1.1 : 1.0;
            if (tempCS->tmpColorSpaceIntraCost[1] > (skipCostRatio*tempCS->tmpColorSpaceIntraCost[0]))
            {
              tempCS->firstColorSpaceTestOnly = bestCS->firstColorSpaceTestOnly = true;
            }
          }
        }
        else
        {
          CHECK(tempCS->tmpColorSpaceIntraCost[1] != MAX_DOUBLE, "the RD test of the second color space should be skipped");
        }
      }
      else
      {
        xCheckRDCostIntra(tempCS, bestCS, partitioner, currTestMode, false);
      }
    }
    else if (currTestMode.type == ETM_PALETTE)
    {
      xCheckPLT( tempCS, bestCS, partitioner, currTestMode );
    }
    else if (currTestMode.type == ETM_IBC)
    {
      xCheckRDCostIBCMode(tempCS, bestCS, partitioner, currTestMode);
    }
    else if (currTestMode.type == ETM_IBC_MERGE)
    {
      xCheckRDCostIBCModeMerge2Nx2N(tempCS, bestCS, partitioner, currTestMode);
    }
    else if( isModeSplit( currTestMode ) )
    {
      if (bestCS->cus.size() != 0)
      {
        splitmode = bestCS->cus[0]->splitSeries;
      }
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
      assert( partitioner.modeType == tempCS->modeType );
      int signalModeConsVal = tempCS->signalModeCons( getPartSplit( currTestMode ), partitioner, modeTypeParent );
      int numRoundRdo = signalModeConsVal == LDT_MODE_TYPE_SIGNAL ? 2 : 1;
      bool skipInterPass = false;
      for( int i = 0; i < numRoundRdo; i++ )
      {
        //change cons modes
        if( signalModeConsVal == LDT_MODE_TYPE_SIGNAL )
        {
          CHECK( numRoundRdo != 2, "numRoundRdo shall be 2 - [LDT_MODE_TYPE_SIGNAL]" );
          tempCS->modeType = partitioner.modeType = (i == 0) ? MODE_TYPE_INTER : MODE_TYPE_INTRA;
        }
        else if( signalModeConsVal == LDT_MODE_TYPE_INFER )
        {
          CHECK( numRoundRdo != 1, "numRoundRdo shall be 1 - [LDT_MODE_TYPE_INFER]" );
          tempCS->modeType = partitioner.modeType = MODE_TYPE_INTRA;
        }
        else if( signalModeConsVal == LDT_MODE_TYPE_INHERIT )
        {
          CHECK( numRoundRdo != 1, "numRoundRdo shall be 1 - [LDT_MODE_TYPE_INHERIT]" );
          tempCS->modeType = partitioner.modeType = modeTypeParent;
        }

        //for lite intra encoding fast algorithm, set the status to save inter coding info
        if( modeTypeParent == MODE_TYPE_ALL && tempCS->modeType == MODE_TYPE_INTER )
        {
          m_pcIntraSearch->setSaveCuCostInSCIPU( true );
          m_pcIntraSearch->setNumCuInSCIPU( 0 );
        }
        else if( modeTypeParent == MODE_TYPE_ALL && tempCS->modeType != MODE_TYPE_INTER )
        {
          m_pcIntraSearch->setSaveCuCostInSCIPU( false );
          if( tempCS->modeType == MODE_TYPE_ALL )
          {
            m_pcIntraSearch->setNumCuInSCIPU( 0 );
          }
        }

        xCheckModeSplit( tempCS, bestCS, partitioner, currTestMode, modeTypeParent, skipInterPass );
#else
      xCheckModeSplit(tempCS, bestCS, partitioner, currTestMode);
#endif
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
        //recover cons modes
        tempCS->modeType = partitioner.modeType = modeTypeParent;
        tempCS->treeType = partitioner.treeType = treeTypeParent;
        partitioner.chType = chTypeParent;
        if( modeTypeParent == MODE_TYPE_ALL )
        {
          m_pcIntraSearch->setSaveCuCostInSCIPU( false );
          if( numRoundRdo == 2 && tempCS->modeType == MODE_TYPE_INTRA )
          {
            m_pcIntraSearch->initCuAreaCostInSCIPU();
          }
        }
        if( skipInterPass )
        {
          break;
        }
      }
#endif
      if (splitmode != bestCS->cus[0]->splitSeries)
      {
        splitmode = bestCS->cus[0]->splitSeries;
        const CodingUnit&     cu = *bestCS->cus.front();
        cu.cs->prevPLT = bestCS->prevPLT;
        for (int i = compBegin; i < (compBegin + numComp); i++)
        {
          ComponentID comID = jointPLT ? (ComponentID)compBegin : ((i > 0) ? COMPONENT_Cb : COMPONENT_Y);
          bestLastPLTSize[comID] = bestCS->cus[0]->cs->prevPLT.curPLTSize[comID];
          memcpy(bestLastPLT[i], bestCS->cus[0]->cs->prevPLT.curPLT[i], bestCS->cus[0]->cs->prevPLT.curPLTSize[comID] * sizeof(Pel));
        }
      }
    }
    else
    {
      THROW( "Don't know how to handle mode: type = " << currTestMode.type << ", options = " << currTestMode.opts );
    }
  } while( m_modeCtrl->nextMode( *tempCS, partitioner ) );


  //////////////////////////////////////////////////////////////////////////
  // Finishing CU
#if ENABLE_SPLIT_PARALLELISM
  if( bestCS->cus.empty() )
  {
    CHECK( bestCS->cost != MAX_DOUBLE, "Cost should be maximal if no encoding found" );
    CHECK( bestCS->picture->scheduler.getSplitJobId() == 0, "Should always get a result in serial case" );

    m_modeCtrl->finishCULevel( partitioner );
    return;
  }

#endif
  if( tempCS->cost == MAX_DOUBLE && bestCS->cost == MAX_DOUBLE )
  {
    //although some coding modes were planned to be tried in RDO, no coding mode actually finished encoding due to early termination
    //thus tempCS->cost and bestCS->cost are both MAX_DOUBLE; in this case, skip the following process for normal case
    m_modeCtrl->finishCULevel( partitioner );
    return;
  }

  // set context states
  m_CABACEstimator->getCtx() = m_CurrCtx->best;

  // QP from last processed CU for further processing
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  //copy the qp of the last non-chroma CU
  int numCUInThisNode = (int)bestCS->cus.size();
  if( numCUInThisNode > 1 && bestCS->cus.back()->chType == CHANNEL_TYPE_CHROMA && !CS::isDualITree( *bestCS ) )
  {
    CHECK( bestCS->cus[numCUInThisNode-2]->chType != CHANNEL_TYPE_LUMA, "wrong chType" );
    bestCS->prevQP[partitioner.chType] = bestCS->cus[numCUInThisNode-2]->qp;
  }
  else
  {
#endif
  bestCS->prevQP[partitioner.chType] = bestCS->cus.back()->qp;
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  }
#endif
  if ((!slice.isIntra() || slice.getSPS()->getIBCFlag())
    && partitioner.chType == CHANNEL_TYPE_LUMA
    && bestCS->cus.size() == 1 && (bestCS->cus.back()->predMode == MODE_INTER || bestCS->cus.back()->predMode == MODE_IBC)
    && bestCS->area.Y() == (*bestCS->cus.back()).Y()
    )
  {
    const CodingUnit&     cu = *bestCS->cus.front();

    bool isIbcSmallBlk = CU::isIBC(cu) && (cu.lwidth() * cu.lheight() <= 16);
    CU::saveMotionInHMVP( cu, isIbcSmallBlk );
  }
  bestCS->picture->getPredBuf(currCsArea).copyFrom(bestCS->getPredBuf(currCsArea));
  bestCS->picture->getRecoBuf( currCsArea ).copyFrom( bestCS->getRecoBuf( currCsArea ) );
  m_modeCtrl->finishCULevel( partitioner );
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  if( m_pcIntraSearch->getSaveCuCostInSCIPU() && bestCS->cus.size() == 1 )
  {
    m_pcIntraSearch->saveCuAreaCostInSCIPU( Area( partitioner.currArea().lumaPos(), partitioner.currArea().lumaSize() ), bestCS->cost );
  }
#endif
#if ENABLE_SPLIT_PARALLELISM
  if( tempCS->picture->scheduler.getSplitJobId() == 0 && m_pcEncCfg->getNumSplitThreads() != 1 )
  {
    tempCS->picture->finishParallelPart( currCsArea );
  }

#endif
  if (bestCS->cus.size() == 1) // no partition
  {
    CHECK(bestCS->cus[0]->tileIdx != bestCS->pps->getTileIdx(bestCS->area.lumaPos()), "Wrong tile index!");
    if (bestCS->cus[0]->predMode == MODE_PLT)
    {
      for (int i = compBegin; i < (compBegin + numComp); i++)
      {
        ComponentID comID = jointPLT ? (ComponentID)compBegin : ((i > 0) ? COMPONENT_Cb : COMPONENT_Y);
        bestCS->prevPLT.curPLTSize[comID] = curLastPLTSize[comID];
        memcpy(bestCS->prevPLT.curPLT[i], curLastPLT[i], curLastPLTSize[comID] * sizeof(Pel));
      }
      bestCS->reorderPrevPLT(bestCS->prevPLT, bestCS->cus[0]->curPLTSize, bestCS->cus[0]->curPLT, bestCS->cus[0]->reuseflag, compBegin, numComp, jointPLT);
    }
    else
    {
      for (int i = compBegin; i<(compBegin + numComp); i++)
      {
        ComponentID comID = jointPLT ? (ComponentID)compBegin : ((i > 0) ? COMPONENT_Cb : COMPONENT_Y);
        bestCS->prevPLT.curPLTSize[comID] = curLastPLTSize[comID];
        memcpy(bestCS->prevPLT.curPLT[i], curLastPLT[i], bestCS->prevPLT.curPLTSize[comID] * sizeof(Pel));
      }
    }
  }
  else
  {
    for (int i = compBegin; i<(compBegin + numComp); i++)
    {
      ComponentID comID = jointPLT ? (ComponentID)compBegin : ((i > 0) ? COMPONENT_Cb : COMPONENT_Y);
      bestCS->prevPLT.curPLTSize[comID] = bestLastPLTSize[comID];
      memcpy(bestCS->prevPLT.curPLT[i], bestLastPLT[i], bestCS->prevPLT.curPLTSize[comID] * sizeof(Pel));
    }
  }
  const CodingUnit&     cu = *bestCS->cus.front();
  cu.cs->prevPLT = bestCS->prevPLT;
  // Assert if Best prediction mode is NONE
  // Selected mode's RD-cost must be not MAX_DOUBLE.
  CHECK( bestCS->cus.empty()                                   , "No possible encoding found" );
  CHECK( bestCS->cus[0]->predMode == NUMBER_OF_PREDICTION_MODES, "No possible encoding found" );
  CHECK( bestCS->cost             == MAX_DOUBLE                , "No possible encoding found" );

 // if (bestCS->cus[0]->tmpFlag)
 //   printf("\n Kalyan TMP is choosen");
 
}

#if SHARP_LUMA_DELTA_QP || ENABLE_QPA_SUB_CTU
void EncCu::updateLambda (Slice* slice, const int dQP,
 #if WCG_EXT && ER_CHROMA_QP_WCG_PPS
                          const bool useWCGChromaControl,
 #endif
                          const bool updateRdCostLambda)
{
#if WCG_EXT && ER_CHROMA_QP_WCG_PPS
  if (useWCGChromaControl)
  {
    const double lambda = m_pcSliceEncoder->initializeLambda (slice, m_pcSliceEncoder->getGopId(), slice->getSliceQp(), (double)dQP);
    const int clippedQP = Clip3 (-slice->getSPS()->getQpBDOffset (CHANNEL_TYPE_LUMA), MAX_QP, dQP);

    m_pcSliceEncoder->setUpLambda (slice, lambda, clippedQP);
    return;
  }
#endif
  int iQP = dQP;
  const double oldQP     = (double)slice->getSliceQpBase();
#if ENABLE_QPA_SUB_CTU
  const double oldLambda = (m_pcEncCfg->getUsePerceptQPA() && !m_pcEncCfg->getUseRateCtrl() && slice->getPPS()->getUseDQP()) ? slice->getLambdas()[0] :
                           m_pcSliceEncoder->calculateLambda (slice, m_pcSliceEncoder->getGopId(), oldQP, oldQP, iQP);
#else
  const double oldLambda = m_pcSliceEncoder->calculateLambda (slice, m_pcSliceEncoder->getGopId(), oldQP, oldQP, iQP);
#endif
  const double newLambda = oldLambda * pow (2.0, ((double)dQP - oldQP) / 3.0);
#if RDOQ_CHROMA_LAMBDA
  const double lambdaArray[MAX_NUM_COMPONENT] = {newLambda / m_pcRdCost->getDistortionWeight (COMPONENT_Y),
                                                 newLambda / m_pcRdCost->getDistortionWeight (COMPONENT_Cb),
                                                 newLambda / m_pcRdCost->getDistortionWeight (COMPONENT_Cr)};
  m_pcTrQuant->setLambdas (lambdaArray);
#else
  m_pcTrQuant->setLambda (newLambda);
#endif
  if (updateRdCostLambda)
  {
    m_pcRdCost->setLambda (newLambda, slice->getSPS()->getBitDepths());
#if WCG_EXT
    if (!m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled())
    {
      m_pcRdCost->saveUnadjustedLambda();
    }
#endif
  }
}
#endif // SHARP_LUMA_DELTA_QP || ENABLE_QPA_SUB_CTU

#if ENABLE_SPLIT_PARALLELISM
//#undef DEBUG_PARALLEL_TIMINGS
//#define DEBUG_PARALLEL_TIMINGS 1
void EncCu::xCompressCUParallel( CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &partitioner )
{
  const unsigned wIdx = gp_sizeIdxInfo->idxFrom( partitioner.currArea().lwidth() );
  const unsigned hIdx = gp_sizeIdxInfo->idxFrom( partitioner.currArea().lheight() );

  Picture* picture = tempCS->picture;

  int numJobs = m_modeCtrl->getNumParallelJobs( *bestCS, partitioner );

  bool    jobUsed                            [NUM_RESERVERD_SPLIT_JOBS];
  std::fill( jobUsed, jobUsed + NUM_RESERVERD_SPLIT_JOBS, false );

  const UnitArea currArea = CS::getArea( *tempCS, partitioner.currArea(), partitioner.chType );
  const bool doParallel   = !m_pcEncCfg->getForceSingleSplitThread();
  omp_set_num_threads( m_pcEncCfg->getNumSplitThreads() );

#pragma omp parallel for schedule(dynamic,1) if(doParallel)
  for( int jId = 1; jId <= numJobs; jId++ )
  {
    // thread start
    picture->scheduler.setSplitThreadId();
    picture->scheduler.setSplitJobId( jId );

    QTBTPartitioner jobPartitioner;
    EncCu*       jobCuEnc       = m_pcEncLib->getCuEncoder( picture->scheduler.getSplitDataId( jId ) );
    auto*        jobBlkCache    = dynamic_cast<CacheBlkInfoCtrl*>( jobCuEnc->m_modeCtrl );
#if REUSE_CU_RESULTS
    auto*        jobBestCache   = dynamic_cast<BestEncInfoCache*>( jobCuEnc->m_modeCtrl );
#endif

    jobPartitioner.copyState( partitioner );
    jobCuEnc      ->copyState( this, jobPartitioner, currArea, true );

    if( jobBlkCache  ) { jobBlkCache ->tick(); }
#if REUSE_CU_RESULTS
    if( jobBestCache ) { jobBestCache->tick(); }

#endif
    CodingStructure *&jobBest = jobCuEnc->m_pBestCS[wIdx][hIdx];
    CodingStructure *&jobTemp = jobCuEnc->m_pTempCS[wIdx][hIdx];

    jobUsed[jId] = true;

    jobCuEnc->xCompressCU( jobTemp, jobBest, jobPartitioner );

    picture->scheduler.setSplitJobId( 0 );
    // thread stop
  }
  picture->scheduler.setSplitThreadId( 0 );

  int    bestJId  = 0;
  double bestCost = bestCS->cost;
  for( int jId = 1; jId <= numJobs; jId++ )
  {
    EncCu* jobCuEnc = m_pcEncLib->getCuEncoder( picture->scheduler.getSplitDataId( jId ) );

    if( jobUsed[jId] && jobCuEnc->m_pBestCS[wIdx][hIdx]->cost < bestCost )
    {
      bestCost = jobCuEnc->m_pBestCS[wIdx][hIdx]->cost;
      bestJId  = jId;
    }
  }

  if( bestJId > 0 )
  {
    copyState( m_pcEncLib->getCuEncoder( picture->scheduler.getSplitDataId( bestJId ) ), partitioner, currArea, false );
    m_CurrCtx->best = m_CABACEstimator->getCtx();

    tempCS = m_pTempCS[wIdx][hIdx];
    bestCS = m_pBestCS[wIdx][hIdx];
  }

  const int      bitDepthY = tempCS->sps->getBitDepth( CH_L );
  const UnitArea clipdArea = clipArea( currArea, *picture );

  CHECK( calcCheckSum( picture->getRecoBuf( clipdArea.Y() ), bitDepthY ) != calcCheckSum( bestCS->getRecoBuf( clipdArea.Y() ), bitDepthY ), "Data copied incorrectly!" );

  picture->finishParallelPart( currArea );

  if( auto *blkCache = dynamic_cast<CacheBlkInfoCtrl*>( m_modeCtrl ) )
  {
    for( int jId = 1; jId <= numJobs; jId++ )
    {
      if( !jobUsed[jId] || jId == bestJId ) continue;

      auto *jobBlkCache = dynamic_cast<CacheBlkInfoCtrl*>( m_pcEncLib->getCuEncoder( picture->scheduler.getSplitDataId( jId ) )->m_modeCtrl );
      CHECK( !jobBlkCache, "If own mode controller has blk info cache capability so should all other mode controllers!" );
      blkCache->CacheBlkInfoCtrl::copyState( *jobBlkCache, partitioner.currArea() );
    }

    blkCache->tick();
  }
#if REUSE_CU_RESULTS

  if( auto *blkCache = dynamic_cast<BestEncInfoCache*>( m_modeCtrl ) )
  {
    for( int jId = 1; jId <= numJobs; jId++ )
    {
      if( !jobUsed[jId] || jId == bestJId ) continue;

      auto *jobBlkCache = dynamic_cast<BestEncInfoCache*>( m_pcEncLib->getCuEncoder( picture->scheduler.getSplitDataId( jId ) )->m_modeCtrl );
      CHECK( !jobBlkCache, "If own mode controller has blk info cache capability so should all other mode controllers!" );
      blkCache->BestEncInfoCache::copyState( *jobBlkCache, partitioner.currArea() );
    }

    blkCache->tick();
  }
#endif
}

void EncCu::copyState( EncCu* other, Partitioner& partitioner, const UnitArea& currArea, const bool isDist )
{
  const unsigned wIdx = gp_sizeIdxInfo->idxFrom( partitioner.currArea().lwidth () );
  const unsigned hIdx = gp_sizeIdxInfo->idxFrom( partitioner.currArea().lheight() );

  if( isDist )
  {
    other->m_pBestCS[wIdx][hIdx]->initSubStructure( *m_pBestCS[wIdx][hIdx], partitioner.chType, partitioner.currArea(), false );
    other->m_pTempCS[wIdx][hIdx]->initSubStructure( *m_pTempCS[wIdx][hIdx], partitioner.chType, partitioner.currArea(), false );
  }
  else
  {
          CodingStructure* dst =        m_pBestCS[wIdx][hIdx];
    const CodingStructure* src = other->m_pBestCS[wIdx][hIdx];
    bool keepResi = KEEP_PRED_AND_RESI_SIGNALS;
    bool keepPred = true;

    dst->useSubStructure( *src, partitioner.chType, currArea, keepPred, true, keepResi, keepResi, true );

    dst->cost           =  src->cost;
    dst->dist           =  src->dist;
    dst->fracBits       =  src->fracBits;
    dst->features       =  src->features;
  }

  if( isDist )
  {
    m_CurrCtx = m_CtxBuffer.data();
  }

  m_pcInterSearch->copyState( *other->m_pcInterSearch );
  m_modeCtrl     ->copyState( *other->m_modeCtrl, partitioner.currArea() );
  m_pcRdCost     ->copyState( *other->m_pcRdCost );
  m_pcTrQuant    ->copyState( *other->m_pcTrQuant );
  if( m_pcEncCfg->getLmcs() )
  {
    EncReshape *encReshapeThis  = dynamic_cast<EncReshape*>(       m_pcReshape);
    EncReshape *encReshapeOther = dynamic_cast<EncReshape*>(other->m_pcReshape);
    encReshapeThis->copyState( *encReshapeOther );
  }

  m_CABACEstimator->getCtx() = other->m_CABACEstimator->getCtx();
}
#endif
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
void EncCu::xCheckModeSplit(CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &partitioner, const EncTestMode& encTestMode )
#else
void EncCu::xCheckModeSplit(CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &partitioner, const EncTestMode& encTestMode, const ModeType modeTypeParent, bool &skipInterPass )
#endif
{
  const int qp                = encTestMode.qp;
  const Slice &slice          = *tempCS->slice;
  const int oldPrevQp         = tempCS->prevQP[partitioner.chType];
  const auto oldMotionLut     = tempCS->motionLut;
#if ENABLE_QPA_SUB_CTU
  const PPS &pps              = *tempCS->pps;
  const uint32_t currDepth    = partitioner.currDepth;
#endif
  const auto oldPLT           = tempCS->prevPLT;

  const PartSplit split = getPartSplit( encTestMode );
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  const ModeType modeTypeChild = partitioner.modeType;
#endif
  CHECK( split == CU_DONT_SPLIT, "No proper split provided!" );

  tempCS->initStructData( qp );

  m_CABACEstimator->getCtx() = m_CurrCtx->start;

  const TempCtx ctxStartSP( m_CtxCache, SubCtx( Ctx::SplitFlag,   m_CABACEstimator->getCtx() ) );
  const TempCtx ctxStartQt( m_CtxCache, SubCtx( Ctx::SplitQtFlag, m_CABACEstimator->getCtx() ) );
  const TempCtx ctxStartHv( m_CtxCache, SubCtx( Ctx::SplitHvFlag, m_CABACEstimator->getCtx() ) );
  const TempCtx ctxStart12( m_CtxCache, SubCtx( Ctx::Split12Flag, m_CABACEstimator->getCtx() ) );
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  const TempCtx ctxStartMC( m_CtxCache, SubCtx( Ctx::ModeConsFlag, m_CABACEstimator->getCtx() ) );
#endif
  m_CABACEstimator->resetBits();

  m_CABACEstimator->split_cu_mode( split, *tempCS, partitioner );
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  m_CABACEstimator->mode_constraint( split, *tempCS, partitioner, modeTypeChild );
#endif
  const double factor = ( tempCS->currQP[partitioner.chType] > 30 ? 1.1 : 1.075 );
  tempCS->useDbCost = m_pcEncCfg->getUseEncDbOpt();
  if (!tempCS->useDbCost)
    CHECK(bestCS->costDbOffset != 0, "error");
  const double cost   = m_pcRdCost->calcRdCost( uint64_t( m_CABACEstimator->getEstFracBits() + ( ( bestCS->fracBits ) / factor ) ), Distortion( bestCS->dist / factor ) ) + bestCS->costDbOffset / factor;

  m_CABACEstimator->getCtx() = SubCtx( Ctx::SplitFlag,   ctxStartSP );
  m_CABACEstimator->getCtx() = SubCtx( Ctx::SplitQtFlag, ctxStartQt );
  m_CABACEstimator->getCtx() = SubCtx( Ctx::SplitHvFlag, ctxStartHv );
  m_CABACEstimator->getCtx() = SubCtx( Ctx::Split12Flag, ctxStart12 );
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  m_CABACEstimator->getCtx() = SubCtx( Ctx::ModeConsFlag, ctxStartMC );
#endif
  if (cost > bestCS->cost + bestCS->costDbOffset
#if ENABLE_QPA_SUB_CTU
    || (m_pcEncCfg->getUsePerceptQPA() && !m_pcEncCfg->getUseRateCtrl() && pps.getUseDQP() && (slice.getCuQpDeltaSubdiv() > 0) && (split == CU_HORZ_SPLIT || split == CU_VERT_SPLIT) &&
        (currDepth == 0)) // force quad-split or no split at CTU level
#endif
    )
  {
    xCheckBestMode( tempCS, bestCS, partitioner, encTestMode );
    return;
  }
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  const bool chromaNotSplit = modeTypeParent == MODE_TYPE_ALL && modeTypeChild == MODE_TYPE_INTRA ? true : false;
  if( partitioner.treeType != TREE_D )
  {
    tempCS->treeType = TREE_L;
  }
  else
  {
    if( chromaNotSplit )
    {
      CHECK( partitioner.chType != CHANNEL_TYPE_LUMA, "chType must be luma" );
      tempCS->treeType = partitioner.treeType = TREE_L;
    }
    else
    {
      tempCS->treeType = partitioner.treeType = TREE_D;
    }
  }
#endif

  partitioner.splitCurrArea( split, *tempCS );
  bool qgEnableChildren = partitioner.currQgEnable(); // QG possible at children level

  m_CurrCtx++;

  tempCS->getRecoBuf().fill( 0 );

  tempCS->getPredBuf().fill(0);
  AffineMVInfo tmpMVInfo;
  bool isAffMVInfoSaved;
  m_pcInterSearch->savePrevAffMVInfo(0, tmpMVInfo, isAffMVInfoSaved);
  BlkUniMvInfo tmpUniMvInfo;
  bool         isUniMvInfoSaved = false;
  if (!tempCS->slice->isIntra())
  {
    m_pcInterSearch->savePrevUniMvInfo(tempCS->area.Y(), tmpUniMvInfo, isUniMvInfoSaved);
  }
#if INTER_LIC
  BlkUniMvInfo tmpUniMvInfoLIC;
  bool         isUniMvInfoSavedLIC = false;
  if (tempCS->slice->getUseLIC() && !tempCS->slice->isIntra())
  {
    m_pcInterSearch->swapUniMvBuffer();
    m_pcInterSearch->savePrevUniMvInfo(tempCS->area.Y(), tmpUniMvInfoLIC, isUniMvInfoSavedLIC);
    m_pcInterSearch->swapUniMvBuffer();
  }
#endif

  do
  {
    const auto &subCUArea  = partitioner.currArea();

    if( tempCS->picture->Y().contains( subCUArea.lumaPos() ) )
    {
      const unsigned wIdx    = gp_sizeIdxInfo->idxFrom( subCUArea.lwidth () );
      const unsigned hIdx    = gp_sizeIdxInfo->idxFrom( subCUArea.lheight() );

      CodingStructure *tempSubCS = m_pTempCS[wIdx][hIdx];
      CodingStructure *bestSubCS = m_pBestCS[wIdx][hIdx];

      tempCS->initSubStructure( *tempSubCS, partitioner.chType, subCUArea, false );
      tempCS->initSubStructure( *bestSubCS, partitioner.chType, subCUArea, false );
      tempSubCS->bestParent = bestSubCS->bestParent = bestCS;
      double newMaxCostAllowed = isLuma(partitioner.chType) ? std::min(encTestMode.maxCostAllowed, bestCS->cost - m_pcRdCost->calcRdCost(tempCS->fracBits, tempCS->dist)) : MAX_DOUBLE;
      newMaxCostAllowed = std::max(0.0, newMaxCostAllowed);
      xCompressCU(tempSubCS, bestSubCS, partitioner, newMaxCostAllowed);
      tempSubCS->bestParent = bestSubCS->bestParent = nullptr;

      if( bestSubCS->cost == MAX_DOUBLE )
      {
        CHECK( split == CU_QUAD_SPLIT, "Split decision reusing cannot skip quad split" );
        tempCS->cost = MAX_DOUBLE;
        tempCS->costDbOffset = 0;
        tempCS->useDbCost = m_pcEncCfg->getUseEncDbOpt();
        m_CurrCtx--;
        partitioner.exitCurrSplit();
        xCheckBestMode( tempCS, bestCS, partitioner, encTestMode );
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
        if( partitioner.chType == CHANNEL_TYPE_LUMA )
        {
          tempCS->motionLut = oldMotionLut;
        }
#endif
        return;
      }

      bool keepResi = KEEP_PRED_AND_RESI_SIGNALS;
      tempCS->useSubStructure( *bestSubCS, partitioner.chType, CS::getArea( *tempCS, subCUArea, partitioner.chType ), KEEP_PRED_AND_RESI_SIGNALS, true, keepResi, keepResi, true );

      if( partitioner.currQgEnable() )
      {
        tempCS->prevQP[partitioner.chType] = bestSubCS->prevQP[partitioner.chType];
      }
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
      if( partitioner.isConsInter() )
      {
        for( int i = 0; i < bestSubCS->cus.size(); i++ )
        {
          CHECK( bestSubCS->cus[i]->predMode != MODE_INTER, "all CUs must be inter mode in an Inter coding region (SCIPU)" );
        }
      }
      else if( partitioner.isConsIntra() )
      {
        for( int i = 0; i < bestSubCS->cus.size(); i++ )
        {
          CHECK( bestSubCS->cus[i]->predMode == MODE_INTER, "all CUs must not be inter mode in an Intra coding region (SCIPU)" );
        }
      }
#endif
      tempSubCS->releaseIntermediateData();
      bestSubCS->releaseIntermediateData();
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
      if( !tempCS->slice->isIntra() && partitioner.isConsIntra() )
      {
        tempCS->cost = m_pcRdCost->calcRdCost( tempCS->fracBits, tempCS->dist );
        if( tempCS->cost > bestCS->cost )
        {
          tempCS->cost = MAX_DOUBLE;
          tempCS->costDbOffset = 0;
          tempCS->useDbCost = m_pcEncCfg->getUseEncDbOpt();
          m_CurrCtx--;
          partitioner.exitCurrSplit();
          if( partitioner.chType == CHANNEL_TYPE_LUMA )
          {
            tempCS->motionLut = oldMotionLut;
          }
          return;
        }
      }
#endif
    }
  } while( partitioner.nextPart( *tempCS ) );

  partitioner.exitCurrSplit();


  m_CurrCtx--;
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  if( chromaNotSplit )
  {
    //Note: In local dual tree region, the chroma CU refers to the central luma CU's QP.
    //If the luma CU QP shall be predQP (no residual in it and before it in the QG), it must be revised to predQP before encoding the chroma CU
    //Otherwise, the chroma CU uses predQP+deltaQP in encoding but is decoded as using predQP, thus causing encoder-decoded mismatch on chroma qp.
    if( tempCS->pps->getUseDQP() )
    {
      //find parent CS that including all coded CUs in the QG before this node
      CodingStructure* qgCS = tempCS;
      bool deltaQpCodedBeforeThisNode = false;
      if( partitioner.currArea().lumaPos() != partitioner.currQgPos )
      {
        int numParentNodeToQgCS = 0;
        while( qgCS->area.lumaPos() != partitioner.currQgPos )
        {
          CHECK( qgCS->parent == nullptr, "parent of qgCS shall exsit" );
          qgCS = qgCS->parent;
          numParentNodeToQgCS++;
        }

        //check whether deltaQP has been coded (in luma CU or luma&chroma CU) before this node
        CodingStructure* parentCS = tempCS->parent;
        for( int i = 0; i < numParentNodeToQgCS; i++ )
        {
          //checking each parent
          CHECK( parentCS == nullptr, "parentCS shall exsit" );
          for( const auto &cu : parentCS->cus )
          {
            if( cu->rootCbf && !isChroma( cu->chType ) )
            {
              deltaQpCodedBeforeThisNode = true;
              break;
            }
          }
          parentCS = parentCS->parent;
        }
      }

      //revise luma CU qp before the first luma CU with residual in the SCIPU to predQP
      if( !deltaQpCodedBeforeThisNode )
      {
        //get pred QP of the QG
        const CodingUnit* cuFirst = qgCS->getCU( CHANNEL_TYPE_LUMA );
        CHECK( cuFirst->lumaPos() != partitioner.currQgPos, "First cu of the Qg is wrong" );
        int predQp = CU::predictQP( *cuFirst, qgCS->prevQP[CHANNEL_TYPE_LUMA] );

        //revise to predQP
        int firstCuHasResidual = (int)tempCS->cus.size();
        for( int i = 0; i < tempCS->cus.size(); i++ )
        {
          if( tempCS->cus[i]->rootCbf )
          {
            firstCuHasResidual = i;
            break;
          }
        }

        for( int i = 0; i < firstCuHasResidual; i++ )
        {
          tempCS->cus[i]->qp = predQp;
        }
      }
    }
    assert( tempCS->treeType == TREE_L );
    uint32_t numCuPuTu[6];
    tempCS->picture->cs->getNumCuPuTuOffset( numCuPuTu );
    tempCS->picture->cs->useSubStructure( *tempCS, partitioner.chType, CS::getArea( *tempCS, partitioner.currArea(), partitioner.chType ), false, true, false, false, false );

    if (isChromaEnabled(tempCS->pcv->chrFormat))
    {
    partitioner.chType = CHANNEL_TYPE_CHROMA;
    tempCS->treeType = partitioner.treeType = TREE_C;

    m_CurrCtx++;

    const unsigned wIdx = gp_sizeIdxInfo->idxFrom( partitioner.currArea().lwidth() );
    const unsigned hIdx = gp_sizeIdxInfo->idxFrom( partitioner.currArea().lheight() );
    CodingStructure *tempCSChroma = m_pTempCS2[wIdx][hIdx];
    CodingStructure *bestCSChroma = m_pBestCS2[wIdx][hIdx];
    tempCS->initSubStructure( *tempCSChroma, partitioner.chType, partitioner.currArea(), false );
    tempCS->initSubStructure( *bestCSChroma, partitioner.chType, partitioner.currArea(), false );
    tempCS->treeType = TREE_D;
    xCompressCU( tempCSChroma, bestCSChroma, partitioner );

    //attach chromaCS to luma CS and update cost
    bool keepResi = KEEP_PRED_AND_RESI_SIGNALS;
    //bestCSChroma->treeType = tempCSChroma->treeType = TREE_C;
    CHECK( bestCSChroma->treeType != TREE_C || tempCSChroma->treeType != TREE_C, "wrong treeType for chroma CS" );
    tempCS->useSubStructure( *bestCSChroma, partitioner.chType, CS::getArea( *bestCSChroma, partitioner.currArea(), partitioner.chType ), KEEP_PRED_AND_RESI_SIGNALS, true, keepResi, true, true );

    //release tmp resource
    tempCSChroma->releaseIntermediateData();
    bestCSChroma->releaseIntermediateData();
    //tempCS->picture->cs->releaseIntermediateData();
      m_CurrCtx--;
    }
    tempCS->picture->cs->clearCuPuTuIdxMap( partitioner.currArea(), numCuPuTu[0], numCuPuTu[1], numCuPuTu[2], numCuPuTu + 3 );


    //recover luma tree status
    partitioner.chType = CHANNEL_TYPE_LUMA;
    partitioner.treeType = TREE_D;
    partitioner.modeType = MODE_TYPE_ALL;
  }
#endif
  // Finally, generate split-signaling bits for RD-cost check
  const PartSplit implicitSplit = partitioner.getImplicitSplit( *tempCS );

  {
    bool enforceQT = implicitSplit == CU_QUAD_SPLIT;

    // LARGE CTU bug
    if( m_pcEncCfg->getUseFastLCTU() )
    {
      unsigned minDepth = 0;
      unsigned maxDepth = floorLog2(tempCS->sps->getCTUSize()) - floorLog2(tempCS->sps->getMinQTSize(slice.getSliceType(), partitioner.chType));

      if( auto ad = dynamic_cast<AdaptiveDepthPartitioner*>( &partitioner ) )
      {
        ad->setMaxMinDepth( minDepth, maxDepth, *tempCS );
      }

      if( minDepth > partitioner.currQtDepth )
      {
        // enforce QT
        enforceQT = true;
      }
    }

    if( !enforceQT )
    {
      m_CABACEstimator->resetBits();

      m_CABACEstimator->split_cu_mode( split, *tempCS, partitioner );
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
      partitioner.modeType = modeTypeParent;
      m_CABACEstimator->mode_constraint( split, *tempCS, partitioner, modeTypeChild );
#endif
      tempCS->fracBits += m_CABACEstimator->getEstFracBits(); // split bits
    }
  }

  tempCS->cost = m_pcRdCost->calcRdCost( tempCS->fracBits, tempCS->dist );

  // Check Delta QP bits for splitted structure
  if( !qgEnableChildren ) // check at deepest QG level only
  xCheckDQP( *tempCS, partitioner, true );

  // If the configuration being tested exceeds the maximum number of bytes for a slice / slice-segment, then
  // a proper RD evaluation cannot be performed. Therefore, termination of the
  // slice/slice-segment must be made prior to this CTU.
  // This can be achieved by forcing the decision to be that of the rpcTempCU.
  // The exception is each slice / slice-segment must have at least one CTU.
  if (bestCS->cost != MAX_DOUBLE)
  {
  }
  else
  {
    bestCS->costDbOffset = 0;
  }
  tempCS->useDbCost = m_pcEncCfg->getUseEncDbOpt();
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  if( tempCS->cus.size() > 0 && modeTypeParent == MODE_TYPE_ALL && modeTypeChild == MODE_TYPE_INTER )
  {
    int areaSizeNoResiCu = 0;
    for( int k = 0; k < tempCS->cus.size(); k++ )
    {
      areaSizeNoResiCu += (tempCS->cus[k]->rootCbf == false) ? tempCS->cus[k]->lumaSize().area() : 0;
    }
    if( areaSizeNoResiCu >= (tempCS->area.lumaSize().area() >> 1) )
    {
      skipInterPass = true;
    }
  }
#endif
  // RD check for sub partitioned coding structure.
  xCheckBestMode( tempCS, bestCS, partitioner, encTestMode );

  if (isAffMVInfoSaved)
    m_pcInterSearch->addAffMVInfo(tmpMVInfo);
  if (!tempCS->slice->isIntra() && isUniMvInfoSaved)
  {
    m_pcInterSearch->addUniMvInfo(tmpUniMvInfo);
  }
#if INTER_LIC
  if (!tempCS->slice->isIntra() && isUniMvInfoSavedLIC)
  {
    m_pcInterSearch->swapUniMvBuffer();
    m_pcInterSearch->addUniMvInfo(tmpUniMvInfoLIC);
    m_pcInterSearch->swapUniMvBuffer();
  }
#endif

  tempCS->motionLut = oldMotionLut;

  tempCS->prevPLT   = oldPLT;

  tempCS->releaseIntermediateData();

  tempCS->prevQP[partitioner.chType] = oldPrevQp;
}

bool EncCu::xCheckRDCostIntra(CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &partitioner, const EncTestMode& encTestMode, bool adaptiveColorTrans)
{
  double          bestInterCost             = m_modeCtrl->getBestInterCost();
  double          costSize2Nx2NmtsFirstPass = m_modeCtrl->getMtsSize2Nx2NFirstPassCost();
  bool            skipSecondMtsPass         = m_modeCtrl->getSkipSecondMTSPass();
  const SPS&      sps                       = *tempCS->sps;
  const int       maxSizeMTS                = MTS_INTRA_MAX_CU_SIZE;
  uint8_t         considerMtsSecondPass     = ( sps.getUseIntraMTS() && isLuma( partitioner.chType ) && partitioner.currArea().lwidth() <= maxSizeMTS && partitioner.currArea().lheight() <= maxSizeMTS ) ? 1 : 0;

  bool   useIntraSubPartitions   = false;
  double maxCostAllowedForChroma = MAX_DOUBLE;
  const  CodingUnit *bestCU      = bestCS->getCU( partitioner.chType );
  Distortion interHad = m_modeCtrl->getInterHad();
#if JVET_W0123_TIMD_FUSION
  int timdMode = 0;
  int timdModeSecondary = 0;
  bool timdIsBlended = false;
  int  timdFusionWeight[2] = { 0 };
#endif


  double dct2Cost                =   MAX_DOUBLE;
  double bestNonDCT2Cost         = MAX_DOUBLE;
  double trGrpBestCost     [ 4 ] = { MAX_DOUBLE, MAX_DOUBLE, MAX_DOUBLE, MAX_DOUBLE };
  double globalBestCost          =   MAX_DOUBLE;
  bool   bestSelFlag       [ 4 ] = { false, false, false, false };
  bool   trGrpCheck        [ 4 ] = { true, true, true, true };
  int    startMTSIdx       [ 4 ] = { 0, 1, 2, 3 };
  int    endMTSIdx         [ 4 ] = { 0, 1, 2, 3 };
#if JVET_W0103_INTRA_MTS
  endMTSIdx[0] = 3; //put all MTS candidates in "Grp 0"
#endif
  double trGrpStopThreshold[ 3 ] = { 1.001, 1.001, 1.001 };
  int    bestMtsFlag             =   0;
  int    bestLfnstIdx            =   0;

#if EXTENDED_LFNST || JVET_W0119_LFNST_EXTENSION
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  const int  maxLfnstIdx         = (CS::isDualITree(*tempCS) && partitioner.chType == CHANNEL_TYPE_CHROMA && (partitioner.currArea().lwidth() < 8 || partitioner.currArea().lheight() < 8))
                                   || (partitioner.currArea().lwidth() > sps.getMaxTbSize() || partitioner.currArea().lheight() > sps.getMaxTbSize()) ? 0 : 3;
#else
  const int  maxLfnstIdx         = ( partitioner.isSepTree( *tempCS ) && partitioner.chType == CHANNEL_TYPE_CHROMA && ( partitioner.currArea().lwidth() < 8 || partitioner.currArea().lheight() < 8 ) )
                                   || ( partitioner.currArea().lwidth() > sps.getMaxTbSize() || partitioner.currArea().lheight() > sps.getMaxTbSize() ) ? 0 : 3;
#endif
#else
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  const int  maxLfnstIdx         = (CS::isDualITree(*tempCS) && partitioner.chType == CHANNEL_TYPE_CHROMA && (partitioner.currArea().lwidth() < 8 || partitioner.currArea().lheight() < 8))
                                   || (partitioner.currArea().lwidth() > sps.getMaxTbSize() || partitioner.currArea().lheight() > sps.getMaxTbSize()) ? 0 : 2;
#else
  const int  maxLfnstIdx         = ( partitioner.isSepTree( *tempCS ) && partitioner.chType == CHANNEL_TYPE_CHROMA && ( partitioner.currArea().lwidth() < 8 || partitioner.currArea().lheight() < 8 ) )
                                   || ( partitioner.currArea().lwidth() > sps.getMaxTbSize() || partitioner.currArea().lheight() > sps.getMaxTbSize() ) ? 0 : 2;
#endif
#endif

  bool       skipOtherLfnst      = false;
  int        startLfnstIdx       = 0;
  int        endLfnstIdx         = sps.getUseLFNST() ? maxLfnstIdx : 0;
#if JVET_W0103_INTRA_MTS
  int grpNumMax = 1;
#else
  int grpNumMax = sps.getUseLFNST() ? m_pcEncCfg->getMTSIntraMaxCand() : 1;
#endif
  m_modeCtrl->setISPWasTested(false);
  m_pcIntraSearch->invalidateBestModeCost();
  if (sps.getUseColorTrans() && !CS::isDualITree(*tempCS))
  {
    if ((m_pcEncCfg->getRGBFormatFlag() && adaptiveColorTrans) || (!m_pcEncCfg->getRGBFormatFlag() && !adaptiveColorTrans))
    {
      m_pcIntraSearch->invalidateBestRdModeFirstColorSpace();
    }
  }

  bool foundZeroRootCbf = false;
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  if (sps.getUseColorTrans())
  {
    CHECK(tempCS->treeType != TREE_D || partitioner.treeType != TREE_D, "localtree should not be applied when adaptive color transform is enabled");
    CHECK(tempCS->modeType != MODE_TYPE_ALL || partitioner.modeType != MODE_TYPE_ALL, "localtree should not be applied when adaptive color transform is enabled");
    CHECK(adaptiveColorTrans && (CS::isDualITree(*tempCS) || partitioner.chType != CHANNEL_TYPE_LUMA), "adaptive color transform cannot be applied to dual-tree");
  }
#endif
#if ENABLE_DIMD
  bool dimd_is_blend = false;
  int  dimdMode = 0;
  int  dimdBlendMode[2] = { 0 };
  int  dimdRelWeight[3] = { 0 };
  bool dimdDerived = false;

  if (isLuma(partitioner.chType))
  {
    CodingUnit cu(tempCS->area);
    cu.cs = tempCS;
    cu.slice = tempCS->slice;
    cu.tileIdx = tempCS->pps->getTileIdx(tempCS->area.lumaPos());
    PredictionUnit pu(tempCS->area);
    pu.cu = &cu;
    pu.cs = tempCS;

    if( cu.slice->getSPS()->getUseDimd() )
    {
      const CompArea &area = cu.Y();
      IntraPrediction::deriveDimdMode(bestCS->picture->getRecoBuf(area), area, cu);

      dimdDerived = true;
      dimd_is_blend = cu.dimd_is_blend;
      dimdMode = cu.dimdMode;
      dimdBlendMode[0] = cu.dimdBlendMode[0];
      dimdBlendMode[1] = cu.dimdBlendMode[1];
      dimdRelWeight[0] = cu.dimdRelWeight[0];
      dimdRelWeight[1] = cu.dimdRelWeight[1];
      dimdRelWeight[2] = cu.dimdRelWeight[2];
    }

#if SECONDARY_MPM
    m_pcIntraSearch->getMpmListSize() = PU::getIntraMPMs(pu, m_pcIntraSearch->getMPMList(), m_pcIntraSearch->getNonMPMList());
#endif
  }
#elif SECONDARY_MPM
  if( isLuma( partitioner.chType ) )
  {
    CodingUnit cu( tempCS->area );
    cu.cs = tempCS;
    cu.slice = tempCS->slice;
    cu.tileIdx = tempCS->pps->getTileIdx( tempCS->area.lumaPos() );
    PredictionUnit pu( tempCS->area );
    pu.cu = &cu;
    pu.cs = tempCS;

    m_pcIntraSearch->getMpmListSize() = PU::getIntraMPMs( pu, m_pcIntraSearch->getMPMList(), m_pcIntraSearch->getNonMPMList() );
  }
#endif

#if JVET_W0123_TIMD_FUSION
  bool timdDerived = false;
#endif
  for( int trGrpIdx = 0; trGrpIdx < grpNumMax; trGrpIdx++ )
  {
    const uint8_t startMtsFlag = trGrpIdx > 0;
    const uint8_t endMtsFlag   = sps.getUseLFNST() ? considerMtsSecondPass : 0;

    if( ( trGrpIdx == 0 || ( !skipSecondMtsPass && considerMtsSecondPass ) ) && trGrpCheck[ trGrpIdx ] )
    {
      for( int lfnstIdx = startLfnstIdx; lfnstIdx <= endLfnstIdx; lfnstIdx++ )
      {
        for( uint8_t mtsFlag = startMtsFlag; mtsFlag <= endMtsFlag; mtsFlag++ )
        {
          if (sps.getUseColorTrans() && !CS::isDualITree(*tempCS))
          {
            m_pcIntraSearch->setSavedRdModeIdx(trGrpIdx*(NUM_LFNST_NUM_PER_SET * 2) + lfnstIdx * 2 + mtsFlag);
          }
          if (mtsFlag > 0 && lfnstIdx > 0)
          {
            continue;
          }
          //3) if interHad is 0, only try further modes if some intra mode was already better than inter
          if( sps.getUseLFNST() && m_pcEncCfg->getUsePbIntraFast() && !tempCS->slice->isIntra() && bestCU && CU::isInter( *bestCS->getCU( partitioner.chType ) ) && interHad == 0 )
          {
            continue;
          }

          tempCS->initStructData( encTestMode.qp );

          CodingUnit &cu      = tempCS->addCU( CS::getArea( *tempCS, tempCS->area, partitioner.chType ), partitioner.chType );

          partitioner.setCUData( cu );
          cu.slice            = tempCS->slice;
          cu.tileIdx          = tempCS->pps->getTileIdx( tempCS->area.lumaPos() );
          cu.skip             = false;
          cu.mmvdSkip = false;
          cu.predMode         = MODE_INTRA;
          cu.chromaQpAdj      = m_cuChromaQpOffsetIdxPlus1;
          cu.qp               = encTestMode.qp;
#if ENABLE_DIMD
          cu.dimd = false;
          if( dimdDerived )
          {
            cu.dimd_is_blend = dimd_is_blend;
            cu.dimdMode = dimdMode;
            cu.dimdBlendMode[0] = dimdBlendMode[0];
            cu.dimdBlendMode[1] = dimdBlendMode[1];
            cu.dimdRelWeight[0] = dimdRelWeight[0];
            cu.dimdRelWeight[1] = dimdRelWeight[1];
            cu.dimdRelWeight[2] = dimdRelWeight[2];
          }
#endif
          cu.lfnstIdx         = lfnstIdx;
          cu.mtsFlag          = mtsFlag;
          cu.ispMode          = NOT_INTRA_SUBPARTITIONS;
          cu.colorTransform = adaptiveColorTrans;

          CU::addPUs( cu );
#if JVET_W0123_TIMD_FUSION
          cu.timd = false;
          if (isLuma(partitioner.chType) && cu.slice->getSPS()->getUseTimd())
          {
            if (cu.lwidth() * cu.lheight() > 1024 && cu.slice->getSliceType() == I_SLICE)
            {
              timdDerived = true;
            }
            if (!timdDerived)
            {
              const CompArea &area = cu.Y();
              cu.timdMode = m_pcIntraSearch->deriveTimdMode(bestCS->picture->getRecoBuf(area), area, cu);
              timdMode = cu.timdMode;
              timdDerived = true;
              timdModeSecondary = cu.timdModeSecondary;
              timdIsBlended     = cu.timdIsBlended;
              timdFusionWeight[0] = cu.timdFusionWeight[0];
              timdFusionWeight[1] = cu.timdFusionWeight[1];
            }
            else
            {
              cu.timdMode = timdMode;
              cu.timdModeSecondary = timdModeSecondary;
              cu.timdIsBlended     = timdIsBlended;
              cu.timdFusionWeight[0] = timdFusionWeight[0];
              cu.timdFusionWeight[1] = timdFusionWeight[1];
            }
          }
#endif

          tempCS->interHad    = interHad;

          m_bestModeUpdated = tempCS->useDbCost = bestCS->useDbCost = false;

          bool validCandRet = false;
          if( isLuma( partitioner.chType ) )
          {
            //ISP uses the value of the best cost so far (luma if it is the fast version) to avoid test non-necessary subpartitions
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
            double bestCostSoFar = CS::isDualITree(*tempCS) ? m_modeCtrl->getBestCostWithoutSplitFlags() : bestCU && bestCU->predMode == MODE_INTRA ? bestCS->lumaCost : bestCS->cost;
            if (CS::isDualITree(*tempCS) && encTestMode.maxCostAllowed < bestCostSoFar)
#else
            double bestCostSoFar = partitioner.isSepTree(*tempCS) ? m_modeCtrl->getBestCostWithoutSplitFlags() : bestCU && bestCU->predMode == MODE_INTRA ? bestCS->lumaCost : bestCS->cost;
            if (partitioner.isSepTree(*tempCS) && encTestMode.maxCostAllowed < bestCostSoFar)
#endif
            {
              bestCostSoFar = encTestMode.maxCostAllowed;
            }
            validCandRet = m_pcIntraSearch->estIntraPredLumaQT(cu, partitioner, bestCostSoFar, mtsFlag, startMTSIdx[trGrpIdx], endMTSIdx[trGrpIdx], (trGrpIdx > 0), !cu.colorTransform ? bestCS : nullptr);
            if ((!validCandRet || (cu.ispMode && cu.firstTU->cbf[COMPONENT_Y] == 0)))
            {
              continue;
            }
#if JVET_W0123_TIMD_FUSION
            PU::spanIpmInfoIntra(*cu.firstPU);
#endif
#if JVET_W0123_TIMD_FUSION
            if (m_pcEncCfg->getUseFastISP() && validCandRet && !mtsFlag && !lfnstIdx && !cu.colorTransform && !cu.timd)
#else
            if (m_pcEncCfg->getUseFastISP() && validCandRet && !mtsFlag && !lfnstIdx && !cu.colorTransform)
#endif
            {
              m_modeCtrl->setISPMode(cu.ispMode);
              m_modeCtrl->setISPLfnstIdx(cu.lfnstIdx);
              m_modeCtrl->setMIPFlagISPPass(cu.mipFlag);
#if JVET_V0130_INTRA_TMP
			        m_modeCtrl->setTPMFlagISPPass(cu.tmpFlag);
#endif
              m_modeCtrl->setBestISPIntraModeRelCU(cu.ispMode ? PU::getFinalIntraMode(*cu.firstPU, CHANNEL_TYPE_LUMA) : UINT8_MAX);
              m_modeCtrl->setBestDCT2NonISPCostRelCU(m_modeCtrl->getMtsFirstPassNoIspCost());
            }

            if (sps.getUseColorTrans() && m_pcEncCfg->getRGBFormatFlag() && !CS::isDualITree(*tempCS) && !cu.colorTransform)
            {
              double curLumaCost = m_pcRdCost->calcRdCost(tempCS->fracBits, tempCS->dist);
              if (curLumaCost > bestCS->cost)
              {
                continue;
              }
            }

            useIntraSubPartitions = cu.ispMode != NOT_INTRA_SUBPARTITIONS;
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
            if (!CS::isDualITree(*tempCS))
#else
            if( !partitioner.isSepTree( *tempCS ) )
#endif
            {
              tempCS->lumaCost = m_pcRdCost->calcRdCost( tempCS->fracBits, tempCS->dist );
              if( useIntraSubPartitions )
              {
                //the difference between the best cost so far and the current luma cost is stored to avoid testing the Cr component if the cost of luma + Cb is larger than the best cost
                maxCostAllowedForChroma = bestCS->cost < MAX_DOUBLE ? bestCS->cost - tempCS->lumaCost : MAX_DOUBLE;
              }
            }

            if (m_pcEncCfg->getUsePbIntraFast() && tempCS->dist == std::numeric_limits<Distortion>::max()
                && tempCS->interHad == 0)
            {
              interHad = 0;
              // JEM assumes only perfect reconstructions can from now on beat the inter mode
              m_modeCtrl->enforceInterHad( 0 );
              continue;
            }
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
            if (!CS::isDualITree(*tempCS))
#else
            if( !partitioner.isSepTree( *tempCS ) )
#endif
            {
              if (!cu.colorTransform)
              {
                cu.cs->picture->getRecoBuf(cu.Y()).copyFrom(cu.cs->getRecoBuf(COMPONENT_Y));
                cu.cs->picture->getPredBuf(cu.Y()).copyFrom(cu.cs->getPredBuf(COMPONENT_Y));
              }
              else
              {
                cu.cs->picture->getRecoBuf(cu).copyFrom(cu.cs->getRecoBuf(cu));
                cu.cs->picture->getPredBuf(cu).copyFrom(cu.cs->getPredBuf(cu));
              }
            }
          }
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
          if (tempCS->area.chromaFormat != CHROMA_400 && (partitioner.chType == CHANNEL_TYPE_CHROMA || !CS::isDualITree(*tempCS)) && !cu.colorTransform)
#else
          if( tempCS->area.chromaFormat != CHROMA_400 && ( partitioner.chType == CHANNEL_TYPE_CHROMA || !cu.isSepTree() ) && !cu.colorTransform )
#endif
          {
            TUIntraSubPartitioner subTuPartitioner( partitioner );
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
            m_pcIntraSearch->estIntraPredChromaQT(cu, (!useIntraSubPartitions || (CS::isDualITree(*tempCS) && !isLuma(CHANNEL_TYPE_CHROMA))) ? partitioner : subTuPartitioner, maxCostAllowedForChroma);
#else
            m_pcIntraSearch->estIntraPredChromaQT( cu, ( !useIntraSubPartitions || ( cu.isSepTree() && !isLuma( CHANNEL_TYPE_CHROMA ) ) ) ? partitioner : subTuPartitioner, maxCostAllowedForChroma );
#endif
            if( useIntraSubPartitions && !cu.ispMode )
            {
              //At this point the temp cost is larger than the best cost. Therefore, we can already skip the remaining calculations
              continue;
            }
          }

          cu.rootCbf = false;

          for( uint32_t t = 0; t < getNumberValidTBlocks( *cu.cs->pcv ); t++ )
          {
            cu.rootCbf |= cu.firstTU->cbf[t] != 0;
          }

          if (!cu.rootCbf)
          {
            cu.colorTransform = false;
            foundZeroRootCbf = true;
          }

          // Get total bits for current mode: encode CU
          m_CABACEstimator->resetBits();

          if ((!cu.cs->slice->isIntra() || cu.cs->slice->getSPS()->getIBCFlag())
            && cu.Y().valid()
            )
          {
            m_CABACEstimator->cu_skip_flag ( cu );
          }
          m_CABACEstimator->pred_mode      ( cu );
#if ENABLE_DIMD
          m_CABACEstimator->cu_dimd_flag   ( cu );
#endif
          m_CABACEstimator->adaptive_color_transform(cu);
          m_CABACEstimator->cu_pred_data   ( cu );

          // Encode Coefficients
          CUCtx cuCtx;
          cuCtx.isDQPCoded = true;
          cuCtx.isChromaQpAdjCoded = true;
          m_CABACEstimator->cu_residual( cu, partitioner, cuCtx );

          tempCS->fracBits = m_CABACEstimator->getEstFracBits();
          tempCS->cost     = m_pcRdCost->calcRdCost(tempCS->fracBits, tempCS->dist);


          double tmpCostWithoutSplitFlags = tempCS->cost;
          xEncodeDontSplit( *tempCS, partitioner );

          xCheckDQP( *tempCS, partitioner );
          xCheckChromaQPOffset( *tempCS, partitioner );

          // Check if low frequency non-separable transform (LFNST) is too expensive
          if( lfnstIdx && !cuCtx.lfnstLastScanPos && !cu.ispMode )
          {
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
            bool cbfAtZeroDepth = CS::isDualITree(*tempCS) ?
#else
            bool cbfAtZeroDepth = cu.isSepTree() ?
#endif
                                       cu.rootCbf
                                     : (tempCS->area.chromaFormat != CHROMA_400 && std::min( cu.firstTU->blocks[ 1 ].width, cu.firstTU->blocks[ 1 ].height ) < 4) ?
                                            TU::getCbfAtDepth( *cu.firstTU, COMPONENT_Y, 0 )
                                          : cu.rootCbf;
            if( cbfAtZeroDepth )
            {
              tempCS->cost = MAX_DOUBLE;
              tmpCostWithoutSplitFlags = MAX_DOUBLE;
            }
          }

          if (isLuma(partitioner.chType) && cu.firstTU->mtsIdx[COMPONENT_Y] > MTS_SKIP)
          {
            CHECK(!cuCtx.mtsLastScanPos, "MTS is disallowed to only contain DC coefficient");
          }

          if( mtsFlag == 0 && lfnstIdx == 0 )
          {
            dct2Cost = tempCS->cost;
          }
          else if (tmpCostWithoutSplitFlags < bestNonDCT2Cost)
          {
            bestNonDCT2Cost = tmpCostWithoutSplitFlags;
          }

          if( tempCS->cost < bestCS->cost )
          {
            m_modeCtrl->setBestCostWithoutSplitFlags( tmpCostWithoutSplitFlags );
          }

          if( !mtsFlag ) static_cast< double& >( costSize2Nx2NmtsFirstPass ) = tempCS->cost;

          if( sps.getUseLFNST() && !tempCS->cus.empty() )
          {
            skipOtherLfnst = m_modeCtrl->checkSkipOtherLfnst( encTestMode, tempCS, partitioner );
          }

          xCalDebCost( *tempCS, partitioner );
          tempCS->useDbCost = m_pcEncCfg->getUseEncDbOpt();

#if WCG_EXT
          DTRACE_MODE_COST( *tempCS, m_pcRdCost->getLambda( true ) );
#else
          DTRACE_MODE_COST( *tempCS, m_pcRdCost->getLambda() );
#endif
          if (sps.getUseColorTrans() && !CS::isDualITree(*tempCS))
          {
            int colorSpaceIdx = ((m_pcEncCfg->getRGBFormatFlag() && adaptiveColorTrans) || (!m_pcEncCfg->getRGBFormatFlag() && !adaptiveColorTrans)) ? 0 : 1;
            if (tempCS->cost < tempCS->tmpColorSpaceIntraCost[colorSpaceIdx])
            {
              tempCS->tmpColorSpaceIntraCost[colorSpaceIdx] = tempCS->cost;
              bestCS->tmpColorSpaceIntraCost[colorSpaceIdx] = tempCS->cost;
            }
          }
         
          if( !sps.getUseLFNST() )
          {
            xCheckBestMode( tempCS, bestCS, partitioner, encTestMode );
          }
          else
          {
            if( xCheckBestMode( tempCS, bestCS, partitioner, encTestMode ) )
            {
              trGrpBestCost[ trGrpIdx ] = globalBestCost = bestCS->cost;
              bestSelFlag  [ trGrpIdx ] = true;
              bestMtsFlag               = mtsFlag;
              bestLfnstIdx              = lfnstIdx;

              if( bestCS->cus.size() == 1 )
              {
                CodingUnit &cu = *bestCS->cus.front();
                if (cu.firstTU->mtsIdx[COMPONENT_Y] == MTS_SKIP)
                {
                  if( ( floorLog2( cu.firstTU->blocks[ COMPONENT_Y ].width ) + floorLog2( cu.firstTU->blocks[ COMPONENT_Y ].height ) ) >= 6 )
                  {
                    endLfnstIdx = 0;
                  }
                }
              }
            }
            
            //we decide to skip the non-DCT-II transforms and LFNST according to the ISP results
            if ((endMtsFlag > 0 || endLfnstIdx > 0) && (cu.ispMode || (bestCS && bestCS->cus[0]->ispMode)) && tempCS->slice->isIntra() && m_pcEncCfg->getUseFastISP())
            {
              double bestCostDct2NoIsp = m_modeCtrl->getMtsFirstPassNoIspCost();
              double bestIspCost       = m_modeCtrl->getIspCost();
              //Kalyan
              //CHECKD( bestCostDct2NoIsp <= bestIspCost, "wrong cost!" );
              double threshold = 1.4;
              
              double lfnstThreshold = 1.01 * threshold;
              if( m_modeCtrl->getStopNonDCT2Transforms() || bestCostDct2NoIsp > bestIspCost*lfnstThreshold )
              {
                endLfnstIdx = lfnstIdx;
              }

              if ( m_modeCtrl->getStopNonDCT2Transforms() || bestCostDct2NoIsp > bestIspCost*threshold )
              {
                skipSecondMtsPass = true;
                m_modeCtrl->setSkipSecondMTSPass( true );
                break;
              }
            }
            //now we check whether the second pass of SIZE_2Nx2N and the whole Intra SIZE_NxN should be skipped or not
            if( !mtsFlag && !tempCS->slice->isIntra() && bestCU && bestCU->predMode != MODE_INTRA )
            {
              const double thEmtInterFastSkipIntra = 1.4; // Skip checking Intra if "2Nx2N using DCT2" is worse than best Inter mode
              if( costSize2Nx2NmtsFirstPass > thEmtInterFastSkipIntra * bestInterCost )
              {
                skipSecondMtsPass = true;
                m_modeCtrl->setSkipSecondMTSPass( true );
                break;
              }
            }
#if JVET_W0103_INTRA_MTS
            if (lfnstIdx && m_modeCtrl->getMtsFirstPassNoIspCost() != MAX_DOUBLE && isLuma(partitioner.chType))
            {
              double threshold = 1.5;
              if (m_modeCtrl->getMtsFirstPassNoIspCost() > threshold * bestCS->cost)
              {
                endLfnstIdx = lfnstIdx;
              }
            }
#endif
          }

        } //for emtCuFlag
        if( skipOtherLfnst )
        {
          startLfnstIdx = lfnstIdx;
          endLfnstIdx   = lfnstIdx;
          break;
        }
      } //for lfnstIdx
    } //if (!skipSecondMtsPass && considerMtsSecondPass && trGrpCheck[iGrpIdx])

    if( sps.getUseLFNST() && trGrpIdx < 3 )
    {
      trGrpCheck[ trGrpIdx + 1 ] = false;

      if( bestSelFlag[ trGrpIdx ] && considerMtsSecondPass )
      {
        double dCostRatio = dct2Cost / trGrpBestCost[ trGrpIdx ];
        trGrpCheck[ trGrpIdx + 1 ] = ( bestMtsFlag != 0 || bestLfnstIdx != 0 ) && dCostRatio < trGrpStopThreshold[ trGrpIdx ];
      }
    }
  } //trGrpIdx
  if(!adaptiveColorTrans)
  m_modeCtrl->setBestNonDCT2Cost(bestNonDCT2Cost);
  return foundZeroRootCbf;
}


void EncCu::xCheckPLT(CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &partitioner, const EncTestMode& encTestMode)
{
  if (((partitioner.currArea().lumaSize().width * partitioner.currArea().lumaSize().height <= 16) && (isLuma(partitioner.chType)) )
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
        || ((partitioner.currArea().chromaSize().width * partitioner.currArea().chromaSize().height <= 16) && (!isLuma(partitioner.chType)) && CS::isDualITree(*tempCS)) )
#else
        || ((partitioner.currArea().chromaSize().width * partitioner.currArea().chromaSize().height <= 16) && (!isLuma(partitioner.chType)) && partitioner.isSepTree(*tempCS) )
        || (partitioner.isLocalSepTree(*tempCS)  && (!isLuma(partitioner.chType))  )  )
#endif
  {
    return;
  }
  tempCS->initStructData(encTestMode.qp);
  CodingUnit &cu = tempCS->addCU(CS::getArea(*tempCS, tempCS->area, partitioner.chType), partitioner.chType);
  partitioner.setCUData(cu);
  cu.slice = tempCS->slice;
  cu.tileIdx = tempCS->pps->getTileIdx(tempCS->area.lumaPos());
  cu.skip = false;
  cu.mmvdSkip = false;
  cu.predMode = MODE_PLT;

  cu.chromaQpAdj = m_cuChromaQpOffsetIdxPlus1;
  cu.qp = encTestMode.qp;
  cu.bdpcmMode = 0;

  tempCS->addPU(CS::getArea(*tempCS, tempCS->area, partitioner.chType), partitioner.chType);
  tempCS->addTU(CS::getArea(*tempCS, tempCS->area, partitioner.chType), partitioner.chType);
  // Search
  tempCS->dist = 0;
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  if (CS::isDualITree(*tempCS))
#else
  if (cu.isSepTree())
#endif
  {
    if (isLuma(partitioner.chType))
    {
      m_pcIntraSearch->PLTSearch(*tempCS, partitioner, COMPONENT_Y, 1);
    }
    if (tempCS->area.chromaFormat != CHROMA_400 && (partitioner.chType == CHANNEL_TYPE_CHROMA))
    {
      m_pcIntraSearch->PLTSearch(*tempCS, partitioner, COMPONENT_Cb, 2);
    }
  }
  else
  {
    if( cu.chromaFormat != CHROMA_400 )
    {
      m_pcIntraSearch->PLTSearch(*tempCS, partitioner, COMPONENT_Y, 3);
    }
    else
    {
      m_pcIntraSearch->PLTSearch(*tempCS, partitioner, COMPONENT_Y, 1);
    }
  }


  m_CABACEstimator->getCtx() = m_CurrCtx->start;
  m_CABACEstimator->resetBits();
  if ((!cu.cs->slice->isIntra() || cu.cs->slice->getSPS()->getIBCFlag())
    && cu.Y().valid())
  {
    m_CABACEstimator->cu_skip_flag(cu);
  }
  m_CABACEstimator->pred_mode(cu);

  // signaling
  CUCtx cuCtx;
  cuCtx.isDQPCoded = true;
  cuCtx.isChromaQpAdjCoded = true;
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  if (CS::isDualITree(*tempCS))
#else
  if (cu.isSepTree())
#endif
  {
    if (isLuma(partitioner.chType))
    {
      m_CABACEstimator->cu_palette_info(cu, COMPONENT_Y, 1, cuCtx);
    }
    if (tempCS->area.chromaFormat != CHROMA_400 && (partitioner.chType == CHANNEL_TYPE_CHROMA))
    {
      m_CABACEstimator->cu_palette_info(cu, COMPONENT_Cb, 2, cuCtx);
    }
  }
  else
  {
    if( cu.chromaFormat != CHROMA_400 )
    {
      m_CABACEstimator->cu_palette_info(cu, COMPONENT_Y, 3, cuCtx);
    }
    else
    {
      m_CABACEstimator->cu_palette_info(cu, COMPONENT_Y, 1, cuCtx);
    }
  }
  tempCS->fracBits = m_CABACEstimator->getEstFracBits();
  tempCS->cost = m_pcRdCost->calcRdCost(tempCS->fracBits, tempCS->dist);

  xEncodeDontSplit(*tempCS, partitioner);
  xCheckDQP(*tempCS, partitioner);
  xCheckChromaQPOffset( *tempCS, partitioner );
  xCalDebCost(*tempCS, partitioner);
  tempCS->useDbCost = m_pcEncCfg->getUseEncDbOpt();

  const Area currCuArea = cu.block(getFirstComponentOfChannel(partitioner.chType));
  cu.slice->m_mapPltCost[isChroma(partitioner.chType)][currCuArea.pos()][currCuArea.size()] = tempCS->cost;
#if WCG_EXT
  DTRACE_MODE_COST(*tempCS, m_pcRdCost->getLambda(true));
#else
  DTRACE_MODE_COST(*tempCS, m_pcRdCost->getLambda());
#endif
  xCheckBestMode(tempCS, bestCS, partitioner, encTestMode);
}

void EncCu::xCheckDQP( CodingStructure& cs, Partitioner& partitioner, bool bKeepCtx )
{
  CHECK( bKeepCtx && cs.cus.size() <= 1 && partitioner.getImplicitSplit( cs ) == CU_DONT_SPLIT, "bKeepCtx should only be set in split case" );
  CHECK( !bKeepCtx && cs.cus.size() > 1, "bKeepCtx should never be set for non-split case" );

  if( !cs.pps->getUseDQP() )
  {
    return;
  }

#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  if (CS::isDualITree(cs) && isChroma(partitioner.chType))
#else
  if (partitioner.isSepTree(cs) && isChroma(partitioner.chType))
#endif
  {
    return;
  }

  if( !partitioner.currQgEnable() ) // do not consider split or leaf/not leaf QG condition (checked by caller)
  {
    return;
  }


  CodingUnit* cuFirst = cs.getCU( partitioner.chType );

  CHECK( !cuFirst, "No CU available" );

  bool hasResidual = false;
  for( const auto &cu : cs.cus )
  {
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
    if ( cu->rootCbf )
#else
    //not include the chroma CU because chroma CU is decided based on corresponding luma QP and deltaQP is not signaled at chroma CU
    if( cu->rootCbf && !isChroma( cu->chType ))
#endif
    {
      hasResidual = true;
      break;
    }
  }

  int predQP = CU::predictQP( *cuFirst, cs.prevQP[partitioner.chType] );

  if( hasResidual )
  {
    TempCtx ctxTemp( m_CtxCache );
    if( !bKeepCtx ) ctxTemp = SubCtx( Ctx::DeltaQP, m_CABACEstimator->getCtx() );

    m_CABACEstimator->resetBits();
    m_CABACEstimator->cu_qp_delta( *cuFirst, predQP, cuFirst->qp );

    cs.fracBits += m_CABACEstimator->getEstFracBits(); // dQP bits
    cs.cost      = m_pcRdCost->calcRdCost(cs.fracBits, cs.dist);


    if( !bKeepCtx ) m_CABACEstimator->getCtx() = SubCtx( Ctx::DeltaQP, ctxTemp );

    // NOTE: reset QPs for CUs without residuals up to first coded CU
    for( const auto &cu : cs.cus )
    {
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
      if ( cu->rootCbf )
#else
      //not include the chroma CU because chroma CU is decided based on corresponding luma QP and deltaQP is not signaled at chroma CU
      if( cu->rootCbf && !isChroma( cu->chType ))
#endif
      {
        break;
      }
      cu->qp = predQP;
    }
  }
  else
  {
    // No residuals: reset CU QP to predicted value
    for( const auto &cu : cs.cus )
    {
      cu->qp = predQP;
    }
  }
}

void EncCu::xCheckChromaQPOffset( CodingStructure& cs, Partitioner& partitioner )
{
  // doesn't apply if CU chroma QP offset is disabled
  if( !cs.slice->getUseChromaQpAdj() )
  {
    return;
  }

  // doesn't apply to luma CUs
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  if (CS::isDualITree(cs) && isLuma(partitioner.chType))
#else
  if( partitioner.isSepTree(cs) && isLuma(partitioner.chType) )
#endif
  {
    return;
  }

  // not needed after the first coded TU in the chroma QG
  if( !partitioner.currQgChromaEnable() )
  {
    return;
  }

  CodingUnit& cu = *cs.getCU( partitioner.chType );

  // check if chroma is coded or not
  bool hasResidual = false;
  for( const TransformUnit &tu : CU::traverseTUs(cu) )
  {
    if( tu.cbf[COMPONENT_Cb] || tu.cbf[COMPONENT_Cr] )
    {
      hasResidual = true;
      break;
    }
  }

  if( hasResidual )
  {
    // estimate cost for coding cu_chroma_qp_offset
    TempCtx ctxTempAdjFlag( m_CtxCache );
    TempCtx ctxTempAdjIdc( m_CtxCache );
    ctxTempAdjFlag = SubCtx( Ctx::ChromaQpAdjFlag, m_CABACEstimator->getCtx() );
    ctxTempAdjIdc = SubCtx( Ctx::ChromaQpAdjIdc,   m_CABACEstimator->getCtx() );
    m_CABACEstimator->resetBits();
    m_CABACEstimator->cu_chroma_qp_offset( cu );
    cs.fracBits += m_CABACEstimator->getEstFracBits();
    cs.cost      = m_pcRdCost->calcRdCost(cs.fracBits, cs.dist);
    m_CABACEstimator->getCtx() = SubCtx( Ctx::ChromaQpAdjFlag, ctxTempAdjFlag );
    m_CABACEstimator->getCtx() = SubCtx( Ctx::ChromaQpAdjIdc,  ctxTempAdjIdc  );
  }
  else
  {
    // reset chroma QP offset to 0 if it will not be coded
    cu.chromaQpAdj = 0;
  }
}

#if!REMOVE_PCM
void EncCu::xFillPCMBuffer( CodingUnit &cu )
{
  const ChromaFormat format        = cu.chromaFormat;
  const uint32_t numberValidComponents = getNumberValidComponents(format);

  for( auto &tu : CU::traverseTUs( cu ) )
  {
    for( uint32_t ch = 0; ch < numberValidComponents; ch++ )
    {
      const ComponentID compID = ComponentID( ch );
      const CompArea &compArea = tu.blocks[ compID ];

      if( tu.cs->slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() && compID == COMPONENT_Y )
      {
        tu.getPcmbuf( compID ).rspSignal( tu.cs->getOrgBuf( compArea ), m_pcReshape->getFwdLUT() );
      }
      else
      {
        tu.getPcmbuf( compID ).copyFrom( tu.cs->getOrgBuf( compArea ) );
      }
    }
  }
}
#endif

void EncCu::xCheckRDCostHashInter( CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &partitioner, const EncTestMode& encTestMode )
{
#if ENABLE_OBMC
  double bestOBMCCost = MAX_DOUBLE;
#endif
  bool isPerfectMatch = false;

  tempCS->initStructData(encTestMode.qp);
  m_pcInterSearch->resetBufferedUniMotions();
  m_pcInterSearch->setAffineModeSelected(false);
  CodingUnit &cu = tempCS->addCU(tempCS->area, partitioner.chType);

  partitioner.setCUData(cu);
  cu.slice = tempCS->slice;
  cu.tileIdx = tempCS->pps->getTileIdx(tempCS->area.lumaPos());
  cu.skip = false;
  cu.predMode = MODE_INTER;
  cu.chromaQpAdj = m_cuChromaQpOffsetIdxPlus1;
  cu.qp = encTestMode.qp;
#if INTER_LIC
  cu.LICFlag = false;
#endif
  CU::addPUs(cu);
  cu.mmvdSkip = false;
  cu.firstPU->mmvdMergeFlag = false;

  if (m_pcInterSearch->predInterHashSearch(cu, partitioner, isPerfectMatch))
  {
    double equBcwCost = MAX_DOUBLE;

    m_bestModeUpdated = tempCS->useDbCost = bestCS->useDbCost = false;
#if ENABLE_OBMC //normal inter
  const unsigned wIdx = gp_sizeIdxInfo->idxFrom(partitioner.currArea().lwidth());
  CodingStructure *prevCS = tempCS;
  PelUnitBuf tempWoOBMCBuf = m_tempWoOBMCBuffer.subBuf(UnitAreaRelative(cu, cu));
  tempWoOBMCBuf.copyFrom(tempCS->getPredBuf(cu));
  cu.isobmcMC = true;
  cu.obmcFlag = true;
  m_pcInterSearch->subBlockOBMC(*cu.firstPU);
  cu.isobmcMC = false;
#endif
    xEncodeInterResidual(tempCS, bestCS, partitioner, encTestMode, 0
      , 0
      , &equBcwCost
    );
#if ENABLE_OBMC // xCheckRDCostInter
  double tempCost = (prevCS == tempCS) ? tempCS->cost : bestCS->cost;
  if (m_pTempCUWoOBMC && tempCost < bestOBMCCost)
  {
    const unsigned hIdx = gp_sizeIdxInfo->idxFrom(prevCS->area.lheight());
    m_pTempCUWoOBMC[wIdx][hIdx]->clearCUs();
    m_pTempCUWoOBMC[wIdx][hIdx]->clearPUs();
    m_pTempCUWoOBMC[wIdx][hIdx]->clearTUs();
    m_pTempCUWoOBMC[wIdx][hIdx]->copyStructure(*prevCS, partitioner.chType);

    m_pPredBufWoOBMC[wIdx][hIdx].copyFrom(tempWoOBMCBuf);
    m_pTempCUWoOBMC[wIdx][hIdx]->getPredBuf(cu).copyFrom(prevCS->getPredBuf(cu));

    bestOBMCCost = tempCost;
  }
#endif
    if ( m_bestModeUpdated && bestCS->cost != MAX_DOUBLE )
    {
      xCalDebCost( *bestCS, partitioner );
    }
  }
  tempCS->initStructData(encTestMode.qp);
  int minSize = min(cu.lwidth(), cu.lheight());
  if (minSize < 64)
  {
    isPerfectMatch = false;
  }
  m_modeCtrl->setIsHashPerfectMatch(isPerfectMatch);
}

void EncCu::xCheckRDCostMerge2Nx2N( CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &partitioner, const EncTestMode& encTestMode )
{
  const Slice &slice = *tempCS->slice;

  CHECK( slice.getSliceType() == I_SLICE, "Merge modes not available for I-slices" );

  tempCS->initStructData( encTestMode.qp );

  MergeCtx mergeCtx;
#if JVET_W0090_ARMC_TM
  MergeCtx mergeCtxtmp;
#endif
  const SPS &sps = *tempCS->sps;

#if MERGE_ENC_OPT
  const bool affineMrgAvail = (sps.getUseAffine() || sps.getSbTMVPEnabledFlag()) && slice.getPicHeader()->getMaxNumAffineMergeCand()
    && !(bestCS->area.lumaSize().width < 8 || bestCS->area.lumaSize().height < 8);

  AffineMergeCtx affineMergeCtx;
#if JVET_W0090_ARMC_TM
  AffineMergeCtx affineMergeCtxTmp;
#endif
  MergeCtx mrgCtx;
#if TM_MRG
  MergeCtx tmMrgCtx;
#endif

  if (sps.getSbTMVPEnabledFlag())
  {
    Size bufSize = g_miScaling.scale(tempCS->area.lumaSize());
    mergeCtx.subPuMvpMiBuf = MotionBuf( m_SubPuMiBuf, bufSize );
    mrgCtx.subPuMvpMiBuf = MotionBuf(m_SubPuMiBuf, bufSize);
    affineMergeCtx.mrgCtx = &mrgCtx;
#if TM_MRG
    tmMrgCtx.subPuMvpMiBuf = MotionBuf(m_SubPuMiBuf, bufSize);
#endif
  }
#endif

#if MULTI_PASS_DMVR
  bool applyBDMVR[MRG_MAX_NUM_CANDS] = { false };
#if TM_MRG && MERGE_ENC_OPT
  bool applyBDMVR4TM[TM_MRG_MAX_NUM_CANDS] = { false };
#endif
#endif
#if !MULTI_PASS_DMVR
  Mv   refinedMvdL0[MAX_NUM_PARTS_IN_CTU][MRG_MAX_NUM_CANDS];
#endif
  setMergeBestSATDCost( MAX_DOUBLE );

  {
    // first get merge candidates
    CodingUnit cu( tempCS->area );
    cu.cs       = tempCS;
    cu.predMode = MODE_INTER;
    cu.slice    = tempCS->slice;
    cu.tileIdx  = tempCS->pps->getTileIdx( tempCS->area.lumaPos() );
#if INTER_LIC
    cu.LICFlag  = false;
#endif

    PredictionUnit pu( tempCS->area );
    pu.cu = &cu;
    pu.cs = tempCS;
#if TM_MRG
    pu.tmMergeFlag = false;
#endif
    PU::getInterMergeCandidates(pu, mergeCtx
      , 0
    );
#if JVET_W0097_GPM_MMVD_TM
    m_mergeCand.copyMergeCtx(mergeCtx);
    m_mergeCandAvail = true;
#endif
    PU::getInterMMVDMergeCandidates(pu, mergeCtx);
#if JVET_W0090_ARMC_TM
    mergeCtxtmp = mergeCtx;
    if (sps.getUseAML())
    {
      m_pcInterSearch->adjustInterMergeCandidates(pu, mergeCtx);
    }
#endif
#if TM_MRG && MERGE_ENC_OPT
    if (cu.cs->sps->getUseDMVDMode())
    {
      cu.firstPU = &pu;
      pu.tmMergeFlag = true;
      PU::getInterMergeCandidates(pu, tmMrgCtx, 0);
#if JVET_W0090_ARMC_TM
      if (sps.getUseAML())
      {
        m_pcInterSearch->adjustInterMergeCandidates(pu, tmMrgCtx);
      }
#endif

      for( uint32_t uiMergeCand = 0; uiMergeCand < tmMrgCtx.numValidMergeCand; uiMergeCand++ )
      {
        tmMrgCtx.setMergeInfo( pu, uiMergeCand );
#if MULTI_PASS_DMVR
        applyBDMVR4TM[uiMergeCand] = PU::checkBDMVRCondition(pu);
        if (applyBDMVR4TM[uiMergeCand])
        {
          pu.bdmvrRefine = true;
          m_pcInterSearch->setBdmvrSubPuMvBuf(m_mvBufBDMVR4TM[uiMergeCand << 1], m_mvBufBDMVR4TM[(uiMergeCand << 1) + 1]);
          applyBDMVR4TM[uiMergeCand] =  m_pcInterSearch->processBDMVR(pu);
        }
        else
        {
          m_pcInterSearch->deriveTMMv(pu);
        }
#else
        m_pcInterSearch->deriveTMMv( pu );
#endif

        // Store refined motion back to tmMrgCtx
        tmMrgCtx.interDirNeighbours[uiMergeCand] = pu.interDir;
        tmMrgCtx.BcwIdx[uiMergeCand] = pu.cu->BcwIdx;  // Bcw may change, because bi may be reduced to uni by deriveTMMv(pu)
        tmMrgCtx.mvFieldNeighbours[2 * uiMergeCand    ].setMvField( pu.mv[0], pu.refIdx[0] );
        tmMrgCtx.mvFieldNeighbours[2 * uiMergeCand + 1].setMvField( pu.mv[1], pu.refIdx[1] );
        if( pu.interDir == 1 )
        {
          tmMrgCtx.mvFieldNeighbours[2 * uiMergeCand + 1].setMvField( Mv(), NOT_VALID );
        }
        if( pu.interDir == 2 )
        {
          tmMrgCtx.mvFieldNeighbours[2 * uiMergeCand    ].setMvField( Mv(), NOT_VALID );
        }
      }

      pu.tmMergeFlag = false;
#if MULTI_PASS_DMVR
      pu.bdmvrRefine = false;
#endif
    }
#endif
#if MERGE_ENC_OPT
    if (affineMrgAvail)
    {
      pu.regularMergeFlag = false;
      cu.affine = true;
      PU::getAffineMergeCand(pu, affineMergeCtx);
#if JVET_W0090_ARMC_TM
      affineMergeCtxTmp = affineMergeCtx;
      if (sps.getUseAML())
      {
        m_pcInterSearch->adjustAffineMergeCandidates(pu, affineMergeCtx);
      }
#endif
      cu.affine = false;
    }
#endif
    pu.regularMergeFlag = true;

#if MULTI_PASS_DMVR
    if (cu.cs->sps->getUseDMVDMode())
    {
      cu.firstPU = &pu;
      for (uint32_t uiMergeCand = 0; uiMergeCand < mergeCtx.numValidMergeCand; uiMergeCand++)
      {
        if( mergeCtx.interDirNeighbours[uiMergeCand] == 3 )
        {
          mergeCtx.setMergeInfo( pu, uiMergeCand );
          applyBDMVR[uiMergeCand] = PU::checkBDMVRCondition(pu);

          if (applyBDMVR[uiMergeCand])
          {
            pu.bdmvrRefine = true;
            m_pcInterSearch->setBdmvrSubPuMvBuf(m_mvBufBDMVR[uiMergeCand << 1], m_mvBufBDMVR[(uiMergeCand << 1) + 1]);

            if (mergeCtx.xCheckSimilarMotion(pu.mergeIdx, PU::getBDMVRMvdThreshold(pu)))
            {
              // span motion to subPU
              for (int subPuIdx = 0; subPuIdx < MAX_NUM_SUBCU_DMVR; subPuIdx++)
              {
                m_mvBufBDMVR[uiMergeCand << 1][subPuIdx] = pu.mv[0];
                m_mvBufBDMVR[(uiMergeCand << 1) + 1][subPuIdx] = pu.mv[1];
              }
            }
            else
            {
              m_pcInterSearch->processBDMVR(pu);
            }
          }
        }
      }
    }
#endif
  }
#if AFFINE_MMVD && MERGE_ENC_OPT
  int  afMmvdBaseIdxToMergeIdxOffset = (int)PU::getMergeIdxFromAfMmvdBaseIdx(affineMergeCtx, 0);
  int  afMmvdBaseCount = std::min<int>((int)AF_MMVD_BASE_NUM, affineMergeCtx.numValidMergeCand - afMmvdBaseIdxToMergeIdxOffset);
  bool affineMmvdAvail = affineMrgAvail && afMmvdBaseCount >= 1 && sps.getUseAffineMmvdMode();
#endif
  bool candHasNoResidual[MRG_MAX_NUM_CANDS + MMVD_ADD_NUM] = { false };
  bool                                        bestIsSkip = false;
  bool                                        bestIsMMVDSkip = true;
#if !MERGE_ENC_OPT
  PelUnitBuf                                  acMergeBuffer[MRG_MAX_NUM_CANDS];
#endif
  PelUnitBuf                                  acMergeTmpBuffer[MRG_MAX_NUM_CANDS];
  PelUnitBuf                                  acMergeRealBuffer[MMVD_MRG_MAX_RD_BUF_NUM];
  PelUnitBuf *                                acMergeTempBuffer[MMVD_MRG_MAX_RD_NUM];
  PelUnitBuf *                                singleMergeTempBuffer;
#if !MERGE_ENC_OPT
  int                                         insertPos;
#endif
  unsigned                                    uiNumMrgSATDCand = mergeCtx.numValidMergeCand + MMVD_ADD_NUM;

#if !MERGE_ENC_OPT
  struct ModeInfo
  {
    uint32_t mergeCand;
    bool     isRegularMerge;
    bool     isMMVD;
    bool     isCIIP;
#if CIIP_PDPC
    bool     isCiipPDPC;
    ModeInfo() : mergeCand(0), isRegularMerge(false), isMMVD(false), isCIIP(false), isCiipPDPC(false) {}
    ModeInfo(const uint32_t mergeCand, const bool isRegularMerge, const bool isMMVD, const bool isCIIP, const bool isCiipPDPC) :
      mergeCand(mergeCand), isRegularMerge(isRegularMerge), isMMVD(isMMVD), isCIIP(isCIIP), isCiipPDPC( isCiipPDPC ) {}
#else
    ModeInfo() : mergeCand(0), isRegularMerge(false), isMMVD(false), isCIIP(false) {}
    ModeInfo(const uint32_t mergeCand, const bool isRegularMerge, const bool isMMVD, const bool isCIIP) :
      mergeCand(mergeCand), isRegularMerge(isRegularMerge), isMMVD(isMMVD), isCIIP(isCIIP) {}
#endif
  };
#endif
  static_vector<ModeInfo, MRG_MAX_NUM_CANDS + MMVD_ADD_NUM>  RdModeList;
  bool                                        mrgTempBufSet = false;
  const int candNum = mergeCtx.numValidMergeCand + (tempCS->sps->getUseMMVD() ? std::min<int>(MMVD_BASE_MV_NUM, mergeCtx.numValidMergeCand) * MMVD_MAX_REFINE_NUM : 0);

  for (int i = 0; i < candNum; i++)
  {
    if (i < mergeCtx.numValidMergeCand)
    {
#if CIIP_PDPC
#if MERGE_ENC_OPT
      RdModeList.push_back(ModeInfo(i, true, false, false, false, false));
#else
      RdModeList.push_back(ModeInfo(i, true, false, false, false));
#endif
#else
#if MERGE_ENC_OPT
      RdModeList.push_back(ModeInfo(i, true, false, false, false));
#else
      RdModeList.push_back(ModeInfo(i, true, false, false));
#endif
#endif
    }
    else
    {
#if MERGE_ENC_OPT
#if CIIP_PDPC
      RdModeList.push_back(ModeInfo(std::min(MMVD_ADD_NUM, i - mergeCtx.numValidMergeCand), false, true, false, false, false));
#else
      RdModeList.push_back(ModeInfo(std::min(MMVD_ADD_NUM, i - mergeCtx.numValidMergeCand), false, true, false, false));
#endif
#else
#if CIIP_PDPC
      RdModeList.push_back(ModeInfo(std::min(MMVD_ADD_NUM, i - mergeCtx.numValidMergeCand), false, true, false, false));
#else
      RdModeList.push_back(ModeInfo(std::min(MMVD_ADD_NUM, i - mergeCtx.numValidMergeCand), false, true, false));
#endif
#endif
    }
  }

  const UnitArea localUnitArea(tempCS->area.chromaFormat, Area(0, 0, tempCS->area.Y().width, tempCS->area.Y().height));
  for (unsigned i = 0; i < MMVD_MRG_MAX_RD_BUF_NUM; i++)
  {
    acMergeRealBuffer[i] = m_acMergeBuffer[i].getBuf(localUnitArea);
    if (i < MMVD_MRG_MAX_RD_NUM)
    {
      acMergeTempBuffer[i] = acMergeRealBuffer + i;
    }
    else
    {
      singleMergeTempBuffer = acMergeRealBuffer + i;
    }
  }

  bool isIntrainterEnabled = sps.getUseCiip();
#if CIIP_RM_BLOCK_SIZE_CONSTRAINTS
#if CTU_256
  const int maxSize = std::min<int>( MAX_TB_SIZEY, MAX_INTRA_SIZE );

  if( bestCS->area.lwidth() * bestCS->area.lheight() < 32 || bestCS->area.lwidth() > maxSize || bestCS->area.lheight() > maxSize )
#else
  if (bestCS->area.lwidth() * bestCS->area.lheight() < 32)
#endif
#else
  if (bestCS->area.lwidth() * bestCS->area.lheight() < 64 || bestCS->area.lwidth() >= MAX_CU_SIZE || bestCS->area.lheight() >= MAX_CU_SIZE)
#endif
  {
    isIntrainterEnabled = false;
  }
  bool isTestSkipMerge[MRG_MAX_NUM_CANDS] = { false };
#if MERGE_ENC_OPT
  if (m_pcEncCfg->getUseFastMerge() || isIntrainterEnabled || affineMrgAvail
    )
#else
  if( m_pcEncCfg->getUseFastMerge() || isIntrainterEnabled)
#endif
  {
#if MERGE_ENC_OPT
    uiNumMrgSATDCand = m_pcEncCfg->getNumFullRDMerge();
#else
    uiNumMrgSATDCand = NUM_MRG_SATD_CAND;
#endif
    if (isIntrainterEnabled)
    {
      uiNumMrgSATDCand += 1;
    }
    bestIsSkip       = false;

    if( auto blkCache = dynamic_cast< CacheBlkInfoCtrl* >( m_modeCtrl ) )
    {
      if (slice.getSPS()->getIBCFlag())
      {
        ComprCUCtx cuECtx = m_modeCtrl->getComprCUCtx();
        bestIsSkip = blkCache->isSkip(tempCS->area) && cuECtx.bestCU;
      }
      else
      {
        bestIsSkip = blkCache->isSkip( tempCS->area );
      }
      bestIsMMVDSkip = blkCache->isMMVDSkip(tempCS->area);
    }

    if (isIntrainterEnabled) // always perform low complexity check
    {
      bestIsSkip = false;
    }
#if MERGE_ENC_OPT
    if (affineMrgAvail)
    {
      bestIsSkip = false;
      uiNumMrgSATDCand += NUM_AFF_MRG_SATD_CAND;
    }

#if TM_MRG
    if (tempCS->sps->getUseDMVDMode())
    {
      bestIsSkip = false;
      uiNumMrgSATDCand += TM_MAX_NUM_SATD_CAND;
    }
#endif
#endif

    static_vector<double, MRG_MAX_NUM_CANDS + MMVD_ADD_NUM> candCostList;

    // 1. Pass: get SATD-cost for selected candidates and reduce their count
    if( !bestIsSkip )
    {
      RdModeList.clear();
      mrgTempBufSet       = true;
      const TempCtx ctxStart(m_CtxCache, m_CABACEstimator->getCtx());

      CodingUnit &cu      = tempCS->addCU( tempCS->area, partitioner.chType );
#if !MERGE_ENC_OPT
      const double sqrtLambdaForFirstPassIntra = m_pcRdCost->getMotionLambda( ) * FRAC_BITS_SCALE;
#endif
      partitioner.setCUData( cu );
      cu.slice            = tempCS->slice;
      cu.tileIdx          = tempCS->pps->getTileIdx( tempCS->area.lumaPos() );
#if INTER_LIC
      cu.LICFlag          = false;
#endif
      cu.skip             = false;
      cu.mmvdSkip = false;
      cu.geoFlag          = false;
    //cu.affine
      cu.predMode         = MODE_INTER;
    //cu.LICFlag
      cu.chromaQpAdj      = m_cuChromaQpOffsetIdxPlus1;
      cu.qp               = encTestMode.qp;
    //cu.emtFlag  is set below

      PredictionUnit &pu  = tempCS->addPU( cu, partitioner.chType );
#if MERGE_ENC_OPT
      cu.affine = false;
      pu.ciipFlag = false;
#if CIIP_PDPC
      pu.ciipPDPC = false;
#endif
      pu.mmvdMergeFlag = false;
#if AFFINE_MMVD
      pu.afMmvdFlag = false;
#endif
#if TM_MRG
      pu.tmMergeFlag = false;
#endif
#if MULTI_PASS_DMVR
      pu.bdmvrRefine = false;
#endif
#endif

      DistParam distParam;
      const bool bUseHadamard = !tempCS->slice->getDisableSATDForRD();
      m_pcRdCost->setDistParam (distParam, tempCS->getOrgBuf().Y(), m_acMergeBuffer[0].Y(), sps.getBitDepth (CHANNEL_TYPE_LUMA), COMPONENT_Y, bUseHadamard);

#if MERGE_ENC_OPT
      xCheckSATDCostRegularMerge(tempCS, cu, pu, mergeCtx, acMergeTempBuffer, singleMergeTempBuffer, acMergeTmpBuffer
#if !MULTI_PASS_DMVR
                               , refinedMvdL0
#endif
                               , uiNumMrgSATDCand, RdModeList, candCostList, distParam, ctxStart
#if MULTI_PASS_DMVR
                               , applyBDMVR
#endif
      );
#else
      const UnitArea localUnitArea( tempCS->area.chromaFormat, Area( 0, 0, tempCS->area.Y().width, tempCS->area.Y().height) );
#if MULTI_HYP_PRED
      const bool testMHP = tempCS->sps->getUseInterMultiHyp()
        && (tempCS->area.lumaSize().area() > MULTI_HYP_PRED_RESTRICT_BLOCK_SIZE 
        && std::min(tempCS->area.lwidth(), tempCS->area.lheight()) >= MULTI_HYP_PRED_RESTRICT_MIN_WH);
#endif
      for( uint32_t uiMergeCand = 0; uiMergeCand < mergeCtx.numValidMergeCand; uiMergeCand++ )
      {
        mergeCtx.setMergeInfo( pu, uiMergeCand );

        PU::spanMotionInfo( pu, mergeCtx );
#if !MULTI_PASS_DMVR
        pu.mvRefine = true;
#endif
        distParam.cur = singleMergeTempBuffer->Y();
        acMergeTmpBuffer[uiMergeCand] = m_acMergeTmpBuffer[uiMergeCand].getBuf(localUnitArea);
#if INTER_LIC
        m_pcInterSearch->m_storeBeforeLIC = mergeCtx.interDirNeighbours[uiMergeCand] == 3 ? false : true;
        if (m_pcInterSearch->m_storeBeforeLIC)
        {
          m_pcInterSearch->m_predictionBeforeLIC = acMergeTmpBuffer[uiMergeCand];
          m_pcInterSearch->motionCompensation(pu, *singleMergeTempBuffer, REF_PIC_LIST_X, true, true);
        }
        else
#endif
#if MULTI_PASS_DMVR
        if (applyBDMVR[uiMergeCand])
        {
          if (isIntrainterEnabled)
          {
#if MULTI_HYP_PRED
            pu.addHypData.clear();
            pu.numMergedAddHyps = 0;
#endif
            pu.mvRefine = false;
            pu.ciipFlag = true;
            m_pcInterSearch->motionCompensation(pu, acMergeTmpBuffer[uiMergeCand]);
            pu.ciipFlag = false;
#if MULTI_HYP_PRED
            mergeCtx.setMergeInfo(pu, uiMergeCand);
#endif
          }
          pu.bdmvrRefine = true;
          m_pcInterSearch->setBdmvrSubPuMvBuf(m_mvBufBDMVR[uiMergeCand << 1], m_mvBufBDMVR[(uiMergeCand << 1) + 1]);
          pu.mvRefine = true;

          m_pcInterSearch->motionCompensation(pu, *singleMergeTempBuffer);

          if( pu.bdmvrRefine )
          {
            ::memcpy( m_mvBufEncBDOF[uiMergeCand], m_pcInterSearch->getBdofSubPuMvOffset(), sizeof( Mv ) * BDOF_SUBPU_MAX_NUM );
          }

          pu.mvRefine = false;
        }
        else
        {
#endif
        m_pcInterSearch->motionCompensation(pu, *singleMergeTempBuffer, REF_PIC_LIST_X, true, true, &(acMergeTmpBuffer[uiMergeCand]));
#if MULTI_PASS_DMVR
        }
#endif
#if INTER_LIC
        m_pcInterSearch->m_storeBeforeLIC = false;
#endif
        acMergeBuffer[uiMergeCand] = m_acRealMergeBuffer[uiMergeCand].getBuf(localUnitArea);
        acMergeBuffer[uiMergeCand].copyFrom(*singleMergeTempBuffer);
#if !MULTI_PASS_DMVR
        pu.mvRefine = false;
        if( mergeCtx.interDirNeighbours[uiMergeCand] == 3 )
        {
          mergeCtx.mvFieldNeighbours[2*uiMergeCand].mv   = pu.mv[0];
          mergeCtx.mvFieldNeighbours[2*uiMergeCand+1].mv = pu.mv[1];
          {
            int dx, dy, i, j, num = 0;
            dy = std::min<int>(pu.lumaSize().height, DMVR_SUBCU_HEIGHT);
            dx = std::min<int>(pu.lumaSize().width, DMVR_SUBCU_WIDTH);
            if (PU::checkDMVRCondition(pu))
            {
              for (i = 0; i < (pu.lumaSize().height); i += dy)
              {
                for (j = 0; j < (pu.lumaSize().width); j += dx)
                {
                  refinedMvdL0[num][uiMergeCand] = pu.mvdL0SubPu[num];
                  num++;
                }
              }
            }
          }
        }
#endif

        Distortion uiSad = distParam.distFunc(distParam);
        m_CABACEstimator->getCtx() = ctxStart;
        uint64_t fracBits = m_pcInterSearch->xCalcPuMeBits(pu);
        double cost = (double)uiSad + (double)fracBits * sqrtLambdaForFirstPassIntra;
#if MULTI_HYP_PRED
        if (testMHP && pu.addHypData.size() < tempCS->sps->getMaxNumAddHyps())
        {
          uint32_t uiBitsCand = uiMergeCand + 1 + 1 + 1; // one bit for merge flag,  one bit for subblock_merge_flag, and one bit for regualr_merge_flag
          MEResult mergeResult;
          mergeResult.cu = cu;
          mergeResult.pu = pu;
          mergeResult.bits = uiBitsCand;
          mergeResult.cost = uiSad + m_pcRdCost->getCost(uiBitsCand);
          m_baseResultsForMH.push_back(mergeResult);
        }
#endif
        insertPos = -1;
#if CIIP_PDPC
        updateCandList(ModeInfo(uiMergeCand, true, false, false, false), cost, RdModeList, candCostList, uiNumMrgSATDCand, &insertPos);
#else
        updateCandList(ModeInfo(uiMergeCand, true, false, false), cost, RdModeList, candCostList, uiNumMrgSATDCand, &insertPos);
#endif
        if (insertPos != -1)
        {
          if (insertPos == RdModeList.size() - 1)
          {
            swap(singleMergeTempBuffer, acMergeTempBuffer[insertPos]);
          }
          else
          {
            for (uint32_t i = uint32_t(RdModeList.size()) - 1; i > insertPos; i--)
            {
              swap(acMergeTempBuffer[i - 1], acMergeTempBuffer[i]);
            }
            swap(singleMergeTempBuffer, acMergeTempBuffer[insertPos]);
          }
        }
        CHECK(std::min(uiMergeCand + 1, uiNumMrgSATDCand) != RdModeList.size(), "");
#if MULTI_PASS_DMVR
        pu.bdmvrRefine = false;
#endif
      }
#endif
      if (isIntrainterEnabled)
      {
#if MERGE_ENC_OPT
        xCheckSATDCostCiipMerge(tempCS, cu, pu, mergeCtx, acMergeTempBuffer, singleMergeTempBuffer, acMergeTmpBuffer, uiNumMrgSATDCand, RdModeList, candCostList, distParam, ctxStart);
#else
        // prepare for Intra bits calculation
        pu.ciipFlag = true;

        // save the to-be-tested merge candidates
        uint32_t CiipMergeCand[NUM_MRG_SATD_CAND];
        for (uint32_t mergeCnt = 0; mergeCnt < std::min(NUM_MRG_SATD_CAND, (const int)mergeCtx.numValidMergeCand); mergeCnt++)
        {
          CiipMergeCand[mergeCnt] = RdModeList[mergeCnt].mergeCand;
        }
        for (uint32_t mergeCnt = 0; mergeCnt < std::min(std::min(NUM_MRG_SATD_CAND, (const int)mergeCtx.numValidMergeCand), 4); mergeCnt++)
        {
          uint32_t mergeCand = CiipMergeCand[mergeCnt];
          acMergeTmpBuffer[mergeCand] = m_acMergeTmpBuffer[mergeCand].getBuf(localUnitArea);

          // estimate merge bits
          mergeCtx.setMergeInfo(pu, mergeCand);

          // first round
          pu.intraDir[0] = PLANAR_IDX;
#if CIIP_PDPC
          for (int intraCnt = 0; intraCnt < 2; intraCnt++)
          {
            pu.ciipPDPC = intraCnt == 1;
#else
          uint32_t intraCnt = 0;
#endif
          PelBuf ciipBuff = m_ciipBuffer[intraCnt].getBuf(localUnitArea.Y());
          // generate intrainter Y prediction
          if (mergeCnt == 0)
          {
            m_pcIntraSearch->initIntraPatternChType(*pu.cu, pu.Y());
            m_pcIntraSearch->predIntraAng(COMPONENT_Y, ciipBuff, pu);
          }
//#if INTER_LIC
//          if( mergeCtx.interDirNeighbours[mergeCand] != 3 )
//          {
//            pu.cs->getPredBuf( pu ).copyFrom( m_acRealMergeBuffer[MRG_MAX_NUM_CANDS + mergeCand].getBuf( localUnitArea ) );
//          }
//          else
//#endif
            if (pu.cs->picHeader->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag())
            {
              m_pcIntraSearch->geneWeightedPred<true>(COMPONENT_Y, singleMergeTempBuffer->Y(), pu, acMergeTmpBuffer[mergeCand].Y(), ciipBuff, m_pcReshape->getFwdLUT().data());
            }
            else
            {
              m_pcIntraSearch->geneWeightedPred<false>(COMPONENT_Y, singleMergeTempBuffer->Y(), pu, acMergeTmpBuffer[mergeCand].Y(), ciipBuff);
            }

          // calculate cost
            if (pu.cs->picHeader->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag())
            {
              PelBuf tmp = m_acGeoWeightedBuffer->getBuf(localUnitArea.Y());
              tmp.rspSignal(singleMergeTempBuffer->Y(), m_pcReshape->getInvLUT());
              distParam.cur = tmp;
            }
            else
            {
              distParam.cur = singleMergeTempBuffer->Y();
            }

          //distParam.cur = pu.cs->getPredBuf(pu).Y();
          Distortion sadValue = distParam.distFunc(distParam);
          if (pu.cs->slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag())
          {
            pu.cs->getPredBuf(pu).Y().rspSignal(m_pcReshape->getFwdLUT());
          }
          m_CABACEstimator->getCtx() = ctxStart;
          pu.regularMergeFlag = false;
          uint64_t fracBits = m_pcInterSearch->xCalcPuMeBits(pu);
          double cost = (double)sadValue + (double)fracBits * sqrtLambdaForFirstPassIntra;
          insertPos = -1;
#if CIIP_PDPC
          updateCandList(ModeInfo(mergeCand, false, false, true, pu.ciipPDPC), cost, RdModeList, candCostList, uiNumMrgSATDCand, &insertPos);
#else
          updateCandList(ModeInfo(mergeCand, false, false, true), cost, RdModeList, candCostList, uiNumMrgSATDCand, &insertPos);
#endif
          if (insertPos != -1)
          {
            for (int i = int(RdModeList.size()) - 1; i > insertPos; i--)
            {
              swap(acMergeTempBuffer[i - 1], acMergeTempBuffer[i]);
            }
            swap(singleMergeTempBuffer, acMergeTempBuffer[insertPos]);
          }
#if CIIP_PDPC
          }
#endif
        }
        pu.ciipFlag = false;
#if CIIP_PDPC
        pu.ciipPDPC = false;
#endif
#endif
      }
      if ( pu.cs->sps->getUseMMVD() )
      {
#if MERGE_ENC_OPT
#if JVET_W0090_ARMC_TM
        xCheckSATDCostMmvdMerge(tempCS, cu, pu, mergeCtxtmp, acMergeTempBuffer, singleMergeTempBuffer, uiNumMrgSATDCand, RdModeList, candCostList, distParam, ctxStart);
#else
        xCheckSATDCostMmvdMerge(tempCS, cu, pu, mergeCtx, acMergeTempBuffer, singleMergeTempBuffer, uiNumMrgSATDCand, RdModeList, candCostList, distParam, ctxStart);
#endif
#else
        cu.mmvdSkip = true;
        pu.regularMergeFlag = true;
        const int tempNum = (mergeCtx.numValidMergeCand > 1) ? MMVD_ADD_NUM : MMVD_ADD_NUM >> 1;
        for (int mmvdMergeCand = 0; mmvdMergeCand < tempNum; mmvdMergeCand++)
        {
          int baseIdx = mmvdMergeCand / MMVD_MAX_REFINE_NUM;
          int refineStep = (mmvdMergeCand - (baseIdx * MMVD_MAX_REFINE_NUM)) / 4;
          if (refineStep >= m_pcEncCfg->getMmvdDisNum())
            continue;
#if JVET_W0090_ARMC_TM
          mergeCtxtmp.setMmvdMergeCandiInfo(pu, mmvdMergeCand);
#else
          mergeCtx.setMmvdMergeCandiInfo(pu, mmvdMergeCand);
#endif

          PU::spanMotionInfo(pu, mergeCtx);
          pu.mvRefine = true;
          distParam.cur = singleMergeTempBuffer->Y();
          pu.mmvdEncOptMode = (refineStep > 2 ? 2 : 1);
          CHECK(!pu.mmvdMergeFlag, "MMVD merge should be set");
          // Don't do chroma MC here
          m_pcInterSearch->motionCompensation(pu, *singleMergeTempBuffer, REF_PIC_LIST_X, true, false);
          pu.mmvdEncOptMode = 0;
          pu.mvRefine = false;
          Distortion uiSad = distParam.distFunc(distParam);

          m_CABACEstimator->getCtx() = ctxStart;
          uint64_t fracBits = m_pcInterSearch->xCalcPuMeBits(pu);
          double cost = (double)uiSad + (double)fracBits * sqrtLambdaForFirstPassIntra;
#if MULTI_HYP_PRED
          if (testMHP && pu.addHypData.size() < tempCS->sps->getMaxNumAddHyps())
          {
            uint32_t uiBitsCand = baseIdx + refineStep + 2 + 1 + 1 + 1; // one bit for merge flag,  one bit for subblock_merge_flag, and one bit for regualr_merge_flag
            MEResult mergeResult;
            mergeResult.cu = cu;
            mergeResult.pu = pu;
            mergeResult.bits = uiBitsCand;
            mergeResult.cost = uiSad + m_pcRdCost->getCost(uiBitsCand);
            m_baseResultsForMH.push_back(mergeResult);
          }
#endif
          insertPos = -1;
#if CIIP_PDPC
          updateCandList(ModeInfo(mmvdMergeCand, false, true, false, false), cost, RdModeList, candCostList, uiNumMrgSATDCand, &insertPos);
#else
          updateCandList(ModeInfo(mmvdMergeCand, false, true, false), cost, RdModeList, candCostList, uiNumMrgSATDCand, &insertPos);
#endif
          if (insertPos != -1)
          {
            for (int i = int(RdModeList.size()) - 1; i > insertPos; i--)
            {
              swap(acMergeTempBuffer[i - 1], acMergeTempBuffer[i]);
            }
            swap(singleMergeTempBuffer, acMergeTempBuffer[insertPos]);
          }
        }
#endif
      }
#if MERGE_ENC_OPT
#if TM_MRG
      if (sps.getUseDMVDMode())
      {
        xCheckSATDCostTMMerge(tempCS, cu, pu, tmMrgCtx, acMergeTempBuffer, singleMergeTempBuffer, uiNumMrgSATDCand, RdModeList, candCostList, distParam, ctxStart
#if MULTI_PASS_DMVR
          , applyBDMVR4TM
#endif
        );
      }
#endif
      if (affineMrgAvail)
      {
        xCheckSATDCostAffineMerge(tempCS, cu, pu, affineMergeCtx, mrgCtx, acMergeTempBuffer, singleMergeTempBuffer, uiNumMrgSATDCand, RdModeList, candCostList, distParam, ctxStart);
      }
#if AFFINE_MMVD
      if (affineMmvdAvail)
      {
#if JVET_W0090_ARMC_TM
        xCheckSATDCostAffineMmvdMerge(tempCS, cu, pu, affineMergeCtxTmp, mrgCtx, acMergeTempBuffer, singleMergeTempBuffer, uiNumMrgSATDCand, RdModeList, candCostList, distParam, ctxStart);
#else
        xCheckSATDCostAffineMmvdMerge(tempCS, cu, pu, affineMergeCtx, mrgCtx, acMergeTempBuffer, singleMergeTempBuffer, uiNumMrgSATDCand, RdModeList, candCostList, distParam, ctxStart);
#endif
      }
#endif
#endif
      // Try to limit number of candidates using SATD-costs
      for( uint32_t i = 1; i < uiNumMrgSATDCand; i++ )
      {
        if( candCostList[i] > MRG_FAST_RATIO * candCostList[0] )
        {
          uiNumMrgSATDCand = i;
          break;
        }
      }

      setMergeBestSATDCost( candCostList[0] );

      if (isIntrainterEnabled && isChromaEnabled(pu.cs->pcv->chrFormat))
      {
        pu.ciipFlag = true;
        bool tag[2] = { false, false };

        for (uint32_t mergeCnt = 0; mergeCnt < uiNumMrgSATDCand; mergeCnt++)
        {
          if (RdModeList[mergeCnt].isCIIP)
          {
            pu.intraDir[0] = PLANAR_IDX;
            pu.intraDir[1] = DM_CHROMA_IDX;
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
            if (pu.chromaSize().width == 2)
              continue;
#endif
#if CIIP_PDPC
            pu.ciipPDPC = RdModeList[mergeCnt].isCiipPDPC;
            uint32_t bufIdx = pu.ciipPDPC ? 1 : 0;
#else
            uint32_t bufIdx = 0;
#endif
            if (!tag[bufIdx])
            {
              m_pcIntraSearch->initIntraPatternChType(*pu.cu, pu.Cb());
              PelBuf ciipBuffCb = m_ciipBuffer[bufIdx].getBuf(localUnitArea.Cb());
              m_pcIntraSearch->predIntraAng(COMPONENT_Cb, ciipBuffCb, pu);

              m_pcIntraSearch->initIntraPatternChType(*pu.cu, pu.Cr());
              PelBuf ciipBuffCr = m_ciipBuffer[bufIdx].getBuf(localUnitArea.Cr());
              m_pcIntraSearch->predIntraAng(COMPONENT_Cr, ciipBuffCr, pu);

              tag[bufIdx] = true;
            }
          }
        }
        pu.ciipFlag = false;
#if CIIP_PDPC
        pu.ciipPDPC = false;
#endif
      }

      tempCS->initStructData( encTestMode.qp );
      m_CABACEstimator->getCtx() = ctxStart;
    }
    else
    {
      if (bestIsMMVDSkip)
      {
        uiNumMrgSATDCand = mergeCtx.numValidMergeCand + ((mergeCtx.numValidMergeCand > 1) ? MMVD_ADD_NUM : MMVD_ADD_NUM >> 1);
      }
      else
      {
        uiNumMrgSATDCand = mergeCtx.numValidMergeCand;
      }
    }
  }
  m_bestModeUpdated = tempCS->useDbCost = bestCS->useDbCost = false;
  uint32_t iteration;
  uint32_t iterationBegin = 0;
  iteration = 2;
  for (uint32_t uiNoResidualPass = iterationBegin; uiNoResidualPass < iteration; ++uiNoResidualPass)
  {
    for( uint32_t uiMrgHADIdx = 0; uiMrgHADIdx < uiNumMrgSATDCand; uiMrgHADIdx++ )
    {
      uint32_t uiMergeCand = RdModeList[uiMrgHADIdx].mergeCand;

      if (uiNoResidualPass != 0 && RdModeList[uiMrgHADIdx].isCIIP) // intrainter does not support skip mode
      {
        if (isTestSkipMerge[uiMergeCand])
        {
          continue;
        }
      }

      if (((uiNoResidualPass != 0) && candHasNoResidual[uiMrgHADIdx])
       || ( (uiNoResidualPass == 0) && bestIsSkip ) )
      {
        continue;
      }
      // first get merge candidates
      CodingUnit &cu      = tempCS->addCU( tempCS->area, partitioner.chType );

      partitioner.setCUData( cu );
      cu.slice            = tempCS->slice;
      cu.tileIdx          = tempCS->pps->getTileIdx( tempCS->area.lumaPos() );
      cu.skip             = false;
      cu.mmvdSkip = false;
#if INTER_LIC
      cu.LICFlag          = false;
#endif
#if MERGE_ENC_OPT
      cu.affine = false;
#endif

      cu.geoFlag          = false;
    //cu.affine
      cu.predMode         = MODE_INTER;
    //cu.LICFlag
      cu.chromaQpAdj      = m_cuChromaQpOffsetIdxPlus1;
      cu.qp               = encTestMode.qp;
      PredictionUnit &pu  = tempCS->addPU( cu, partitioner.chType );
#if AFFINE_MMVD && MERGE_ENC_OPT
      pu.afMmvdFlag       = false;
#endif
#if TM_MRG
      pu.tmMergeFlag      = false;
#endif
#if MULTI_PASS_DMVR
      bool isDMVR         = false;
#endif
#if ENABLE_OBMC
      pu.ciipFlag = false;
#endif
      if (uiNoResidualPass == 0 && RdModeList[uiMrgHADIdx].isCIIP)
      {
        cu.mmvdSkip = false;
#if MULTI_HYP_PRED
        pu.ciipFlag = true;
        mergeCtx.setMergeInfo(pu, uiMergeCand);
#else
        mergeCtx.setMergeInfo(pu, uiMergeCand);
        pu.ciipFlag = true;
#endif
#if CIIP_PDPC
        pu.ciipPDPC = RdModeList[uiMrgHADIdx].isCiipPDPC;
#endif
        pu.regularMergeFlag = false;
        pu.intraDir[0] = PLANAR_IDX;
        CHECK(pu.intraDir[0]<0 || pu.intraDir[0]>(NUM_LUMA_MODE - 1), "out of intra mode");
        pu.intraDir[1] = DM_CHROMA_IDX;
      }
      else if (RdModeList[uiMrgHADIdx].isMMVD)
      {
        cu.mmvdSkip = true;
        pu.regularMergeFlag = true;
#if JVET_W0090_ARMC_TM
        mergeCtxtmp.setMmvdMergeCandiInfo(pu, uiMergeCand);
#else
        mergeCtx.setMmvdMergeCandiInfo(pu, uiMergeCand);
#endif
      }
#if MERGE_ENC_OPT
#if AFFINE_MMVD
      else if (RdModeList[uiMrgHADIdx].isAffineMmvd)
      {
        int baseIdx = (int)uiMergeCand / AF_MMVD_MAX_REFINE_NUM;
        int stepIdx = (int)uiMergeCand - baseIdx * AF_MMVD_MAX_REFINE_NUM;
        int dirIdx  = stepIdx % AF_MMVD_OFFSET_DIR;
            stepIdx = stepIdx / AF_MMVD_OFFSET_DIR;

        cu.affine           = true;
        cu.imv              = IMV_OFF;
        cu.mmvdSkip         = false;
        pu.regularMergeFlag = false;
        pu.mmvdMergeFlag    = false;

        pu.mergeFlag      = true;
        pu.afMmvdFlag     = true;
        pu.afMmvdBaseIdx  = (uint8_t)baseIdx;
        pu.afMmvdDir      = (uint8_t)dirIdx;
        pu.afMmvdStep     = (uint8_t)stepIdx;
        pu.mergeIdx       = (uint8_t)(baseIdx + afMmvdBaseIdxToMergeIdxOffset);
#if JVET_W0090_ARMC_TM
        pu.mergeType = affineMergeCtxTmp.mergeType[pu.mergeIdx];
#if INTER_LIC
        pu.cu->LICFlag = affineMergeCtxTmp.LICFlags[pu.mergeIdx];
#endif
        pu.interDir = affineMergeCtxTmp.interDirNeighbours[pu.mergeIdx];
        pu.cu->affineType = affineMergeCtxTmp.affineType[pu.mergeIdx];
        pu.cu->BcwIdx = affineMergeCtxTmp.BcwIdx[pu.mergeIdx];
        pu.mmvdMergeFlag = false;
        pu.ciipFlag = false;

        CHECK(pu.mergeIdx >= affineMergeCtxTmp.numValidMergeCand, "Invalid merge index for AffineMMVD");

        MvField mvfMmvd[2][3];
        PU::getAfMmvdMvf(pu, affineMergeCtxTmp, mvfMmvd, pu.mergeIdx, pu.afMmvdStep, pu.afMmvdDir);
#else
        pu.mergeType      = affineMergeCtx.mergeType         [pu.mergeIdx];
#if INTER_LIC
        pu.cu->LICFlag    = affineMergeCtx.LICFlags          [pu.mergeIdx];
#endif
        pu.interDir       = affineMergeCtx.interDirNeighbours[pu.mergeIdx];
        pu.cu->affineType = affineMergeCtx.affineType        [pu.mergeIdx];
        pu.cu->BcwIdx     = affineMergeCtx.BcwIdx            [pu.mergeIdx];
        pu.mmvdMergeFlag  = false;
        pu.ciipFlag       = false;

        CHECK(pu.mergeIdx >= affineMergeCtx.numValidMergeCand, "Invalid merge index for AffineMMVD");

        MvField mvfMmvd[2][3];
        PU::getAfMmvdMvf(pu, affineMergeCtx, mvfMmvd, pu.mergeIdx, pu.afMmvdStep, pu.afMmvdDir);
#endif

        for (int i = 0; i < 2; i++)
        {
          pu.refIdx[i] = mvfMmvd[i][0].refIdx;
          pu.mvAffi[i][0] = mvfMmvd[i][0].mv;
          pu.mvAffi[i][1] = mvfMmvd[i][1].mv;
          pu.mvAffi[i][2] = mvfMmvd[i][2].mv;
        }

        PU::spanMotionInfo(pu);
      }
#endif
      else if (RdModeList[uiMrgHADIdx].isAffine)
      {
        CHECK(uiMergeCand >= affineMergeCtx.numValidMergeCand, "");
        cu.mmvdSkip = false;
        cu.affine = true;
        cu.imv = 0;
        pu.regularMergeFlag = false;
        pu.mergeFlag = true;
        pu.mergeIdx = uiMergeCand;
        pu.mmvdMergeFlag = false;
        pu.interDir = affineMergeCtx.interDirNeighbours[uiMergeCand];
        cu.affineType = affineMergeCtx.affineType[uiMergeCand];
        cu.BcwIdx = affineMergeCtx.BcwIdx[uiMergeCand];
#if INTER_LIC
        cu.LICFlag = affineMergeCtx.LICFlags[uiMergeCand];
#endif
        pu.mv[0].setZero();
        pu.mv[1].setZero();
        pu.mvd[REF_PIC_LIST_0] = Mv();
        pu.mvd[REF_PIC_LIST_1] = Mv();
        pu.mvpIdx[REF_PIC_LIST_0] = NOT_VALID;
        pu.mvpIdx[REF_PIC_LIST_1] = NOT_VALID;
        pu.mvpNum[REF_PIC_LIST_0] = NOT_VALID;
        pu.mvpNum[REF_PIC_LIST_1] = NOT_VALID;
        pu.mergeType = affineMergeCtx.mergeType[uiMergeCand];
        if (pu.mergeType == MRG_TYPE_SUBPU_ATMVP)
        {
          pu.refIdx[0] = affineMergeCtx.mvFieldNeighbours[(uiMergeCand << 1) + 0][0].refIdx;
          pu.refIdx[1] = affineMergeCtx.mvFieldNeighbours[(uiMergeCand << 1) + 1][0].refIdx;
          PU::spanMotionInfo(pu, mrgCtx);
        }
        else
        {
          for (int i = 0; i < 2; i++)
          {
            pu.refIdx[i] = affineMergeCtx.mvFieldNeighbours[(uiMergeCand << 1) + i][0].refIdx;
            pu.mvAffi[i][0] = affineMergeCtx.mvFieldNeighbours[(uiMergeCand << 1) + i][0].mv;
            pu.mvAffi[i][1] = affineMergeCtx.mvFieldNeighbours[(uiMergeCand << 1) + i][1].mv;
            pu.mvAffi[i][2] = affineMergeCtx.mvFieldNeighbours[(uiMergeCand << 1) + i][2].mv;
          }

          PU::spanMotionInfo(pu);
        }
      }
#if TM_MRG && MERGE_ENC_OPT
      else if (RdModeList[uiMrgHADIdx].isTMMrg)
      {
        cu.mmvdSkip         = false;
        pu.regularMergeFlag = true;
        pu.tmMergeFlag      = true;
        tmMrgCtx.setMergeInfo(pu, uiMergeCand);
#if MULTI_PASS_DMVR
        if (applyBDMVR4TM[uiMergeCand])
        {
          isDMVR = true;
          pu.bdmvrRefine = true;
          m_pcInterSearch->setBdmvrSubPuMvBuf(m_mvBufBDMVR4TM[uiMergeCand << 1], m_mvBufBDMVR4TM[(uiMergeCand << 1) + 1]);
        }
#endif
      }
#endif
#endif
      else
      {
        cu.mmvdSkip = false;
        pu.regularMergeFlag = true;
        mergeCtx.setMergeInfo(pu, uiMergeCand);
#if MULTI_PASS_DMVR
        if (applyBDMVR[uiMergeCand])
        {
          isDMVR = true;
          pu.bdmvrRefine = true;
          m_pcInterSearch->setBdmvrSubPuMvBuf(m_mvBufBDMVR[uiMergeCand << 1], m_mvBufBDMVR[(uiMergeCand << 1) + 1]);
        }
#endif
      }
#if MERGE_ENC_OPT
      if (!RdModeList[uiMrgHADIdx].isAffine && !RdModeList[uiMrgHADIdx].isGeo)
#endif
#if MULTI_PASS_DMVR
      if( !pu.bdmvrRefine )
      {
        PU::spanMotionInfo( pu, mergeCtx );
      }
#else
      PU::spanMotionInfo( pu, mergeCtx );
#endif

      if( m_pcEncCfg->getMCTSEncConstraint() )
      {
#if !MULTI_PASS_DMVR
        bool isDMVR = PU::checkDMVRCondition( pu );
#endif
        if( ( isDMVR && MCTSHelper::isRefBlockAtRestrictedTileBoundary( pu ) ) || ( !isDMVR && !( MCTSHelper::checkMvBufferForMCTSConstraint( pu ) ) ) )
        {
          // Do not use this mode
          tempCS->initStructData( encTestMode.qp );
          continue;
        }
      }
#if MERGE_ENC_OPT
      if (mrgTempBufSet && uiMrgHADIdx < MMVD_MRG_MAX_RD_NUM)
#else
      if( mrgTempBufSet )
#endif
      {
#if !MULTI_PASS_DMVR
        {
          int dx, dy, i, j, num = 0;
          dy = std::min<int>(pu.lumaSize().height, DMVR_SUBCU_HEIGHT);
          dx = std::min<int>(pu.lumaSize().width, DMVR_SUBCU_WIDTH);
          if (PU::checkDMVRCondition(pu))
          {
            for (i = 0; i < (pu.lumaSize().height); i += dy)
            {
              for (j = 0; j < (pu.lumaSize().width); j += dx)
              {
                pu.mvdL0SubPu[num] = refinedMvdL0[num][uiMergeCand];
                num++;
              }
            }
          }
        }
#endif
        if (pu.ciipFlag)
        {
#if CIIP_PDPC
          uint32_t bufIdx = pu.ciipPDPC ? 1 : 0;
#else
          uint32_t bufIdx = 0;
#endif
          // Luma CIIP was already done in SATD check stage and stored
          tempCS->getPredBuf().Y().copyFrom( acMergeTempBuffer[uiMrgHADIdx]->Y() );

#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
          if (isChromaEnabled(pu.chromaFormat))
#else
          if( isChromaEnabled( pu.chromaFormat ) && pu.chromaSize().width > 2 )
#endif
          {
            m_pcIntraSearch->geneWeightedPred<false>( COMPONENT_Cb, tempCS->getPredBuf( pu ).Cb(), pu, acMergeTmpBuffer[uiMergeCand].Cb(), m_ciipBuffer[bufIdx].getBuf( localUnitArea.Cb() ) );
            m_pcIntraSearch->geneWeightedPred<false>( COMPONENT_Cr, tempCS->getPredBuf( pu ).Cr(), pu, acMergeTmpBuffer[uiMergeCand].Cr(), m_ciipBuffer[bufIdx].getBuf( localUnitArea.Cr() ) );
          }
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
          else if (isChromaEnabled(pu.chromaFormat))
          {
            tempCS->getPredBuf().Cb().copyFrom(acMergeTmpBuffer[uiMergeCand].Cb());
            tempCS->getPredBuf().Cr().copyFrom(acMergeTmpBuffer[uiMergeCand].Cr());
          }
#endif
        }
        else
        {
          if (RdModeList[uiMrgHADIdx].isMMVD
#if AFFINE_MMVD && MERGE_ENC_OPT
            || RdModeList[uiMrgHADIdx].isAffineMmvd
#endif
            )
          {
            pu.mmvdEncOptMode = 0;
            m_pcInterSearch->motionCompensation(pu);
          }
#if MERGE_ENC_OPT
          else if (uiNoResidualPass != 0 && RdModeList[uiMrgHADIdx].isCIIP)
          {
            // perform regular MC instead, i.e. test skip mode
            pu.mvRefine = true;
            m_pcInterSearch->motionCompensation(pu);
            pu.mvRefine = false;
#if MULTI_PASS_DMVR
            if (!RdModeList[uiMrgHADIdx].isAffine && !RdModeList[uiMrgHADIdx].isGeo && pu.bdmvrRefine)
            {
#if TM_MRG
              if( pu.tmMergeFlag )
              {
                PU::spanMotionInfo( pu, mergeCtx, m_mvBufBDMVR4TM[uiMergeCand << 1], m_mvBufBDMVR4TM[( uiMergeCand << 1 ) + 1], m_pcInterSearch->getBdofSubPuMvOffset() );
              }
              else
#endif
              PU::spanMotionInfo(pu, mergeCtx, m_mvBufBDMVR[uiMergeCand << 1], m_mvBufBDMVR[(uiMergeCand << 1) + 1], m_pcInterSearch->getBdofSubPuMvOffset());
            }
#endif
          }
#else
          else if (uiNoResidualPass != 0 && RdModeList[uiMrgHADIdx].isCIIP)
          {
            tempCS->getPredBuf().copyFrom(acMergeBuffer[uiMergeCand]);
#if MULTI_PASS_DMVR
            if (pu.bdmvrRefine)
            {
#if TM_MRG
              if( pu.tmMergeFlag )
              {
                PU::spanMotionInfo( pu, mergeCtx, m_mvBufBDMVR4TM[uiMergeCand << 1], m_mvBufBDMVR4TM[( uiMergeCand << 1 ) + 1], m_mvBufEncBDOF4TM[uiMergeCand] );
              }
              else
#endif
                PU::spanMotionInfo(pu, mergeCtx, m_mvBufBDMVR[uiMergeCand << 1], m_mvBufBDMVR[(uiMergeCand << 1) + 1], m_mvBufEncBDOF[uiMergeCand]);
            }
#endif
          }
#endif
#if MERGE_ENC_OPT
          else if (RdModeList[uiMrgHADIdx].isAffine)
          {
            tempCS->getPredBuf().copyFrom(*acMergeTempBuffer[uiMrgHADIdx], true);
            m_pcInterSearch->motionCompensation(pu, REF_PIC_LIST_X, false, true);
          }
#endif
          else
          {
            tempCS->getPredBuf().copyFrom(*acMergeTempBuffer[uiMrgHADIdx]);
#if MULTI_PASS_DMVR
#if MERGE_ENC_OPT
            if (!RdModeList[uiMrgHADIdx].isAffine && !RdModeList[uiMrgHADIdx].isGeo && pu.bdmvrRefine)
#else
            if(pu.bdmvrRefine)
#endif
            {
#if TM_MRG
              if( pu.tmMergeFlag )
              {
                PU::spanMotionInfo( pu, mergeCtx, m_mvBufBDMVR4TM[uiMergeCand << 1], m_mvBufBDMVR4TM[( uiMergeCand << 1 ) + 1], m_mvBufEncBDOF4TM[uiMergeCand] );
              }
              else
#endif
              PU::spanMotionInfo(pu, mergeCtx, m_mvBufBDMVR[uiMergeCand << 1], m_mvBufBDMVR[(uiMergeCand << 1) + 1], m_mvBufEncBDOF[uiMergeCand]);
            }
#endif
          }
        }
      }
      else
      {
        pu.mvRefine = true;
        m_pcInterSearch->motionCompensation( pu );
        pu.mvRefine = false;
#if MULTI_PASS_DMVR
        if (pu.bdmvrRefine)
        {
#if TM_MRG
          if( pu.tmMergeFlag )
          {
            PU::spanMotionInfo( pu, mergeCtx, m_mvBufBDMVR4TM[uiMergeCand << 1], m_mvBufBDMVR4TM[( uiMergeCand << 1 ) + 1], m_pcInterSearch->getBdofSubPuMvOffset() );
          }
          else
#endif
            PU::spanMotionInfo(pu, mergeCtx, m_mvBufBDMVR[uiMergeCand << 1], m_mvBufBDMVR[(uiMergeCand << 1) + 1], m_pcInterSearch->getBdofSubPuMvOffset());
        }
#endif
      }
      if (!cu.mmvdSkip && !pu.ciipFlag && uiNoResidualPass != 0 && !cu.affine 
#if TM_MRG
          && !pu.tmMergeFlag
#endif
#if AFFINE_MMVD && MERGE_ENC_OPT
          && !pu.afMmvdFlag
#endif
        )
      {
        CHECK(uiMergeCand >= mergeCtx.numValidMergeCand, "out of normal merge");
        isTestSkipMerge[uiMergeCand] = true;
      }

#if ENABLE_OBMC
      cu.isobmcMC = true;
      cu.obmcFlag = true;
      m_pcInterSearch->subBlockOBMC( pu );
      cu.isobmcMC = false;
#endif

      xEncodeInterResidual( tempCS, bestCS, partitioner, encTestMode, uiNoResidualPass, uiNoResidualPass == 0 ? &candHasNoResidual[uiMrgHADIdx] : NULL );

      if( m_pcEncCfg->getUseFastDecisionForMerge() && !bestIsSkip && !pu.ciipFlag)
      {
        bestIsSkip = !bestCS->cus.empty() && bestCS->getCU( partitioner.chType )->rootCbf == 0;
      }
      tempCS->initStructData( encTestMode.qp );
    }// end loop uiMrgHADIdx

    if( uiNoResidualPass == 0 && m_pcEncCfg->getUseEarlySkipDetection() )
    {
      const CodingUnit     &bestCU = *bestCS->getCU( partitioner.chType );
      const PredictionUnit &bestPU = *bestCS->getPU( partitioner.chType );

      if( bestCU.rootCbf == 0 )
      {
        if( bestPU.mergeFlag )
        {
          m_modeCtrl->setEarlySkipDetected();
        }
        else if( m_pcEncCfg->getMotionEstimationSearchMethod() != MESEARCH_SELECTIVE )
        {
          int absolute_MV = 0;

          for( uint32_t uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++ )
          {
            if( slice.getNumRefIdx( RefPicList( uiRefListIdx ) ) > 0 )
            {
              absolute_MV += bestPU.mvd[uiRefListIdx].getAbsHor() + bestPU.mvd[uiRefListIdx].getAbsVer();
            }
          }

          if( absolute_MV == 0 )
          {
            m_modeCtrl->setEarlySkipDetected();
          }
        }
      }
    }
  }
  if ( m_bestModeUpdated && bestCS->cost != MAX_DOUBLE )
  {
    xCalDebCost( *bestCS, partitioner );
  }
}

#if JVET_W0097_GPM_MMVD_TM
void EncCu::xCheckRDCostMergeGeoComb2Nx2N(CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &pm, const EncTestMode& encTestMode, bool isSecondPass)
{
  int numSATDCands = ((fastGpmMmvdSearch && isSecondPass) ? 60 : 70);

  tempCS->initStructData(encTestMode.qp);
#if TM_MRG
  MergeCtx mergeCtx[GEO_NUM_TM_MV_CAND];
#else
  MergeCtx mergeCtx;
#endif
  const SPS &sps = *tempCS->sps;
  CodedCUInfo& relatedCU = ((EncModeCtrlMTnoRQT *)m_modeCtrl)->getBlkInfo(pm.currArea());
  bool extMMVD = tempCS->picHeader->getGPMMMVDTableFlag();

  if (sps.getSbTMVPEnabledFlag())
  {
    Size bufSize = g_miScaling.scale(tempCS->area.lumaSize());
#if TM_MRG
    for (int i = 0; i < GEO_NUM_TM_MV_CAND; i++)
    {
      mergeCtx[i].subPuMvpMiBuf = MotionBuf(m_SubPuMiBuf, bufSize);
    }
#else
    mergeCtx.subPuMvpMiBuf = MotionBuf(m_SubPuMiBuf, bufSize);
#endif
  }

  // 1. bit estimation
  const double sqrtLambdaFracBits = m_pcRdCost->getMotionLambda() * FRAC_BITS_SCALE;
  uint8_t maxNumMergeCandidates = tempCS->sps->getMaxNumGeoCand();
  const TempCtx ctxStart(m_CtxCache, m_CABACEstimator->getCtx());

  double geoModeCost[GEO_NUM_PARTITION_MODE], geoMergeIdxCost[MRG_MAX_NUM_CANDS], geoMMVDFlagCost[2], geoMMVDIdxCost[GPM_EXT_MMVD_MAX_REFINE_NUM];
#if TM_MRG
  double geoTMFlagCost[2];
#endif
  for (int idx = 0; idx < GEO_NUM_PARTITION_MODE; idx++)
  {
    uint64_t fracBits = m_CABACEstimator->geo_mode_est(ctxStart, idx);
    geoModeCost[idx] = (double)fracBits * sqrtLambdaFracBits;
  }
  for (int idx = 0; idx < maxNumMergeCandidates; idx++)
  {
    uint64_t fracBits = m_CABACEstimator->geo_mergeIdx_est(ctxStart, idx, maxNumMergeCandidates);
    geoMergeIdxCost[idx] = (double)fracBits * sqrtLambdaFracBits;
  }
  for (int idx = 0; idx < 2; idx++)
  {
    uint64_t fracBits = m_CABACEstimator->geo_mmvdFlag_est(ctxStart, idx);
    geoMMVDFlagCost[idx] = (double)fracBits * sqrtLambdaFracBits;
  }
  for (int idx = 0; idx < (extMMVD ? GPM_EXT_MMVD_MAX_REFINE_NUM : GPM_MMVD_MAX_REFINE_NUM); idx++)
  {
    uint64_t fractBits = m_CABACEstimator->geo_mmvdIdx_est(ctxStart, idx, extMMVD);
    geoMMVDIdxCost[idx] = (double)fractBits * sqrtLambdaFracBits;
  }
#if TM_MRG
  for (int idx = 0; idx < 2 && sps.getUseDMVDMode(); idx++)
  {
    uint64_t fracBits = m_CABACEstimator->geo_tmFlag_est(ctxStart, idx);
    geoTMFlagCost[idx] = (double)fracBits * sqrtLambdaFracBits;
  }
#endif
  m_CABACEstimator->getCtx() = ctxStart;

  // 2. get SAD for all candidates
  CodingUnit &cu = tempCS->addCU(tempCS->area, pm.chType);
  pm.setCUData(cu);
  cu.predMode = MODE_INTER;
  cu.slice = tempCS->slice;
  cu.tileIdx = tempCS->pps->getTileIdx(tempCS->area.lumaPos());
  cu.qp = encTestMode.qp;
  cu.affine = false;
  cu.mtsFlag = false;
#if INTER_LIC
  cu.LICFlag = false;
#endif
  cu.BcwIdx = BCW_DEFAULT;
  cu.geoFlag = true;
  cu.imv = 0;
  cu.mmvdSkip = false;
  cu.skip = false;
  cu.mipFlag = false;
#if JVET_V0130_INTRA_TMP
  cu.tmpFlag = false;
#endif
  cu.bdpcmMode = 0;

  PredictionUnit &pu = tempCS->addPU(cu, pm.chType);
  pu.mergeFlag = true;
  pu.regularMergeFlag = false;
#if TM_MRG
  pu.tmMergeFlag = false;
#endif
  CHECK(!m_mergeCandAvail, "merge candidates are not available");
#if TM_MRG
  PU::getGeoMergeCandidates(pu, mergeCtx[GEO_TM_OFF], &m_mergeCand);
  maxNumMergeCandidates = min((int)maxNumMergeCandidates, mergeCtx[GEO_TM_OFF].numValidMergeCand);
#else
  PU::getGeoMergeCandidates(pu, mergeCtx, &m_mergeCand);
  maxNumMergeCandidates = min((int)maxNumMergeCandidates, mergeCtx.numValidMergeCand);
#endif

  PelUnitBuf geoBuffer[GEO_MAX_NUM_UNI_CANDS];
  PelUnitBuf geoTempBuf[GEO_MAX_NUM_UNI_CANDS];
  PelUnitBuf geoMMVDBuf[GEO_MAX_NUM_UNI_CANDS][GPM_EXT_MMVD_MAX_REFINE_NUM];
  PelUnitBuf geoMMVDTempBuf[GEO_MAX_NUM_UNI_CANDS][GPM_EXT_MMVD_MAX_REFINE_NUM];
  PelUnitBuf geoCombinations[GEO_MAX_TRY_WEIGHTED_SAD];
  DistParam  distParam;

  const UnitArea localUnitArea(tempCS->area.chromaFormat, Area(0, 0, tempCS->area.Y().width, tempCS->area.Y().height));
  DistParam distParamWholeBlk;
  m_pcRdCost->setDistParam(distParamWholeBlk, tempCS->getOrgBuf().Y(), m_acMergeBuffer[0].Y().buf, m_acMergeBuffer[0].Y().stride, sps.getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y);
  Distortion sadWholeBlk[GEO_MAX_NUM_UNI_CANDS], sadMMVDWholeBlk[GEO_MAX_NUM_UNI_CANDS][GPM_EXT_MMVD_MAX_REFINE_NUM];

  int  pocMrg[GEO_MAX_NUM_UNI_CANDS];
  Mv   mrgMv[GEO_MAX_NUM_UNI_CANDS];
  bool mrgDuplicated[GEO_MAX_NUM_UNI_CANDS];

  double bestMrgCost = MAX_DOUBLE;
  double bestNormalMrgCost = MAX_DOUBLE;

  bool isGeoChromaAvail[GEO_MAX_NUM_UNI_CANDS];
  bool isGeoMMVDChromaAvail[GEO_MAX_NUM_UNI_CANDS][GPM_EXT_MMVD_MAX_REFINE_NUM];
  memset(isGeoChromaAvail, false, sizeof(bool) * GEO_MAX_NUM_UNI_CANDS);
  memset(isGeoMMVDChromaAvail, false, sizeof(bool) * GEO_MAX_NUM_UNI_CANDS * GPM_EXT_MMVD_MAX_REFINE_NUM);
#if TM_MRG
  bool isGeoTmChromaAvail[GEO_TM_MAX_NUM_CANDS];
  memset(isGeoTmChromaAvail, false, sizeof(bool) * GEO_TM_MAX_NUM_CANDS);
#endif

  for (uint8_t mergeCand = 0; mergeCand < maxNumMergeCandidates; mergeCand++)
  {
#if TM_MRG
    int mrgList = mergeCtx[GEO_TM_OFF].mvFieldNeighbours[(mergeCand << 1) + 0].refIdx == -1 ? 1 : 0;
    int mrgRefIdx = mergeCtx[GEO_TM_OFF].mvFieldNeighbours[(mergeCand << 1) + mrgList].refIdx;
#else
    int mrgList = mergeCtx.mvFieldNeighbours[(mergeCand << 1) + 0].refIdx == -1 ? 1 : 0;
    int mrgRefIdx = mergeCtx.mvFieldNeighbours[(mergeCand << 1) + mrgList].refIdx;
#endif
    pocMrg[mergeCand] = tempCS->slice->getRefPic((RefPicList)mrgList, mrgRefIdx)->getPOC();
#if TM_MRG
    mrgMv[mergeCand] = mergeCtx[GEO_TM_OFF].mvFieldNeighbours[(mergeCand << 1) + mrgList].mv;
#else
    mrgMv[mergeCand] = mergeCtx.mvFieldNeighbours[(mergeCand << 1) + mrgList].mv;
#endif
    mrgDuplicated[mergeCand] = false;
    if (mergeCand)
    {
      for (int i = 0; i < mergeCand; i++)
      {
        if (pocMrg[mergeCand] == pocMrg[i] && mrgMv[mergeCand] == mrgMv[i])
        {
          mrgDuplicated[mergeCand] = true;
          break;
        }
      }
    }
#if !MULTI_HYP_PRED
    if (mrgDuplicated[mergeCand])
    {
      continue;
    }
#endif
    geoBuffer[mergeCand] = m_acMergeBuffer[mergeCand].getBuf(localUnitArea);
#if TM_MRG
    mergeCtx[GEO_TM_OFF].setMergeInfo(pu, mergeCand);
#else
    mergeCtx.setMergeInfo(pu, mergeCand);
#endif
    if (m_pcEncCfg->getMCTSEncConstraint() && (!(MCTSHelper::checkMvBufferForMCTSConstraint(pu))))
    {
      tempCS->initStructData(encTestMode.qp);
      return;
    }
    m_pcInterSearch->motionCompensation(pu, geoBuffer[mergeCand], REF_PIC_LIST_X, true, false);
#if MULTI_HYP_PRED
    geoTempBuf[mergeCand] = m_acRealMergeBuffer[MRG_MAX_NUM_CANDS + mergeCand].getBuf(localUnitArea);
#else
    geoTempBuf[mergeCand] = m_acMergeTmpBuffer[mergeCand].getBuf(localUnitArea);
#endif
    geoTempBuf[mergeCand].Y().copyFrom(geoBuffer[mergeCand].Y());
    geoTempBuf[mergeCand].Y().roundToOutputBitdepth(geoTempBuf[mergeCand].Y(), cu.slice->clpRng(COMPONENT_Y));
    distParamWholeBlk.cur.buf = geoTempBuf[mergeCand].Y().buf;
    distParamWholeBlk.cur.stride = geoTempBuf[mergeCand].Y().stride;
    sadWholeBlk[mergeCand] = distParamWholeBlk.distFunc(distParamWholeBlk);
    double curCost = sadWholeBlk[mergeCand] + geoMergeIdxCost[mergeCand];
    if (curCost < bestNormalMrgCost)
    {
      bestNormalMrgCost = curCost;
    }
    curCost += geoMMVDFlagCost[0];
    if (curCost < bestMrgCost)
    {
      bestMrgCost = curCost;
    }
  }
#if MULTI_HYP_PRED
#if TM_MRG
  m_pcInterSearch->setGeoTmpBuffer(mergeCtx[GEO_TM_OFF]);
#else
  m_pcInterSearch->setGeoTmpBuffer(mergeCtx);
#endif
#endif

  int wIdx = floorLog2(cu.lwidth()) - GEO_MIN_CU_LOG2;
  int hIdx = floorLog2(cu.lheight()) - GEO_MIN_CU_LOG2;
  Distortion sadSmall = 0, sadLarge = 0;
  int maskStride = 0, maskStride2 = 0, stepX = 1;
  Pel* SADmask;
  static_vector<int, GEO_NUM_PARTITION_MODE> selGeoModeList;
  static_vector<double, GEO_NUM_PARTITION_MODE> selGeoModeRDList;
  static_vector<int, 5> mergeCandList0[GEO_NUM_PARTITION_MODE];
  static_vector<int, 5> mergeCandList1[GEO_NUM_PARTITION_MODE];
  static_vector<int, 5> mmvdCandList0[GEO_NUM_PARTITION_MODE];
  static_vector<int, 5> mmvdCandList1[GEO_NUM_PARTITION_MODE];
  static_vector<double, 5> sadCostList0[GEO_NUM_PARTITION_MODE];
  static_vector<double, 5> sadCostList1[GEO_NUM_PARTITION_MODE];

  for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; splitDir++)
  {
    int16_t angle = g_GeoParams[splitDir][0];
    if (g_angle2mirror[angle] == 2)
    {
      stepX = 1;
      maskStride = -GEO_WEIGHT_MASK_SIZE;
      maskStride2 = -(int)cu.lwidth();
      SADmask = &g_globalGeoEncSADmask[g_angle2mask[g_GeoParams[splitDir][0]]][(GEO_WEIGHT_MASK_SIZE - 1 - g_weightOffset[splitDir][hIdx][wIdx][1]) * GEO_WEIGHT_MASK_SIZE + g_weightOffset[splitDir][hIdx][wIdx][0]];
    }
    else if (g_angle2mirror[angle] == 1)
    {
      stepX = -1;
      maskStride2 = cu.lwidth();
      maskStride = GEO_WEIGHT_MASK_SIZE;
      SADmask = &g_globalGeoEncSADmask[g_angle2mask[g_GeoParams[splitDir][0]]][g_weightOffset[splitDir][hIdx][wIdx][1] * GEO_WEIGHT_MASK_SIZE + (GEO_WEIGHT_MASK_SIZE - 1 - g_weightOffset[splitDir][hIdx][wIdx][0])];
    }
    else
    {
      stepX = 1;
      maskStride = GEO_WEIGHT_MASK_SIZE;
      maskStride2 = -(int)cu.lwidth();
      SADmask = &g_globalGeoEncSADmask[g_angle2mask[g_GeoParams[splitDir][0]]][g_weightOffset[splitDir][hIdx][wIdx][1] * GEO_WEIGHT_MASK_SIZE + g_weightOffset[splitDir][hIdx][wIdx][0]];
    }
    for (uint8_t mergeCand = 0; mergeCand < maxNumMergeCandidates; mergeCand++)
    {
      if (mrgDuplicated[mergeCand])
      {
        continue;
      }
      m_pcRdCost->setDistParam(distParam, tempCS->getOrgBuf().Y(), geoTempBuf[mergeCand].Y().buf, geoTempBuf[mergeCand].Y().stride, SADmask, maskStride, stepX, maskStride2, sps.getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y);
      sadLarge = distParam.distFunc(distParam);
      double tempCost = (double)sadLarge + geoMergeIdxCost[mergeCand] + geoMMVDFlagCost[0];
      m_GeoMMVDCostList.insert(splitDir, 0, mergeCand, 0, tempCost);
      sortCandList(tempCost, mergeCand, 0, sadCostList0[splitDir], mergeCandList0[splitDir], mmvdCandList0[splitDir], numCandPerPar);
      sadSmall = sadWholeBlk[mergeCand] - sadLarge;
      tempCost = (double)sadSmall + geoMergeIdxCost[mergeCand] + geoMMVDFlagCost[0];
      m_GeoMMVDCostList.insert(splitDir, 1, mergeCand, 0, tempCost);
      sortCandList(tempCost, mergeCand, 0, sadCostList1[splitDir], mergeCandList1[splitDir], mmvdCandList1[splitDir], numCandPerPar);
    }
    updateCandList(splitDir, (sadCostList0[splitDir][0] + sadCostList1[splitDir][0]), selGeoModeList, selGeoModeRDList, GEO_NUM_PARTITION_MODE);
  }

  static_vector<int, GEO_MAX_TRY_WEIGHTED_SAD> geoSplitDirList;
  static_vector<int, GEO_MAX_TRY_WEIGHTED_SAD> geoMergeCand0;
  static_vector<int, GEO_MAX_TRY_WEIGHTED_SAD> geoMergeCand1;
  static_vector<int, GEO_MAX_TRY_WEIGHTED_SAD> geoMmvdCand0;
  static_vector<int, GEO_MAX_TRY_WEIGHTED_SAD> geoMmvdCand1;
  static_vector<double, GEO_MAX_TRY_WEIGHTED_SAD> geoSADCostList;

  for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; splitDir++)
  {
    int numCandPart0 = min(numCandPerPar, (int)mergeCandList0[splitDir].size());
    int numCandPart1 = min(numCandPerPar, (int)mergeCandList1[splitDir].size());
    for (int candIdx0 = 0; candIdx0 < numCandPart0; candIdx0++)
    {
      for (int candIdx1 = 0; candIdx1 < numCandPart1; candIdx1++)
      {
        int mergeCand0 = mergeCandList0[splitDir][candIdx0];
        int mergeCand1 = mergeCandList1[splitDir][candIdx1];

        if (mergeCand0 == mergeCand1)
        {
          continue;
        }

        double tempCost = m_GeoMMVDCostList.singleDistList[0][splitDir][mergeCand0][0].cost + m_GeoMMVDCostList.singleDistList[1][splitDir][mergeCand1][0].cost;
        tempCost = tempCost + geoModeCost[splitDir];
#if TM_MRG
        if (sps.getUseDMVDMode())
        {
          tempCost += geoTMFlagCost[0];
        }
#endif
        updateGeoMMVDCandList(tempCost, splitDir, mergeCand0, mergeCand1, 0, 0, geoSADCostList, geoSplitDirList, geoMergeCand0, geoMergeCand1, geoMmvdCand0, geoMmvdCand1, numSATDCands);
      }
    }
  }

  static_vector<uint8_t, GEO_MAX_TRY_WEIGHTED_SAD>  geoRdModeList;
  static_vector<bool, GEO_MAX_TRY_WEIGHTED_SAD>  isNonMMVDListIdx;
  static_vector<int, GEO_MAX_TRY_WEIGHTED_SAD>  geoPartitionModeList;
  static_vector<double, GEO_MAX_TRY_WEIGHTED_SAD>  geocandCostList;

  DistParam distParamSAD2;
  const bool useHadamard = !tempCS->slice->getDisableSATDForRD();
  m_pcRdCost->setDistParam(distParamSAD2, tempCS->getOrgBuf().Y(), m_acMergeBuffer[0].Y(), sps.getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, useHadamard);

  int numberGeoCandChecked = (int)geoSADCostList.size();
  int geoNumMrgSATDCand = min(GEO_MAX_TRY_WEIGHTED_SATD, numberGeoCandChecked);
  int numStoredCands = geoNumMrgSATDCand;
  for (uint8_t candidateIdx = 0; candidateIdx < numberGeoCandChecked; candidateIdx++)
  {
    int splitDir = geoSplitDirList[candidateIdx];
    int mergeCand0 = geoMergeCand0[candidateIdx];
    int mergeCand1 = geoMergeCand1[candidateIdx];
    bool mmvdFlag0 = false;
    bool mmvdFlag1 = false;

    geoCombinations[candidateIdx] = m_acGeoWeightedBuffer[candidateIdx].getBuf(localUnitArea);
    m_pcInterSearch->weightedGeoBlk(pu, splitDir, CHANNEL_TYPE_LUMA, geoCombinations[candidateIdx], geoBuffer[mergeCand0], geoBuffer[mergeCand1]);
    distParamSAD2.cur = geoCombinations[candidateIdx].Y();
    Distortion sad = distParamSAD2.distFunc(distParamSAD2);

    double updateCost = geoModeCost[splitDir] + geoMergeIdxCost[mergeCand0] + (fastGpmMmvdSearch ? geoMergeIdxCost[mergeCand1 > mergeCand0 ? (mergeCand1 - 1) : mergeCand1] : geoMergeIdxCost[mergeCand1]) + geoMMVDFlagCost[mmvdFlag0] + geoMMVDFlagCost[mmvdFlag1];
#if TM_MRG
    if (sps.getUseDMVDMode())
    {
      updateCost += geoTMFlagCost[0];
    }
#endif
    updateCost += (double)sad;
    orderCandList(candidateIdx, true, splitDir, updateCost, geoRdModeList, isNonMMVDListIdx, geoPartitionModeList, geocandCostList, numStoredCands);
  }

  for (uint8_t i = 1; i < geoNumMrgSATDCand; i++)
  {
#if MERGE_ENC_OPT
    if (geocandCostList[i] > MRG_FAST_RATIO * geocandCostList[0] || geocandCostList[i] > getMergeBestSATDCost())
#else
    if (geocandCostList[i] > MRG_FAST_RATIO * geocandCostList[0] || geocandCostList[i] > getMergeBestSATDCost() || geocandCostList[i] > getAFFBestSATDCost())
#endif
    {
      geoNumMrgSATDCand = i;
      break;
    }
  }

  for (uint8_t i = 0; i < geoNumMrgSATDCand && isChromaEnabled(pu.chromaFormat); i++)
  {
    uint8_t candidateIdx = geoRdModeList[i];
    int splitDir = geoSplitDirList[candidateIdx];
    int mergeCand0 = geoMergeCand0[candidateIdx];
    int mergeCand1 = geoMergeCand1[candidateIdx];

    if (!isGeoChromaAvail[mergeCand0])
    {
#if TM_MRG
      mergeCtx[GEO_TM_OFF].setMergeInfo(pu, mergeCand0);
#else
      mergeCtx.setMergeInfo(pu, mergeCand0);
#endif
      m_pcInterSearch->motionCompensation(pu, geoBuffer[mergeCand0], REF_PIC_LIST_X, false, true);
      isGeoChromaAvail[mergeCand0] = true;
    }

    if (!isGeoChromaAvail[mergeCand1])
    {
#if TM_MRG
      mergeCtx[GEO_TM_OFF].setMergeInfo(pu, mergeCand1);
#else
      mergeCtx.setMergeInfo(pu, mergeCand1);
#endif
      m_pcInterSearch->motionCompensation(pu, geoBuffer[mergeCand1], REF_PIC_LIST_X, false, true);
      isGeoChromaAvail[mergeCand1] = true;
    }

    geoCombinations[candidateIdx] = m_acGeoWeightedBuffer[candidateIdx].getBuf(localUnitArea);
    m_pcInterSearch->weightedGeoBlk(pu, splitDir, CHANNEL_TYPE_CHROMA, geoCombinations[candidateIdx], geoBuffer[mergeCand0], geoBuffer[mergeCand1]);
  }

  bool geocandHasNoResidual[GEO_MAX_TRY_WEIGHTED_SAD];
  bool bestIsSkip = false;
  std::memset(geocandHasNoResidual, false, GEO_MAX_TRY_WEIGHTED_SAD * sizeof(bool));

  m_bestModeUpdated = tempCS->useDbCost = bestCS->useDbCost = false;
  tempCS->initStructData(encTestMode.qp);
  uint8_t iteration = 2, iterationBegin = 0;
  for (uint8_t noResidualPass = iterationBegin; noResidualPass < iteration; ++noResidualPass)
  {
    for (uint8_t mrgHADIdx = 0; mrgHADIdx < geoNumMrgSATDCand; mrgHADIdx++)
    {
      uint8_t candidateIdx = geoRdModeList[mrgHADIdx];
      if (((noResidualPass != 0) && geocandHasNoResidual[candidateIdx])
        || ((noResidualPass == 0) && bestIsSkip))
      {
        continue;
      }
      CodingUnit &cu = tempCS->addCU(tempCS->area, pm.chType);
      pm.setCUData(cu);
      cu.predMode = MODE_INTER;
      cu.slice = tempCS->slice;
      cu.tileIdx = tempCS->pps->getTileIdx(tempCS->area.lumaPos());
      cu.qp = encTestMode.qp;
      cu.affine = false;
      cu.mtsFlag = false;
#if INTER_LIC
      cu.LICFlag = false;
#endif
      cu.BcwIdx = BCW_DEFAULT;
      cu.geoFlag = true;
      cu.imv = 0;
      cu.mmvdSkip = false;
      cu.skip = false;
      cu.mipFlag = false;
      cu.bdpcmMode = 0;
      PredictionUnit &pu = tempCS->addPU(cu, pm.chType);
      pu.mergeFlag = true;
      pu.regularMergeFlag = false;
      pu.geoSplitDir = geoSplitDirList[candidateIdx];
      pu.geoMergeIdx0 = geoMergeCand0[candidateIdx];
      pu.geoMergeIdx1 = geoMergeCand1[candidateIdx];
#if TM_MRG
      pu.tmMergeFlag = false;
      pu.geoTmFlag0 = false;
      pu.geoTmFlag1 = false;
#endif
      pu.geoMMVDFlag0 = false;
      pu.geoMMVDFlag1 = false;

      pu.mmvdMergeFlag = false;
      pu.mmvdMergeIdx = MAX_UCHAR;
#if TM_MRG
      MergeCtx *mergeTmCtx0 = nullptr;
      MergeCtx *mergeTmCtx1 = nullptr;
      PU::spanGeoMMVDMotionInfo(pu, mergeCtx[GEO_TM_OFF], *mergeTmCtx0, *mergeTmCtx1, pu.geoSplitDir, pu.geoMergeIdx0, pu.geoMergeIdx1, pu.geoTmFlag0, pu.geoMMVDFlag0, pu.geoMMVDIdx0, pu.geoTmFlag1, pu.geoMMVDFlag1, pu.geoMMVDIdx1);
#else
      PU::spanGeoMMVDMotionInfo(pu, mergeCtx, pu.geoSplitDir, pu.geoMergeIdx0, pu.geoMergeIdx1, pu.geoMMVDFlag0, pu.geoMMVDIdx0, pu.geoMMVDFlag1, pu.geoMMVDIdx1);
#endif
      tempCS->getPredBuf().copyFrom(geoCombinations[candidateIdx]);
#if ENABLE_OBMC
      cu.isobmcMC = true;
      cu.obmcFlag = true;
      m_pcInterSearch->subBlockOBMC(pu);
      cu.isobmcMC = false;
#endif

      xEncodeInterResidual(tempCS, bestCS, pm, encTestMode, noResidualPass, (noResidualPass == 0 ? &geocandHasNoResidual[candidateIdx] : NULL));

      if (m_pcEncCfg->getUseFastDecisionForMerge() && !bestIsSkip)
      {
        bestIsSkip = bestCS->getCU(pm.chType)->rootCbf == 0;
      }
      tempCS->initStructData(encTestMode.qp);
    }
  }

  CodingUnit *bestCU = bestCS->getCU(CHANNEL_TYPE_LUMA);
  bool skipGPMMMVD = false;
  if (geoNumMrgSATDCand > 0)
  {
    if (bestCU->skip && !bestCU->geoFlag && !bestCU->affine && !bestCU->mmvdSkip && !bestCU->firstPU->mmvdMergeFlag)
    {
      skipGPMMMVD = true;
    }
    else if (bestCU->affine && bestCU->skip && (bestCU->lwidth() >= 16 || bestCU->lheight() >= 16))
    {
      skipGPMMMVD = true;
    }
  }

  bool isBaseMergeCandIncluded[GEO_MAX_NUM_UNI_CANDS];
  std::memset(isBaseMergeCandIncluded, false, GEO_MAX_NUM_UNI_CANDS * sizeof(bool));
  bool isGPMModeIncludedForMMVD[GEO_NUM_PARTITION_MODE];
  std::memset(isGPMModeIncludedForMMVD, false, GEO_NUM_PARTITION_MODE * sizeof(bool));

  if (!skipGPMMMVD)
  {
    skipGPMMMVD = (selGeoModeRDList[0] > (bestNormalMrgCost * 1.1));
  }

  if (!skipGPMMMVD)
  {
    if (isSecondPass)
    {
      for (int i = 0; i < relatedCU.numGeoDirCand; i++)
      {
        isGPMModeIncludedForMMVD[relatedCU.geoDirCandList[i]] = true;
        if (fastGpmMmvdRelatedCU)
        {
          isBaseMergeCandIncluded[relatedCU.geoMrgIdx0List[i]] = true;
          isBaseMergeCandIncluded[relatedCU.geoMrgIdx1List[i]] = true;
        }
      }
      if (!fastGpmMmvdRelatedCU)
      {
        std::memset(isBaseMergeCandIncluded, true, GEO_MAX_NUM_UNI_CANDS * sizeof(bool));
      }
    }
    else
    {
      double dirCostThresh = (selGeoModeRDList[0] * 1.2);
      isGPMModeIncludedForMMVD[selGeoModeList[0]] = true;
      isBaseMergeCandIncluded[mergeCandList0[selGeoModeList[0]][0]] = true;
      isBaseMergeCandIncluded[mergeCandList1[selGeoModeList[0]][0]] = true;

      for (int i = 1; i < maxNumGPMDirFirstPass; i++)
      {
        if (selGeoModeRDList[i] > dirCostThresh)
        {
          break;
        }
        else
        {
          isGPMModeIncludedForMMVD[selGeoModeList[i]] = true;
          isBaseMergeCandIncluded[mergeCandList0[selGeoModeList[i]][0]] = true;
          isBaseMergeCandIncluded[mergeCandList1[selGeoModeList[i]][0]] = true;
        }
      }
      if (includeMoreMMVDCandFirstPass)
      {
        int num = 0;
        // add more cands from best combo results obtained in weighted blended nonmmvd combo
        num = min((int)geocandCostList.size(), GEO_MAX_TRY_WEIGHTED_SATD);
        for (int i = 0; i < num; i++)
        {
          if (fastGpmMmvdSearch && (geocandCostList[i] > dirCostThresh))
          {
            break;
          }
          else
          {
            isGPMModeIncludedForMMVD[geoPartitionModeList[i]] = true;
            isBaseMergeCandIncluded[geoMergeCand0[geoRdModeList[i]]] = true;
            isBaseMergeCandIncluded[geoMergeCand1[geoRdModeList[i]]] = true;
          }
        }
      }
    }
  }

  if (!skipGPMMMVD)
  {
    CodingUnit &cu = tempCS->addCU(tempCS->area, pm.chType);
    pm.setCUData(cu);
    cu.predMode = MODE_INTER;
    cu.slice = tempCS->slice;
    cu.tileIdx = tempCS->pps->getTileIdx(tempCS->area.lumaPos());
    cu.qp = encTestMode.qp;
    cu.affine = false;
    cu.mtsFlag = false;
#if INTER_LIC
    cu.LICFlag = false;
#endif
    cu.BcwIdx = BCW_DEFAULT;
    cu.geoFlag = true;
    cu.imv = 0;
    cu.mmvdSkip = false;
    cu.skip = false;
    cu.mipFlag = false;
    cu.bdpcmMode = 0;

    PredictionUnit &pu = tempCS->addPU(cu, pm.chType);
    pu.mergeFlag = true;
    pu.regularMergeFlag = false;
#if TM_MRG
    pu.tmMergeFlag = false;
#endif
    bool simpleGPMMMVDStep = (m_pcEncCfg->getIntraPeriod() == -1);
    double mmvdMrgCost[GEO_MAX_NUM_UNI_CANDS][GPM_EXT_MMVD_MAX_REFINE_NUM];
    for (uint8_t mergeCand = 0; mergeCand < maxNumMergeCandidates; mergeCand++)
    {
      if (mrgDuplicated[mergeCand])
      {
        continue;
      }
      if (!isBaseMergeCandIncluded[mergeCand])
      {
        continue;
      }
      for (uint8_t mmvdCand = 0; mmvdCand < (extMMVD ? GPM_EXT_MMVD_MAX_REFINE_NUM : GPM_MMVD_MAX_REFINE_NUM); mmvdCand++)
      {
        if (simpleGPMMMVDStep)
        {
          int mmvdStep = (extMMVD ? (mmvdCand >> 3) : (mmvdCand >> 2));
          if (mmvdStep >= 5 && (!fastGpmMmvdSearch || (fastGpmMmvdSearch && !isSecondPass)))
          {
            continue;
          }
        }
        geoMMVDBuf[mergeCand][mmvdCand] = m_acGeoMMVDBuffer[mergeCand][mmvdCand].getBuf(localUnitArea);
#if TM_MRG
        mergeCtx[GEO_TM_OFF].setGeoMmvdMergeInfo(pu, mergeCand, mmvdCand);
#else
        mergeCtx.setGeoMmvdMergeInfo(pu, mergeCand, mmvdCand);
#endif
        if (m_pcEncCfg->getMCTSEncConstraint() && (!(MCTSHelper::checkMvBufferForMCTSConstraint(pu))))
        {
          tempCS->initStructData(encTestMode.qp);
          return;
        }
        m_pcInterSearch->motionCompensation(pu, geoMMVDBuf[mergeCand][mmvdCand], REF_PIC_LIST_X, true, false);
        geoMMVDTempBuf[mergeCand][mmvdCand] = m_acGeoMMVDTmpBuffer[mergeCand][mmvdCand].getBuf(localUnitArea);
        geoMMVDTempBuf[mergeCand][mmvdCand].Y().copyFrom(geoMMVDBuf[mergeCand][mmvdCand].Y());
        geoMMVDTempBuf[mergeCand][mmvdCand].Y().roundToOutputBitdepth(geoMMVDTempBuf[mergeCand][mmvdCand].Y(), cu.slice->clpRng(COMPONENT_Y));
        distParamWholeBlk.cur.buf = geoMMVDTempBuf[mergeCand][mmvdCand].Y().buf;
        distParamWholeBlk.cur.stride = geoMMVDTempBuf[mergeCand][mmvdCand].Y().stride;
        sadMMVDWholeBlk[mergeCand][mmvdCand] = distParamWholeBlk.distFunc(distParamWholeBlk);
        mmvdMrgCost[mergeCand][mmvdCand] = sadMMVDWholeBlk[mergeCand][mmvdCand] + geoMergeIdxCost[mergeCand] + geoMMVDFlagCost[1] + geoMMVDIdxCost[mmvdCand];
        if (mmvdMrgCost[mergeCand][mmvdCand] < bestMrgCost)
        {
          bestMrgCost = mmvdMrgCost[mergeCand][mmvdCand];
        }
      }
    }

    double mrgCostThres = (bestMrgCost * 3.0);
    for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; splitDir++)
    {
      if (!isGPMModeIncludedForMMVD[splitDir])
      {
        continue;
      }
      int16_t angle = g_GeoParams[splitDir][0];
      if (g_angle2mirror[angle] == 2)
      {
        stepX = 1;
        maskStride = -GEO_WEIGHT_MASK_SIZE;
        maskStride2 = -(int)cu.lwidth();
        SADmask = &g_globalGeoEncSADmask[g_angle2mask[g_GeoParams[splitDir][0]]][(GEO_WEIGHT_MASK_SIZE - 1 - g_weightOffset[splitDir][hIdx][wIdx][1]) * GEO_WEIGHT_MASK_SIZE + g_weightOffset[splitDir][hIdx][wIdx][0]];
      }
      else if (g_angle2mirror[angle] == 1)
      {
        stepX = -1;
        maskStride2 = cu.lwidth();
        maskStride = GEO_WEIGHT_MASK_SIZE;
        SADmask = &g_globalGeoEncSADmask[g_angle2mask[g_GeoParams[splitDir][0]]][g_weightOffset[splitDir][hIdx][wIdx][1] * GEO_WEIGHT_MASK_SIZE + (GEO_WEIGHT_MASK_SIZE - 1 - g_weightOffset[splitDir][hIdx][wIdx][0])];
      }
      else
      {
        stepX = 1;
        maskStride = GEO_WEIGHT_MASK_SIZE;
        maskStride2 = -(int)cu.lwidth();
        SADmask = &g_globalGeoEncSADmask[g_angle2mask[g_GeoParams[splitDir][0]]][g_weightOffset[splitDir][hIdx][wIdx][1] * GEO_WEIGHT_MASK_SIZE + g_weightOffset[splitDir][hIdx][wIdx][0]];
      }
      for (uint8_t mergeCand = 0; mergeCand < maxNumMergeCandidates; mergeCand++)
      {
        if (mrgDuplicated[mergeCand])
        {
          continue;
        }
        if (!isBaseMergeCandIncluded[mergeCand])
        {
          continue;
        }
        for (uint8_t mmvdCand = 0; mmvdCand < (extMMVD ? GPM_EXT_MMVD_MAX_REFINE_NUM : GPM_MMVD_MAX_REFINE_NUM); mmvdCand++)
        {
          if (simpleGPMMMVDStep)
          {
            int mmvdStep = (extMMVD ? (mmvdCand >> 3) : (mmvdCand >> 2));
            if (mmvdStep >= 5 && (!fastGpmMmvdSearch || (fastGpmMmvdSearch && !isSecondPass)))
            {
              continue;
            }
          }
          if (mmvdMrgCost[mergeCand][mmvdCand] > mrgCostThres)
          {
            continue;
          }
          m_pcRdCost->setDistParam(distParam, tempCS->getOrgBuf().Y(), geoMMVDTempBuf[mergeCand][mmvdCand].Y().buf, geoMMVDTempBuf[mergeCand][mmvdCand].Y().stride, SADmask, maskStride, stepX, maskStride2, sps.getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y);
          sadLarge = distParam.distFunc(distParam);
          double tempCost = (double)sadLarge + geoMergeIdxCost[mergeCand] + geoMMVDFlagCost[1] + geoMMVDIdxCost[mmvdCand];
          m_GeoMMVDCostList.insert(splitDir, 0, mergeCand, (mmvdCand + 1), tempCost);
          sortCandList(tempCost, mergeCand, (mmvdCand + 1), sadCostList0[splitDir], mergeCandList0[splitDir], mmvdCandList0[splitDir], numCandPerPar);

          sadSmall = sadMMVDWholeBlk[mergeCand][mmvdCand] - sadLarge;
          tempCost = (double)sadSmall + geoMergeIdxCost[mergeCand] + geoMMVDFlagCost[1] + geoMMVDIdxCost[mmvdCand];
          m_GeoMMVDCostList.insert(splitDir, 1, mergeCand, (mmvdCand + 1), tempCost);
          sortCandList(tempCost, mergeCand, (mmvdCand + 1), sadCostList1[splitDir], mergeCandList1[splitDir], mmvdCandList1[splitDir], numCandPerPar);
        }
      }
    }

    for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; splitDir++)
    {
      int numCandPart0 = min(numCandPerPar, (int)mergeCandList0[splitDir].size());
      int numCandPart1 = min(numCandPerPar, (int)mergeCandList1[splitDir].size());
      for (int candIdx0 = 0; candIdx0 < numCandPart0; candIdx0++)
      {
        for (int candIdx1 = 0; candIdx1 < numCandPart1; candIdx1++)
        {
          int mergeCand0 = mergeCandList0[splitDir][candIdx0];
          int mergeCand1 = mergeCandList1[splitDir][candIdx1];
          int mmvdCand0 = mmvdCandList0[splitDir][candIdx0];
          int mmvdCand1 = mmvdCandList1[splitDir][candIdx1];
#if TM_MRG
          bool geoTmFlag0 = (mmvdCand0 == (GPM_EXT_MMVD_MAX_REFINE_NUM + 1));
          bool geoTmFlag1 = (mmvdCand1 == (GPM_EXT_MMVD_MAX_REFINE_NUM + 1));
          CHECK(geoTmFlag0 || geoTmFlag1, "GPM TM has not been tested by far");
#endif
          if ((mmvdCand0 == 0) && (mmvdCand1 == 0))
          {
            continue;
          }
          if ((mmvdCand0 == mmvdCand1) && (mmvdCand0 > 0))
          {
            if (mergeCand0 == mergeCand1)
            {
              continue;
            }
          }

          double tempCost = m_GeoMMVDCostList.singleDistList[0][splitDir][mergeCand0][mmvdCand0].cost + m_GeoMMVDCostList.singleDistList[1][splitDir][mergeCand1][mmvdCand1].cost;
          tempCost = tempCost + geoModeCost[splitDir];
#if TM_MRG
          if (sps.getUseDMVDMode())
          {
            tempCost += geoTMFlagCost[0];
          }
#endif
          updateGeoMMVDCandList(tempCost, splitDir, mergeCand0, mergeCand1, mmvdCand0, mmvdCand1,
            geoSADCostList, geoSplitDirList, geoMergeCand0, geoMergeCand1, geoMmvdCand0, geoMmvdCand1, numSATDCands);
        }
      }
    }

#if TM_MRG
    uint8_t maxNumTmMrgCand = maxNumMergeCandidates;
    PelUnitBuf geoTmBuffer[GEO_TM_MAX_NUM_CANDS];
    PelUnitBuf geoTmTempBuf[GEO_TM_MAX_NUM_CANDS];
    if (sps.getUseDMVDMode())
    {
      for (int i = GEO_TM_SHAPE_AL; i < GEO_NUM_TM_MV_CAND; i++)
      {
        mergeCtx[i].numValidMergeCand = maxNumTmMrgCand;
        for (int idx = 0; idx < maxNumTmMrgCand; idx++)
        {
          if (mrgDuplicated[idx])
          {
            continue;
          }
          mergeCtx[i].BcwIdx[idx] = BCW_DEFAULT;
          mergeCtx[i].useAltHpelIf[idx] = false;
#if INTER_LIC
          mergeCtx[i].LICFlags[idx] = false;
#endif
          mergeCtx[i].interDirNeighbours[idx] = mergeCtx[GEO_TM_OFF].interDirNeighbours[idx];
          mergeCtx[i].mvFieldNeighbours[(idx << 1)].mv = mergeCtx[GEO_TM_OFF].mvFieldNeighbours[(idx << 1)].mv;
          mergeCtx[i].mvFieldNeighbours[(idx << 1) + 1].mv = mergeCtx[GEO_TM_OFF].mvFieldNeighbours[(idx << 1) + 1].mv;
          mergeCtx[i].mvFieldNeighbours[(idx << 1)].refIdx = mergeCtx[GEO_TM_OFF].mvFieldNeighbours[(idx << 1)].refIdx;
          mergeCtx[i].mvFieldNeighbours[(idx << 1) + 1].refIdx = mergeCtx[GEO_TM_OFF].mvFieldNeighbours[(idx << 1) + 1].refIdx;
        }
      }

      pu.tmMergeFlag = true;
      Distortion sadTmWholeBlk[GEO_TM_MAX_NUM_CANDS];
      for (uint8_t tmType = GEO_TM_SHAPE_AL; tmType < GEO_NUM_TM_MV_CAND; tmType++)
      {
        pu.geoTmType = tmType;
        for (uint8_t mrgIdx = 0; mrgIdx < maxNumTmMrgCand; mrgIdx++)
        {
          if (mrgDuplicated[mrgIdx])
          {
            continue;
          }
          uint8_t mergeCand = mrgIdx + (tmType - 1) * GEO_MAX_NUM_UNI_CANDS;
          mergeCtx[tmType].setMergeInfo(pu, mrgIdx);
          m_pcInterSearch->deriveTMMv(pu);
          mergeCtx[tmType].mvFieldNeighbours[(mrgIdx << 1)].mv = pu.mv[0];
          mergeCtx[tmType].mvFieldNeighbours[(mrgIdx << 1) + 1].mv = pu.mv[1];

          geoTmBuffer[mergeCand] = m_acGeoMergeTmpBuffer[mergeCand].getBuf(localUnitArea);
          m_pcInterSearch->motionCompensation(pu, geoTmBuffer[mergeCand]);

          // calculate SAD for each candidate
          geoTmTempBuf[mergeCand] = m_acGeoSADTmpBuffer[mergeCand].getBuf(localUnitArea);
          geoTmTempBuf[mergeCand].Y().copyFrom(geoTmBuffer[mergeCand].Y());
          geoTmTempBuf[mergeCand].Y().roundToOutputBitdepth(geoTmTempBuf[mergeCand].Y(), cu.slice->clpRng(COMPONENT_Y));
          distParamWholeBlk.cur.buf = geoTmTempBuf[mergeCand].Y().buf;
          distParamWholeBlk.cur.stride = geoTmTempBuf[mergeCand].Y().stride;
          sadTmWholeBlk[mergeCand] = distParamWholeBlk.distFunc(distParamWholeBlk);
        }
      }
      pu.tmMergeFlag = false;

      for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; splitDir++)
      {
        int16_t angle = g_GeoParams[splitDir][0];
        if (g_angle2mirror[angle] == 2)
        {
          stepX = 1;
          maskStride = -GEO_WEIGHT_MASK_SIZE;
          maskStride2 = -(int)cu.lwidth();
          SADmask = &g_globalGeoEncSADmask[g_angle2mask[g_GeoParams[splitDir][0]]][(GEO_WEIGHT_MASK_SIZE - 1 - g_weightOffset[splitDir][hIdx][wIdx][1]) * GEO_WEIGHT_MASK_SIZE + g_weightOffset[splitDir][hIdx][wIdx][0]];
        }
        else if (g_angle2mirror[angle] == 1)
        {
          stepX = -1;
          maskStride2 = cu.lwidth();
          maskStride = GEO_WEIGHT_MASK_SIZE;
          SADmask = &g_globalGeoEncSADmask[g_angle2mask[g_GeoParams[splitDir][0]]][g_weightOffset[splitDir][hIdx][wIdx][1] * GEO_WEIGHT_MASK_SIZE + (GEO_WEIGHT_MASK_SIZE - 1 - g_weightOffset[splitDir][hIdx][wIdx][0])];
        }
        else
        {
          stepX = 1;
          maskStride = GEO_WEIGHT_MASK_SIZE;
          maskStride2 = -(int)cu.lwidth();
          SADmask = &g_globalGeoEncSADmask[g_angle2mask[g_GeoParams[splitDir][0]]][g_weightOffset[splitDir][hIdx][wIdx][1] * GEO_WEIGHT_MASK_SIZE + g_weightOffset[splitDir][hIdx][wIdx][0]];
        }
        for (uint8_t mergeCand = 0; mergeCand < maxNumTmMrgCand; mergeCand++)
        {
          if (mrgDuplicated[mergeCand])
          {
            continue;
          }

          uint8_t mergeCand0 = mergeCand + (g_geoTmShape[0][g_GeoParams[splitDir][0]] - 1) * GEO_MAX_NUM_UNI_CANDS;
          m_pcRdCost->setDistParam(distParam, tempCS->getOrgBuf().Y(), geoTmTempBuf[mergeCand0].Y().buf, geoTmTempBuf[mergeCand0].Y().stride, SADmask, maskStride, stepX, maskStride2, sps.getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y);
          sadLarge = distParam.distFunc(distParam);
          double tempCost = (double)sadLarge + geoMergeIdxCost[mergeCand] + geoMMVDFlagCost[0];
          m_GeoMMVDCostList.insert(splitDir, 0, mergeCand, (GPM_EXT_MMVD_MAX_REFINE_NUM + 1), tempCost);

          uint8_t mergeCand1 = mergeCand + (g_geoTmShape[1][g_GeoParams[splitDir][0]] - 1) * GEO_MAX_NUM_UNI_CANDS;
          m_pcRdCost->setDistParam(distParam, tempCS->getOrgBuf().Y(), geoTmTempBuf[mergeCand1].Y().buf, geoTmTempBuf[mergeCand1].Y().stride, SADmask, maskStride, stepX, maskStride2, sps.getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y);
          sadSmall = sadTmWholeBlk[mergeCand1] - distParam.distFunc(distParam);
          tempCost = (double)sadSmall + geoMergeIdxCost[mergeCand] + geoMMVDFlagCost[0];
          m_GeoMMVDCostList.insert(splitDir, 1, mergeCand, (GPM_EXT_MMVD_MAX_REFINE_NUM + 1), tempCost);
        }
      }

      for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; splitDir++)
      {
        for (int mergeCand0 = 0; mergeCand0 < maxNumTmMrgCand; mergeCand0++)
        {
          if (mrgDuplicated[mergeCand0])
          {
            continue;
          }
          for (int mergeCand1 = 0; mergeCand1 < maxNumTmMrgCand; mergeCand1++)
          {
            if (mrgDuplicated[mergeCand1])
            {
              continue;
            }
            if (mergeCand0 == mergeCand1)
            {
              continue;
            }
            double tempCost = m_GeoMMVDCostList.singleDistList[0][splitDir][mergeCand0][GPM_EXT_MMVD_MAX_REFINE_NUM + 1].cost + m_GeoMMVDCostList.singleDistList[1][splitDir][mergeCand1][GPM_EXT_MMVD_MAX_REFINE_NUM + 1].cost;
            tempCost = tempCost + geoModeCost[splitDir] + geoTMFlagCost[1];
            updateGeoMMVDCandList(tempCost, splitDir, mergeCand0, mergeCand1, (GPM_EXT_MMVD_MAX_REFINE_NUM + 1), (GPM_EXT_MMVD_MAX_REFINE_NUM + 1),
              geoSADCostList, geoSplitDirList, geoMergeCand0, geoMergeCand1, geoMmvdCand0, geoMmvdCand1, numSATDCands);
          }
        }
      }
    }
#endif

    int numberGeoCandChecked = (int)geoSADCostList.size();
    if (numberGeoCandChecked == 0)
    {
      return;
    }

    geoNumMrgSATDCand = min(GEO_MAX_TRY_WEIGHTED_SATD, numberGeoCandChecked);
    if (geoRdModeList.size() > geoNumMrgSATDCand)
    {
      geoRdModeList.resize(geoNumMrgSATDCand);
      isNonMMVDListIdx.resize(geoNumMrgSATDCand);
      geoPartitionModeList.resize(geoNumMrgSATDCand);
      geocandCostList.resize(geoNumMrgSATDCand);
    }
    for (uint8_t candidateIdx = 0; candidateIdx < numberGeoCandChecked; candidateIdx++)
    {
      int splitDir = geoSplitDirList[candidateIdx];
      int mergeCand0 = geoMergeCand0[candidateIdx];
      int mergeCand1 = geoMergeCand1[candidateIdx];
#if TM_MRG
      bool tmFlag0 = (geoMmvdCand0[candidateIdx] == (GPM_EXT_MMVD_MAX_REFINE_NUM + 1));
      bool tmFlag1 = (geoMmvdCand1[candidateIdx] == (GPM_EXT_MMVD_MAX_REFINE_NUM + 1));
      bool mmvdFlag0 = (geoMmvdCand0[candidateIdx] >= 1 && geoMmvdCand0[candidateIdx] <= GPM_EXT_MMVD_MAX_REFINE_NUM);
      bool mmvdFlag1 = (geoMmvdCand1[candidateIdx] >= 1 && geoMmvdCand1[candidateIdx] <= GPM_EXT_MMVD_MAX_REFINE_NUM);
      int  mmvdCand0 = (mmvdFlag0 ? (geoMmvdCand0[candidateIdx] - 1) : MAX_INT);
      int  mmvdCand1 = (mmvdFlag1 ? (geoMmvdCand1[candidateIdx] - 1) : MAX_INT);
      CHECK(tmFlag0 != tmFlag1, "TM flag cannot be enabled/disabled for two partitions separately");

      if (!tmFlag0 && !tmFlag1 && !mmvdFlag0 && !mmvdFlag1)
      {
        continue;
      }
      if (tmFlag0 && mergeCand0 == mergeCand1)
      {
        continue;
      }
#else
      int mmvdCand0 = geoMmvdCand0[candidateIdx] - 1;
      int mmvdCand1 = geoMmvdCand1[candidateIdx] - 1;
      bool mmvdFlag0 = (mmvdCand0 >= 0);
      bool mmvdFlag1 = (mmvdCand1 >= 0);

      if (!mmvdFlag0 && !mmvdFlag1)
      {
        continue;
      }
#endif

      geoCombinations[candidateIdx] = m_acGeoWeightedBuffer[candidateIdx].getBuf(localUnitArea);
#if TM_MRG
      PelUnitBuf predSrc0, predSrc1;
      if (tmFlag0)
      {
        int mrgTmCand0 = mergeCand0 + (g_geoTmShape[0][g_GeoParams[splitDir][0]] - 1) * GEO_MAX_NUM_UNI_CANDS;
        predSrc0 = geoTmBuffer[mrgTmCand0];
      }
      else if (mmvdFlag0)
      {
        predSrc0 = geoMMVDBuf[mergeCand0][mmvdCand0];
      }
      else
      {
        predSrc0 = geoBuffer[mergeCand0];
      }

      if (tmFlag1)
      {
        int mrgTmCand1 = mergeCand1 + (g_geoTmShape[1][g_GeoParams[splitDir][0]] - 1) * GEO_MAX_NUM_UNI_CANDS;
        predSrc1 = geoTmBuffer[mrgTmCand1];
      }
      else if (mmvdFlag1)
      {
        predSrc1 = geoMMVDBuf[mergeCand1][mmvdCand1];
      }
      else
      {
        predSrc1 = geoBuffer[mergeCand1];
      }

      m_pcInterSearch->weightedGeoBlk(pu, splitDir, CHANNEL_TYPE_LUMA, geoCombinations[candidateIdx], predSrc0, predSrc1);
#else
      m_pcInterSearch->weightedGeoBlk(pu, splitDir, CHANNEL_TYPE_LUMA, geoCombinations[candidateIdx], (mmvdFlag0 ? geoMMVDBuf[mergeCand0][mmvdCand0] : geoBuffer[mergeCand0]), (mmvdFlag1 ? geoMMVDBuf[mergeCand1][mmvdCand1] : geoBuffer[mergeCand1]));
#endif
      distParamSAD2.cur = geoCombinations[candidateIdx].Y();
      Distortion sad = distParamSAD2.distFunc(distParamSAD2);

      double updateCost = geoModeCost[splitDir] + geoMergeIdxCost[mergeCand0] + geoMergeIdxCost[mergeCand1] + geoMMVDFlagCost[mmvdFlag0] + geoMMVDFlagCost[mmvdFlag1];
      if (mmvdFlag0)
      {
        updateCost += geoMMVDIdxCost[mmvdCand0];
      }
      if (mmvdFlag1)
      {
        updateCost += geoMMVDIdxCost[mmvdCand1];
      }
#if TM_MRG
      if (sps.getUseDMVDMode())
      {
        if (!mmvdFlag0 && !mmvdFlag1)
        {
          updateCost += geoTMFlagCost[tmFlag0];
        }
      }
#endif
      updateCost += (double)sad;
      orderCandList(candidateIdx, false, splitDir, updateCost, geoRdModeList, isNonMMVDListIdx, geoPartitionModeList, geocandCostList, geoNumMrgSATDCand);
    }

    if (fastGpmMmvdRelatedCU)
    {
      int cnt = 0;
      for (uint8_t i = 0; i < geoNumMrgSATDCand; i++)
      {
        if (isNonMMVDListIdx[i] == false)
        {
          relatedCU.geoDirCandList[cnt] = geoPartitionModeList[i];
          relatedCU.geoMrgIdx0List[cnt] = geoMergeCand0[geoRdModeList[i]];
          relatedCU.geoMrgIdx1List[cnt] = geoMergeCand1[geoRdModeList[i]];
          cnt++;
        }
      }
      relatedCU.numGeoDirCand = cnt;
    }
    else
    {
      relatedCU.numGeoDirCand = geoNumMrgSATDCand;
      for (uint8_t i = 0; i < geoNumMrgSATDCand; i++)
      {
        relatedCU.geoDirCandList[i] = geoPartitionModeList[i];
      }
    }

    for (uint8_t i = 1; i < geoNumMrgSATDCand; i++)
    {
#if MERGE_ENC_OPT
      if (geocandCostList[i] > MRG_FAST_RATIO * geocandCostList[0] || geocandCostList[i] > getMergeBestSATDCost())
#else
      if (geocandCostList[i] > MRG_FAST_RATIO * geocandCostList[0] || geocandCostList[i] > getMergeBestSATDCost() || geocandCostList[i] > getAFFBestSATDCost())
#endif
      {
        geoNumMrgSATDCand = i;
        break;
      }
    }
    for (uint8_t i = 0; i < geoNumMrgSATDCand && isChromaEnabled(pu.chromaFormat); i++)
    {
      if (isNonMMVDListIdx[i])
      {
        continue;
      }
      uint8_t candidateIdx = geoRdModeList[i];
      int splitDir = geoSplitDirList[candidateIdx];
      int mergeCand0 = geoMergeCand0[candidateIdx];
      int mergeCand1 = geoMergeCand1[candidateIdx];
#if TM_MRG
      bool tmFlag0 = (geoMmvdCand0[candidateIdx] == (GPM_EXT_MMVD_MAX_REFINE_NUM + 1));
      bool tmFlag1 = (geoMmvdCand1[candidateIdx] == (GPM_EXT_MMVD_MAX_REFINE_NUM + 1));
      bool mmvdFlag0 = (geoMmvdCand0[candidateIdx] >= 1 && geoMmvdCand0[candidateIdx] <= GPM_EXT_MMVD_MAX_REFINE_NUM);
      bool mmvdFlag1 = (geoMmvdCand1[candidateIdx] >= 1 && geoMmvdCand1[candidateIdx] <= GPM_EXT_MMVD_MAX_REFINE_NUM);
      int  mmvdCand0 = (mmvdFlag0 ? (geoMmvdCand0[candidateIdx] - 1) : MAX_INT);
      int  mmvdCand1 = (mmvdFlag1 ? (geoMmvdCand1[candidateIdx] - 1) : MAX_INT);
#else
      int mmvdCand0 = geoMmvdCand0[candidateIdx] - 1;
      int mmvdCand1 = geoMmvdCand1[candidateIdx] - 1;

      bool mmvdFlag0 = (mmvdCand0 >= 0);
      bool mmvdFlag1 = (mmvdCand1 >= 0);
#endif

#if TM_MRG
      int mrgTmCand0 = MAX_INT, mrgTmCand1 = MAX_INT;
      if (tmFlag0)
      {
        int geoTmType = g_geoTmShape[0][g_GeoParams[splitDir][0]];
        mrgTmCand0 = mergeCand0 + (geoTmType - 1) * GEO_MAX_NUM_UNI_CANDS;
        if (!isGeoTmChromaAvail[mrgTmCand0])
        {
          mergeCtx[geoTmType].setMergeInfo(pu, mergeCand0);
          m_pcInterSearch->motionCompensation(pu, geoTmBuffer[mrgTmCand0], REF_PIC_LIST_X, false, true);
          isGeoTmChromaAvail[mrgTmCand0] = true;
        }
      }
      else
#endif
        if (mmvdFlag0)
        {
          if (!isGeoMMVDChromaAvail[mergeCand0][mmvdCand0])
          {
#if TM_MRG
            mergeCtx[GEO_TM_OFF].setGeoMmvdMergeInfo(pu, mergeCand0, mmvdCand0);
#else
            mergeCtx.setGeoMmvdMergeInfo(pu, mergeCand0, mmvdCand0);
#endif
            m_pcInterSearch->motionCompensation(pu, geoMMVDBuf[mergeCand0][mmvdCand0], REF_PIC_LIST_X, false, true);
            isGeoMMVDChromaAvail[mergeCand0][mmvdCand0] = true;
          }
        }
        else
        {
          if (!isGeoChromaAvail[mergeCand0])
          {
#if TM_MRG
            mergeCtx[GEO_TM_OFF].setMergeInfo(pu, mergeCand0);
#else
            mergeCtx.setMergeInfo(pu, mergeCand0);
#endif
            m_pcInterSearch->motionCompensation(pu, geoBuffer[mergeCand0], REF_PIC_LIST_X, false, true);
            isGeoChromaAvail[mergeCand0] = true;
          }
        }

#if TM_MRG
      if (tmFlag1)
      {
        int geoTmType = g_geoTmShape[1][g_GeoParams[splitDir][0]];
        mrgTmCand1 = mergeCand1 + (geoTmType - 1) * GEO_MAX_NUM_UNI_CANDS;
        if (!isGeoTmChromaAvail[mrgTmCand1])
        {
          mergeCtx[geoTmType].setMergeInfo(pu, mergeCand1);
          m_pcInterSearch->motionCompensation(pu, geoTmBuffer[mrgTmCand1], REF_PIC_LIST_X, false, true);
          isGeoTmChromaAvail[mrgTmCand1] = true;
        }
      }
      else
#endif
        if (mmvdFlag1)
        {
          if (!isGeoMMVDChromaAvail[mergeCand1][mmvdCand1])
          {
#if TM_MRG
            mergeCtx[GEO_TM_OFF].setGeoMmvdMergeInfo(pu, mergeCand1, mmvdCand1);
#else
            mergeCtx.setGeoMmvdMergeInfo(pu, mergeCand1, mmvdCand1);
#endif
            m_pcInterSearch->motionCompensation(pu, geoMMVDBuf[mergeCand1][mmvdCand1], REF_PIC_LIST_X, false, true);
            isGeoMMVDChromaAvail[mergeCand1][mmvdCand1] = true;
          }
        }
        else
        {
          if (!isGeoChromaAvail[mergeCand1])
          {
#if TM_MRG
            mergeCtx[GEO_TM_OFF].setMergeInfo(pu, mergeCand1);
#else
            mergeCtx.setMergeInfo(pu, mergeCand1);
#endif
            m_pcInterSearch->motionCompensation(pu, geoBuffer[mergeCand1], REF_PIC_LIST_X, false, true);
            isGeoChromaAvail[mergeCand1] = true;
          }
        }

      geoCombinations[candidateIdx] = m_acGeoWeightedBuffer[candidateIdx].getBuf(localUnitArea);
#if TM_MRG
      PelUnitBuf predSrc0, predSrc1;
      if (tmFlag0)
      {
        predSrc0 = geoTmBuffer[mrgTmCand0];
      }
      else if (mmvdFlag0)
      {
        predSrc0 = geoMMVDBuf[mergeCand0][mmvdCand0];
      }
      else
      {
        predSrc0 = geoBuffer[mergeCand0];
      }

      if (tmFlag1)
      {
        predSrc1 = geoTmBuffer[mrgTmCand1];
      }
      else if (mmvdFlag1)
      {
        predSrc1 = geoMMVDBuf[mergeCand1][mmvdCand1];
      }
      else
      {
        predSrc1 = geoBuffer[mergeCand1];
      }

      m_pcInterSearch->weightedGeoBlk(pu, splitDir, CHANNEL_TYPE_CHROMA, geoCombinations[candidateIdx], predSrc0, predSrc1);
#else
      m_pcInterSearch->weightedGeoBlk(pu, splitDir, CHANNEL_TYPE_CHROMA, geoCombinations[candidateIdx], (mmvdFlag0 ? geoMMVDBuf[mergeCand0][mmvdCand0] : geoBuffer[mergeCand0]), (mmvdFlag1 ? geoMMVDBuf[mergeCand1][mmvdCand1] : geoBuffer[mergeCand1]));
#endif
    }

    std::memset(geocandHasNoResidual, false, GEO_MAX_TRY_WEIGHTED_SAD * sizeof(bool));
    tempCS->initStructData(encTestMode.qp);
    uint8_t iteration = 2, iterationBegin = 0;
    for (uint8_t noResidualPass = iterationBegin; noResidualPass < iteration; ++noResidualPass)
    {
      for (uint8_t mrgHADIdx = 0; mrgHADIdx < geoNumMrgSATDCand; mrgHADIdx++)
      {
        if (isNonMMVDListIdx[mrgHADIdx])
        {
          continue;
        }
        uint8_t candidateIdx = geoRdModeList[mrgHADIdx];
        if (((noResidualPass != 0) && geocandHasNoResidual[candidateIdx])
          || ((noResidualPass == 0) && bestIsSkip))
        {
          continue;
        }
        CodingUnit &cu = tempCS->addCU(tempCS->area, pm.chType);
        pm.setCUData(cu);
        cu.predMode = MODE_INTER;
        cu.slice = tempCS->slice;
        cu.tileIdx = tempCS->pps->getTileIdx(tempCS->area.lumaPos());
        cu.qp = encTestMode.qp;
        cu.affine = false;
        cu.mtsFlag = false;
#if INTER_LIC
        cu.LICFlag = false;
#endif
        cu.BcwIdx = BCW_DEFAULT;
        cu.geoFlag = true;
        cu.imv = 0;
        cu.mmvdSkip = false;
        cu.skip = false;
        cu.mipFlag = false;
        cu.bdpcmMode = 0;
        PredictionUnit &pu = tempCS->addPU(cu, pm.chType);
        pu.mergeFlag = true;
        pu.regularMergeFlag = false;
        pu.geoSplitDir = geoSplitDirList[candidateIdx];
        pu.geoMergeIdx0 = geoMergeCand0[candidateIdx];
        pu.geoMergeIdx1 = geoMergeCand1[candidateIdx];
#if TM_MRG
        pu.geoTmFlag0 = (geoMmvdCand0[candidateIdx] == (GPM_EXT_MMVD_MAX_REFINE_NUM + 1));
        pu.geoTmFlag1 = (geoMmvdCand1[candidateIdx] == (GPM_EXT_MMVD_MAX_REFINE_NUM + 1));
        pu.geoMMVDFlag0 = (geoMmvdCand0[candidateIdx] >= 1 && geoMmvdCand0[candidateIdx] <= GPM_EXT_MMVD_MAX_REFINE_NUM);
        pu.geoMMVDFlag1 = (geoMmvdCand1[candidateIdx] >= 1 && geoMmvdCand1[candidateIdx] <= GPM_EXT_MMVD_MAX_REFINE_NUM);
        CHECK(pu.geoTmFlag0 != pu.geoTmFlag1, "TM flag cannot be enabled/disabled for two partitions separately");
        pu.tmMergeFlag = pu.geoTmFlag0;
#else
        pu.geoMMVDFlag0 = (geoMmvdCand0[candidateIdx] > 0);
        pu.geoMMVDFlag1 = (geoMmvdCand1[candidateIdx] > 0);
#endif

        if (pu.geoMMVDFlag0)
        {
          pu.geoMMVDIdx0 = (geoMmvdCand0[candidateIdx] - 1);
        }
        if (pu.geoMMVDFlag1)
        {
          pu.geoMMVDIdx1 = (geoMmvdCand1[candidateIdx] - 1);
        }
        pu.mmvdMergeFlag = false;
        pu.mmvdMergeIdx = MAX_UCHAR;
#if TM_MRG
        MergeCtx* mrgTmCtx0 = (pu.geoTmFlag0 == 0 ? nullptr : &mergeCtx[g_geoTmShape[0][g_GeoParams[pu.geoSplitDir][0]]]);
        MergeCtx* mrgTmCtx1 = (pu.geoTmFlag1 == 0 ? nullptr : &mergeCtx[g_geoTmShape[1][g_GeoParams[pu.geoSplitDir][0]]]);
        PU::spanGeoMMVDMotionInfo(pu, mergeCtx[GEO_TM_OFF], *mrgTmCtx0, *mrgTmCtx1, pu.geoSplitDir, pu.geoMergeIdx0, pu.geoMergeIdx1, pu.geoTmFlag0, pu.geoMMVDFlag0, pu.geoMMVDIdx0, pu.geoTmFlag1, pu.geoMMVDFlag1, pu.geoMMVDIdx1);
#else
        PU::spanGeoMMVDMotionInfo(pu, mergeCtx, pu.geoSplitDir, pu.geoMergeIdx0, pu.geoMergeIdx1, pu.geoMMVDFlag0, pu.geoMMVDIdx0, pu.geoMMVDFlag1, pu.geoMMVDIdx1);
#endif
        tempCS->getPredBuf().copyFrom(geoCombinations[candidateIdx]);
#if ENABLE_OBMC
        cu.isobmcMC = true;
        cu.obmcFlag = true;
        m_pcInterSearch->subBlockOBMC(pu);
        cu.isobmcMC = false;
#endif

        xEncodeInterResidual(tempCS, bestCS, pm, encTestMode, noResidualPass, (noResidualPass == 0 ? &geocandHasNoResidual[candidateIdx] : NULL));

        if (m_pcEncCfg->getUseFastDecisionForMerge() && !bestIsSkip)
        {
          bestIsSkip = bestCS->getCU(pm.chType)->rootCbf == 0;
        }
        tempCS->initStructData(encTestMode.qp);
      }
    }
  }
  if (m_bestModeUpdated && bestCS->cost != MAX_DOUBLE)
  {
    xCalDebCost(*bestCS, pm);
  }
}
#else
void EncCu::xCheckRDCostMergeGeo2Nx2N(CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &pm, const EncTestMode& encTestMode)
{
  const Slice &slice = *tempCS->slice;
  CHECK(slice.getSliceType() == I_SLICE, "Merge modes not available for I-slices");

  tempCS->initStructData(encTestMode.qp);

  MergeCtx mergeCtx;
  const SPS &sps = *tempCS->sps;

  if (sps.getSbTMVPEnabledFlag())
  {
    Size bufSize = g_miScaling.scale(tempCS->area.lumaSize());
    mergeCtx.subPuMvpMiBuf = MotionBuf(m_SubPuMiBuf, bufSize);
  }

  CodingUnit &cu = tempCS->addCU(tempCS->area, pm.chType);
  pm.setCUData(cu);
  cu.predMode = MODE_INTER;
  cu.slice = tempCS->slice;
  cu.tileIdx = tempCS->pps->getTileIdx(tempCS->area.lumaPos());
  cu.qp = encTestMode.qp;
  cu.affine = false;
  cu.mtsFlag = false;
#if INTER_LIC
  cu.LICFlag = false;
#endif
  cu.BcwIdx = BCW_DEFAULT;
  cu.geoFlag = true;
  cu.imv = 0;
  cu.mmvdSkip = false;
  cu.skip = false;
  cu.mipFlag = false;
#if JVET_V0130_INTRA_TMP
  cu.tmpFlag = false;
#endif
  cu.bdpcmMode = 0;

  PredictionUnit &pu = tempCS->addPU(cu, pm.chType);
  pu.mergeFlag = true;
  pu.regularMergeFlag = false;
#if TM_MRG
  pu.tmMergeFlag = false;
#endif
  PU::getGeoMergeCandidates(pu, mergeCtx);

  GeoComboCostList comboList;
  int bitsCandTB = floorLog2(GEO_NUM_PARTITION_MODE);
  PelUnitBuf geoBuffer[GEO_MAX_NUM_UNI_CANDS];
  PelUnitBuf geoTempBuf[GEO_MAX_NUM_UNI_CANDS];
  PelUnitBuf geoCombinations[GEO_MAX_TRY_WEIGHTED_SAD];
  DistParam distParam;

  const UnitArea localUnitArea(tempCS->area.chromaFormat, Area(0, 0, tempCS->area.Y().width, tempCS->area.Y().height));
  const double sqrtLambdaForFirstPass = m_pcRdCost->getMotionLambda();
  uint8_t maxNumMergeCandidates = cu.cs->sps->getMaxNumGeoCand();
  DistParam distParamWholeBlk;
  m_pcRdCost->setDistParam(distParamWholeBlk, tempCS->getOrgBuf().Y(), m_acMergeBuffer[0].Y().buf, m_acMergeBuffer[0].Y().stride, sps.getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y);
  Distortion bestWholeBlkSad = MAX_UINT64;
  double bestWholeBlkCost = MAX_DOUBLE;

  Distortion sadWholeBlk[GEO_MAX_NUM_UNI_CANDS];
  int pocMrg[GEO_MAX_NUM_UNI_CANDS];
  Mv MrgMv[GEO_MAX_NUM_UNI_CANDS];
  bool isSkipThisCand[GEO_MAX_NUM_UNI_CANDS] = { false };

  for (uint8_t mergeCand = 0; mergeCand < maxNumMergeCandidates; mergeCand++)
  {
    geoBuffer[mergeCand] = m_acMergeBuffer[mergeCand].getBuf(localUnitArea);
    mergeCtx.setMergeInfo(pu, mergeCand);
    int MrgList = mergeCtx.mvFieldNeighbours[(mergeCand << 1) + 0].refIdx == -1 ? 1 : 0;
    RefPicList MrgeRefPicList = (MrgList ? REF_PIC_LIST_1 : REF_PIC_LIST_0);
    int MrgrefIdx = mergeCtx.mvFieldNeighbours[(mergeCand << 1) + MrgList].refIdx;
    pocMrg[mergeCand] = tempCS->slice->getRefPic(MrgeRefPicList, MrgrefIdx)->getPOC();
    MrgMv[mergeCand] = mergeCtx.mvFieldNeighbours[(mergeCand << 1) + MrgList].mv;

    for( int i = 0; i < mergeCand; i++ )
    {
      if( pocMrg[mergeCand] == pocMrg[i] && MrgMv[mergeCand] == MrgMv[i] )
      {
        isSkipThisCand[mergeCand] = true;
        break;
      }
    }

    if (m_pcEncCfg->getMCTSEncConstraint() && (!(MCTSHelper::checkMvBufferForMCTSConstraint(pu))))
    {
      tempCS->initStructData(encTestMode.qp);
      return;
    }
    m_pcInterSearch->motionCompensation(pu, geoBuffer[mergeCand]);
#if MULTI_HYP_PRED
    geoTempBuf[mergeCand] = m_acRealMergeBuffer[MRG_MAX_NUM_CANDS + mergeCand].getBuf(localUnitArea);
#else
    geoTempBuf[mergeCand] = m_acMergeTmpBuffer[mergeCand].getBuf(localUnitArea);
#endif
    geoTempBuf[mergeCand].Y().copyFrom(geoBuffer[mergeCand].Y());
    geoTempBuf[mergeCand].Y().roundToOutputBitdepth(geoTempBuf[mergeCand].Y(), cu.slice->clpRng(COMPONENT_Y));
    distParamWholeBlk.cur.buf = geoTempBuf[mergeCand].Y().buf;
    distParamWholeBlk.cur.stride = geoTempBuf[mergeCand].Y().stride;
    sadWholeBlk[mergeCand] = distParamWholeBlk.distFunc(distParamWholeBlk);
    if (sadWholeBlk[mergeCand] < bestWholeBlkSad)
    {
      bestWholeBlkSad = sadWholeBlk[mergeCand];
      bestWholeBlkCost = ( double ) bestWholeBlkSad + ( mergeCand + 1 ) * sqrtLambdaForFirstPass;
    }
  }
#if MULTI_HYP_PRED
  m_pcInterSearch->setGeoTmpBuffer(mergeCtx);
#endif
  bool isGeo = true;
  for (uint8_t mergeCand = 1; mergeCand < maxNumMergeCandidates; mergeCand++)
  {
    isGeo &= isSkipThisCand[mergeCand];
  }
  if (isGeo)
  {
    return;
  }

  int wIdx = floorLog2(cu.lwidth()) - GEO_MIN_CU_LOG2;
  int hIdx = floorLog2(cu.lheight()) - GEO_MIN_CU_LOG2;
  for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; splitDir++)
  {
    int maskStride = 0, maskStride2 = 0;
    int stepX = 1;
    Pel* SADmask;
    int16_t angle = g_GeoParams[splitDir][0];
    if (g_angle2mirror[angle] == 2)
    {
      maskStride = -GEO_WEIGHT_MASK_SIZE;
      maskStride2 = -(int)cu.lwidth();
      SADmask = &g_globalGeoEncSADmask[g_angle2mask[g_GeoParams[splitDir][0]]][(GEO_WEIGHT_MASK_SIZE - 1 - g_weightOffset[splitDir][hIdx][wIdx][1]) * GEO_WEIGHT_MASK_SIZE + g_weightOffset[splitDir][hIdx][wIdx][0]];
    }
    else if (g_angle2mirror[angle] == 1)
    {
      stepX = -1;
      maskStride2 = cu.lwidth();
      maskStride = GEO_WEIGHT_MASK_SIZE;
      SADmask = &g_globalGeoEncSADmask[g_angle2mask[g_GeoParams[splitDir][0]]][g_weightOffset[splitDir][hIdx][wIdx][1] * GEO_WEIGHT_MASK_SIZE + (GEO_WEIGHT_MASK_SIZE - 1 - g_weightOffset[splitDir][hIdx][wIdx][0])];
    }
    else
    {
      maskStride = GEO_WEIGHT_MASK_SIZE;
      maskStride2 = -(int)cu.lwidth();
      SADmask = &g_globalGeoEncSADmask[g_angle2mask[g_GeoParams[splitDir][0]]][g_weightOffset[splitDir][hIdx][wIdx][1] * GEO_WEIGHT_MASK_SIZE + g_weightOffset[splitDir][hIdx][wIdx][0]];
    }
    Distortion sadSmall = 0, sadLarge = 0;
    for (uint8_t mergeCand = 0; mergeCand < maxNumMergeCandidates; mergeCand++)
    {
      int bitsCand = mergeCand + 1;

      m_pcRdCost->setDistParam(distParam, tempCS->getOrgBuf().Y(), geoTempBuf[mergeCand].Y().buf, geoTempBuf[mergeCand].Y().stride, SADmask, maskStride, stepX, maskStride2, sps.getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y);
      sadLarge = distParam.distFunc(distParam);
      m_GeoCostList.insert(splitDir, 0, mergeCand, (double)sadLarge + (double)bitsCand * sqrtLambdaForFirstPass);
      sadSmall = sadWholeBlk[mergeCand] - sadLarge;
      m_GeoCostList.insert(splitDir, 1, mergeCand, (double)sadSmall + (double)bitsCand * sqrtLambdaForFirstPass);
    }
  }

  for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; splitDir++)
  {
    for (int GeoMotionIdx = 0; GeoMotionIdx < maxNumMergeCandidates * (maxNumMergeCandidates - 1); GeoMotionIdx++)
    {
      unsigned int mergeCand0 = m_GeoModeTest[GeoMotionIdx].m_candIdx0;
      unsigned int mergeCand1 = m_GeoModeTest[GeoMotionIdx].m_candIdx1;
      double tempCost = m_GeoCostList.singleDistList[0][splitDir][mergeCand0].cost + m_GeoCostList.singleDistList[1][splitDir][mergeCand1].cost;
      if( tempCost > bestWholeBlkCost )
      {
        continue;
      }
      tempCost = tempCost + (double)bitsCandTB * sqrtLambdaForFirstPass;
      comboList.list.push_back(GeoMergeCombo(splitDir, mergeCand0, mergeCand1, tempCost));
    }
  }
  if( comboList.list.empty() )
  {
    return;
  }

  comboList.sortByCost();

  bool geocandHasNoResidual[GEO_MAX_TRY_WEIGHTED_SAD] = { false };
  bool bestIsSkip = false;
  int geoNumCobo = (int)comboList.list.size();
  static_vector<uint8_t, GEO_MAX_TRY_WEIGHTED_SAD> geoRdModeList;
  static_vector<double, GEO_MAX_TRY_WEIGHTED_SAD> geocandCostList;

  DistParam distParamSAD2;
  const bool useHadamard = !tempCS->slice->getDisableSATDForRD();
  m_pcRdCost->setDistParam(distParamSAD2, tempCS->getOrgBuf().Y(), m_acMergeBuffer[0].Y(), sps.getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, useHadamard);
  int geoNumMrgSATDCand = min(GEO_MAX_TRY_WEIGHTED_SATD, geoNumCobo);

  for (uint8_t candidateIdx = 0; candidateIdx < min(geoNumCobo, GEO_MAX_TRY_WEIGHTED_SAD); candidateIdx++)
  {
    int splitDir = comboList.list[candidateIdx].splitDir;
    int mergeCand0 = comboList.list[candidateIdx].mergeIdx0;
    int mergeCand1 = comboList.list[candidateIdx].mergeIdx1;

    geoCombinations[candidateIdx] = m_acGeoWeightedBuffer[candidateIdx].getBuf(localUnitArea);
    m_pcInterSearch->weightedGeoBlk(pu, splitDir, CHANNEL_TYPE_LUMA, geoCombinations[candidateIdx], geoBuffer[mergeCand0], geoBuffer[mergeCand1]);
    distParamSAD2.cur = geoCombinations[candidateIdx].Y();
    Distortion sad = distParamSAD2.distFunc(distParamSAD2);
    int mvBits = 2;
    mergeCand1 -= mergeCand1 < mergeCand0 ? 0 : 1;
    mvBits += mergeCand0;
    mvBits += mergeCand1;
    double updateCost = (double)sad + (double)(bitsCandTB + mvBits) * sqrtLambdaForFirstPass;
    comboList.list[candidateIdx].cost = updateCost;
    updateCandList(candidateIdx, updateCost, geoRdModeList, geocandCostList, geoNumMrgSATDCand);
  }
  for (uint8_t i = 0; i < geoNumMrgSATDCand; i++)
  {
    if (geocandCostList[i] > MRG_FAST_RATIO * geocandCostList[0] || geocandCostList[i] > getMergeBestSATDCost() || geocandCostList[i] > getAFFBestSATDCost())
    {
      geoNumMrgSATDCand = i;
      break;
    }
  }
  for (uint8_t i = 0; i < geoNumMrgSATDCand && isChromaEnabled(pu.chromaFormat); i++)
  {
    uint8_t candidateIdx = geoRdModeList[i];
    int splitDir = comboList.list[candidateIdx].splitDir;
    int mergeCand0 = comboList.list[candidateIdx].mergeIdx0;
    int mergeCand1 = comboList.list[candidateIdx].mergeIdx1;
    geoCombinations[candidateIdx] = m_acGeoWeightedBuffer[candidateIdx].getBuf(localUnitArea);
    m_pcInterSearch->weightedGeoBlk(pu, splitDir, CHANNEL_TYPE_CHROMA, geoCombinations[candidateIdx], geoBuffer[mergeCand0], geoBuffer[mergeCand1]);
  }

  m_bestModeUpdated = tempCS->useDbCost = bestCS->useDbCost = false;
  tempCS->initStructData(encTestMode.qp);
  uint8_t iteration;
  uint8_t iterationBegin = 0;
  iteration = 2;
  for (uint8_t noResidualPass = iterationBegin; noResidualPass < iteration; ++noResidualPass)
  {
    for (uint8_t mrgHADIdx = 0; mrgHADIdx < geoNumMrgSATDCand; mrgHADIdx++)
    {
      uint8_t candidateIdx = geoRdModeList[mrgHADIdx];
      if (((noResidualPass != 0) && geocandHasNoResidual[candidateIdx])
        || ((noResidualPass == 0) && bestIsSkip))
      {
        continue;
      }
      CodingUnit &cu = tempCS->addCU(tempCS->area, pm.chType);
      pm.setCUData(cu);
      cu.predMode = MODE_INTER;
      cu.slice = tempCS->slice;
      cu.tileIdx = tempCS->pps->getTileIdx(tempCS->area.lumaPos());
      cu.qp = encTestMode.qp;
      cu.affine = false;
      cu.mtsFlag = false;
#if INTER_LIC
      cu.LICFlag = false;
#endif
      cu.BcwIdx = BCW_DEFAULT;
      cu.geoFlag = true;
      cu.imv = 0;
      cu.mmvdSkip = false;
      cu.skip = false;
      cu.mipFlag = false;
#if JVET_V0130_INTRA_TMP
	    cu.tmpFlag = false;
#endif
      cu.bdpcmMode = 0;
      PredictionUnit &pu = tempCS->addPU(cu, pm.chType);
      pu.mergeFlag = true;
      pu.regularMergeFlag = false;
      pu.geoSplitDir = comboList.list[candidateIdx].splitDir;
      pu.geoMergeIdx0 = comboList.list[candidateIdx].mergeIdx0;
      pu.geoMergeIdx1 = comboList.list[candidateIdx].mergeIdx1;
      pu.mmvdMergeFlag = false;
      pu.mmvdMergeIdx = MAX_UCHAR;

      PU::spanGeoMotionInfo(pu, mergeCtx, pu.geoSplitDir, pu.geoMergeIdx0, pu.geoMergeIdx1);
      tempCS->getPredBuf().copyFrom(geoCombinations[candidateIdx]);
#if ENABLE_OBMC
      cu.isobmcMC = true;
      cu.obmcFlag = true;
      m_pcInterSearch->subBlockOBMC(pu);
      cu.isobmcMC = false;
#endif
      xEncodeInterResidual(tempCS, bestCS, pm, encTestMode, noResidualPass, (noResidualPass == 0 ? &geocandHasNoResidual[candidateIdx] : NULL));

      if (m_pcEncCfg->getUseFastDecisionForMerge() && !bestIsSkip)
      {
        bestIsSkip = bestCS->getCU(pm.chType)->rootCbf == 0;
      }
      tempCS->initStructData(encTestMode.qp);
    }
  }
  if (m_bestModeUpdated && bestCS->cost != MAX_DOUBLE)
  {
    xCalDebCost(*bestCS, pm);
  }
}
#endif
#if MERGE_ENC_OPT
void EncCu::xCheckSATDCostRegularMerge(CodingStructure *&tempCS, CodingUnit &cu, PredictionUnit &pu, MergeCtx mergeCtx, PelUnitBuf *acMergeTempBuffer[MMVD_MRG_MAX_RD_NUM], PelUnitBuf *&singleMergeTempBuffer, PelUnitBuf  acMergeTmpBuffer[MRG_MAX_NUM_CANDS]
#if !MULTI_PASS_DMVR
  , Mv   refinedMvdL0[MAX_NUM_PARTS_IN_CTU][MRG_MAX_NUM_CANDS]
#endif
  , unsigned& uiNumMrgSATDCand, static_vector<ModeInfo, MRG_MAX_NUM_CANDS + MMVD_ADD_NUM>  &RdModeList, static_vector<double, MRG_MAX_NUM_CANDS + MMVD_ADD_NUM> &candCostList, DistParam distParam, const TempCtx &ctxStart
#if MULTI_PASS_DMVR
  , bool* applyBDMVR
#endif
)
{
#if INTER_LIC
  cu.LICFlag = false;
#endif
  cu.mmvdSkip = false;
  cu.geoFlag = false;
  cu.affine = false;
  cu.imv = 0;

  pu.ciipFlag = false;
#if CIIP_PDPC
  pu.ciipPDPC = false;
#endif
  pu.mmvdMergeFlag = false;
  pu.regularMergeFlag = true;

  const double sqrtLambdaForFirstPassIntra = m_pcRdCost->getMotionLambda() * FRAC_BITS_SCALE;
  const UnitArea localUnitArea(tempCS->area.chromaFormat, Area(0, 0, tempCS->area.Y().width, tempCS->area.Y().height));
#if MULTI_HYP_PRED
  const bool testMHP = tempCS->sps->getUseInterMultiHyp()
    && (tempCS->area.lumaSize().area() > MULTI_HYP_PRED_RESTRICT_BLOCK_SIZE 
    && std::min(tempCS->area.lwidth(), tempCS->area.lheight()) >= MULTI_HYP_PRED_RESTRICT_MIN_WH);
#endif

  int insertPos = -1;
  for (uint32_t uiMergeCand = 0; uiMergeCand < mergeCtx.numValidMergeCand; uiMergeCand++)
  {
    mergeCtx.setMergeInfo(pu, uiMergeCand);
#if MULTI_PASS_DMVR
    pu.bdmvrRefine = false; // init as false
#endif
    pu.mvRefine = true;
    distParam.cur = singleMergeTempBuffer->Y();
    acMergeTmpBuffer[uiMergeCand] = m_acMergeTmpBuffer[uiMergeCand].getBuf(localUnitArea);
#if INTER_LIC
    m_pcInterSearch->m_storeBeforeLIC = mergeCtx.interDirNeighbours[uiMergeCand] == 3 ? false : true;
    if (m_pcInterSearch->m_storeBeforeLIC)
    {
      m_pcInterSearch->m_predictionBeforeLIC = acMergeTmpBuffer[uiMergeCand];
      m_pcInterSearch->motionCompensation(pu, *singleMergeTempBuffer, REF_PIC_LIST_X, true, true);
    }
    else
#endif
#if MULTI_PASS_DMVR
    if (applyBDMVR[uiMergeCand])
    {
      if (pu.cu->cs->sps->getUseCiip())
      {
#if MULTI_HYP_PRED
        pu.addHypData.clear();
        pu.numMergedAddHyps = 0;
#endif
        pu.mvRefine = false;
        pu.ciipFlag = true;
        m_pcInterSearch->motionCompensation(pu, acMergeTmpBuffer[uiMergeCand]);
        pu.ciipFlag = false;
#if MULTI_HYP_PRED
        mergeCtx.setMergeInfo(pu, uiMergeCand);
#endif
      }
      pu.bdmvrRefine = true;
      m_pcInterSearch->setBdmvrSubPuMvBuf(m_mvBufBDMVR[uiMergeCand << 1], m_mvBufBDMVR[(uiMergeCand << 1) + 1]);

      pu.mvRefine = true;
      m_pcInterSearch->motionCompensation(pu, *singleMergeTempBuffer);

      if( pu.bdmvrRefine )
      {
        ::memcpy( m_mvBufEncBDOF[uiMergeCand], m_pcInterSearch->getBdofSubPuMvOffset(), sizeof( Mv ) * BDOF_SUBPU_MAX_NUM );
      }

      pu.mvRefine = false;
    }
    else
#endif
    {
      m_pcInterSearch->motionCompensation(pu, *singleMergeTempBuffer, REF_PIC_LIST_X, true, true, &(acMergeTmpBuffer[uiMergeCand]));
    }
#if INTER_LIC
    m_pcInterSearch->m_storeBeforeLIC = false;
#endif
    pu.mvRefine = false;
#if !MULTI_PASS_DMVR
    if (mergeCtx.interDirNeighbours[uiMergeCand] == 3 )
    {
      mergeCtx.mvFieldNeighbours[2 * uiMergeCand].mv = pu.mv[0];
      mergeCtx.mvFieldNeighbours[2 * uiMergeCand + 1].mv = pu.mv[1];
      {
        int dx, dy, i, j, num = 0;
        dy = std::min<int>(pu.lumaSize().height, DMVR_SUBCU_HEIGHT);
        dx = std::min<int>(pu.lumaSize().width, DMVR_SUBCU_WIDTH);
        if (PU::checkDMVRCondition(pu))
        {
          for (i = 0; i < (pu.lumaSize().height); i += dy)
          {
            for (j = 0; j < (pu.lumaSize().width); j += dx)
            {
              refinedMvdL0[num][uiMergeCand] = pu.mvdL0SubPu[num];
              num++;
            }
          }
        }
      }
    }
#endif

    Distortion uiSad = distParam.distFunc(distParam);
    m_CABACEstimator->getCtx() = ctxStart;
    uint64_t fracBits = m_pcInterSearch->xCalcPuMeBits(pu);
    double cost = (double)uiSad + (double)fracBits * sqrtLambdaForFirstPassIntra;
#if MULTI_HYP_PRED
    if (testMHP && pu.addHypData.size() < tempCS->sps->getMaxNumAddHyps())
    {
      CHECK(pu.bdmvrRefine && !applyBDMVR[uiMergeCand], "");
      uint32_t uiBitsCand = uiMergeCand + 1 + 1 + 1; // one bit for merge flag,  one bit for subblock_merge_flag, and one bit for regualr_merge_flag
      MEResult mergeResult;
      mergeResult.cu = cu;
      mergeResult.pu = pu;
      mergeResult.bits = uiBitsCand;
      mergeResult.cost = uiSad + m_pcRdCost->getCost(uiBitsCand);
      m_baseResultsForMH.push_back(mergeResult);
    }
#endif
    insertPos = -1;
    updateCandList(ModeInfo(cu, pu), cost, RdModeList, candCostList, uiNumMrgSATDCand, &insertPos);
    if (insertPos != -1 && insertPos < MMVD_MRG_MAX_RD_NUM)
    {
      if (insertPos == RdModeList.size() - 1)
      {
        swap(singleMergeTempBuffer, acMergeTempBuffer[insertPos]);
      }
      else
      {
        for (uint32_t i = uint32_t(RdModeList.size()) - 1; i > insertPos; i--)
        {
          swap(acMergeTempBuffer[i - 1], acMergeTempBuffer[i]);
        }
        swap(singleMergeTempBuffer, acMergeTempBuffer[insertPos]);
      }
    }
    CHECK(std::min(uiMergeCand + 1, uiNumMrgSATDCand) != RdModeList.size(), "");
  }
#if MULTI_PASS_DMVR
  pu.bdmvrRefine = false;
#endif
}

void EncCu::xCheckSATDCostCiipMerge(CodingStructure *&tempCS, CodingUnit &cu, PredictionUnit &pu, MergeCtx mergeCtx, PelUnitBuf *acMergeTempBuffer[MMVD_MRG_MAX_RD_NUM], PelUnitBuf *&singleMergeTempBuffer, PelUnitBuf  acMergeTmpBuffer[MRG_MAX_NUM_CANDS]
  , unsigned& uiNumMrgSATDCand, static_vector<ModeInfo, MRG_MAX_NUM_CANDS + MMVD_ADD_NUM>  &RdModeList, static_vector<double, MRG_MAX_NUM_CANDS + MMVD_ADD_NUM> &candCostList, DistParam distParam, const TempCtx &ctxStart)
{
#if INTER_LIC
  cu.LICFlag = false;
#endif
  cu.mmvdSkip = false;
  cu.geoFlag = false;
  cu.affine = false;
  //cu.imv = 0;

  pu.ciipFlag = true;
#if CIIP_PDPC
  pu.ciipPDPC = false;
#endif
  pu.mmvdMergeFlag = false;
  pu.regularMergeFlag = false;

  const double sqrtLambdaForFirstPassIntra = m_pcRdCost->getMotionLambda() * FRAC_BITS_SCALE;
  const UnitArea localUnitArea(tempCS->area.chromaFormat, Area(0, 0, tempCS->area.Y().width, tempCS->area.Y().height));

  int insertPos = -1;
  for (uint32_t mergeCand = 0; mergeCand < mergeCtx.numValidMergeCand; mergeCand++)
  {
    //acMergeTmpBuffer[mergeCand] = m_acMergeTmpBuffer[mergeCand].getBuf(localUnitArea);

    // estimate merge bits
    mergeCtx.setMergeInfo(pu, mergeCand);

    // first round
    pu.intraDir[0] = PLANAR_IDX;
#if CIIP_PDPC
    for (int intraCnt = 0; intraCnt < 2; intraCnt++)
    {
      pu.ciipPDPC = intraCnt == 1;
#else
    uint32_t intraCnt = 0;
#endif
    PelBuf ciipBuff = m_ciipBuffer[intraCnt].getBuf(localUnitArea.Y());

    // generate intrainter Y prediction
    if (mergeCand == 0)
    {
      m_pcIntraSearch->initIntraPatternChType(*pu.cu, pu.Y());
      m_pcIntraSearch->predIntraAng(COMPONENT_Y, ciipBuff, pu);
    }

    if( pu.cs->picHeader->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() )
    {
      m_pcIntraSearch->geneWeightedPred<true>( COMPONENT_Y, singleMergeTempBuffer->Y(), pu, acMergeTmpBuffer[mergeCand].Y(), ciipBuff, m_pcReshape->getFwdLUT().data() );
    }
    else
    {
      m_pcIntraSearch->geneWeightedPred<false>( COMPONENT_Y, singleMergeTempBuffer->Y(), pu, acMergeTmpBuffer[mergeCand].Y(), ciipBuff );
    }

    // calculate cost
    if (pu.cs->picHeader->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag())
    {
      PelBuf tmp = m_acGeoWeightedBuffer->getBuf(localUnitArea.Y());
      tmp.rspSignal(singleMergeTempBuffer->Y(), m_pcReshape->getInvLUT());
      distParam.cur = tmp;
    }
    else
    {
      distParam.cur = singleMergeTempBuffer->Y();
    }

    Distortion sadValue = distParam.distFunc(distParam);
    m_CABACEstimator->getCtx() = ctxStart;
    pu.regularMergeFlag = false;
    uint64_t fracBits = m_pcInterSearch->xCalcPuMeBits(pu);
    double cost = (double)sadValue + (double)fracBits * sqrtLambdaForFirstPassIntra; // need to check the cost calculation again???
    insertPos = -1;
    updateCandList(ModeInfo(cu, pu), cost, RdModeList, candCostList, uiNumMrgSATDCand, &insertPos);
    if (insertPos != -1 && insertPos < MMVD_MRG_MAX_RD_NUM)
    {
      for (int i = int(RdModeList.size()) - 1; i > insertPos; i--)
      {
        swap(acMergeTempBuffer[i - 1], acMergeTempBuffer[i]);
      }
      swap(singleMergeTempBuffer, acMergeTempBuffer[insertPos]);
    }
#if CIIP_PDPC
    }
#endif
  }
pu.ciipFlag = false;
#if CIIP_PDPC
pu.ciipPDPC = false;
#endif
}

void EncCu::xCheckSATDCostMmvdMerge(CodingStructure *&tempCS, CodingUnit &cu, PredictionUnit &pu, MergeCtx mergeCtx, PelUnitBuf *acMergeTempBuffer[MMVD_MRG_MAX_RD_NUM], PelUnitBuf *&singleMergeTempBuffer
  , unsigned& uiNumMrgSATDCand, static_vector<ModeInfo, MRG_MAX_NUM_CANDS + MMVD_ADD_NUM>  &RdModeList, static_vector<double, MRG_MAX_NUM_CANDS + MMVD_ADD_NUM> &candCostList, DistParam distParam, const TempCtx &ctxStart)
{
#if INTER_LIC
  cu.LICFlag = false;
#endif
  cu.mmvdSkip = true;
  cu.geoFlag = false;
  cu.affine = false;
  cu.imv = 0;

  pu.ciipFlag = false;
#if CIIP_PDPC
  pu.ciipPDPC = false;
#endif
  pu.mmvdMergeFlag = true;
  pu.regularMergeFlag = true;

  const double sqrtLambdaForFirstPassIntra = m_pcRdCost->getMotionLambda() * FRAC_BITS_SCALE;
  int insertPos = -1;
#if MULTI_HYP_PRED
  const bool testMHP = tempCS->sps->getUseInterMultiHyp()
    && (tempCS->area.lumaSize().area() > MULTI_HYP_PRED_RESTRICT_BLOCK_SIZE
      && std::min(tempCS->area.lwidth(), tempCS->area.lheight()) >= MULTI_HYP_PRED_RESTRICT_MIN_WH);
#endif

  const int tempNum = (mergeCtx.numValidMergeCand > 1) ? MMVD_ADD_NUM : MMVD_ADD_NUM >> 1;
  for (int mmvdMergeCand = 0; mmvdMergeCand < tempNum; mmvdMergeCand++)
  {
    int baseIdx = mmvdMergeCand / MMVD_MAX_REFINE_NUM;
    int refineStep = (mmvdMergeCand - (baseIdx * MMVD_MAX_REFINE_NUM)) / 4;
    if( refineStep >= m_pcEncCfg->getMmvdDisNum() )
    {
      continue;
    }

    mergeCtx.setMmvdMergeCandiInfo(pu, mmvdMergeCand);
    pu.mvRefine = true;
    distParam.cur = singleMergeTempBuffer->Y();
    pu.mmvdEncOptMode = (refineStep > 2 ? 2 : 1);
    CHECK(!pu.mmvdMergeFlag, "MMVD merge should be set");
    // Don't do chroma MC here
    m_pcInterSearch->motionCompensation(pu, *singleMergeTempBuffer, REF_PIC_LIST_X, true, false);
    pu.mmvdEncOptMode = 0;
    pu.mvRefine = false;
    Distortion uiSad = distParam.distFunc(distParam);

    m_CABACEstimator->getCtx() = ctxStart;
    uint64_t fracBits = m_pcInterSearch->xCalcPuMeBits(pu);
    double cost = (double)uiSad + (double)fracBits * sqrtLambdaForFirstPassIntra;
#if MULTI_HYP_PRED
    if (testMHP && pu.addHypData.size() < tempCS->sps->getMaxNumAddHyps())
    {
      uint32_t uiBitsCand = baseIdx + refineStep + 2 + 1 + 1 + 1; // one bit for merge flag,  one bit for subblock_merge_flag, and one bit for regualr_merge_flag
      MEResult mergeResult;
      mergeResult.cu = cu;
      mergeResult.pu = pu;
      mergeResult.bits = uiBitsCand;
      mergeResult.cost = uiSad + m_pcRdCost->getCost(uiBitsCand);

      m_baseResultsForMH.push_back(mergeResult);
    }
#endif
    insertPos = -1;
    updateCandList(ModeInfo(cu, pu), cost, RdModeList, candCostList, uiNumMrgSATDCand, &insertPos);
    if (insertPos != -1 && insertPos < MMVD_MRG_MAX_RD_NUM)
    {
      for (int i = int(RdModeList.size()) - 1; i > insertPos; i--)
      {
        swap(acMergeTempBuffer[i - 1], acMergeTempBuffer[i]);
      }
      swap(singleMergeTempBuffer, acMergeTempBuffer[insertPos]);
    }
  }
}

void EncCu::xCheckSATDCostAffineMerge(CodingStructure *&tempCS, CodingUnit &cu, PredictionUnit &pu, AffineMergeCtx affineMergeCtx, MergeCtx& mrgCtx, PelUnitBuf *acMergeTempBuffer[MMVD_MRG_MAX_RD_NUM], PelUnitBuf *&singleMergeTempBuffer
  , unsigned& uiNumMrgSATDCand, static_vector<ModeInfo, MRG_MAX_NUM_CANDS + MMVD_ADD_NUM>  &RdModeList, static_vector<double, MRG_MAX_NUM_CANDS + MMVD_ADD_NUM> &candCostList, DistParam distParam, const TempCtx &ctxStart)
{
  cu.mmvdSkip = false;
  cu.geoFlag = false;
  cu.affine = true;
#if INTER_LIC
  cu.LICFlag = false;
#endif

  pu.mergeFlag = true;
  pu.ciipFlag = false;
#if CIIP_PDPC
  pu.ciipPDPC = false;
#endif
  pu.mmvdMergeFlag = false;
  pu.regularMergeFlag = false;
#if MULTI_HYP_PRED
  pu.addHypData.clear();
  pu.numMergedAddHyps = 0;
  const bool testMHP = tempCS->sps->getUseInterMultiHyp()
    && (tempCS->area.lumaSize().area() > MULTI_HYP_PRED_RESTRICT_BLOCK_SIZE
      && std::min(tempCS->area.lwidth(), tempCS->area.lheight()) >= MULTI_HYP_PRED_RESTRICT_MIN_WH);
#endif

  const double sqrtLambdaForFirstPassIntra = m_pcRdCost->getMotionLambda() * FRAC_BITS_SCALE;
  int insertPos = -1;
  for (uint32_t uiAffMergeCand = 0; uiAffMergeCand < affineMergeCtx.numValidMergeCand; uiAffMergeCand++)
  {
    // set merge information
    pu.interDir = affineMergeCtx.interDirNeighbours[uiAffMergeCand];
    pu.mergeIdx = uiAffMergeCand;
    cu.affineType = affineMergeCtx.affineType[uiAffMergeCand];
    cu.BcwIdx = affineMergeCtx.BcwIdx[uiAffMergeCand];
#if INTER_LIC
    cu.LICFlag = affineMergeCtx.LICFlags[uiAffMergeCand];
#endif
    pu.mv[0].setZero();
    pu.mv[1].setZero();
    cu.imv = 0;

    pu.mergeType = affineMergeCtx.mergeType[uiAffMergeCand];
    if (pu.mergeType == MRG_TYPE_SUBPU_ATMVP)
    {
      pu.refIdx[0] = affineMergeCtx.mvFieldNeighbours[(uiAffMergeCand << 1) + 0][0].refIdx;
      pu.refIdx[1] = affineMergeCtx.mvFieldNeighbours[(uiAffMergeCand << 1) + 1][0].refIdx;
      // the SbTmvp use xSubPuMC which will need to access the motion buffer for subblock MV
      PU::spanMotionInfo(pu, mrgCtx);
    }
    else
    {
      PU::setAllAffineMvField(pu, affineMergeCtx.mvFieldNeighbours[(uiAffMergeCand << 1) + 0], REF_PIC_LIST_0);
      PU::setAllAffineMvField(pu, affineMergeCtx.mvFieldNeighbours[(uiAffMergeCand << 1) + 1], REF_PIC_LIST_1);
    }

    distParam.cur = singleMergeTempBuffer->Y();

    m_pcInterSearch->motionCompensation(pu, *singleMergeTempBuffer, REF_PIC_LIST_X, true, false);
    Distortion uiSad = distParam.distFunc(distParam);

    m_CABACEstimator->getCtx() = ctxStart;
    uint64_t fracBits = m_pcInterSearch->xCalcPuMeBits(pu);
    double cost = (double)uiSad + (double)fracBits * sqrtLambdaForFirstPassIntra;
#if MULTI_HYP_PRED
    if (testMHP && pu.addHypData.size() < tempCS->sps->getMaxNumAddHyps())
    {
      uint32_t   uiBitsCand = uiAffMergeCand + 1;
      if (uiAffMergeCand == tempCS->picHeader->getMaxNumAffineMergeCand() - 1)
      {
        uiBitsCand--;
      }
      uiBitsCand = uiBitsCand + 1 + 1; // one bit for merge flag, and one bit for subblock_merge_flag
      MEResult mergeResult;
      mergeResult.cu = cu;
      mergeResult.pu = pu;
      mergeResult.bits = uiBitsCand;
      mergeResult.cost = uiSad + m_pcRdCost->getCost(uiBitsCand);
      m_baseResultsForMH.push_back(mergeResult);
    }
#endif
    insertPos = -1;
    updateCandList(ModeInfo(cu, pu), cost, RdModeList, candCostList, uiNumMrgSATDCand, &insertPos);
#if MERGE_ENC_OPT
    if (insertPos != -1 && insertPos < MMVD_MRG_MAX_RD_NUM)
#else
    if (insertPos != -1)
#endif
    {
      for (int i = int(RdModeList.size()) - 1; i > insertPos; i--)
      {
        swap(acMergeTempBuffer[i - 1], acMergeTempBuffer[i]);
      }
      swap(singleMergeTempBuffer, acMergeTempBuffer[insertPos]);
    }
  }
  pu.regularMergeFlag = true;
  cu.affine = false;
}

#if TM_MRG && MERGE_ENC_OPT
void EncCu::xCheckSATDCostTMMerge(       CodingStructure*& tempCS,
                                         CodingUnit&       cu,
                                         PredictionUnit&   pu,
                                         MergeCtx&         mrgCtx,
                                         PelUnitBuf*       acMergeTempBuffer[MMVD_MRG_MAX_RD_NUM],
                                         PelUnitBuf*&      singleMergeTempBuffer,
                                         unsigned&         uiNumMrgSATDCand,
                                         static_vector<ModeInfo, MRG_MAX_NUM_CANDS + MMVD_ADD_NUM>  &RdModeList,
                                         static_vector<double, MRG_MAX_NUM_CANDS + MMVD_ADD_NUM>    &candCostList,
                                         DistParam         distParam,
                                   const TempCtx&          ctxStart
#if MULTI_PASS_DMVR
                                       , bool*             applyBDMVR
#endif
)
{
#if MULTI_PASS_DMVR
  CHECK(applyBDMVR == nullptr, "Unexpected error");
#endif
  pu.mergeFlag        = true;
  cu.mmvdSkip         = false;
  cu.geoFlag          = false;
  cu.affine           = false;
  cu.imv              = IMV_OFF;
  pu.ciipFlag         = false;
#if CIIP_PDPC
  pu.ciipPDPC         = false;
#endif
  pu.mmvdMergeFlag    = false;
  pu.regularMergeFlag = false;
  pu.tmMergeFlag      = true;

  const double sqrtLambdaForFirstPassIntra = m_pcRdCost->getMotionLambda() * FRAC_BITS_SCALE;
  int insertPos = -1;
  for (uint32_t uiMergeCand = 0; uiMergeCand < mrgCtx.numValidMergeCand; uiMergeCand++)
  {
    mrgCtx.setMergeInfo( pu, uiMergeCand );
#if MULTI_PASS_DMVR
    if (applyBDMVR[uiMergeCand])
    {
      pu.bdmvrRefine = true;
      m_pcInterSearch->setBdmvrSubPuMvBuf(m_mvBufBDMVR4TM[uiMergeCand << 1], m_mvBufBDMVR4TM[(uiMergeCand << 1) + 1]);
    }
#endif

    pu.mvRefine = false;
#if INTER_LIC
    m_pcInterSearch->m_storeBeforeLIC = false;
#endif
    m_pcInterSearch->motionCompensation(pu, *singleMergeTempBuffer);
#if MULTI_PASS_DMVR
    if( pu.bdmvrRefine )
    {
      ::memcpy( m_mvBufEncBDOF4TM[uiMergeCand], m_pcInterSearch->getBdofSubPuMvOffset(), sizeof( Mv ) * BDOF_SUBPU_MAX_NUM );
    }
#endif
    distParam.cur = singleMergeTempBuffer->Y();
    Distortion uiSad = distParam.distFunc(distParam);

    m_CABACEstimator->getCtx() = ctxStart;
    uint64_t fracBits = m_pcInterSearch->xCalcPuMeBits(pu);
    double cost = (double)uiSad + (double)fracBits * sqrtLambdaForFirstPassIntra;
    insertPos = -1;
    updateCandList(ModeInfo(cu, pu), cost, RdModeList, candCostList, uiNumMrgSATDCand, &insertPos);

    if (insertPos != -1 && insertPos < MMVD_MRG_MAX_RD_NUM)
    {
      for (int i = int(RdModeList.size()) - 1; i > insertPos; i--)
      {
        swap(acMergeTempBuffer[i - 1], acMergeTempBuffer[i]);
      }
      swap(singleMergeTempBuffer, acMergeTempBuffer[insertPos]);
    }
  }
  pu.regularMergeFlag = true;
  cu.affine           = false;
#if AFFINE_MMVD
  pu.afMmvdFlag       = false;
#endif
  pu.tmMergeFlag      = false;
#if MULTI_PASS_DMVR
  pu.bdmvrRefine      = false;
#endif
}
#endif

#if AFFINE_MMVD && MERGE_ENC_OPT
void EncCu::xCheckSATDCostAffineMmvdMerge(       CodingStructure*& tempCS,
                                                 CodingUnit&       cu,
                                                 PredictionUnit&   pu,
                                                 AffineMergeCtx    affineMergeCtx,
                                                 MergeCtx&         mrgCtx,
                                                 PelUnitBuf*       acMergeTempBuffer[MMVD_MRG_MAX_RD_NUM],
                                                 PelUnitBuf*&      singleMergeTempBuffer,
                                                 unsigned&         uiNumMrgSATDCand,
                                                 static_vector<ModeInfo, MRG_MAX_NUM_CANDS + MMVD_ADD_NUM>  &RdModeList,
                                                 static_vector<double, MRG_MAX_NUM_CANDS + MMVD_ADD_NUM>    &candCostList,
                                                 DistParam         distParam,
                                           const TempCtx&          ctxStart
)
{
  cu.mmvdSkip         = false;
  cu.geoFlag          = false;
  cu.affine           = true;
  cu.imv              = IMV_OFF;
#if INTER_LIC
  cu.LICFlag = false;
#endif

  pu.mergeFlag        = true;
  pu.ciipFlag         = false;
#if CIIP_PDPC
  pu.ciipPDPC         = false;
#endif
  pu.mmvdMergeFlag    = false;
  pu.regularMergeFlag = false;
#if MULTI_HYP_PRED
  pu.addHypData.clear();
  pu.numMergedAddHyps = 0;
#endif

  int baseIdxToMergeIdxOffset = (int)PU::getMergeIdxFromAfMmvdBaseIdx(affineMergeCtx, 0);
  int baseCount               = std::min<int>((int)AF_MMVD_BASE_NUM, affineMergeCtx.numValidMergeCand - baseIdxToMergeIdxOffset);
  int afMmvdCandCount         = baseCount * AF_MMVD_MAX_REFINE_NUM;
  if (baseCount < 1)
  {
    return;
  }

#if MULTI_HYP_PRED
  const bool testMHP = tempCS->sps->getUseInterMultiHyp()
    && (tempCS->area.lumaSize().area() > MULTI_HYP_PRED_RESTRICT_BLOCK_SIZE
      && std::min(tempCS->area.lwidth(), tempCS->area.lheight()) >= MULTI_HYP_PRED_RESTRICT_MIN_WH);
#endif
  const double sqrtLambdaForFirstPassIntra = m_pcRdCost->getMotionLambda() * FRAC_BITS_SCALE;
  int insertPos = -1;
  for (uint32_t uiMergeCand = 0; uiMergeCand < afMmvdCandCount; uiMergeCand++)
  {
    int baseIdx = (int)uiMergeCand / AF_MMVD_MAX_REFINE_NUM;
    int stepIdx = (int)uiMergeCand - baseIdx * AF_MMVD_MAX_REFINE_NUM;
    int dirIdx  = stepIdx % AF_MMVD_OFFSET_DIR;
        stepIdx = stepIdx / AF_MMVD_OFFSET_DIR;

    // Pass Affine MMVD parameters from candidate to PU
    {
      pu.afMmvdFlag     = true;
      pu.afMmvdBaseIdx  = (uint8_t)baseIdx;
      pu.afMmvdDir      = (uint8_t)dirIdx;
      pu.afMmvdStep     = (uint8_t)stepIdx;
      pu.mergeIdx       = (uint8_t)(baseIdx + baseIdxToMergeIdxOffset);
      pu.mergeType      = affineMergeCtx.mergeType         [pu.mergeIdx];
      pu.interDir       = affineMergeCtx.interDirNeighbours[pu.mergeIdx];
      pu.cu->affineType = affineMergeCtx.affineType        [pu.mergeIdx];
#if INTER_LIC
      pu.cu->LICFlag    = affineMergeCtx.LICFlags          [pu.mergeIdx];
#endif
      pu.cu->BcwIdx     = affineMergeCtx.BcwIdx            [pu.mergeIdx];
      CHECK(pu.afMmvdDir >= AF_MMVD_OFFSET_DIR || pu.afMmvdStep >= AF_MMVD_STEP_NUM, "Affine MMVD dir or Affine MMVD step is out of range ");
      CHECK(pu.mergeType != MRG_TYPE_DEFAULT_N, "Affine MMVD must have non-SbTMVP base!");
    }

    MvField mvfMmvd[2][3];
    PU::getAfMmvdMvf(pu, affineMergeCtx, mvfMmvd, pu.mergeIdx, pu.afMmvdStep, pu.afMmvdDir);
    PU::setAllAffineMvField(pu, mvfMmvd[0], REF_PIC_LIST_0);
    PU::setAllAffineMvField(pu, mvfMmvd[1], REF_PIC_LIST_1);
    distParam.cur = singleMergeTempBuffer->Y();
    pu.mmvdEncOptMode = (stepIdx > 2 ? 3 : 0);
    m_pcInterSearch->motionCompensation(pu, *singleMergeTempBuffer, REF_PIC_LIST_X, true, false);
    pu.mmvdEncOptMode = 0;
    Distortion uiSad = distParam.distFunc(distParam);

    m_CABACEstimator->getCtx() = ctxStart;
    uint64_t fracBits = m_pcInterSearch->xCalcPuMeBits(pu);
    double cost = (double)uiSad + (double)fracBits * sqrtLambdaForFirstPassIntra;
#if MULTI_HYP_PRED
    if (testMHP && pu.addHypData.size() < tempCS->sps->getMaxNumAddHyps())
    {
      uint32_t uiBitsCand = baseIdx + stepIdx + 2 + 1 + 1 + 1; // one bit for merge flag,  one bit for subblock_merge_flag, and one bit for afMmvdFlag
      MEResult mergeResult;
      mergeResult.cu = cu;
      mergeResult.pu = pu;
      mergeResult.bits = uiBitsCand;
      mergeResult.cost = uiSad + m_pcRdCost->getCost(uiBitsCand);

      m_baseResultsForMH.push_back(mergeResult);
    }
#endif
    insertPos = -1;
    updateCandList(ModeInfo(cu, pu), cost, RdModeList, candCostList, uiNumMrgSATDCand, &insertPos);

    if (insertPos != -1 && insertPos < MMVD_MRG_MAX_RD_NUM)
    {
      for (int i = int(RdModeList.size()) - 1; i > insertPos; i--)
      {
        swap(acMergeTempBuffer[i - 1], acMergeTempBuffer[i]);
      }
      swap(singleMergeTempBuffer, acMergeTempBuffer[insertPos]);
    }
  }
  pu.regularMergeFlag = true;
  cu.affine = false;
  pu.afMmvdFlag = false;
}
#endif
#if !JVET_W0097_GPM_MMVD_TM
void EncCu::xCheckSATDCostGeoMerge(CodingStructure *&tempCS, CodingUnit &cu, PredictionUnit &pu, MergeCtx geoMergeCtx, PelUnitBuf *acMergeTempBuffer[MMVD_MRG_MAX_RD_NUM], PelUnitBuf *&singleMergeTempBuffer
  , unsigned& uiNumMrgSATDCand, static_vector<ModeInfo, MRG_MAX_NUM_CANDS + MMVD_ADD_NUM>  &RdModeList, static_vector<double, MRG_MAX_NUM_CANDS + MMVD_ADD_NUM> &candCostList, DistParam distParam, const TempCtx &ctxStart)
{
  const SPS &sps = *tempCS->sps;
  int numGeoChecked = 0;
  GeoComboCostList comboList;
  int bitsCandTB = floorLog2(GEO_NUM_PARTITION_MODE);
  PelUnitBuf geoBuffer[GEO_MAX_NUM_UNI_CANDS];
  PelUnitBuf geoTempBuf[GEO_MAX_NUM_UNI_CANDS];

  uint8_t maxNumMergeCandidates = cu.cs->sps->getMaxNumGeoCand();
  DistParam distParamWholeBlk;
  m_pcRdCost->setDistParam(distParamWholeBlk, tempCS->getOrgBuf().Y(), m_acMergeBuffer[0].Y().buf, m_acMergeBuffer[0].Y().stride, sps.getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y);
  Distortion bestWholeBlkSad = MAX_UINT64;
  double bestWholeBlkCost = MAX_DOUBLE;

  Distortion sadWholeBlk[GEO_MAX_NUM_UNI_CANDS];
  int pocMrg[GEO_MAX_NUM_UNI_CANDS];
  Mv MrgMv[GEO_MAX_NUM_UNI_CANDS];
  bool isSkipThisCand[GEO_MAX_NUM_UNI_CANDS] = { false };

  cu.affine = false;
  cu.mtsFlag = false;
#if INTER_LIC
  cu.LICFlag = false;
#endif
  cu.BcwIdx = BCW_DEFAULT;
  cu.geoFlag = true;
  cu.imv = 0;
  cu.mmvdSkip = false;
  cu.skip = false;
  cu.mipFlag = false;
#if JVET_V0130_INTRA_TMP
  cu.tmpFlag = false;
#endif
  cu.bdpcmMode = 0;
  pu.mergeFlag = true;
  pu.regularMergeFlag = false;
  pu.mmvdMergeFlag = false;
  pu.mmvdMergeIdx = MAX_UCHAR;

  const UnitArea localUnitArea(tempCS->area.chromaFormat, Area(0, 0, tempCS->area.Y().width, tempCS->area.Y().height));
  const double sqrtLambdaForFirstPass = m_pcRdCost->getMotionLambda();
  const double sqrtLambdaForFirstPassIntra = m_pcRdCost->getMotionLambda() * FRAC_BITS_SCALE;
  int insertPos = -1;

  for (uint8_t mergeCand = 0; mergeCand < maxNumMergeCandidates; mergeCand++)
  {
    geoBuffer[mergeCand] = m_acRealMergeBuffer[mergeCand].getBuf(localUnitArea);
    geoMergeCtx.setMergeInfo(pu, mergeCand);
    int MrgList = geoMergeCtx.mvFieldNeighbours[(mergeCand << 1) + 0].refIdx == -1 ? 1 : 0;
    RefPicList MrgeRefPicList = (MrgList ? REF_PIC_LIST_1 : REF_PIC_LIST_0);
    int MrgrefIdx = geoMergeCtx.mvFieldNeighbours[(mergeCand << 1) + MrgList].refIdx;
    pocMrg[mergeCand] = tempCS->slice->getRefPic(MrgeRefPicList, MrgrefIdx)->getPOC();
    MrgMv[mergeCand] = geoMergeCtx.mvFieldNeighbours[(mergeCand << 1) + MrgList].mv;

    for( int i = 0; i < mergeCand; i++ )
    {
      if( pocMrg[mergeCand] == pocMrg[i] && MrgMv[mergeCand] == MrgMv[i] )
      {
        isSkipThisCand[mergeCand] = true;
        break;
      }
    }

    m_pcInterSearch->motionCompensation(pu, geoBuffer[mergeCand]);
    geoTempBuf[mergeCand] = m_acRealMergeBuffer[MRG_MAX_NUM_CANDS + mergeCand].getBuf(localUnitArea);
    geoTempBuf[mergeCand].Y().copyFrom(geoBuffer[mergeCand].Y());
    geoTempBuf[mergeCand].Y().roundToOutputBitdepth(geoTempBuf[mergeCand].Y(), cu.slice->clpRng(COMPONENT_Y));
    distParamWholeBlk.cur.buf = geoTempBuf[mergeCand].Y().buf;
    distParamWholeBlk.cur.stride = geoTempBuf[mergeCand].Y().stride;
    sadWholeBlk[mergeCand] = distParamWholeBlk.distFunc(distParamWholeBlk);
    if (sadWholeBlk[mergeCand] < bestWholeBlkSad)
    {
      bestWholeBlkSad = sadWholeBlk[mergeCand];
      bestWholeBlkCost = (double)bestWholeBlkSad + (double)( mergeCand + 1 ) * sqrtLambdaForFirstPass;
    }
  }

  bool skipGeo = true;
  for (uint8_t mergeCand = 1; mergeCand < maxNumMergeCandidates; mergeCand++)
  {
    if( !isSkipThisCand[mergeCand] )
    {
      skipGeo = false;
      break;
    }
  }

  if( !skipGeo )
  {
    DistParam distParamGeo;
    int wIdx = floorLog2(cu.lwidth()) - GEO_MIN_CU_LOG2;
    int hIdx = floorLog2(cu.lheight()) - GEO_MIN_CU_LOG2;
    for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; splitDir++)
    {
      int maskStride = 0, maskStride2 = 0;
      int stepX = 1;
      Pel* SADmask;
      int16_t angle = g_GeoParams[splitDir][0];
      if (g_angle2mirror[angle] == 2)
      {
        maskStride = -GEO_WEIGHT_MASK_SIZE;
        maskStride2 = -(int)cu.lwidth();
        SADmask = &g_globalGeoEncSADmask[g_angle2mask[g_GeoParams[splitDir][0]]][(GEO_WEIGHT_MASK_SIZE - 1 - g_weightOffset[splitDir][hIdx][wIdx][1]) * GEO_WEIGHT_MASK_SIZE + g_weightOffset[splitDir][hIdx][wIdx][0]];
      }
      else if (g_angle2mirror[angle] == 1)
      {
        stepX = -1;
        maskStride2 = cu.lwidth();
        maskStride = GEO_WEIGHT_MASK_SIZE;
        SADmask = &g_globalGeoEncSADmask[g_angle2mask[g_GeoParams[splitDir][0]]][g_weightOffset[splitDir][hIdx][wIdx][1] * GEO_WEIGHT_MASK_SIZE + (GEO_WEIGHT_MASK_SIZE - 1 - g_weightOffset[splitDir][hIdx][wIdx][0])];
      }
      else
      {
        maskStride = GEO_WEIGHT_MASK_SIZE;
        maskStride2 = -(int)cu.lwidth();
        SADmask = &g_globalGeoEncSADmask[g_angle2mask[g_GeoParams[splitDir][0]]][g_weightOffset[splitDir][hIdx][wIdx][1] * GEO_WEIGHT_MASK_SIZE + g_weightOffset[splitDir][hIdx][wIdx][0]];
      }
      Distortion sadSmall = 0, sadLarge = 0;
      for (uint8_t mergeCand = 0; mergeCand < maxNumMergeCandidates; mergeCand++)
      {
        int bitsCand = mergeCand + 1;

        m_pcRdCost->setDistParam(distParamGeo, tempCS->getOrgBuf().Y(), geoTempBuf[mergeCand].Y().buf, geoTempBuf[mergeCand].Y().stride, SADmask, maskStride, stepX, maskStride2, sps.getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y);
        sadLarge = distParamGeo.distFunc(distParamGeo);
        m_GeoCostList.insert(splitDir, 0, mergeCand, (double)sadLarge + (double)bitsCand * sqrtLambdaForFirstPass);
        sadSmall = sadWholeBlk[mergeCand] - sadLarge;
        m_GeoCostList.insert(splitDir, 1, mergeCand, (double)sadSmall + (double)bitsCand * sqrtLambdaForFirstPass);
      }
    }

    for (int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; splitDir++)
    {
      for (int GeoMotionIdx = 0; GeoMotionIdx < maxNumMergeCandidates * (maxNumMergeCandidates - 1); GeoMotionIdx++)
      {
        unsigned int mergeCand0 = m_GeoModeTest[GeoMotionIdx].m_candIdx0;
        unsigned int mergeCand1 = m_GeoModeTest[GeoMotionIdx].m_candIdx1;
        double tempCost = m_GeoCostList.singleDistList[0][splitDir][mergeCand0].cost + m_GeoCostList.singleDistList[1][splitDir][mergeCand1].cost;
        if( tempCost > bestWholeBlkCost )
        {
          continue;
        }
        tempCost = tempCost + (double)bitsCandTB * sqrtLambdaForFirstPass;
        comboList.list.push_back(GeoMergeCombo(splitDir, mergeCand0, mergeCand1, tempCost));
      }
    }
    if (!comboList.list.empty())
    {
      comboList.sortByCost();
      int geoNumCobo = (int)comboList.list.size();
      const int numGeoSATD = min(geoNumCobo, GEO_MAX_TRY_WEIGHTED_SAD);
      double bestPrevCost = candCostList.size() > 0 ? candCostList[0] : MAX_DOUBLE;
      for (uint8_t candidateIdx = 0; candidateIdx < numGeoSATD; candidateIdx++)
      {
        int splitDir = comboList.list[candidateIdx].splitDir;
        int mergeCand0 = comboList.list[candidateIdx].mergeIdx0;
        int mergeCand1 = comboList.list[candidateIdx].mergeIdx1;

        pu.geoSplitDir = comboList.list[candidateIdx].splitDir;
        pu.geoMergeIdx0 = comboList.list[candidateIdx].mergeIdx0;
        pu.geoMergeIdx1 = comboList.list[candidateIdx].mergeIdx1;

        m_pcInterSearch->weightedGeoBlk(pu, splitDir, CHANNEL_TYPE_LUMA, *singleMergeTempBuffer, geoBuffer[mergeCand0], geoBuffer[mergeCand1]);

        distParam.cur = singleMergeTempBuffer->Y();
        Distortion uiSad = distParam.distFunc(distParam);

        m_CABACEstimator->getCtx() = ctxStart;
        uint64_t fracBits = m_pcInterSearch->xCalcPuMeBits(pu);
        double cost = (double)uiSad + (double)fracBits * sqrtLambdaForFirstPassIntra;

        if (cost > bestPrevCost) // skip GEO candidate if cost larger than the best from regular and affine as in original design
        {
          continue;
        }
        numGeoChecked++;
        insertPos = -1;
        updateCandList(ModeInfo(cu, pu), cost, RdModeList, candCostList, uiNumMrgSATDCand, &insertPos);
#if MERGE_ENC_OPT
        if (insertPos != -1 && insertPos < MMVD_MRG_MAX_RD_NUM)
#else
        if (insertPos != -1)
#endif
        {
          m_pcInterSearch->weightedGeoBlk(pu, splitDir, CHANNEL_TYPE_CHROMA, *singleMergeTempBuffer, geoBuffer[pu.geoMergeIdx0], geoBuffer[pu.geoMergeIdx1]); //have to use pu.geoMergeIdx1 since  mergeCand1 maybe changed
          for (int i = int(RdModeList.size()) - 1; i > insertPos; i--)
          {
            swap(acMergeTempBuffer[i - 1], acMergeTempBuffer[i]);
          }
          swap(singleMergeTempBuffer, acMergeTempBuffer[insertPos]);
        }
      }
    }
  }

  cu.geoFlag = false;
  if( numGeoChecked < GEO_MAX_TRY_WEIGHTED_SATD ) // try to match the original number of full RD
  {
    uiNumMrgSATDCand = uiNumMrgSATDCand - GEO_MAX_TRY_WEIGHTED_SATD + numGeoChecked;
  }

  if (uiNumMrgSATDCand > RdModeList.size()) //to make sure we have engough candidates in the list
  {
    uiNumMrgSATDCand = (unsigned int)RdModeList.size();
  }
}
#endif
#else
void EncCu::xCheckRDCostAffineMerge2Nx2N( CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &partitioner, const EncTestMode& encTestMode )
{
  if( m_modeCtrl->getFastDeltaQp() )
  {
    return;
  }

  if ( bestCS->area.lumaSize().width < 8 || bestCS->area.lumaSize().height < 8 )
  {
    return;
  }
  m_bestModeUpdated = tempCS->useDbCost = bestCS->useDbCost = false;
  const Slice &slice = *tempCS->slice;

  CHECK( slice.getSliceType() == I_SLICE, "Affine Merge modes not available for I-slices" );

  tempCS->initStructData( encTestMode.qp );

  AffineMergeCtx affineMergeCtx;
  const SPS &sps = *tempCS->sps;
  if (sps.getMaxNumAffineMergeCand() == 0)
  {
    return;
  }

  setAFFBestSATDCost(MAX_DOUBLE);

  MergeCtx mrgCtx;
  if (sps.getSbTMVPEnabledFlag())
  {
    Size bufSize = g_miScaling.scale( tempCS->area.lumaSize() );
    mrgCtx.subPuMvpMiBuf = MotionBuf( m_SubPuMiBuf, bufSize );
    affineMergeCtx.mrgCtx = &mrgCtx;
  }

  {
    // first get merge candidates
    CodingUnit cu( tempCS->area );
    cu.cs = tempCS;
    cu.predMode = MODE_INTER;
    cu.slice = tempCS->slice;
    cu.tileIdx          = tempCS->pps->getTileIdx( tempCS->area.lumaPos() );
    cu.mmvdSkip = false;
#if INTER_LIC
    cu.LICFlag = false;
#endif

    PredictionUnit pu( tempCS->area );
    pu.cu = &cu;
    pu.cs = tempCS;
    pu.regularMergeFlag = false;
    PU::getAffineMergeCand( pu, affineMergeCtx );
#if JVET_W0090_ARMC_TM
    if (sps.getUseAML())
    {
      m_pcInterSearch->adjustAffineMergeCandidates(pu, affineMergeCtx);
    }
#endif

    if ( affineMergeCtx.numValidMergeCand <= 0 )
    {
      return;
    }
  }

  bool candHasNoResidual[AFFINE_MRG_MAX_NUM_CANDS];
  for ( uint32_t ui = 0; ui < affineMergeCtx.numValidMergeCand; ui++ )
  {
    candHasNoResidual[ui] = false;
  }

  bool                                        bestIsSkip = false;
  uint32_t                                    uiNumMrgSATDCand = affineMergeCtx.numValidMergeCand;
  PelUnitBuf                                  acMergeBuffer[AFFINE_MRG_MAX_NUM_CANDS];
  static_vector<uint32_t, AFFINE_MRG_MAX_NUM_CANDS>  RdModeList;
  bool                                        mrgTempBufSet = false;

  for ( uint32_t i = 0; i < AFFINE_MRG_MAX_NUM_CANDS; i++ )
  {
    RdModeList.push_back( i );
  }

  if ( m_pcEncCfg->getUseFastMerge() )
  {
    uiNumMrgSATDCand = std::min( NUM_AFF_MRG_SATD_CAND, affineMergeCtx.numValidMergeCand );
    bestIsSkip = false;

    if ( auto blkCache = dynamic_cast<CacheBlkInfoCtrl*>(m_modeCtrl) )
    {
      bestIsSkip = blkCache->isSkip( tempCS->area );
    }

    static_vector<double, AFFINE_MRG_MAX_NUM_CANDS> candCostList;

    // 1. Pass: get SATD-cost for selected candidates and reduce their count
    if ( !bestIsSkip )
    {
      RdModeList.clear();
      mrgTempBufSet = true;
#if JVET_W0097_GPM_MMVD_TM
      const double sqrtLambdaForFirstPassIntra = m_pcRdCost->getMotionLambda() * FRAC_BITS_SCALE;
      const TempCtx ctxStart(m_CtxCache, m_CABACEstimator->getCtx());
#else
      const double sqrtLambdaForFirstPass = m_pcRdCost->getMotionLambda( );
#endif

      CodingUnit &cu = tempCS->addCU( tempCS->area, partitioner.chType );

      partitioner.setCUData( cu );
      cu.slice = tempCS->slice;
      cu.tileIdx          = tempCS->pps->getTileIdx( tempCS->area.lumaPos() );
      cu.skip = false;
      cu.affine = true;
#if INTER_LIC
      cu.LICFlag = false;
#endif
      cu.predMode = MODE_INTER;
      cu.chromaQpAdj = m_cuChromaQpOffsetIdxPlus1;
      cu.qp = encTestMode.qp;

      PredictionUnit &pu = tempCS->addPU( cu, partitioner.chType );

      DistParam distParam;
      const bool bUseHadamard = !tempCS->slice->getDisableSATDForRD();
      m_pcRdCost->setDistParam( distParam, tempCS->getOrgBuf().Y(), m_acMergeBuffer[0].Y(), sps.getBitDepth( CHANNEL_TYPE_LUMA ), COMPONENT_Y, bUseHadamard );

      const UnitArea localUnitArea( tempCS->area.chromaFormat, Area( 0, 0, tempCS->area.Y().width, tempCS->area.Y().height ) );
#if MULTI_HYP_PRED
      const bool testMHP = tempCS->sps->getUseInterMultiHyp()
        && (tempCS->area.lumaSize().area() > MULTI_HYP_PRED_RESTRICT_BLOCK_SIZE 
        && std::min(tempCS->area.lwidth(), tempCS->area.lheight()) >= MULTI_HYP_PRED_RESTRICT_MIN_WH);
#endif

      for ( uint32_t uiMergeCand = 0; uiMergeCand < affineMergeCtx.numValidMergeCand; uiMergeCand++ )
      {
        acMergeBuffer[uiMergeCand] = m_acMergeBuffer[uiMergeCand].getBuf( localUnitArea );

        // set merge information
        pu.interDir = affineMergeCtx.interDirNeighbours[uiMergeCand];
        pu.mergeFlag = true;
        pu.regularMergeFlag = false;
        pu.mergeIdx = uiMergeCand;
        cu.affineType = affineMergeCtx.affineType[uiMergeCand];
#if AFFINE_MMVD
        pu.afMmvdFlag = false;
#endif
        cu.BcwIdx = affineMergeCtx.BcwIdx[uiMergeCand];
#if INTER_LIC
        cu.LICFlag = affineMergeCtx.LICFlags[uiMergeCand];
#endif

        pu.mergeType = affineMergeCtx.mergeType[uiMergeCand];
        if ( pu.mergeType == MRG_TYPE_SUBPU_ATMVP )
        {
          pu.refIdx[0] = affineMergeCtx.mvFieldNeighbours[(uiMergeCand << 1) + 0][0].refIdx;
          pu.refIdx[1] = affineMergeCtx.mvFieldNeighbours[(uiMergeCand << 1) + 1][0].refIdx;
          PU::spanMotionInfo( pu, mrgCtx );
        }
        else
        {
          PU::setAllAffineMvField( pu, affineMergeCtx.mvFieldNeighbours[(uiMergeCand << 1) + 0], REF_PIC_LIST_0 );
          PU::setAllAffineMvField( pu, affineMergeCtx.mvFieldNeighbours[(uiMergeCand << 1) + 1], REF_PIC_LIST_1 );

          PU::spanMotionInfo( pu );
        }

        distParam.cur = acMergeBuffer[uiMergeCand].Y();

        m_pcInterSearch->motionCompensation( pu, acMergeBuffer[uiMergeCand], REF_PIC_LIST_X, true, false );

        Distortion uiSad = distParam.distFunc( distParam );
#if JVET_W0097_GPM_MMVD_TM
        m_CABACEstimator->getCtx() = ctxStart;
        uint64_t fracBits = m_pcInterSearch->xCalcPuMeBits(pu);
        double cost = (double)uiSad + (double)fracBits * sqrtLambdaForFirstPassIntra;
#else
        uint32_t   uiBitsCand = uiMergeCand + 1;
        if ( uiMergeCand == tempCS->picHeader->getMaxNumAffineMergeCand() - 1 )
        {
          uiBitsCand--;
        }
        double cost = (double)uiSad + (double)uiBitsCand * sqrtLambdaForFirstPass;
#endif
#if MULTI_HYP_PRED
        if (testMHP && pu.addHypData.size() < tempCS->sps->getMaxNumAddHyps())
        {
          uiBitsCand = uiBitsCand + 1 + 1; // one bit for merge flag, and one bit for subblock_merge_flag
          MEResult mergeResult;
          mergeResult.cu = cu;
          mergeResult.pu = pu;
          mergeResult.bits = uiBitsCand;
          mergeResult.cost = uiSad + m_pcRdCost->getCost(uiBitsCand);
          m_baseResultsForMH.push_back(mergeResult);
        }
#endif
        updateCandList( uiMergeCand, cost, RdModeList, candCostList
          , uiNumMrgSATDCand );

        CHECK( std::min( uiMergeCand + 1, uiNumMrgSATDCand ) != RdModeList.size(), "" );
      }

      // Try to limit number of candidates using SATD-costs
      for ( uint32_t i = 1; i < uiNumMrgSATDCand; i++ )
      {
        if ( candCostList[i] > MRG_FAST_RATIO * candCostList[0] )
        {
          uiNumMrgSATDCand = i;
          break;
        }
      }

      tempCS->initStructData( encTestMode.qp );
      setAFFBestSATDCost(candCostList[0]);
#if JVET_W0097_GPM_MMVD_TM
      m_CABACEstimator->getCtx() = ctxStart;
#endif
    }
    else
    {
      uiNumMrgSATDCand = affineMergeCtx.numValidMergeCand;
    }
  }

  uint32_t iteration;
  uint32_t iterationBegin = 0;
  iteration = 2;
  for (uint32_t uiNoResidualPass = iterationBegin; uiNoResidualPass < iteration; ++uiNoResidualPass)
  {
    for ( uint32_t uiMrgHADIdx = 0; uiMrgHADIdx < uiNumMrgSATDCand; uiMrgHADIdx++ )
    {
      uint32_t uiMergeCand = RdModeList[uiMrgHADIdx];

      if ( ((uiNoResidualPass != 0) && candHasNoResidual[uiMergeCand])
        || ((uiNoResidualPass == 0) && bestIsSkip) )
      {
        continue;
      }

      // first get merge candidates
      CodingUnit &cu = tempCS->addCU( tempCS->area, partitioner.chType );

      partitioner.setCUData( cu );
      cu.slice = tempCS->slice;
      cu.tileIdx          = tempCS->pps->getTileIdx( tempCS->area.lumaPos() );
#if INTER_LIC
      cu.LICFlag = false;
#endif
      cu.skip = false;
      cu.affine = true;
      cu.predMode = MODE_INTER;
      cu.chromaQpAdj = m_cuChromaQpOffsetIdxPlus1;
      cu.qp = encTestMode.qp;
      PredictionUnit &pu = tempCS->addPU( cu, partitioner.chType );

      // set merge information
      pu.mergeFlag = true;
      pu.mergeIdx = uiMergeCand;
      pu.interDir = affineMergeCtx.interDirNeighbours[uiMergeCand];
      cu.affineType = affineMergeCtx.affineType[uiMergeCand];
#if AFFINE_MMVD
      pu.afMmvdFlag = false;
#endif
      cu.BcwIdx = affineMergeCtx.BcwIdx[uiMergeCand];
#if INTER_LIC
      cu.LICFlag = affineMergeCtx.LICFlags[uiMergeCand];
#endif

      pu.mergeType = affineMergeCtx.mergeType[uiMergeCand];
      if ( pu.mergeType == MRG_TYPE_SUBPU_ATMVP )
      {
        pu.refIdx[0] = affineMergeCtx.mvFieldNeighbours[(uiMergeCand << 1) + 0][0].refIdx;
        pu.refIdx[1] = affineMergeCtx.mvFieldNeighbours[(uiMergeCand << 1) + 1][0].refIdx;
        PU::spanMotionInfo( pu, mrgCtx );
      }
      else
      {
        PU::setAllAffineMvField( pu, affineMergeCtx.mvFieldNeighbours[(uiMergeCand << 1) + 0], REF_PIC_LIST_0 );
        PU::setAllAffineMvField( pu, affineMergeCtx.mvFieldNeighbours[(uiMergeCand << 1) + 1], REF_PIC_LIST_1 );

        PU::spanMotionInfo( pu );
      }

      if( m_pcEncCfg->getMCTSEncConstraint() && ( !( MCTSHelper::checkMvBufferForMCTSConstraint( *cu.firstPU ) ) ) )
      {
        // Do not use this mode
        tempCS->initStructData( encTestMode.qp );
        return;
      }
      if ( mrgTempBufSet )
      {
        tempCS->getPredBuf().copyFrom(acMergeBuffer[uiMergeCand], true, false);   // Copy Luma Only
        m_pcInterSearch->motionCompensation(pu, REF_PIC_LIST_X, false, true);
      }
      else
      {
        m_pcInterSearch->motionCompensation( pu );
      }
#if ENABLE_OBMC
      cu.isobmcMC = true;
      m_pcInterSearch->subBlockOBMC(*cu.firstPU);
      cu.isobmcMC = false;
#endif
      xEncodeInterResidual( tempCS, bestCS, partitioner, encTestMode, uiNoResidualPass, ( uiNoResidualPass == 0 ? &candHasNoResidual[uiMergeCand] : NULL ) );

      if ( m_pcEncCfg->getUseFastDecisionForMerge() && !bestIsSkip )
      {
        bestIsSkip = bestCS->getCU( partitioner.chType )->rootCbf == 0;
      }
      tempCS->initStructData( encTestMode.qp );
    }// end loop uiMrgHADIdx

    if ( uiNoResidualPass == 0 && m_pcEncCfg->getUseEarlySkipDetection() )
    {
      const CodingUnit     &bestCU = *bestCS->getCU( partitioner.chType );
      const PredictionUnit &bestPU = *bestCS->getPU( partitioner.chType );

      if ( bestCU.rootCbf == 0 )
      {
        if ( bestPU.mergeFlag )
        {
          m_modeCtrl->setEarlySkipDetected();
        }
        else if ( m_pcEncCfg->getMotionEstimationSearchMethod() != MESEARCH_SELECTIVE )
        {
          int absolute_MV = 0;

          for ( uint32_t uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++ )
          {
            if ( slice.getNumRefIdx( RefPicList( uiRefListIdx ) ) > 0 )
            {
              absolute_MV += bestPU.mvd[uiRefListIdx].getAbsHor() + bestPU.mvd[uiRefListIdx].getAbsVer();
            }
          }

          if ( absolute_MV == 0 )
          {
            m_modeCtrl->setEarlySkipDetected();
          }
        }
      }
    }
  }
  if ( m_bestModeUpdated && bestCS->cost != MAX_DOUBLE )
  {
    xCalDebCost( *bestCS, partitioner );
  }
}
#endif

#if AFFINE_MMVD && !MERGE_ENC_OPT
void EncCu::xCheckRDCostAffineMmvd2Nx2N(CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &partitioner, const EncTestMode& encTestMode)
{
  if (m_modeCtrl->getFastDeltaQp())
  {
    return;
  }
  if (bestCS->area.lumaSize().width < 8 || bestCS->area.lumaSize().height < 8)
  {
    return;
  }
  m_bestModeUpdated = tempCS->useDbCost = bestCS->useDbCost = false;
  const Slice &slice = *tempCS->slice;

  CHECK(slice.getSliceType() == I_SLICE, "Affine MMVD modes not available for I-slices");

  tempCS->initStructData(encTestMode.qp);

  AffineMergeCtx affineMergeCtx;
  const SPS &sps = *tempCS->sps;
  MergeCtx mrgCtx;
  if (sps.getSbTMVPEnabledFlag())
  {
    Size bufSize = g_miScaling.scale(tempCS->area.lumaSize());
    mrgCtx.subPuMvpMiBuf = MotionBuf(m_SubPuMiBuf, bufSize);
    affineMergeCtx.mrgCtx = &mrgCtx;
  }

  {
    // first get merge candidates
    CodingUnit cu(tempCS->area);
    cu.cs       = tempCS;
    cu.predMode = MODE_INTER;
    cu.slice    = tempCS->slice;
    cu.tileIdx  = tempCS->pps->getTileIdx( tempCS->area.lumaPos() );
    cu.mmvdSkip = false;
    cu.geoFlag  = false;
#if INTER_LIC
    cu.LICFlag  = false;
#endif

    PredictionUnit pu(tempCS->area);
    pu.cu = &cu;
    pu.cs = tempCS;

    PU::getAffineMergeCand(pu, affineMergeCtx);
  }
  int baseIdxToMergeIdxOffset = (int)PU::getMergeIdxFromAfMmvdBaseIdx(affineMergeCtx, 0);
  int baseCount = std::min<int>((int)AF_MMVD_BASE_NUM, affineMergeCtx.numValidMergeCand - baseIdxToMergeIdxOffset);
  if (baseCount < 1)
  {
    return;
  }

  bool candHasNoResidual[AF_MMVD_NUM];
  for (uint32_t ui = 0; ui < AF_MMVD_NUM; ui++)
  {
    candHasNoResidual[ui] = false;
  }

  const UnitArea localUnitArea(tempCS->area.chromaFormat, Area(0, 0, tempCS->area.Y().width, tempCS->area.Y().height));

  bool                                    bestIsSkip       = false;
  int                                     afMmvdCandCount  = baseCount * AF_MMVD_MAX_REFINE_NUM;
  uint32_t                                uiNumMrgSATDCand = std::min( AF_MMVD_NUM, afMmvdCandCount );
  static_vector<uint32_t, AF_MMVD_NUM>    RdModeList;

  for (uint32_t i = 0; i < AF_MMVD_NUM; i++)
  {
    RdModeList.push_back(i);
  }

  if (m_pcEncCfg->getUseFastMerge())
  {
    uiNumMrgSATDCand = std::min((uint32_t)NUM_AF_MMVD_SATD_CAND, uiNumMrgSATDCand);
    bestIsSkip = false;
    if (auto blkCache = dynamic_cast<CacheBlkInfoCtrl*>(m_modeCtrl))
    {
      bestIsSkip = blkCache->isSkip(tempCS->area);
    }

    static_vector<double, AF_MMVD_NUM> candCostList;

    // 1. Pass: get SATD-cost for selected candidates and reduce their count
    if (!bestIsSkip)
    {
      RdModeList.clear();
      const double sqrtLambdaForFirstPass = m_pcRdCost->getMotionLambda();

      CodingUnit &cu = tempCS->addCU(tempCS->area, partitioner.chType);

      partitioner.setCUData(cu);
      cu.slice       = tempCS->slice;
      cu.tileIdx     = tempCS->pps->getTileIdx( tempCS->area.lumaPos() );
      cu.skip        = false;
      cu.affine      = true;
#if INTER_LIC
      cu.LICFlag     = false;
#endif
      cu.predMode    = MODE_INTER;
      cu.chromaQpAdj = m_cuChromaQpOffsetIdxPlus1;
      cu.qp          = encTestMode.qp;

      PredictionUnit &pu = tempCS->addPU(cu, partitioner.chType);

#if MULTI_HYP_PRED
      const bool testMHP = tempCS->sps->getUseInterMultiHyp()
        && (tempCS->area.lumaSize().area() > MULTI_HYP_PRED_RESTRICT_BLOCK_SIZE
          && std::min(tempCS->area.lwidth(), tempCS->area.lheight()) >= MULTI_HYP_PRED_RESTRICT_MIN_WH);
#endif
      DistParam distParam;
      const bool bUseHadamard = !tempCS->slice->getDisableSATDForRD();
      m_pcRdCost->setDistParam(distParam, tempCS->getOrgBuf().Y(), m_acMergeBuffer[0].Y(), sps.getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, bUseHadamard);

      {
        for (uint32_t uiMergeCand = 0; uiMergeCand < afMmvdCandCount; uiMergeCand++)
        {
          int baseIdx = (int)uiMergeCand / AF_MMVD_MAX_REFINE_NUM;
          int stepIdx = (int)uiMergeCand - baseIdx * AF_MMVD_MAX_REFINE_NUM;
          int dirIdx  = stepIdx % AF_MMVD_OFFSET_DIR;
              stepIdx = stepIdx / AF_MMVD_OFFSET_DIR;

          PelUnitBuf predTempBuf = m_acMergeBuffer[0].getBuf(localUnitArea);

          // Pass Affine MMVD parameters from candidate to PU
          {
            pu.mergeFlag      = true;
            pu.afMmvdFlag     = true;
            pu.afMmvdBaseIdx  = (uint8_t)baseIdx;
            pu.afMmvdDir      = (uint8_t)dirIdx;
            pu.afMmvdStep     = (uint8_t)stepIdx;
            pu.mergeIdx       = (uint8_t)(baseIdx + baseIdxToMergeIdxOffset);
            pu.mergeType      = affineMergeCtx.mergeType         [pu.mergeIdx];
            pu.interDir       = affineMergeCtx.interDirNeighbours[pu.mergeIdx];
            pu.cu->affineType = affineMergeCtx.affineType        [pu.mergeIdx];
#if INTER_LIC
            pu.cu->LICFlag    = affineMergeCtx.LICFlags          [pu.mergeIdx];
#endif
            pu.cu->BcwIdx     = affineMergeCtx.BcwIdx            [pu.mergeIdx];
            pu.mmvdMergeFlag  = false;
            pu.ciipFlag       = false;

            CHECK(pu.afMmvdDir >= AF_MMVD_OFFSET_DIR || pu.afMmvdStep >= AF_MMVD_STEP_NUM, "Affine MMVD dir or Affine MMVD step is out of range ");
            CHECK(pu.mergeType != MRG_TYPE_DEFAULT_N, "Affine MMVD must have non-SbTMVP base!");
          }

          MvField mvfMmvd[2][3];
          PU::getAfMmvdMvf(pu, affineMergeCtx, mvfMmvd, pu.mergeIdx, pu.afMmvdStep, pu.afMmvdDir);
          PU::setAllAffineMvField(pu, mvfMmvd[0], REF_PIC_LIST_0);
          PU::setAllAffineMvField(pu, mvfMmvd[1], REF_PIC_LIST_1);
          PU::spanMotionInfo(pu);
          pu.mmvdEncOptMode = (stepIdx > 2 ? 3 : 0);
          m_pcInterSearch->motionCompensation(pu, predTempBuf, REF_PIC_LIST_X, true, false);
          pu.mmvdEncOptMode = 0;

          distParam.cur = predTempBuf.Y();
          Distortion uiSad = distParam.distFunc(distParam);
          uint32_t   uiBitsCand = PU::getAfMmvdEstBits(pu);
          double cost = (double)uiSad + (double)uiBitsCand * sqrtLambdaForFirstPass;
#if MULTI_HYP_PRED
          if (testMHP && pu.addHypData.size() < tempCS->sps->getMaxNumAddHyps())
          {
            uint32_t uiBitsCand = baseIdx + stepIdx + 2 + 1 + 1 + 1; // one bit for merge flag,  one bit for subblock_merge_flag, and one bit for afMmvdFlag
            uiBitsCand++; // for mmvd_flag
            MEResult mergeResult;
            mergeResult.cu = cu;
            mergeResult.pu = pu;
            mergeResult.bits = uiBitsCand;
            mergeResult.cost = uiSad + m_pcRdCost->getCost(uiBitsCand);

            m_baseResultsForMH.push_back(mergeResult);
          }
#endif
          updateCandList(uiMergeCand, cost, RdModeList, candCostList, uiNumMrgSATDCand);

          CHECK(std::min(uiMergeCand + 1, uiNumMrgSATDCand) != RdModeList.size(), "");
        }
      }

      // Try to limit number of candidates using SATD-costs
      for (uint32_t i = 1; i < uiNumMrgSATDCand; i++)
      {
        if (candCostList[i] > MRG_FAST_RATIO * candCostList[0])
        {
          uiNumMrgSATDCand = i;
          break;
        }
      }

      tempCS->initStructData(encTestMode.qp);
#if !MERGE_ENC_OPT
      setAFFBestSATDCost(std::min<double>(getAFFBestSATDCost(), candCostList[0]));
#endif
    }
    else
    {
      uiNumMrgSATDCand = afMmvdCandCount;
    }
  }

  bool bBufferMcFromNoResidualPass = (uiNumMrgSATDCand <= NUM_AF_MMVD_SATD_CAND);

  // 2. Pass: check candidates using full RD test
  for (uint32_t uiNoResidualPass = 0; uiNoResidualPass < 2; ++uiNoResidualPass)
  {
    for (uint32_t uiMrgHADIdx = 0; uiMrgHADIdx < uiNumMrgSATDCand; uiMrgHADIdx++)
    {
      uint32_t uiMergeCand = RdModeList[uiMrgHADIdx];

      if (((uiNoResidualPass != 0) && candHasNoResidual[uiMergeCand])
        || ((uiNoResidualPass == 0) && bestIsSkip))
      {
        continue;
      }
      // first get merge candidates
      CodingUnit &cu = tempCS->addCU(tempCS->area, partitioner.chType);

      partitioner.setCUData(cu);
      cu.slice       = tempCS->slice;
      cu.tileIdx     = tempCS->pps->getTileIdx( tempCS->area.lumaPos() );
      cu.skip        = false;
      cu.affine      = true;
#if INTER_LIC
      cu.LICFlag     = false;
#endif
      cu.predMode    = MODE_INTER;
      cu.chromaQpAdj = m_cuChromaQpOffsetIdxPlus1;
      cu.qp          = encTestMode.qp;
      PredictionUnit &pu = tempCS->addPU(cu, partitioner.chType);

      // Pass Affine MMVD parameters from candidate to PU
      int baseIdx = (int)uiMergeCand / AF_MMVD_MAX_REFINE_NUM;
      int stepIdx = (int)uiMergeCand - baseIdx * AF_MMVD_MAX_REFINE_NUM;
      int dirIdx  = stepIdx % AF_MMVD_OFFSET_DIR;
          stepIdx = stepIdx / AF_MMVD_OFFSET_DIR;
      {
        // set merge information
        pu.mergeFlag      = true;
        pu.afMmvdFlag     = true;
        pu.afMmvdBaseIdx  = (uint8_t)baseIdx;
        pu.afMmvdDir      = (uint8_t)dirIdx;
        pu.afMmvdStep     = (uint8_t)stepIdx;
        pu.mergeIdx       = (uint8_t)(baseIdx + baseIdxToMergeIdxOffset);
        pu.mergeType      = affineMergeCtx.mergeType         [pu.mergeIdx];
#if INTER_LIC
        pu.cu->LICFlag    = affineMergeCtx.LICFlags          [pu.mergeIdx];
#endif
        pu.interDir       = affineMergeCtx.interDirNeighbours[pu.mergeIdx];
        pu.cu->affineType = affineMergeCtx.affineType        [pu.mergeIdx];
        pu.cu->BcwIdx     = affineMergeCtx.BcwIdx            [pu.mergeIdx];
        pu.mmvdMergeFlag  = false;
        pu.ciipFlag       = false;
      }

      MvField mvfMmvd[2][3];
      PU::getAfMmvdMvf(pu, affineMergeCtx, mvfMmvd, pu.mergeIdx, pu.afMmvdStep, pu.afMmvdDir);
      PU::setAllAffineMvField(pu, mvfMmvd[0], REF_PIC_LIST_0);
      PU::setAllAffineMvField(pu, mvfMmvd[1], REF_PIC_LIST_1);
      PU::spanMotionInfo(pu);
      if (m_pcEncCfg->getMCTSEncConstraint() && (!(MCTSHelper::checkMvBufferForMCTSConstraint(*cu.firstPU))))
      {
        // Do not use this mode
        tempCS->initStructData(encTestMode.qp);
        return;
      }

      pu.mmvdEncOptMode = 0;
      if (bBufferMcFromNoResidualPass)
      {
        PelUnitBuf predTempBuf = m_acMergeBuffer[uiMrgHADIdx].getBuf(localUnitArea);
        if (uiNoResidualPass == 0)
        {
          m_pcInterSearch->motionCompensation(pu, predTempBuf, REF_PIC_LIST_X);
        }
        pu.cs->getPredBuf(pu).copyFrom(predTempBuf);
      }
      else
      {
        m_pcInterSearch->motionCompensation(pu);
      }

      xEncodeInterResidual( tempCS, bestCS, partitioner, encTestMode, uiNoResidualPass, ( uiNoResidualPass == 0 ? &candHasNoResidual[uiMergeCand] : NULL ) );

      if (m_pcEncCfg->getUseFastDecisionForMerge() && !bestIsSkip)
      {
        bestIsSkip = bestCS->getCU(partitioner.chType)->rootCbf == 0;
      }
      tempCS->initStructData(encTestMode.qp);
    }// end loop uiMrgHADIdx

    if (uiNoResidualPass == 0 && m_pcEncCfg->getUseEarlySkipDetection())
    {
      const CodingUnit     &bestCU = *bestCS->getCU(partitioner.chType);
      const PredictionUnit &bestPU = *bestCS->getPU(partitioner.chType);

      if (bestCU.rootCbf == 0)
      {
        if (bestPU.mergeFlag)
        {
          m_modeCtrl->setEarlySkipDetected();
        }
        else if (m_pcEncCfg->getMotionEstimationSearchMethod() != MESEARCH_SELECTIVE)
        {
          int absolute_MV = 0;

          for (uint32_t uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++)
          {
            if (slice.getNumRefIdx(RefPicList(uiRefListIdx)) > 0)
            {
              absolute_MV += bestPU.mvd[uiRefListIdx].getAbsHor() + bestPU.mvd[uiRefListIdx].getAbsVer();
            }
          }

          if (absolute_MV == 0)
          {
            m_modeCtrl->setEarlySkipDetected();
          }
        }
      }
    }
  }
}
#endif

#if TM_MRG && !MERGE_ENC_OPT
void EncCu::xCheckRDCostTMMerge2Nx2N(CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &partitioner, const EncTestMode& encTestMode)
{
  const Slice &slice = *tempCS->slice;
  CHECK(slice.getSliceType() == I_SLICE, "Merge modes not available for I-slices");
  tempCS->initStructData(encTestMode.qp);

  MergeCtx mergeCtx;
  const SPS &sps = *tempCS->sps;
  if (sps.getSbTMVPEnabledFlag())
  {
    Size bufSize = g_miScaling.scale(tempCS->area.lumaSize());
    mergeCtx.subPuMvpMiBuf = MotionBuf(m_SubPuMiBuf, bufSize);
  }
#if MULTI_PASS_DMVR
  bool applyBDMVR4TM[TM_MRG_MAX_NUM_CANDS] = { false };
#endif

  {
    // first get merge candidates
    CodingUnit cu( tempCS->area );
    cu.cs       = tempCS;
    cu.predMode = MODE_INTER;
    cu.slice    = tempCS->slice;
    cu.tileIdx  = tempCS->pps->getTileIdx( tempCS->area.lumaPos() );
#if INTER_LIC
    cu.LICFlag  = false;
#endif

    PredictionUnit pu( tempCS->area );
    pu.cu = &cu;
    pu.cs = tempCS;
    cu.firstPU = &pu;

    pu.tmMergeFlag = true;
    PU::getInterMergeCandidates(pu, mergeCtx, 0);
#if JVET_W0090_ARMC_TM
    if (sps.getUseAML())
    {
      m_pcInterSearch->adjustInterMergeCandidates(pu, mergeCtx);
    }
#endif

    for( uint32_t uiMergeCand = 0; uiMergeCand < mergeCtx.numValidMergeCand; uiMergeCand++ )
    {
      mergeCtx.setMergeInfo( pu, uiMergeCand );

#if MULTI_PASS_DMVR
      applyBDMVR4TM[uiMergeCand] = PU::checkBDMVRCondition(pu);
      if (applyBDMVR4TM[uiMergeCand])
      {
        pu.bdmvrRefine = true;
        m_pcInterSearch->setBdmvrSubPuMvBuf(m_mvBufBDMVR4TM[uiMergeCand << 1], m_mvBufBDMVR4TM[(uiMergeCand << 1) + 1]);
        applyBDMVR4TM[uiMergeCand] = m_pcInterSearch->processBDMVR(pu);
      }
      else
      {
        m_pcInterSearch->deriveTMMv(pu);
      }
#else
      m_pcInterSearch->deriveTMMv( pu );
#endif

      // Store refined motion back to mergeCtx
      mergeCtx.interDirNeighbours[uiMergeCand] = pu.interDir;
      mergeCtx.BcwIdx[uiMergeCand] = pu.cu->BcwIdx;  // Bcw may change, because bi may be reduced to uni by deriveTMMv(pu)
      mergeCtx.mvFieldNeighbours[2 * uiMergeCand].setMvField( pu.mv[0], pu.refIdx[0] );
      mergeCtx.mvFieldNeighbours[2 * uiMergeCand + 1].setMvField( pu.mv[1], pu.refIdx[1] );
      if( pu.interDir == 1 )
      {
        mergeCtx.mvFieldNeighbours[2 * uiMergeCand + 1].setMvField( Mv(), NOT_VALID );
      }
      if( pu.interDir == 2 )
      {
        mergeCtx.mvFieldNeighbours[2 * uiMergeCand].setMvField( Mv(), NOT_VALID );
      }
    }
    pu.regularMergeFlag = true;
  }

  bool candHasNoResidual[TM_MRG_MAX_NUM_CANDS];
  for (uint32_t ui = 0; ui < TM_MRG_MAX_NUM_CANDS; ui++)
  {
    candHasNoResidual[ui] = false;
  }

  bool    bestIsSkip     = false;
  int32_t candNum        = std::min(TM_MRG_MAX_NUM_CANDS, mergeCtx.numValidMergeCand);
  int32_t numMrgSATDCand = candNum;
  bool    mrgTempBufSet = false;

  static_vector<uint32_t, TM_MRG_MAX_NUM_CANDS> RdModeList;
  for (uint32_t i = 0; i < TM_MRG_MAX_NUM_CANDS; i++)
  {
    RdModeList.push_back(i);
  }

  const UnitArea localUnitArea(tempCS->area.chromaFormat, Area(0, 0, tempCS->area.Y().width, tempCS->area.Y().height));
  PelUnitBuf  acMergeRealBuffer[TM_MRG_MAX_NUM_CANDS];  CHECK(TM_MRG_MAX_NUM_CANDS > MMVD_MRG_MAX_RD_BUF_NUM, "TM cannot buffer a larger number of merge candidates than that of regular merge");
  for (unsigned i = 0; i < TM_MRG_MAX_NUM_CANDS; i++)
  {
    acMergeRealBuffer[i] = m_acMergeBuffer[i].getBuf(localUnitArea);
  }

  if (m_pcEncCfg->getUseFastMerge())
  {
    numMrgSATDCand = std::min(candNum, TM_MAX_NUM_SATD_CAND);
    bestIsSkip     = false;
    if (auto blkCache = dynamic_cast<CacheBlkInfoCtrl*>(m_modeCtrl))
    {
      bestIsSkip = blkCache->isSkip(tempCS->area);
      if (slice.getSPS()->getIBCFlag())
      {
        ComprCUCtx cuECtx = m_modeCtrl->getComprCUCtx();
        bestIsSkip = bestIsSkip && cuECtx.bestCU;
      }
    }

    static_vector<double, TM_MRG_MAX_NUM_CANDS> candCostList;

    // 1. Pass: get SATD-cost for selected candidates and reduce their count
    if( !bestIsSkip )
    {
      RdModeList.clear();
      mrgTempBufSet = true;
      const TempCtx ctxStart(m_CtxCache, m_CABACEstimator->getCtx());

      CodingUnit &cu      = tempCS->addCU( tempCS->area, partitioner.chType );
      const double sqrtLambdaForFirstPassIntra = m_pcRdCost->getMotionLambda( ) * FRAC_BITS_SCALE;
      partitioner.setCUData( cu );
      cu.slice            = tempCS->slice;
      cu.tileIdx          = tempCS->pps->getTileIdx( tempCS->area.lumaPos() );
#if INTER_LIC
      cu.LICFlag          = false;
#endif
      cu.skip             = false;
      cu.mmvdSkip         = false;
      cu.geoFlag          = false;
      cu.affine           = false;
      cu.predMode         = MODE_INTER;
      cu.chromaQpAdj      = m_cuChromaQpOffsetIdxPlus1;
      cu.qp               = encTestMode.qp;
    //cu.emtFlag  is set below

      PredictionUnit &pu  = tempCS->addPU( cu, partitioner.chType );
#if AFFINE_MMVD
      pu.afMmvdFlag       = false;
#endif
      pu.tmMergeFlag      = true;

      DistParam distParam;
      const bool bUseHadamard = !tempCS->slice->getDisableSATDForRD();
      m_pcRdCost->setDistParam (distParam, tempCS->getOrgBuf().Y(), m_acMergeBuffer[0].Y(), sps.getBitDepth (CHANNEL_TYPE_LUMA), COMPONENT_Y, bUseHadamard);

      for( uint32_t uiMergeCand = 0; uiMergeCand < candNum; uiMergeCand++ )
      {
        mergeCtx.setMergeInfo( pu, uiMergeCand );
#if MULTI_PASS_DMVR
        if (applyBDMVR4TM[uiMergeCand])
        {
          pu.bdmvrRefine = true;
          m_pcInterSearch->setBdmvrSubPuMvBuf(m_mvBufBDMVR4TM[uiMergeCand << 1], m_mvBufBDMVR4TM[(uiMergeCand << 1) + 1]);
        }
        else
        {
          PU::spanMotionInfo(pu, mergeCtx);
        }
#else
        PU::spanMotionInfo(pu, mergeCtx);
#endif

        pu.mvRefine = false;
#if INTER_LIC
        m_pcInterSearch->m_storeBeforeLIC = false;
#endif

        m_pcInterSearch->motionCompensation(pu, acMergeRealBuffer[uiMergeCand], REF_PIC_LIST_X, true, true);

#if MULTI_PASS_DMVR
        if( pu.bdmvrRefine )
        {
          ::memcpy( m_mvBufEncBDOF4TM[uiMergeCand], m_pcInterSearch->getBdofSubPuMvOffset(), sizeof( Mv ) * BDOF_SUBPU_MAX_NUM );
          PU::spanMotionInfo( pu, mergeCtx, m_mvBufBDMVR4TM[uiMergeCand << 1], m_mvBufBDMVR4TM[( uiMergeCand << 1 ) + 1], m_mvBufEncBDOF4TM[uiMergeCand] );
        }
#endif
        distParam.cur = acMergeRealBuffer[uiMergeCand].Y();
        Distortion uiSad = distParam.distFunc(distParam);
        m_CABACEstimator->getCtx() = ctxStart;
        uint64_t fracBits = m_pcInterSearch->xCalcPuMeBits(pu);
        double cost = (double)uiSad + (double)fracBits * sqrtLambdaForFirstPassIntra;

        updateCandList(uiMergeCand, cost, RdModeList, candCostList, numMrgSATDCand);

        CHECK(std::min(uiMergeCand + 1, (uint32_t)numMrgSATDCand) != RdModeList.size(), "");
      }

      // Try to limit number of candidates using SATD-costs
      for( uint32_t i = 1; i < numMrgSATDCand; i++ )
      {
        if( candCostList[i] > MRG_FAST_RATIO * candCostList[0] )
        {
          numMrgSATDCand = i;
          break;
        }
      }

      tempCS->initStructData( encTestMode.qp );
      m_CABACEstimator->getCtx() = ctxStart;
    }
    else
    {
      numMrgSATDCand = candNum;
    }
  }

  // 2. Pass: check candidates using full RD test
  m_bestModeUpdated = tempCS->useDbCost = bestCS->useDbCost = false;

  for (uint32_t uiNoResidualPass = 0; uiNoResidualPass < 2; ++uiNoResidualPass)
  {
    for( uint32_t uiMrgHADIdx = 0; uiMrgHADIdx < numMrgSATDCand; uiMrgHADIdx++ )
    {
      uint32_t uiMergeCand = RdModeList[uiMrgHADIdx];

      if (((uiNoResidualPass != 0) && candHasNoResidual[uiMrgHADIdx])
       || ( (uiNoResidualPass == 0) && bestIsSkip ) )
      {
        continue;
      }

      // first get merge candidates
      CodingUnit &cu      = tempCS->addCU( tempCS->area, partitioner.chType );

      partitioner.setCUData( cu );
      cu.slice            = tempCS->slice;
      cu.tileIdx          = tempCS->pps->getTileIdx( tempCS->area.lumaPos() );
      cu.skip             = false;
      cu.mmvdSkip         = false;
#if INTER_LIC
      cu.LICFlag          = false;
#endif
      cu.affine           = false;
      cu.geoFlag          = false;
      cu.predMode         = MODE_INTER;
      cu.chromaQpAdj      = m_cuChromaQpOffsetIdxPlus1;
      cu.qp               = encTestMode.qp;

      PredictionUnit &pu  = tempCS->addPU( cu, partitioner.chType );
#if AFFINE_MMVD
      pu.afMmvdFlag       = false;
#endif
      pu.tmMergeFlag      = true;

      {
        mergeCtx.setMergeInfo(pu, uiMergeCand);
#if MULTI_PASS_DMVR
        if (applyBDMVR4TM[uiMergeCand])
        {
          pu.bdmvrRefine = true;
          m_pcInterSearch->setBdmvrSubPuMvBuf(m_mvBufBDMVR4TM[uiMergeCand << 1], m_mvBufBDMVR4TM[(uiMergeCand << 1) + 1]);
        }
#endif
      }
#if MULTI_PASS_DMVR
      if (!pu.bdmvrRefine)
      {
        PU::spanMotionInfo(pu, mergeCtx);
      }
#else
      PU::spanMotionInfo(pu, mergeCtx);
#endif

      if( mrgTempBufSet )
      {
        tempCS->getPredBuf().copyFrom(acMergeRealBuffer[uiMergeCand]);
#if MULTI_PASS_DMVR
        if( pu.bdmvrRefine )
        {
          PU::spanMotionInfo( pu, mergeCtx, m_mvBufBDMVR4TM[uiMergeCand << 1], m_mvBufBDMVR4TM[( uiMergeCand << 1 ) + 1], m_mvBufEncBDOF4TM[uiMergeCand] );
        }
#endif
      }
      else
      {
        pu.mvRefine = false;
#if INTER_LIC
        m_pcInterSearch->m_storeBeforeLIC = false;
#endif
        m_pcInterSearch->motionCompensation( pu );
#if MULTI_PASS_DMVR
        if( pu.bdmvrRefine )
        {
          ::memcpy( m_mvBufEncBDOF4TM[uiMergeCand], m_pcInterSearch->getBdofSubPuMvOffset(), sizeof( Mv ) * BDOF_SUBPU_MAX_NUM );
          PU::spanMotionInfo( pu, mergeCtx, m_mvBufBDMVR4TM[uiMergeCand << 1], m_mvBufBDMVR4TM[( uiMergeCand << 1 ) + 1], m_mvBufEncBDOF4TM[uiMergeCand] );
        }
#endif
      }

      xEncodeInterResidual( tempCS, bestCS, partitioner, encTestMode, uiNoResidualPass, uiNoResidualPass == 0 ? &candHasNoResidual[uiMrgHADIdx] : NULL );

      if( m_pcEncCfg->getUseFastDecisionForMerge() && !bestIsSkip )
      {
        bestIsSkip = !bestCS->cus.empty() && bestCS->getCU( partitioner.chType )->rootCbf == 0;
      }
      tempCS->initStructData( encTestMode.qp );
    }// end loop uiMrgHADIdx

    if( uiNoResidualPass == 0 && m_pcEncCfg->getUseEarlySkipDetection() )
    {
      const CodingUnit     &bestCU = *bestCS->getCU( partitioner.chType );
      const PredictionUnit &bestPU = *bestCS->getPU( partitioner.chType );

      if( bestCU.rootCbf == 0 )
      {
        if( bestPU.mergeFlag )
        {
          m_modeCtrl->setEarlySkipDetected();
        }
        else if( m_pcEncCfg->getMotionEstimationSearchMethod() != MESEARCH_SELECTIVE )
        {
          int absolute_MV = 0;

          for( uint32_t uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++ )
          {
            if( slice.getNumRefIdx( RefPicList( uiRefListIdx ) ) > 0 )
            {
              absolute_MV += bestPU.mvd[uiRefListIdx].getAbsHor() + bestPU.mvd[uiRefListIdx].getAbsVer();
            }
          }

          if( absolute_MV == 0 )
          {
            m_modeCtrl->setEarlySkipDetected();
          }
        }
      }
    }
  }
  if ( m_bestModeUpdated && bestCS->cost != MAX_DOUBLE )
  {
    xCalDebCost( *bestCS, partitioner );
  }
}
#endif

//////////////////////////////////////////////////////////////////////////////////////////////
// ibc merge/skip mode check
void EncCu::xCheckRDCostIBCModeMerge2Nx2N(CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &partitioner, const EncTestMode& encTestMode)
{
  assert(partitioner.chType != CHANNEL_TYPE_CHROMA); // chroma IBC is derived
#if CTU_256
  if( tempCS->area.lwidth() >= 128 || tempCS->area.lheight() >= 128 ) // disable IBC mode larger than 64x64
#else
  if (tempCS->area.lwidth() == 128 || tempCS->area.lheight() == 128) // disable IBC mode larger than 64x64
#endif
  {
    return;
  }
  const SPS &sps = *tempCS->sps;

  tempCS->initStructData(encTestMode.qp);
  MergeCtx mergeCtx;

  if (sps.getSbTMVPEnabledFlag())
  {
    Size bufSize = g_miScaling.scale(tempCS->area.lumaSize());
    mergeCtx.subPuMvpMiBuf = MotionBuf(m_SubPuMiBuf, bufSize);
  }

  {
    // first get merge candidates
    CodingUnit cu(tempCS->area);
    cu.cs = tempCS;
    cu.predMode = MODE_IBC;
    cu.slice = tempCS->slice;
    cu.tileIdx          = tempCS->pps->getTileIdx( tempCS->area.lumaPos() );
    PredictionUnit pu(tempCS->area);
    pu.cu = &cu;
    pu.cs = tempCS;
    cu.mmvdSkip = false;
    pu.mmvdMergeFlag = false;
    pu.regularMergeFlag = false;
#if INTER_LIC
    cu.LICFlag = false;
#endif

    cu.geoFlag = false;
    PU::getIBCMergeCandidates(pu, mergeCtx);
  }

  int candHasNoResidual[MRG_MAX_NUM_CANDS];
  for (unsigned int ui = 0; ui < mergeCtx.numValidMergeCand; ui++)
  {
    candHasNoResidual[ui] = 0;
  }

  bool                                        bestIsSkip = false;
  unsigned                                    numMrgSATDCand = mergeCtx.numValidMergeCand;
  static_vector<unsigned, MRG_MAX_NUM_CANDS>  RdModeList(MRG_MAX_NUM_CANDS);
  for (unsigned i = 0; i < MRG_MAX_NUM_CANDS; i++)
  {
    RdModeList[i] = i;
  }

  //{
    static_vector<double, MRG_MAX_NUM_CANDS>  candCostList(MRG_MAX_NUM_CANDS, MAX_DOUBLE);
    // 1. Pass: get SATD-cost for selected candidates and reduce their count
    {
      const double sqrtLambdaForFirstPass = m_pcRdCost->getMotionLambda( );

      CodingUnit &cu = tempCS->addCU(CS::getArea(*tempCS, tempCS->area, (const ChannelType)partitioner.chType), (const ChannelType)partitioner.chType);

      partitioner.setCUData(cu);
      cu.slice = tempCS->slice;
      cu.tileIdx = tempCS->pps->getTileIdx( tempCS->area.lumaPos() );
      cu.skip = false;
      cu.predMode = MODE_IBC;
      cu.chromaQpAdj = m_cuChromaQpOffsetIdxPlus1;
      cu.qp = encTestMode.qp;
      cu.mmvdSkip = false;
#if INTER_LIC
      cu.LICFlag = false;
#endif

      cu.geoFlag = false;
      DistParam distParam;
      const bool bUseHadamard = !cu.slice->getDisableSATDForRD();
      PredictionUnit &pu = tempCS->addPU(cu, partitioner.chType); //tempCS->addPU(cu);
      pu.mmvdMergeFlag = false;
      pu.regularMergeFlag = false;
      Picture* refPic = pu.cu->slice->getPic();
      const CPelBuf refBuf = refPic->getRecoBuf(pu.blocks[COMPONENT_Y]);
      const Pel*        piRefSrch = refBuf.buf;
      if (tempCS->slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag())
      {
        const CompArea &area = cu.blocks[COMPONENT_Y];
        CompArea    tmpArea(COMPONENT_Y, area.chromaFormat, Position(0, 0), area.size());
        PelBuf tmpLuma = m_tmpStorageLCU->getBuf(tmpArea);
        tmpLuma.rspSignal( tempCS->getOrgBuf().Y(), m_pcReshape->getFwdLUT() );
        m_pcRdCost->setDistParam(distParam, tmpLuma, refBuf, sps.getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, bUseHadamard);
      }
      else
      m_pcRdCost->setDistParam(distParam, tempCS->getOrgBuf().Y(), refBuf, sps.getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, bUseHadamard);
      int refStride = refBuf.stride;
      const UnitArea localUnitArea(tempCS->area.chromaFormat, Area(0, 0, tempCS->area.Y().width, tempCS->area.Y().height));
      int numValidBv = mergeCtx.numValidMergeCand;
      for (unsigned int mergeCand = 0; mergeCand < mergeCtx.numValidMergeCand; mergeCand++)
      {
        mergeCtx.setMergeInfo(pu, mergeCand); // set bv info in merge mode
        const int cuPelX = pu.Y().x;
        const int cuPelY = pu.Y().y;
        int roiWidth = pu.lwidth();
        int roiHeight = pu.lheight();
        const int picWidth = pu.cs->slice->getPPS()->getPicWidthInLumaSamples();
        const int picHeight = pu.cs->slice->getPPS()->getPicHeightInLumaSamples();
        const unsigned int  lcuWidth = pu.cs->slice->getSPS()->getMaxCUWidth();
        int xPred = pu.bv.getHor();
        int yPred = pu.bv.getVer();

        if (!m_pcInterSearch->searchBv(pu, cuPelX, cuPelY, roiWidth, roiHeight, picWidth, picHeight, xPred, yPred, lcuWidth)) // not valid bv derived
        {
          numValidBv--;
          continue;
        }
        PU::spanMotionInfo(pu, mergeCtx);

        distParam.cur.buf = piRefSrch + refStride * yPred + xPred;

        Distortion sad = distParam.distFunc(distParam);
        unsigned int bitsCand = mergeCand + 1;
        if (mergeCand == tempCS->sps->getMaxNumMergeCand() - 1)
        {
          bitsCand--;
        }
        double cost = (double)sad + (double)bitsCand * sqrtLambdaForFirstPass;

        updateCandList(mergeCand, cost, RdModeList, candCostList
         , numMrgSATDCand);
      }

      // Try to limit number of candidates using SATD-costs
      if (numValidBv)
      {
        numMrgSATDCand = numValidBv;
        for (unsigned int i = 1; i < numValidBv; i++)
        {
          if (candCostList[i] > MRG_FAST_RATIO*candCostList[0])
          {
            numMrgSATDCand = i;
            break;
          }
        }
      }
      else
      {
        tempCS->dist = 0;
        tempCS->fracBits = 0;
        tempCS->cost = MAX_DOUBLE;
        tempCS->costDbOffset = 0;
        tempCS->initStructData(encTestMode.qp);
        return;
      }

      tempCS->initStructData(encTestMode.qp);
    }
  //}


  const unsigned int iteration = 2;
  m_bestModeUpdated = tempCS->useDbCost = bestCS->useDbCost = false;
  // 2. Pass: check candidates using full RD test
  for (unsigned int numResidualPass = 0; numResidualPass < iteration; numResidualPass++)
  {
    for (unsigned int mrgHADIdx = 0; mrgHADIdx < numMrgSATDCand; mrgHADIdx++)
    {
      unsigned int mergeCand = RdModeList[mrgHADIdx];
      if (!(numResidualPass == 1 && candHasNoResidual[mergeCand] == 1))
      {
        if (!(bestIsSkip && (numResidualPass == 0)))
        {
          {

            // first get merge candidates
            CodingUnit &cu = tempCS->addCU(CS::getArea(*tempCS, tempCS->area, (const ChannelType)partitioner.chType), (const ChannelType)partitioner.chType);

            partitioner.setCUData(cu);
            cu.slice = tempCS->slice;
            cu.tileIdx = tempCS->pps->getTileIdx( tempCS->area.lumaPos() );
            cu.skip = false;
            cu.predMode = MODE_IBC;
            cu.chromaQpAdj = m_cuChromaQpOffsetIdxPlus1;
            cu.qp = encTestMode.qp;
            cu.sbtInfo = 0;
#if INTER_LIC
            cu.LICFlag = false;
#endif

            PredictionUnit &pu = tempCS->addPU(cu, partitioner.chType);// tempCS->addPU(cu);
            pu.intraDir[0] = DC_IDX; // set intra pred for ibc block
            pu.intraDir[1] = PLANAR_IDX; // set intra pred for ibc block
            cu.mmvdSkip = false;
            pu.mmvdMergeFlag = false;
            pu.regularMergeFlag = false;
            cu.geoFlag = false;
            mergeCtx.setMergeInfo(pu, mergeCand);
            PU::spanMotionInfo(pu, mergeCtx);
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
            const bool chroma = !(CS::isDualITree(*tempCS));
#else
            const bool chroma = !pu.cu->isSepTree();
#endif
            //  MC
            m_pcInterSearch->motionCompensation(pu,REF_PIC_LIST_0, true, chroma);
            m_CABACEstimator->getCtx() = m_CurrCtx->start;

            m_pcInterSearch->encodeResAndCalcRdInterCU(*tempCS, partitioner, (numResidualPass != 0), true, chroma);
            if (tempCS->slice->getSPS()->getUseColorTrans())
            {
              bestCS->tmpColorSpaceCost = tempCS->tmpColorSpaceCost;
              bestCS->firstColorSpaceSelected = tempCS->firstColorSpaceSelected;
            }
            xEncodeDontSplit(*tempCS, partitioner);

#if ENABLE_QPA_SUB_CTU
            xCheckDQP (*tempCS, partitioner);
#else
            // this if-check is redundant
            if (tempCS->pps->getUseDQP() && partitioner.currQgEnable())
            {
              xCheckDQP(*tempCS, partitioner);
            }
#endif
            xCheckChromaQPOffset( *tempCS, partitioner );


            DTRACE_MODE_COST(*tempCS, m_pcRdCost->getLambda());
            xCheckBestMode(tempCS, bestCS, partitioner, encTestMode);

            tempCS->initStructData(encTestMode.qp);
          }

            if (m_pcEncCfg->getUseFastDecisionForMerge() && !bestIsSkip)
            {
              if (bestCS->getCU(partitioner.chType) == NULL)
                bestIsSkip = 0;
              else
              bestIsSkip = bestCS->getCU(partitioner.chType)->rootCbf == 0;
            }
        }
      }
    }
  }
  if ( m_bestModeUpdated && bestCS->cost != MAX_DOUBLE )
  {
    xCalDebCost( *bestCS, partitioner );
  }
}

void EncCu::xCheckRDCostIBCMode(CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &partitioner, const EncTestMode& encTestMode)
{
#if CTU_256
  if( tempCS->area.lwidth() >= 128 || tempCS->area.lheight() >= 128 ) // disable IBC mode larger than 64x64
#else
  if (tempCS->area.lwidth() == 128 || tempCS->area.lheight() == 128) // disable IBC mode larger than 64x64
#endif
  {
    return;
  }

    tempCS->initStructData(encTestMode.qp);

    m_bestModeUpdated = tempCS->useDbCost = bestCS->useDbCost = false;

    CodingUnit &cu = tempCS->addCU(CS::getArea(*tempCS, tempCS->area, partitioner.chType), partitioner.chType);

    partitioner.setCUData(cu);
    cu.slice = tempCS->slice;
    cu.tileIdx = tempCS->pps->getTileIdx( tempCS->area.lumaPos() );
    cu.skip = false;
    cu.predMode = MODE_IBC;
    cu.chromaQpAdj = m_cuChromaQpOffsetIdxPlus1;
    cu.qp = encTestMode.qp;
    cu.imv = 0;
    cu.sbtInfo = 0;

    CU::addPUs(cu);

    m_bestModeUpdated = tempCS->useDbCost = bestCS->useDbCost = false;

    PredictionUnit& pu = *cu.firstPU;
    cu.mmvdSkip = false;
    pu.mmvdMergeFlag = false;
    pu.regularMergeFlag = false;
#if INTER_LIC
    cu.LICFlag = false;
#endif

    pu.intraDir[0] = DC_IDX; // set intra pred for ibc block
    pu.intraDir[1] = PLANAR_IDX; // set intra pred for ibc block

    pu.interDir = 1; // use list 0 for IBC mode
    pu.refIdx[REF_PIC_LIST_0] = MAX_NUM_REF; // last idx in the list
      bool bValid = m_pcInterSearch->predIBCSearch(cu, partitioner, m_ctuIbcSearchRangeX, m_ctuIbcSearchRangeY, m_ibcHashMap);

      if (bValid)
      {
        PU::spanMotionInfo(pu);
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
        const bool chroma = !(CS::isDualITree(*tempCS));
#else
        const bool chroma = !pu.cu->isSepTree();
#endif
        //  MC
        m_pcInterSearch->motionCompensation(pu, REF_PIC_LIST_0, true, chroma);

        {

          m_pcInterSearch->encodeResAndCalcRdInterCU(*tempCS, partitioner, false, true, chroma);
          if (tempCS->slice->getSPS()->getUseColorTrans())
          {
            bestCS->tmpColorSpaceCost = tempCS->tmpColorSpaceCost;
            bestCS->firstColorSpaceSelected = tempCS->firstColorSpaceSelected;
          }

          xEncodeDontSplit(*tempCS, partitioner);

#if ENABLE_QPA_SUB_CTU
          xCheckDQP (*tempCS, partitioner);
#else
          // this if-check is redundant
          if (tempCS->pps->getUseDQP() && partitioner.currQgEnable())
          {
            xCheckDQP(*tempCS, partitioner);
          }
#endif
          xCheckChromaQPOffset( *tempCS, partitioner );

          tempCS->useDbCost = m_pcEncCfg->getUseEncDbOpt();
          if ( m_bestModeUpdated )
          {
            xCalDebCost( *tempCS, partitioner );
          }

          DTRACE_MODE_COST(*tempCS, m_pcRdCost->getLambda());
          xCheckBestMode(tempCS, bestCS, partitioner, encTestMode);

        }

      } // bValid
      else
      {
        tempCS->dist = 0;
        tempCS->fracBits = 0;
        tempCS->cost = MAX_DOUBLE;
        tempCS->costDbOffset = 0;
      }
}
  // check ibc mode in encoder RD
  //////////////////////////////////////////////////////////////////////////////////////////////

void EncCu::xCheckRDCostInter( CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &partitioner, const EncTestMode& encTestMode )
{
#if INTER_LIC
  if (m_pcInterSearch->m_fastLicCtrl.skipRDCheckForLIC((encTestMode.opts & ETO_LIC) > 0, IMV_OFF, bestCS->cost, tempCS->area.Y().area()))
  {
    return;
  }
#endif
#if ENABLE_OBMC
  double bestOBMCCost = MAX_DOUBLE;
#endif
  tempCS->initStructData( encTestMode.qp );
  m_pcInterSearch->setAffineModeSelected(false);

#if INTER_LIC
  bool lic = encTestMode.opts & ETO_LIC;
  if( tempCS->slice->getCheckLDC() && !lic )
#else
  if( tempCS->slice->getCheckLDC() )
#endif
  {
    m_bestBcwCost[0] = m_bestBcwCost[1] = std::numeric_limits<double>::max();
    m_bestBcwIdx[0] = m_bestBcwIdx[1] = -1;
  }

  m_pcInterSearch->resetBufferedUniMotions();
  int bcwLoopNum = (tempCS->slice->isInterB() ? BCW_NUM : 1);
  bcwLoopNum = (tempCS->sps->getUseBcw() ? bcwLoopNum : 1);
#if INTER_LIC
  bcwLoopNum = lic ? 1 : bcwLoopNum;
#endif

  if( tempCS->area.lwidth() * tempCS->area.lheight() < BCW_SIZE_CONSTRAINT )
  {
    bcwLoopNum = 1;
  }

  double curBestCost = bestCS->cost;
  double equBcwCost = MAX_DOUBLE;

  m_bestModeUpdated = tempCS->useDbCost = bestCS->useDbCost = false;

  for( int bcwLoopIdx = 0; bcwLoopIdx < bcwLoopNum; bcwLoopIdx++ )
  {
    if( m_pcEncCfg->getUseBcwFast() )
    {
      auto blkCache = dynamic_cast< CacheBlkInfoCtrl* >(m_modeCtrl);

      if( blkCache )
      {
        bool isBestInter = blkCache->getInter(bestCS->area);
        uint8_t bestBcwIdx = blkCache->getBcwIdx(bestCS->area);

        if( isBestInter && g_BcwSearchOrder[bcwLoopIdx] != BCW_DEFAULT && g_BcwSearchOrder[bcwLoopIdx] != bestBcwIdx )
        {
          continue;
        }
      }
    }
    if( !tempCS->slice->getCheckLDC() )
    {
      if( bcwLoopIdx != 0 && bcwLoopIdx != 3 && bcwLoopIdx != 4 )
      {
        continue;
      }
    }

  CodingUnit &cu      = tempCS->addCU( tempCS->area, partitioner.chType );

  partitioner.setCUData( cu );
  cu.slice            = tempCS->slice;
  cu.tileIdx          = tempCS->pps->getTileIdx( tempCS->area.lumaPos() );
  cu.skip             = false;
  cu.mmvdSkip = false;
//cu.affine
  cu.predMode         = MODE_INTER;
#if INTER_LIC
  cu.LICFlag          = lic;
#endif
  cu.chromaQpAdj      = m_cuChromaQpOffsetIdxPlus1;
  cu.qp               = encTestMode.qp;
  CU::addPUs( cu );

  cu.BcwIdx = g_BcwSearchOrder[bcwLoopIdx];
  uint8_t bcwIdx = cu.BcwIdx;
  bool  testBcw = (bcwIdx != BCW_DEFAULT);

#if INTER_LIC
  if (cu.slice->getUseLIC() && lic) { m_pcInterSearch->swapUniMvBuffer(); }
#endif
  m_pcInterSearch->predInterSearch( cu, partitioner );
#if INTER_LIC
  if (cu.slice->getUseLIC() && lic) { m_pcInterSearch->swapUniMvBuffer(); }
#endif

  bcwIdx = CU::getValidBcwIdx(cu);
  if( testBcw && bcwIdx == BCW_DEFAULT ) // Enabled Bcw but the search results is uni.
  {
    tempCS->initStructData(encTestMode.qp);
    continue;
  }
  CHECK(!(testBcw || (!testBcw && bcwIdx == BCW_DEFAULT)), " !( bTestBcw || (!bTestBcw && bcwIdx == BCW_DEFAULT ) )");

  bool isEqualUni = false;
  if( m_pcEncCfg->getUseBcwFast() )
  {
    if( cu.firstPU->interDir != 3 && testBcw == 0 )
    {
      isEqualUni = true;
    }
  }
#if ENABLE_OBMC //normal inter
  const unsigned wIdx = gp_sizeIdxInfo->idxFrom(partitioner.currArea().lwidth());
  CodingStructure *prevCS = tempCS;
  PelUnitBuf tempWoOBMCBuf = m_tempWoOBMCBuffer.subBuf(UnitAreaRelative(cu, cu));
  tempWoOBMCBuf.copyFrom(tempCS->getPredBuf(cu));
  cu.isobmcMC = true;
  cu.obmcFlag = true;
  m_pcInterSearch->subBlockOBMC(*cu.firstPU);
  cu.isobmcMC = false;
#endif
  xEncodeInterResidual( tempCS, bestCS, partitioner, encTestMode, 0
                        , 0
                        , &equBcwCost
  );
#if ENABLE_OBMC // xCheckRDCostInter
  double tempCost = (prevCS == tempCS) ? tempCS->cost : bestCS->cost;
  if (m_pTempCUWoOBMC && tempCost < bestOBMCCost)
  {
    const unsigned hIdx = gp_sizeIdxInfo->idxFrom(prevCS->area.lheight());
    m_pTempCUWoOBMC[wIdx][hIdx]->clearCUs();
    m_pTempCUWoOBMC[wIdx][hIdx]->clearPUs();
    m_pTempCUWoOBMC[wIdx][hIdx]->clearTUs();
    m_pTempCUWoOBMC[wIdx][hIdx]->copyStructure(*prevCS, partitioner.chType);

    m_pPredBufWoOBMC[wIdx][hIdx].copyFrom(tempWoOBMCBuf);
    m_pTempCUWoOBMC[wIdx][hIdx]->getPredBuf(cu).copyFrom(prevCS->getPredBuf(cu));

    bestOBMCCost = tempCost;
  }
#endif
  if( g_BcwSearchOrder[bcwLoopIdx] == BCW_DEFAULT )
    m_pcInterSearch->setAffineModeSelected((bestCS->cus.front()->affine && !(bestCS->cus.front()->firstPU->mergeFlag)));

  tempCS->initStructData(encTestMode.qp);

  double skipTH = MAX_DOUBLE;
  skipTH = (m_pcEncCfg->getUseBcwFast() ? 1.05 : MAX_DOUBLE);
  if( equBcwCost > curBestCost * skipTH )
  {
    break;
  }

  if( m_pcEncCfg->getUseBcwFast() )
  {
    if( isEqualUni == true && m_pcEncCfg->getIntraPeriod() == -1 )
    {
      break;
    }
  }
  if( g_BcwSearchOrder[bcwLoopIdx] == BCW_DEFAULT && xIsBcwSkip(cu) && m_pcEncCfg->getUseBcwFast() )
  {
    break;
  }
 }  // for( UChar bcwLoopIdx = 0; bcwLoopIdx < bcwLoopNum; bcwLoopIdx++ )
  if ( m_bestModeUpdated && bestCS->cost != MAX_DOUBLE )
  {
    xCalDebCost( *bestCS, partitioner );
  }
}




bool EncCu::xCheckRDCostInterIMV(CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &partitioner, const EncTestMode& encTestMode, double &bestIntPelCost)
{
#if ENABLE_OBMC
  double bestOBMCCost = MAX_DOUBLE;
#endif
  int iIMV = int( ( encTestMode.opts & ETO_IMV ) >> ETO_IMV_SHIFT );
  m_pcInterSearch->setAffineModeSelected(false);
  // Only Half-Pel, int-Pel, 4-Pel and fast 4-Pel allowed
  CHECK(iIMV < 1 || iIMV > 4, "Unsupported IMV Mode");
  const bool testAltHpelFilter = iIMV == 4;
  // Fast 4-Pel Mode

#if INTER_LIC
  if (m_pcInterSearch->m_fastLicCtrl.skipRDCheckForLIC((encTestMode.opts & ETO_LIC) > 0, (iIMV <= 2 ? iIMV : iIMV - 1), bestCS->cost, tempCS->area.Y().area()))
  {
    return false;
  }
#endif

  m_bestModeUpdated = tempCS->useDbCost = bestCS->useDbCost = false;

  EncTestMode encTestModeBase = encTestMode;                                        // copy for clearing non-IMV options
  encTestModeBase.opts        = EncTestModeOpts( encTestModeBase.opts & ETO_IMV );  // clear non-IMV options (is that intended?)

  tempCS->initStructData( encTestMode.qp );

#if INTER_LIC
  bool lic = encTestMode.opts & ETO_LIC;
#endif

  m_pcInterSearch->resetBufferedUniMotions();
  int bcwLoopNum = (tempCS->slice->isInterB() ? BCW_NUM : 1);
  bcwLoopNum = (tempCS->slice->getSPS()->getUseBcw() ? bcwLoopNum : 1);
#if INTER_LIC
  bcwLoopNum = lic ? 1 : bcwLoopNum;
#endif

  if( tempCS->area.lwidth() * tempCS->area.lheight() < BCW_SIZE_CONSTRAINT )
  {
    bcwLoopNum = 1;
  }

  bool validMode = false;
  double curBestCost = bestCS->cost;
  double equBcwCost = MAX_DOUBLE;

  for( int bcwLoopIdx = 0; bcwLoopIdx < bcwLoopNum; bcwLoopIdx++ )
  {
    if( m_pcEncCfg->getUseBcwFast() )
    {
      auto blkCache = dynamic_cast< CacheBlkInfoCtrl* >(m_modeCtrl);

      if( blkCache )
      {
        bool isBestInter = blkCache->getInter(bestCS->area);
        uint8_t bestBcwIdx = blkCache->getBcwIdx(bestCS->area);

        if( isBestInter && g_BcwSearchOrder[bcwLoopIdx] != BCW_DEFAULT && g_BcwSearchOrder[bcwLoopIdx] != bestBcwIdx )
        {
          continue;
        }
      }
    }

    if( !tempCS->slice->getCheckLDC() )
    {
      if( bcwLoopIdx != 0 && bcwLoopIdx != 3 && bcwLoopIdx != 4 )
      {
        continue;
      }
    }

    if( m_pcEncCfg->getUseBcwFast() && tempCS->slice->getCheckLDC() && g_BcwSearchOrder[bcwLoopIdx] != BCW_DEFAULT
      && (m_bestBcwIdx[0] >= 0 && g_BcwSearchOrder[bcwLoopIdx] != m_bestBcwIdx[0])
      && (m_bestBcwIdx[1] >= 0 && g_BcwSearchOrder[bcwLoopIdx] != m_bestBcwIdx[1]))
    {
      continue;
    }

  CodingUnit &cu = tempCS->addCU( tempCS->area, partitioner.chType );

  partitioner.setCUData( cu );
  cu.slice            = tempCS->slice;
  cu.tileIdx          = tempCS->pps->getTileIdx( tempCS->area.lumaPos() );
  cu.skip             = false;
  cu.mmvdSkip = false;
//cu.affine
  cu.predMode         = MODE_INTER;
#if INTER_LIC
  cu.LICFlag          = lic;
#endif
  cu.chromaQpAdj      = m_cuChromaQpOffsetIdxPlus1;
  cu.qp               = encTestMode.qp;

  CU::addPUs( cu );

  if (testAltHpelFilter)
  {
    cu.imv = IMV_HPEL;
  }
  else
  {
    cu.imv = iIMV == 1 ? IMV_FPEL : IMV_4PEL;
  }

  bool testBcw;
  uint8_t bcwIdx;
  bool affineAmvrEanbledFlag = !testAltHpelFilter && cu.slice->getSPS()->getAffineAmvrEnabledFlag();

  cu.BcwIdx = g_BcwSearchOrder[bcwLoopIdx];
  bcwIdx = cu.BcwIdx;
  testBcw = (bcwIdx != BCW_DEFAULT);

  cu.firstPU->interDir = 10;

#if INTER_LIC
  if (cu.slice->getUseLIC() && lic) { m_pcInterSearch->swapUniMvBuffer(); }
#endif
  m_pcInterSearch->predInterSearch( cu, partitioner );
#if INTER_LIC
  if (cu.slice->getUseLIC() && lic) { m_pcInterSearch->swapUniMvBuffer(); }
#endif

  if ( cu.firstPU->interDir <= 3 )
  {
    bcwIdx = CU::getValidBcwIdx(cu);
  }
  else
  {
    return false;
  }

  if( m_pcEncCfg->getMCTSEncConstraint() && ( ( cu.firstPU->refIdx[L0] < 0 && cu.firstPU->refIdx[L1] < 0 ) || ( !( MCTSHelper::checkMvBufferForMCTSConstraint( *cu.firstPU ) ) ) ) )
  {
    // Do not use this mode
    tempCS->initStructData( encTestMode.qp );
    continue;
  }
  if( testBcw && bcwIdx == BCW_DEFAULT ) // Enabled Bcw but the search results is uni.
  {
    tempCS->initStructData(encTestMode.qp);
    continue;
  }
  CHECK(!(testBcw || (!testBcw && bcwIdx == BCW_DEFAULT)), " !( bTestBcw || (!bTestBcw && bcwIdx == BCW_DEFAULT ) )");

  bool isEqualUni = false;
  if( m_pcEncCfg->getUseBcwFast() )
  {
    if( cu.firstPU->interDir != 3 && testBcw == 0 )
    {
      isEqualUni = true;
    }
  }

  if ( !CU::hasSubCUNonZeroMVd( cu ) && !CU::hasSubCUNonZeroAffineMVd( cu ) )
  {
    if (m_modeCtrl->useModeResult(encTestModeBase, tempCS, partitioner))
    {
      std::swap(tempCS, bestCS);
      // store temp best CI for next CU coding
      m_CurrCtx->best = m_CABACEstimator->getCtx();
    }
    if ( affineAmvrEanbledFlag )
    {
      tempCS->initStructData( encTestMode.qp );
      continue;
    }
    else
    {
      return false;
    }
  }
#if ENABLE_OBMC //normal inter IMV
  CodingStructure *prevCS = tempCS;
  PelUnitBuf tempWoOBMCBuf = m_tempWoOBMCBuffer.subBuf(UnitAreaRelative(cu, cu));
  tempWoOBMCBuf.copyFrom(tempCS->getPredBuf(cu));
  cu.isobmcMC = true;
  cu.obmcFlag = true;
  m_pcInterSearch->subBlockOBMC(*cu.firstPU);
  cu.isobmcMC = false;
#endif
  xEncodeInterResidual( tempCS, bestCS, partitioner, encTestModeBase, 0
                        , 0
                        , &equBcwCost
  );
#if ENABLE_OBMC
  double tempCost = (prevCS == tempCS) ? tempCS->cost : bestCS->cost;
  if (m_pTempCUWoOBMC && tempCost < bestOBMCCost)
  {
    const unsigned wIdx = gp_sizeIdxInfo->idxFrom(tempCS->area.lwidth());
    const unsigned hIdx = gp_sizeIdxInfo->idxFrom(tempCS->area.lheight());

    m_pTempCUWoOBMC[wIdx][hIdx]->clearCUs();
    m_pTempCUWoOBMC[wIdx][hIdx]->clearPUs();
    m_pTempCUWoOBMC[wIdx][hIdx]->clearTUs();
    m_pTempCUWoOBMC[wIdx][hIdx]->copyStructure(*prevCS, partitioner.chType);

    m_pPredBufWoOBMC[wIdx][hIdx].copyFrom(tempWoOBMCBuf);
    m_pTempCUWoOBMC[wIdx][hIdx]->getPredBuf(cu).copyFrom(prevCS->getPredBuf(cu));

    bestOBMCCost = tempCost;
  }
#endif
  if( cu.imv == IMV_FPEL && tempCS->cost < bestIntPelCost )
  {
    bestIntPelCost = tempCS->cost;
  }
  tempCS->initStructData(encTestMode.qp);

  double skipTH = MAX_DOUBLE;
  skipTH = (m_pcEncCfg->getUseBcwFast() ? 1.05 : MAX_DOUBLE);
  if( equBcwCost > curBestCost * skipTH )
  {
    break;
  }

  if( m_pcEncCfg->getUseBcwFast() )
  {
    if( isEqualUni == true && m_pcEncCfg->getIntraPeriod() == -1 )
    {
      break;
    }
  }
  if( g_BcwSearchOrder[bcwLoopIdx] == BCW_DEFAULT && xIsBcwSkip(cu) && m_pcEncCfg->getUseBcwFast() )
  {
    break;
  }
  validMode = true;
 } // for( UChar bcwLoopIdx = 0; bcwLoopIdx < bcwLoopNum; bcwLoopIdx++ )

  if ( m_bestModeUpdated && bestCS->cost != MAX_DOUBLE )
  {
    xCalDebCost( *bestCS, partitioner );
  }

  return tempCS->slice->getSPS()->getAffineAmvrEnabledFlag() ? validMode : true;
}

void EncCu::xCalDebCost( CodingStructure &cs, Partitioner &partitioner, bool calDist )
{
  if ( cs.cost == MAX_DOUBLE )
  {
    cs.costDbOffset = 0;
  }

  if ( cs.slice->getDeblockingFilterDisable() || ( !m_pcEncCfg->getUseEncDbOpt() && !calDist ) )
  {
    return;
  }

  m_pcLoopFilter->setEnc(true);
  const ChromaFormat format = cs.area.chromaFormat;
  CodingUnit*                cu = cs.getCU(partitioner.chType);
  const Position lumaPos = cu->Y().valid() ? cu->Y().pos() : recalcPosition( format, cu->chType, CHANNEL_TYPE_LUMA, cu->blocks[cu->chType].pos() );
  bool topEdgeAvai = lumaPos.y > 0 && ((lumaPos.y % 4) == 0);
  bool leftEdgeAvai = lumaPos.x > 0 && ((lumaPos.x % 4) == 0);
  bool anyEdgeAvai = topEdgeAvai || leftEdgeAvai;
  cs.costDbOffset = 0;

  if ( calDist )
  {
    const UnitArea currCsArea = clipArea( CS::getArea( cs, cs.area, partitioner.chType ), *cs.picture );
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
    ComponentID compStr = (CS::isDualITree(cs) && !isLuma(partitioner.chType)) ? COMPONENT_Cb : COMPONENT_Y;
    ComponentID compEnd = (CS::isDualITree(cs) && isLuma(partitioner.chType)) ? COMPONENT_Y : COMPONENT_Cr;
#else
    ComponentID compStr = ( cu->isSepTree() && !isLuma( partitioner.chType ) ) ? COMPONENT_Cb : COMPONENT_Y;
    ComponentID compEnd = ( cu->isSepTree() && isLuma( partitioner.chType ) ) ? COMPONENT_Y : COMPONENT_Cr;
#endif
    Distortion finalDistortion = 0;
    for ( int comp = compStr; comp <= compEnd; comp++ )
    {
      const ComponentID compID = ComponentID( comp );
      CPelBuf org = cs.getOrgBuf( compID );
      CPelBuf reco = cs.getRecoBuf( compID );
      finalDistortion += getDistortionDb( cs, org, reco, compID, currCsArea.block( compID ), false );
    }
    //updated distortion
    cs.dist = finalDistortion;
  }

  if ( anyEdgeAvai && m_pcEncCfg->getUseEncDbOpt() )
  {
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
    ComponentID compStr = (CS::isDualITree(cs) && !isLuma(partitioner.chType)) ? COMPONENT_Cb : COMPONENT_Y;
    ComponentID compEnd = (CS::isDualITree(cs) && isLuma(partitioner.chType)) ? COMPONENT_Y : COMPONENT_Cr;
#else
    ComponentID compStr = ( cu->isSepTree() && !isLuma( partitioner.chType ) ) ? COMPONENT_Cb : COMPONENT_Y;
    ComponentID compEnd = ( cu->isSepTree() &&  isLuma( partitioner.chType ) ) ? COMPONENT_Y : COMPONENT_Cr;
#endif
    const UnitArea currCsArea = clipArea( CS::getArea( cs, cs.area, partitioner.chType ), *cs.picture );

    PelStorage&          picDbBuf = m_pcLoopFilter->getDbEncPicYuvBuffer();

    //deblock neighbour pixels
    const Size     lumaSize = cu->Y().valid() ? cu->Y().size() : recalcSize( format, cu->chType, CHANNEL_TYPE_LUMA, cu->blocks[cu->chType].size() );

    const int verOffset = lumaPos.y > 7 ? 8 : 4;
    const int horOffset = lumaPos.x > 7 ? 8 : 4;
    const UnitArea areaTop(  format, Area( lumaPos.x, lumaPos.y - verOffset, lumaSize.width, verOffset  ) );
    const UnitArea areaLeft( format, Area( lumaPos.x - horOffset, lumaPos.y, horOffset, lumaSize.height ) );
    for ( int compIdx = compStr; compIdx <= compEnd; compIdx++ )
    {
      ComponentID compId = (ComponentID)compIdx;

      //Copy current CU's reco to Deblock Pic Buffer
      const CompArea&  curCompArea = currCsArea.block( compId );

      if( cs.slice->getLmcsEnabledFlag() && m_pcReshape->getSliceReshaperInfo().getUseSliceReshaper() && isLuma( compId ) )
      {
        picDbBuf.getBuf( curCompArea ).rspSignal( cs.getRecoBuf( curCompArea ), m_pcReshape->getInvLUT() );
      }
      else
      {
        picDbBuf.getBuf( curCompArea ).copyFrom( cs.getRecoBuf( curCompArea ) );
      }

      //left neighbour
      if ( leftEdgeAvai )
      {
        const CompArea&  compArea = areaLeft.block(compId);

        if( cs.slice->getLmcsEnabledFlag() && m_pcReshape->getSliceReshaperInfo().getUseSliceReshaper() && isLuma( compId ) )
        {
          picDbBuf.getBuf( compArea ).rspSignal( cs.picture->getRecoBuf( compArea ), m_pcReshape->getInvLUT() );
        }
        else
        {
          picDbBuf.getBuf( compArea ).copyFrom( cs.picture->getRecoBuf( compArea ) );
        }
      }
      //top neighbour
      if ( topEdgeAvai )
      {
        const CompArea&  compArea = areaTop.block( compId );

        if( cs.slice->getLmcsEnabledFlag() && m_pcReshape->getSliceReshaperInfo().getUseSliceReshaper() && isLuma( compId ) )
        {
          picDbBuf.getBuf( compArea ).rspSignal( cs.picture->getRecoBuf( compArea ), m_pcReshape->getInvLUT() );
        }
        else
        {
          picDbBuf.getBuf( compArea ).copyFrom( cs.picture->getRecoBuf( compArea ) );
        }
      }
    }
    
#if JVET_V0094_BILATERAL_FILTER
    // Bilateral:
    // The CU itself, the above area and the area to the left have been copied into
    //     PelStorage&          picDbBuf = m_pcLoopFilter->getDbEncPicYuvBuffer();
    //  It is now possible to insert the code for bilateral filtering here.
    
    if(cs.pps->getUseBIF())
    {
      for (auto &currTU : CU::traverseTUs(*cu))
      {
        bool isInter = (cu->predMode == MODE_INTER) ? true : false;
        if ((TU::getCbf(currTU, COMPONENT_Y) || isInter == false) && (currTU.cu->qp > 17) && (128 > std::max(currTU.lumaSize().width, currTU.lumaSize().height)))
        {
          if ((isInter == false) || (32 > std::min(currTU.lumaSize().width, currTU.lumaSize().height)))
          {
            CompArea &compArea = currTU.block(COMPONENT_Y);
            PelBuf    recBuf = picDbBuf.getBuf(compArea);
            PelBuf recIPredBuf = recBuf;
            std::vector<Pel>        my_invLUT;
            m_bilateralFilter->bilateralFilterRDOdiamond5x5(recBuf, recBuf, recBuf, currTU.cu->qp, recIPredBuf, cs.slice->clpRng(COMPONENT_Y), currTU, true, false, my_invLUT);
          }
        }
      }
    }
#endif
    
    //deblock
    if ( leftEdgeAvai )
    {
      m_pcLoopFilter->resetFilterLengths();
      m_pcLoopFilter->xDeblockCU( *cu, EDGE_VER );
    }

    if (topEdgeAvai)
    {
      m_pcLoopFilter->resetFilterLengths();
      m_pcLoopFilter->xDeblockCU( *cu, EDGE_HOR );
    }

    //update current CU SSE
    Distortion distCur = 0;
    for ( int compIdx = compStr; compIdx <= compEnd; compIdx++ )
    {
      ComponentID compId = (ComponentID)compIdx;
      CPelBuf reco = picDbBuf.getBuf( currCsArea.block( compId ) );
      CPelBuf org = cs.getOrgBuf( compId );
      distCur += getDistortionDb( cs, org, reco, compId, currCsArea.block( compId ), true );
    }

    //calculate difference between DB_before_SSE and DB_after_SSE for neighbouring CUs
    Distortion distBeforeDb = 0, distAfterDb = 0;
    for (int compIdx = compStr; compIdx <= compEnd; compIdx++)
    {
      ComponentID compId = (ComponentID)compIdx;
      if ( leftEdgeAvai )
      {
        const CompArea&  compArea = areaLeft.block( compId );
        CPelBuf org = cs.picture->getOrigBuf( compArea );
        CPelBuf reco = cs.picture->getRecoBuf( compArea );
        CPelBuf recoDb = picDbBuf.getBuf( compArea );
        distBeforeDb += getDistortionDb( cs, org, reco,   compId, compArea, false );
        distAfterDb  += getDistortionDb( cs, org, recoDb, compId, compArea, true  );
      }
      if ( topEdgeAvai )
      {
        const CompArea&  compArea = areaTop.block( compId );
        CPelBuf org = cs.picture->getOrigBuf( compArea );
        CPelBuf reco = cs.picture->getRecoBuf( compArea );
        CPelBuf recoDb = picDbBuf.getBuf( compArea );
        distBeforeDb += getDistortionDb( cs, org, reco,   compId, compArea, false );
        distAfterDb  += getDistortionDb( cs, org, recoDb, compId, compArea, true  );
      }
    }

    //updated cost
    int64_t distTmp = distCur - cs.dist + distAfterDb - distBeforeDb;
    int sign = distTmp < 0 ? -1 : 1;
    distTmp = distTmp < 0 ? -distTmp : distTmp;
    cs.costDbOffset = sign * m_pcRdCost->calcRdCost( 0, distTmp );
  }

  m_pcLoopFilter->setEnc( false );
}

Distortion EncCu::getDistortionDb( CodingStructure &cs, CPelBuf org, CPelBuf reco, ComponentID compID, const CompArea& compArea, bool afterDb )
{
  Distortion dist = 0;
#if WCG_EXT
  m_pcRdCost->setChromaFormat(cs.sps->getChromaFormatIdc());
  CPelBuf orgLuma = cs.picture->getOrigBuf( cs.area.blocks[COMPONENT_Y] );
  if (m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled() || (
    m_pcEncCfg->getLmcs() && (cs.slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag())))
  {
    if ( compID == COMPONENT_Y && !afterDb && !m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled())
    {
      CompArea    tmpArea( COMPONENT_Y, cs.area.chromaFormat, Position( 0, 0 ), compArea.size() );
      PelBuf tmpRecLuma = m_tmpStorageLCU->getBuf( tmpArea );
      tmpRecLuma.rspSignal( reco, m_pcReshape->getInvLUT() );
      dist += m_pcRdCost->getDistPart( org, tmpRecLuma, cs.sps->getBitDepth( toChannelType( compID ) ), compID, DF_SSE_WTD, &orgLuma );
    }
    else
    {
      dist += m_pcRdCost->getDistPart( org, reco, cs.sps->getBitDepth( toChannelType( compID ) ), compID, DF_SSE_WTD, &orgLuma );
    }
  }
  else if (m_pcEncCfg->getLmcs() && cs.slice->getLmcsEnabledFlag() && cs.slice->isIntra()) //intra slice
  {
    if ( compID == COMPONENT_Y && afterDb )
    {
      CompArea    tmpArea( COMPONENT_Y, cs.area.chromaFormat, Position( 0, 0 ), compArea.size() );
      PelBuf tmpRecLuma = m_tmpStorageLCU->getBuf( tmpArea );
      tmpRecLuma.rspSignal( reco, m_pcReshape->getFwdLUT() );
      dist += m_pcRdCost->getDistPart( org, tmpRecLuma, cs.sps->getBitDepth( toChannelType( compID ) ), compID, DF_SSE );
    }
    else
    {
      dist += m_pcRdCost->getDistPart( org, reco, cs.sps->getBitDepth(toChannelType( compID ) ), compID, DF_SSE );
    }
  }
  else
#endif
  {
    dist = m_pcRdCost->getDistPart( org, reco, cs.sps->getBitDepth( toChannelType( compID ) ), compID, DF_SSE );
  }
  return dist;
}

void EncCu::xEncodeInterResidual(   CodingStructure *&tempCS
                                  , CodingStructure *&bestCS
                                  , Partitioner &partitioner
                                  , const EncTestMode& encTestMode
                                  , int residualPass
                                  , bool* bestHasNonResi
                                  , double* equBcwCost
  )
{

  CodingUnit*            cu        = tempCS->getCU( partitioner.chType );
  double   bestCostInternal        = MAX_DOUBLE;
  double           bestCost        = bestCS->cost;
  double           bestCostBegin   = bestCS->cost;
  CodingUnit*      prevBestCU      = bestCS->getCU( partitioner.chType );
  uint8_t          prevBestSbt     = ( prevBestCU == nullptr ) ? 0 : prevBestCU->sbtInfo;
  bool              swapped        = false; // avoid unwanted data copy
  bool             reloadCU        = false;

  const PredictionUnit& pu = *cu->firstPU;

  // clang-format off
  const int affineShiftTab[3] =
  {
    MV_PRECISION_INTERNAL - MV_PRECISION_QUARTER,
    MV_PRECISION_INTERNAL - MV_PRECISION_SIXTEENTH,
    MV_PRECISION_INTERNAL - MV_PRECISION_INT
  };

  const int normalShiftTab[NUM_IMV_MODES] =
  {
    MV_PRECISION_INTERNAL - MV_PRECISION_QUARTER,
    MV_PRECISION_INTERNAL - MV_PRECISION_INT,
    MV_PRECISION_INTERNAL - MV_PRECISION_4PEL,
    MV_PRECISION_INTERNAL - MV_PRECISION_HALF,
  };
  // clang-format on

  int mvShift;

  for (int refList = 0; refList < NUM_REF_PIC_LIST_01; refList++)
  {
    if (pu.refIdx[refList] >= 0)
    {
      if (!cu->affine)
      {
        mvShift = normalShiftTab[cu->imv];
        Mv signaledmvd(pu.mvd[refList].getHor() >> mvShift, pu.mvd[refList].getVer() >> mvShift);
        if (!((signaledmvd.getHor() >= MVD_MIN) && (signaledmvd.getHor() <= MVD_MAX)) || !((signaledmvd.getVer() >= MVD_MIN) && (signaledmvd.getVer() <= MVD_MAX)))
          return;
      }
      else
      {
        for (int ctrlP = 1 + (cu->affineType == AFFINEMODEL_6PARAM); ctrlP >= 0; ctrlP--)
        {
          mvShift = affineShiftTab[cu->imv];
          Mv signaledmvd(pu.mvdAffi[refList][ctrlP].getHor() >> mvShift, pu.mvdAffi[refList][ctrlP].getVer() >> mvShift);
          if (!((signaledmvd.getHor() >= MVD_MIN) && (signaledmvd.getHor() <= MVD_MAX)) || !((signaledmvd.getVer() >= MVD_MIN) && (signaledmvd.getVer() <= MVD_MAX)))
            return;
        }
      }
    }
  }
  // avoid MV exceeding 18-bit dynamic range
  const int maxMv = 1 << 17;
  if (!cu->affine && !pu.mergeFlag)
  {
    if ( (pu.refIdx[0] >= 0 && (pu.mv[0].getAbsHor() >= maxMv || pu.mv[0].getAbsVer() >= maxMv))
      || (pu.refIdx[1] >= 0 && (pu.mv[1].getAbsHor() >= maxMv || pu.mv[1].getAbsVer() >= maxMv)))
    {
      return;
    }
  }
  if (cu->affine && !pu.mergeFlag)
  {
    for (int refList = 0; refList < NUM_REF_PIC_LIST_01; refList++)
    {
      if (pu.refIdx[refList] >= 0)
      {
        for (int ctrlP = 1 + (cu->affineType == AFFINEMODEL_6PARAM); ctrlP >= 0; ctrlP--)
        {
          if (pu.mvAffi[refList][ctrlP].getAbsHor() >= maxMv || pu.mvAffi[refList][ctrlP].getAbsVer() >= maxMv)
          {
            return;
          }
        }
      }
    }
  }
  const bool mtsAllowed = tempCS->sps->getUseInterMTS() && CU::isInter( *cu ) && partitioner.currArea().lwidth() <= MTS_INTER_MAX_CU_SIZE && partitioner.currArea().lheight() <= MTS_INTER_MAX_CU_SIZE;
  uint8_t sbtAllowed = cu->checkAllowedSbt();
  //SBT resolution-dependent fast algorithm: not try size-64 SBT in RDO for low-resolution sequences (now resolution below HD)
  if( tempCS->pps->getPicWidthInLumaSamples() < (uint32_t)m_pcEncCfg->getSBTFast64WidthTh() )
  {
    sbtAllowed = ((cu->lwidth() > 32 || cu->lheight() > 32)) ? 0 : sbtAllowed;
  }
  uint8_t numRDOTried = 0;
  Distortion sbtOffDist = 0;
  bool    sbtOffRootCbf = 0;
  double  sbtOffCost      = MAX_DOUBLE;
  double  currBestCost = MAX_DOUBLE;
  bool    doPreAnalyzeResi = ( sbtAllowed || mtsAllowed ) && residualPass == 0;

  m_pcInterSearch->initTuAnalyzer();
  if( doPreAnalyzeResi )
  {
    m_pcInterSearch->calcMinDistSbt( *tempCS, *cu, sbtAllowed );
  }

  auto    slsSbt = dynamic_cast<SaveLoadEncInfoSbt*>( m_modeCtrl );
  int     slShift = 4 + std::min( (int)gp_sizeIdxInfo->idxFrom( cu->lwidth() ) + (int)gp_sizeIdxInfo->idxFrom( cu->lheight() ), 9 );
  Distortion curPuSse = m_pcInterSearch->getEstDistSbt( NUMBER_SBT_MODE );
  uint8_t currBestSbt = 0;
  uint8_t currBestTrs = MAX_UCHAR;
  uint8_t histBestSbt = MAX_UCHAR;
  uint8_t histBestTrs = MAX_UCHAR;
  m_pcInterSearch->setHistBestTrs( MAX_UCHAR, MAX_UCHAR );
  if( doPreAnalyzeResi )
  {
    if( m_pcInterSearch->getSkipSbtAll() && !mtsAllowed ) //emt is off
    {
      histBestSbt = 0; //try DCT2
      m_pcInterSearch->setHistBestTrs( histBestSbt, histBestTrs );
    }
    else
    {
      assert( curPuSse != std::numeric_limits<uint64_t>::max() );
      uint16_t compositeSbtTrs = slsSbt->findBestSbt( cu->cs->area, (uint32_t)( curPuSse >> slShift ) );
      histBestSbt = ( compositeSbtTrs >> 0 ) & 0xff;
      histBestTrs = ( compositeSbtTrs >> 8 ) & 0xff;
      if( m_pcInterSearch->getSkipSbtAll() && CU::isSbtMode( histBestSbt ) ) //special case, skip SBT when loading SBT
      {
        histBestSbt = 0; //try DCT2
      }
      m_pcInterSearch->setHistBestTrs( histBestSbt, histBestTrs );
    }
  }

  {
    if( reloadCU )
    {
      if( bestCost == bestCS->cost ) //The first EMT pass didn't become the bestCS, so we clear the TUs generated
      {
        tempCS->clearTUs();
      }
      else if( false == swapped )
      {
        tempCS->initStructData( encTestMode.qp );
        tempCS->copyStructure( *bestCS, partitioner.chType );
        tempCS->getPredBuf().copyFrom( bestCS->getPredBuf() );
        bestCost = bestCS->cost;
        cu       = tempCS->getCU( partitioner.chType );
        swapped = true;
      }
      else
      {
        tempCS->clearTUs();
        bestCost = bestCS->cost;
        cu       = tempCS->getCU( partitioner.chType );
      }

      //we need to restart the distortion for the new tempCS, the bit count and the cost
      tempCS->dist     = 0;
      tempCS->fracBits = 0;
      tempCS->cost     = MAX_DOUBLE;
      tempCS->costDbOffset = 0;
    }

    reloadCU    = true; // enable cu reloading
    cu->skip    = false;
    cu->sbtInfo = 0;

    const bool skipResidual = residualPass == 1;
    if( skipResidual || histBestSbt == MAX_UCHAR || !CU::isSbtMode( histBestSbt ) )
    {
    m_pcInterSearch->encodeResAndCalcRdInterCU( *tempCS, partitioner, skipResidual );
    if (tempCS->slice->getSPS()->getUseColorTrans())
    {
      bestCS->tmpColorSpaceCost = tempCS->tmpColorSpaceCost;
      bestCS->firstColorSpaceSelected = tempCS->firstColorSpaceSelected;
    }
    numRDOTried += mtsAllowed ? 2 : 1;
    xEncodeDontSplit( *tempCS, partitioner );

    xCheckDQP( *tempCS, partitioner );
    xCheckChromaQPOffset( *tempCS, partitioner );


    if( NULL != bestHasNonResi && (bestCostInternal > tempCS->cost) )
    {
      bestCostInternal = tempCS->cost;
      if (!(tempCS->getPU(partitioner.chType)->ciipFlag))
      *bestHasNonResi  = !cu->rootCbf;
    }

    if (cu->rootCbf == false)
    {
      if (tempCS->getPU(partitioner.chType)->ciipFlag)
      {
        tempCS->cost = MAX_DOUBLE;
        tempCS->costDbOffset = 0;
        return;
      }
    }
    currBestCost = tempCS->cost;
    sbtOffCost = tempCS->cost;
    sbtOffDist = tempCS->dist;
    sbtOffRootCbf = cu->rootCbf;
    currBestSbt = CU::getSbtInfo(cu->firstTU->mtsIdx[COMPONENT_Y] > MTS_SKIP ? SBT_OFF_MTS : SBT_OFF_DCT, 0);
    currBestTrs = cu->firstTU->mtsIdx[COMPONENT_Y];

#if WCG_EXT
    DTRACE_MODE_COST( *tempCS, m_pcRdCost->getLambda( true ) );
#else
    DTRACE_MODE_COST( *tempCS, m_pcRdCost->getLambda() );
#endif
    xCheckBestMode( tempCS, bestCS, partitioner, encTestMode );

    }

    uint8_t numSbtRdo = CU::numSbtModeRdo( sbtAllowed );
    //early termination if all SBT modes are not allowed
    //normative
    if( !sbtAllowed || skipResidual )
    {
      numSbtRdo = 0;
    }
    //fast algorithm
    if( ( histBestSbt != MAX_UCHAR && !CU::isSbtMode( histBestSbt ) ) || m_pcInterSearch->getSkipSbtAll() )
    {
      numSbtRdo = 0;
    }
    if( bestCost != MAX_DOUBLE && sbtOffCost != MAX_DOUBLE )
    {
      double th = 1.07;
      if( !( prevBestSbt == 0 || m_sbtCostSave[0] == MAX_DOUBLE ) )
      {
        assert( m_sbtCostSave[1] <= m_sbtCostSave[0] );
        th *= ( m_sbtCostSave[0] / m_sbtCostSave[1] );
      }
      if( sbtOffCost > bestCost * th )
      {
        numSbtRdo = 0;
      }
    }
    if( !sbtOffRootCbf && sbtOffCost != MAX_DOUBLE )
    {
      double th = Clip3( 0.05, 0.55, ( 27 - cu->qp ) * 0.02 + 0.35 );
      if( sbtOffCost < m_pcRdCost->calcRdCost( ( cu->lwidth() * cu->lheight() ) << SCALE_BITS, 0 ) * th )
      {
        numSbtRdo = 0;
      }
    }

    if( histBestSbt != MAX_UCHAR && numSbtRdo != 0 )
    {
      numSbtRdo = 1;
      m_pcInterSearch->initSbtRdoOrder( CU::getSbtMode( CU::getSbtIdx( histBestSbt ), CU::getSbtPos( histBestSbt ) ) );
    }

    for( int sbtModeIdx = 0; sbtModeIdx < numSbtRdo; sbtModeIdx++ )
    {
      uint8_t sbtMode = m_pcInterSearch->getSbtRdoOrder( sbtModeIdx );
      uint8_t sbtIdx = CU::getSbtIdxFromSbtMode( sbtMode );
      uint8_t sbtPos = CU::getSbtPosFromSbtMode( sbtMode );

      //fast algorithm (early skip, save & load)
      if( histBestSbt == MAX_UCHAR )
      {
        uint8_t skipCode = m_pcInterSearch->skipSbtByRDCost( cu->lwidth(), cu->lheight(), cu->mtDepth, sbtIdx, sbtPos, bestCS->cost, sbtOffDist, sbtOffCost, sbtOffRootCbf );
        if( skipCode != MAX_UCHAR )
        {
          continue;
        }

        if( sbtModeIdx > 0 )
        {
          uint8_t prevSbtMode = m_pcInterSearch->getSbtRdoOrder( sbtModeIdx - 1 );
          //make sure the prevSbtMode is the same size as the current SBT mode (otherwise the estimated dist may not be comparable)
          if( CU::isSameSbtSize( prevSbtMode, sbtMode ) )
          {
            Distortion currEstDist = m_pcInterSearch->getEstDistSbt( sbtMode );
            Distortion prevEstDist = m_pcInterSearch->getEstDistSbt( prevSbtMode );
            if( currEstDist > prevEstDist * 1.15 )
            {
              continue;
            }
          }
        }
      }

      //init tempCS and TU
      if( bestCost == bestCS->cost ) //The first EMT pass didn't become the bestCS, so we clear the TUs generated
      {
        tempCS->clearTUs();
      }
      else if( false == swapped )
      {
        tempCS->initStructData( encTestMode.qp );
        tempCS->copyStructure( *bestCS, partitioner.chType );
        tempCS->getPredBuf().copyFrom( bestCS->getPredBuf() );
        bestCost = bestCS->cost;
        cu = tempCS->getCU( partitioner.chType );
        swapped = true;
      }
      else
      {
        tempCS->clearTUs();
        bestCost = bestCS->cost;
        cu = tempCS->getCU( partitioner.chType );
      }

      //we need to restart the distortion for the new tempCS, the bit count and the cost
      tempCS->dist = 0;
      tempCS->fracBits = 0;
      tempCS->cost = MAX_DOUBLE;
      cu->skip = false;

      //set SBT info
      cu->setSbtIdx( sbtIdx );
      cu->setSbtPos( sbtPos );

      //try residual coding
      m_pcInterSearch->encodeResAndCalcRdInterCU( *tempCS, partitioner, skipResidual );
      if (tempCS->slice->getSPS()->getUseColorTrans())
      {
        bestCS->tmpColorSpaceCost = tempCS->tmpColorSpaceCost;
        bestCS->firstColorSpaceSelected = tempCS->firstColorSpaceSelected;
      }
      numRDOTried++;

      xEncodeDontSplit( *tempCS, partitioner );

      xCheckDQP( *tempCS, partitioner );
      xCheckChromaQPOffset( *tempCS, partitioner );

      if( NULL != bestHasNonResi && ( bestCostInternal > tempCS->cost ) )
      {
        bestCostInternal = tempCS->cost;
        if( !( tempCS->getPU( partitioner.chType )->ciipFlag ) )
          *bestHasNonResi = !cu->rootCbf;
      }

      if( tempCS->cost < currBestCost )
      {
        currBestSbt = cu->sbtInfo;
        currBestTrs = tempCS->tus[cu->sbtInfo ? cu->getSbtPos() : 0]->mtsIdx[COMPONENT_Y];
        assert( currBestTrs == 0 || currBestTrs == 1 );
        currBestCost = tempCS->cost;
      }

#if WCG_EXT
      DTRACE_MODE_COST( *tempCS, m_pcRdCost->getLambda( true ) );
#else
      DTRACE_MODE_COST( *tempCS, m_pcRdCost->getLambda() );
#endif
      xCheckBestMode( tempCS, bestCS, partitioner, encTestMode );
    }

    if( bestCostBegin != bestCS->cost )
    {
      m_sbtCostSave[0] = sbtOffCost;
      m_sbtCostSave[1] = currBestCost;
    }
  } //end emt loop

  if( histBestSbt == MAX_UCHAR && doPreAnalyzeResi && numRDOTried > 1 )
  {
    slsSbt->saveBestSbt( cu->cs->area, (uint32_t)( curPuSse >> slShift ), currBestSbt, currBestTrs );
  }
  tempCS->cost = currBestCost;
  if( ETM_INTER_ME == encTestMode.type )
  {
    if( equBcwCost != NULL )
    {
      if( tempCS->cost < ( *equBcwCost ) && cu->BcwIdx == BCW_DEFAULT )
      {
        ( *equBcwCost ) = tempCS->cost;
      }
    }
    else
    {
#if ENABLE_OBMC
      if (cu->obmcFlag)
#endif
      CHECK( equBcwCost == NULL, "equBcwCost == NULL" );
    }
    if( tempCS->slice->getCheckLDC() && !cu->imv && cu->BcwIdx != BCW_DEFAULT && tempCS->cost < m_bestBcwCost[1] )
    {
      if( tempCS->cost < m_bestBcwCost[0] )
      {
        m_bestBcwCost[1] = m_bestBcwCost[0];
        m_bestBcwCost[0] = tempCS->cost;
        m_bestBcwIdx[1] = m_bestBcwIdx[0];
        m_bestBcwIdx[0] = cu->BcwIdx;
      }
      else
      {
        m_bestBcwCost[1] = tempCS->cost;
        m_bestBcwIdx[1] = cu->BcwIdx;
      }
    }
#if INTER_LIC
    m_pcInterSearch->m_fastLicCtrl.setBestAmvpRDBeforeLIC(*cu, currBestCost);
#endif
  }
}


void EncCu::xEncodeDontSplit( CodingStructure &cs, Partitioner &partitioner )
{
  m_CABACEstimator->resetBits();

  m_CABACEstimator->split_cu_mode( CU_DONT_SPLIT, cs, partitioner );
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  if( partitioner.treeType == TREE_C )
    CHECK( m_CABACEstimator->getEstFracBits() != 0, "must be 0 bit" );
#endif
  cs.fracBits += m_CABACEstimator->getEstFracBits(); // split bits
  cs.cost      = m_pcRdCost->calcRdCost( cs.fracBits, cs.dist );

}

#if REUSE_CU_RESULTS
void EncCu::xReuseCachedResult( CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &partitioner )
{
  m_pcRdCost->setChromaFormat(tempCS->sps->getChromaFormatIdc());
  BestEncInfoCache* bestEncCache = dynamic_cast<BestEncInfoCache*>( m_modeCtrl );
  CHECK( !bestEncCache, "If this mode is chosen, mode controller has to implement the mode caching capabilities" );
  EncTestMode cachedMode;

  if( bestEncCache->setCsFrom( *tempCS, cachedMode, partitioner ) )
  {
    CodingUnit& cu = *tempCS->cus.front();
    partitioner.setCUData( cu );
#if MULTI_PASS_DMVR
    for( auto &pu : CU::traversePUs( cu ) )
    {
      pu.bdmvrRefine = false;
    }
#endif
    if( CU::isIntra( cu )
    || CU::isPLT(cu)
    )
    {
      xReconIntraQT( cu );
    }
    else
    {
      xDeriveCUMV( cu );
      xReconInter( cu );
    }

    Distortion finalDistortion = 0;
    tempCS->useDbCost = m_pcEncCfg->getUseEncDbOpt();
    if ( m_pcEncCfg->getUseEncDbOpt() )
    {
      xCalDebCost( *tempCS, partitioner, true );
      finalDistortion = tempCS->dist;
    }
    else
    {
    const SPS &sps = *tempCS->sps;
    const int  numValidComponents = getNumberValidComponents( tempCS->area.chromaFormat );

    for( int comp = 0; comp < numValidComponents; comp++ )
    {
      const ComponentID compID = ComponentID( comp );
#if INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
      if( CS::isDualITree(*tempCS) && toChannelType(compID) != partitioner.chType )
#else
      if( partitioner.isSepTree( *tempCS ) && toChannelType( compID ) != partitioner.chType )
#endif
      {
        continue;
      }

      CPelBuf reco = tempCS->getRecoBuf( compID );
      CPelBuf org  = tempCS->getOrgBuf ( compID );

      
#if JVET_V0094_BILATERAL_FILTER
        const CompArea &area = cu.blocks[COMPONENT_Y];
        CompArea    tmpArea(COMPONENT_Y, area.chromaFormat, Position(0, 0), area.size());
        PelBuf tmpRecLuma = m_tmpStorageLCU->getBuf(tmpArea);
        if(isLuma(compID))
        {
          tmpRecLuma.copyFrom(reco);

          if (m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled() || (
            m_pcEncCfg->getLmcs() && (tempCS->slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag())))
          {
            tmpRecLuma.rspSignal(m_pcReshape->getInvLUT());
          }
        }
        if(tempCS->pps->getUseBIF() && isLuma(compID) && (cu.qp > 17))
        {
          for (auto &currTU : CU::traverseTUs(cu))
          {
            Position tuPosInCu = currTU.lumaPos() - cu.lumaPos();
            PelBuf tmpSubBuf = tmpRecLuma.subBuf(tuPosInCu, currTU.lumaSize());

            bool isInter = (cu.predMode == MODE_INTRA) ? false : true;
            
            if ((TU::getCbf(currTU, COMPONENT_Y) || isInter == false) && (currTU.cu->qp > 17) && (128 > std::max(currTU.lumaSize().width, currTU.lumaSize().height)) && ((isInter == false) || (32 > std::min(currTU.lumaSize().width, currTU.lumaSize().height))))
            {
              CompArea compArea = currTU.blocks[compID];
              PelBuf recIPredBuf = tempCS->slice->getPic()->getRecoBuf(compArea);
              // Do we need to use clipArea?
              
              // Only reshape surrounding samples if reshaping is on
               if(m_pcEncCfg->getLmcs() && (tempCS->slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag() ) && !(m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled()))
               {
                 m_bilateralFilter->bilateralFilterRDOdiamond5x5(tmpSubBuf, tmpSubBuf, tmpSubBuf, currTU.cu->qp, recIPredBuf, tempCS->slice->clpRng(compID), currTU, true, true, m_pcReshape->getInvLUT());
               }
               else
               {
                 std::vector<Pel> dummy_invLUT;
                 m_bilateralFilter->bilateralFilterRDOdiamond5x5(tmpSubBuf, tmpSubBuf, tmpSubBuf, currTU.cu->qp, recIPredBuf, tempCS->slice->clpRng(compID), currTU, true, false, dummy_invLUT);
               }
            }
          }
        }

#if WCG_EXT
        if (m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled() || (
          m_pcEncCfg->getLmcs() && (tempCS->slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag())))
        {
          const CPelBuf orgLuma = tempCS->getOrgBuf(tempCS->area.blocks[COMPONENT_Y]);
          if (compID == COMPONENT_Y && !(m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled()))
          {
            finalDistortion += m_pcRdCost->getDistPart(org, tmpRecLuma, sps.getBitDepth(toChannelType(compID)), compID, DF_SSE_WTD, &orgLuma);
          }
          else
          {
            finalDistortion += m_pcRdCost->getDistPart( org, reco, sps.getBitDepth( toChannelType( compID ) ), compID, DF_SSE_WTD, &orgLuma );
          }
        }
        else
#endif
        {
          finalDistortion += m_pcRdCost->getDistPart( org, reco, sps.getBitDepth( toChannelType( compID ) ), compID, DF_SSE );
        }
#else
#if WCG_EXT
      if (m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled() || (m_pcEncCfg->getLmcs() && (tempCS->slice->getLmcsEnabledFlag() && m_pcReshape->getCTUFlag())))
      {
        const CPelBuf orgLuma = tempCS->getOrgBuf(tempCS->area.blocks[COMPONENT_Y]);
        if (compID == COMPONENT_Y && !(m_pcEncCfg->getLumaLevelToDeltaQPMapping().isEnabled()))
        {
          const CompArea &area = cu.blocks[COMPONENT_Y];
          CompArea    tmpArea(COMPONENT_Y, area.chromaFormat, Position(0, 0), area.size());
          PelBuf tmpRecLuma = m_tmpStorageLCU->getBuf(tmpArea);
          tmpRecLuma.rspSignal( reco, m_pcReshape->getInvLUT() );
          finalDistortion += m_pcRdCost->getDistPart(org, tmpRecLuma, sps.getBitDepth(toChannelType(compID)), compID, DF_SSE_WTD, &orgLuma);
        }
        else
        {
          finalDistortion += m_pcRdCost->getDistPart( org, reco, sps.getBitDepth( toChannelType( compID ) ), compID, DF_SSE_WTD, &orgLuma );
        }
      }
      else
#endif
      {
        finalDistortion += m_pcRdCost->getDistPart( org, reco, sps.getBitDepth( toChannelType( compID ) ), compID, DF_SSE );
      }
#endif
    }
    }

    m_CABACEstimator->getCtx() = m_CurrCtx->start;
    m_CABACEstimator->resetBits();

    CUCtx cuCtx;
    cuCtx.isDQPCoded = true;
    cuCtx.isChromaQpAdjCoded = true;
    m_CABACEstimator->coding_unit( cu, partitioner, cuCtx );


    tempCS->dist     = finalDistortion;
    tempCS->fracBits = m_CABACEstimator->getEstFracBits();
    tempCS->cost     = m_pcRdCost->calcRdCost( tempCS->fracBits, tempCS->dist );

    xEncodeDontSplit( *tempCS,         partitioner );
    xCheckDQP       ( *tempCS,         partitioner );
    xCheckChromaQPOffset( *tempCS,     partitioner );
    xCheckBestMode  (  tempCS, bestCS, partitioner, cachedMode );
  }
  else
  {
    THROW( "Should never happen!" );
  }
}
#endif

#if MULTI_HYP_PRED
void EncCu::predInterSearchAdditionalHypothesisMulti(const MEResultVec& in, MEResultVec& out, PredictionUnit& pu, const MergeCtx &mrgCtx)
{
  for (const auto &x : in)
  {
    *pu.cu = x.cu;
    pu = x.pu;

    if (pu.mergeType == MRG_TYPE_SUBPU_ATMVP)
    {
      // the SbTmvp use xSubPuMC which will need to access the motion buffer for subblock MV
      PU::spanMotionInfo(pu, mrgCtx);
    }
    else if (x.cu.affine)
    {
      PU::setAllAffineMv(pu, pu.mvAffi[0][0], pu.mvAffi[0][1], pu.mvAffi[0][2], REF_PIC_LIST_0);
      PU::setAllAffineMv(pu, pu.mvAffi[1][0], pu.mvAffi[1][1], pu.mvAffi[1][2], REF_PIC_LIST_1);
    }
#if MULTI_PASS_DMVR
    else if( pu.bdmvrRefine )
    {
#if TM_MRG
      if( pu.tmMergeFlag )
      {
        m_pcInterSearch->setBdmvrSubPuMvBuf( m_mvBufBDMVR4TM[pu.mergeIdx << 1], m_mvBufBDMVR4TM[( pu.mergeIdx << 1 ) + 1] );
      }
      else
#endif
      {
        m_pcInterSearch->setBdmvrSubPuMvBuf( m_mvBufBDMVR[pu.mergeIdx << 1], m_mvBufBDMVR[( pu.mergeIdx << 1 ) + 1] );
      }
    }
#endif
    m_pcInterSearch->predInterSearchAdditionalHypothesis(pu, x, out);
  }
}

void EncCu::xCheckRDCostInterMultiHyp2Nx2N(CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &partitioner, const EncTestMode& encTestMode)
{
#if ENABLE_OBMC
  double bestOBMCCost = MAX_DOUBLE;
#endif
  if (tempCS->area.Y().area() <= MULTI_HYP_PRED_RESTRICT_BLOCK_SIZE || std::min(tempCS->area.Y().width, tempCS->area.Y().height) < MULTI_HYP_PRED_RESTRICT_MIN_WH)
  {
    return;
  }
  const SPS &sps = *tempCS->sps;
  CHECK(!sps.getUseInterMultiHyp(), "Multi Hyp is not active");
  CHECK(!tempCS->slice->isInterB(), "Multi Hyp only allowed in B slices");
  CHECK(encTestMode.opts != ETO_STANDARD, "unknown encoding option to EncCu::xCheckRDCostInterMultiHyp2Nx2N()");

  if ((m_pcEncCfg->getBaseQP() > 32) || (m_pcEncCfg->getBaseQP() > 27 && tempCS->slice->getTLayer() >= 4)) // KRNOTE: explicit QP tests
    return;

  if (m_modeCtrl->getFastDeltaQp())
  {
    if (tempCS->area.lumaSize().width > tempCS->pcv->fastDeltaQPCuMaxSize)
    {
      return; // only check necessary 2Nx2N Inter in fast deltaqp mode
    }
  }

  MEResultVec mhResults;
  const auto RDCostComp = [](const MEResult &x, const MEResult &y) { return x.cost < y.cost; };

  MergeCtx mrgCtx;
  if (sps.getSbTMVPEnabledFlag())
  {
    Size bufSize = g_miScaling.scale(tempCS->area.lumaSize());
    mrgCtx.subPuMvpMiBuf = MotionBuf(m_SubPuMiBuf, bufSize);
  }

  // Hadamard-based pre-search
  {
    MEResultVec base = m_baseResultsForMH;

    if (base.empty())
      return;

    tempCS->initStructData(encTestMode.qp);
    CodingUnit &cu = tempCS->addCU(tempCS->area, partitioner.chType);
    PredictionUnit &pu = tempCS->addPU(cu, partitioner.chType);

    int iter = 0;
    do
    {
      MEResultVec out;
      const auto survivors = iter > 0 ? 2 : m_pcEncCfg->getNumMHPCandsToTest();
      iter++;

      std::stable_sort(base.begin(), base.end(), RDCostComp);
      if (base.size() > survivors)
        base.resize(survivors);
      predInterSearchAdditionalHypothesisMulti(base, out, pu, mrgCtx);
      mhResults.insert(mhResults.end(), out.begin(), out.end());
      base = out;
    } while (!base.empty());
  }

  std::stable_sort(mhResults.begin(), mhResults.end(), RDCostComp);

  // actual testing with "true" RD costs
  for (int i = 0; i < std::min((int)mhResults.size(), m_pcEncCfg->getAddHypTries()); ++i)
  {
    tempCS->initStructData(encTestMode.qp);
    CodingUnit &cu = tempCS->addCU(tempCS->area, partitioner.chType);
    PredictionUnit &pu = tempCS->addPU(cu, partitioner.chType);

    pu = mhResults[i].pu;
    cu = mhResults[i].cu;

#if MULTI_PASS_DMVR
    if (pu.bdmvrRefine)
    {
#if TM_MRG
      if( pu.tmMergeFlag )
      {
        m_pcInterSearch->setBdmvrSubPuMvBuf( m_mvBufBDMVR4TM[pu.mergeIdx << 1], m_mvBufBDMVR4TM[( pu.mergeIdx << 1 ) + 1] );
      }
      else
#endif
      m_pcInterSearch->setBdmvrSubPuMvBuf(m_mvBufBDMVR[pu.mergeIdx << 1], m_mvBufBDMVR[(pu.mergeIdx << 1) + 1]);
    }
    else
    {
      PU::spanMotionInfo(pu, mrgCtx);
    }
#else
    PU::spanMotionInfo(pu, mrgCtx);
#endif
    cu.skip = false;
    cu.mmvdSkip = false;

    CHECK(cu.qtDepth != partitioner.currQtDepth, "Mismatch");
    CHECK(cu.btDepth != partitioner.currBtDepth, "Mismatch");
    CHECK(cu.mtDepth != partitioner.currMtDepth, "Mismatch");
    pu.mvRefine = true;
    m_pcInterSearch->motionCompensation(pu);

#if MULTI_PASS_DMVR
    if (pu.bdmvrRefine)
    {
#if TM_MRG
      if( pu.tmMergeFlag )
      {
        PU::spanMotionInfo( pu, mrgCtx, m_mvBufBDMVR4TM[pu.mergeIdx << 1], m_mvBufBDMVR4TM[( pu.mergeIdx << 1 ) + 1], m_pcInterSearch->getBdofSubPuMvOffset() );
      }
      else
#endif
      PU::spanMotionInfo(pu, mrgCtx, m_mvBufBDMVR[pu.mergeIdx << 1], m_mvBufBDMVR[(pu.mergeIdx << 1) + 1], m_pcInterSearch->getBdofSubPuMvOffset());
    }
#endif
    pu.mvRefine = false;

#if ENABLE_OBMC //multi hyp inter IMV
    CodingStructure *prevCS = tempCS;
    PelUnitBuf tempWoOBMCBuf = m_tempWoOBMCBuffer.subBuf(UnitAreaRelative(cu, cu));
    tempWoOBMCBuf.copyFrom(tempCS->getPredBuf(cu));
    cu.isobmcMC = true;
    cu.obmcFlag = true;
    m_pcInterSearch->subBlockOBMC(*cu.firstPU);
    cu.isobmcMC = false;
#endif
    xEncodeInterResidual(tempCS, bestCS, partitioner, encTestMode);
#if ENABLE_OBMC
    double tempCost = (prevCS == tempCS) ? tempCS->cost : bestCS->cost;
    if (m_pTempCUWoOBMC && tempCost < bestOBMCCost)
    {
      const unsigned wIdx = gp_sizeIdxInfo->idxFrom(tempCS->area.lwidth());
      const unsigned hIdx = gp_sizeIdxInfo->idxFrom(tempCS->area.lheight());

      m_pTempCUWoOBMC[wIdx][hIdx]->clearCUs();
      m_pTempCUWoOBMC[wIdx][hIdx]->clearPUs();
      m_pTempCUWoOBMC[wIdx][hIdx]->clearTUs();
      m_pTempCUWoOBMC[wIdx][hIdx]->copyStructure(*prevCS, partitioner.chType);

      m_pPredBufWoOBMC[wIdx][hIdx].copyFrom(tempWoOBMCBuf);
      m_pTempCUWoOBMC[wIdx][hIdx]->getPredBuf(cu).copyFrom(prevCS->getPredBuf(cu));

      bestOBMCCost = tempCost;
    }
#endif
  }
}
#endif
#if ENABLE_OBMC
void EncCu::xCheckRDCostInterWoOBMC(CodingStructure *&tempCS, CodingStructure *&bestCS, Partitioner &partitioner, const EncTestMode& encTestMode)
{
  if (!tempCS->sps->getUseOBMC())
  {
    return;
  }

  if (m_modeCtrl->getFastDeltaQp())
  {
    if (tempCS->area.lumaSize().width > tempCS->pcv->fastDeltaQPCuMaxSize)
    {
      return; // only check necessary 2Nx2N Inter in fast deltaqp mode
    }
  }

  tempCS->initStructData(encTestMode.qp);

  const SPS &sps = *tempCS->sps;
  const unsigned wIdx = gp_sizeIdxInfo->idxFrom(tempCS->area.lwidth());
  const unsigned hIdx = gp_sizeIdxInfo->idxFrom(tempCS->area.lheight());

  CodingStructure* CSWoOBMC = m_pTempCUWoOBMC[wIdx][hIdx];

  if(CSWoOBMC->cus.size() == 0)
    return;

  CodingUnit *cu = CSWoOBMC->getCU(partitioner.chType);

  if (
    !cu->obmcFlag
    || cu->predMode == MODE_INTRA
    || cu->firstPU->mergeFlag
    || CU::isIBC(*cu)
    || cu->geoFlag
#if INTER_LIC
    || cu->LICFlag
#endif
    )
  {
    return;
  }

  const Distortion uiSADOBMCOff = m_pcRdCost->getDistPart(tempCS->getOrgBuf(cu->Y()), m_pPredBufWoOBMC[wIdx][hIdx].Y(),
    sps.getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, DF_SAD_FULL_NBIT);
  const Distortion uiSADOBMCOn = m_pcRdCost->getDistPart(tempCS->getOrgBuf(cu->Y()), CSWoOBMC->getPredBuf(cu->Y()),
    sps.getBitDepth(CHANNEL_TYPE_LUMA), COMPONENT_Y, DF_SAD_FULL_NBIT);

  const double    dOBMCThOff = 1.0;
  const bool   bCheckOBMCOff = uiSADOBMCOff * dOBMCThOff < uiSADOBMCOn;

  if (!bCheckOBMCOff)
  {
    return;
  }

  tempCS->copyStructure(*CSWoOBMC, partitioner.chType);
  tempCS->getPredBuf(*cu).copyFrom(m_pPredBufWoOBMC[wIdx][hIdx]);
  cu = tempCS->getCU(partitioner.chType);
  cu->obmcFlag = false;
  //
  xEncodeInterResidual(tempCS, bestCS, partitioner, encTestMode, 0);
}
#endif
//! \}
