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

#include "DepQuant.h"
#include "TrQuant.h"
#include "CodingStructure.h"
#include "UnitTools.h"

#include <bitset>






namespace DQIntern
{
  /*================================================================================*/
  /*=====                                                                      =====*/
  /*=====   R A T E   E S T I M A T O R                                        =====*/
  /*=====                                                                      =====*/
  /*================================================================================*/

  struct NbInfoSbb
  {
    uint8_t   num;
    uint8_t   inPos[5];
  };
  struct NbInfoOut
  {
    uint16_t  maxDist;
    uint16_t  num;
    uint16_t  outPos[5];
  };
  struct CoeffFracBits
  {
    int32_t   bits[6];
  };


  enum ScanPosType { SCAN_ISCSBB = 0, SCAN_SOCSBB = 1, SCAN_EOCSBB = 2 };

  struct ScanInfo
  {
    ScanInfo() {}
    int           sbbSize;
    int           numSbb;
    int           scanIdx;
    int           rasterPos;
    int           sbbPos;
    int           insidePos;
    bool          eosbb;
    ScanPosType   spt;
    unsigned      sigCtxOffsetNext;
    unsigned      gtxCtxOffsetNext;
    int           nextInsidePos;
    NbInfoSbb     nextNbInfoSbb;
    int           nextSbbRight;
    int           nextSbbBelow;
    int           posX;
    int           posY;
    ChannelType   chType;
    int           sbtInfo;
    int           tuWidth;
    int           tuHeight;
  };

  class Rom;
  struct TUParameters
  {
    TUParameters ( const Rom& rom, const unsigned width, const unsigned height, const ChannelType chType );
    ~TUParameters()
    {
      delete [] m_scanInfo;
    }

    ChannelType       m_chType;
    unsigned          m_width;
    unsigned          m_height;
    unsigned          m_numCoeff;
    unsigned          m_numSbb;
    unsigned          m_log2SbbWidth;
    unsigned          m_log2SbbHeight;
    unsigned          m_log2SbbSize;
    unsigned          m_sbbSize;
    unsigned          m_sbbMask;
    unsigned          m_widthInSbb;
    unsigned          m_heightInSbb;
    CoeffScanType     m_scanType;
    const ScanElement *m_scanSbbId2SbbPos;
    const ScanElement *m_scanId2BlkPos;
    const NbInfoSbb*  m_scanId2NbInfoSbb;
    const NbInfoOut*  m_scanId2NbInfoOut;
    ScanInfo*         m_scanInfo;
  private:
    void xSetScanInfo( ScanInfo& scanInfo, int scanIdx );
  };

  class Rom
  {
  public:
    Rom() : m_scansInitialized(false) {}
    ~Rom() { xUninitScanArrays(); }
    void                init        ()                       { xInitScanArrays(); }
    const NbInfoSbb*    getNbInfoSbb( int hd, int vd ) const { return m_scanId2NbInfoSbbArray[hd][vd]; }
    const NbInfoOut*    getNbInfoOut( int hd, int vd ) const { return m_scanId2NbInfoOutArray[hd][vd]; }
    const TUParameters* getTUPars   ( const CompArea& area, const ComponentID compID ) const
    {
      return m_tuParameters[floorLog2(area.width)][floorLog2(area.height)][toChannelType(compID)];
    }
  private:
    void  xInitScanArrays   ();
    void  xUninitScanArrays ();
  private:
    bool          m_scansInitialized;
    NbInfoSbb*    m_scanId2NbInfoSbbArray[ MAX_CU_DEPTH+1 ][ MAX_CU_DEPTH+1 ];
    NbInfoOut*    m_scanId2NbInfoOutArray[ MAX_CU_DEPTH+1 ][ MAX_CU_DEPTH+1 ];
    TUParameters* m_tuParameters         [ MAX_CU_DEPTH+1 ][ MAX_CU_DEPTH+1 ][ MAX_NUM_CHANNEL_TYPE ];
  };

  void Rom::xInitScanArrays()
  {
    if( m_scansInitialized )
    {
      return;
    }
    ::memset( m_scanId2NbInfoSbbArray, 0, sizeof(m_scanId2NbInfoSbbArray) );
    ::memset( m_scanId2NbInfoOutArray, 0, sizeof(m_scanId2NbInfoOutArray) );
    ::memset( m_tuParameters,          0, sizeof(m_tuParameters) );

    uint32_t raster2id[ MAX_CU_SIZE * MAX_CU_SIZE ];
    ::memset(raster2id, 0, sizeof(raster2id));

    for( int hd = 0; hd <= MAX_CU_DEPTH; hd++ )
    {
      for( int vd = 0; vd <= MAX_CU_DEPTH; vd++ )
      {
        if( (hd == 0 && vd <= 1) || (hd <= 1 && vd == 0) )
        {
          continue;
        }
        const uint32_t      blockWidth    = (1 << hd);
        const uint32_t      blockHeight   = (1 << vd);
        const uint32_t      log2CGWidth   = g_log2SbbSize[hd][vd][0];
        const uint32_t      log2CGHeight  = g_log2SbbSize[hd][vd][1];
        const uint32_t      groupWidth    = 1 << log2CGWidth;
        const uint32_t      groupHeight   = 1 << log2CGHeight;
        const uint32_t      groupSize     = groupWidth * groupHeight;
        const CoeffScanType scanType      = SCAN_DIAG;
        const SizeType      blkWidthIdx   = gp_sizeIdxInfo->idxFrom( blockWidth  );
        const SizeType      blkHeightIdx  = gp_sizeIdxInfo->idxFrom( blockHeight );
        const ScanElement * scanId2RP     = g_scanOrder[SCAN_GROUPED_4x4][scanType][blkWidthIdx][blkHeightIdx];
        NbInfoSbb*&         sId2NbSbb     = m_scanId2NbInfoSbbArray[hd][vd];
        NbInfoOut*&         sId2NbOut     = m_scanId2NbInfoOutArray[hd][vd];
        // consider only non-zero-out region
        const uint32_t      blkWidthNZOut = std::min<unsigned>( JVET_C0024_ZERO_OUT_TH, blockWidth  );
        const uint32_t      blkHeightNZOut= std::min<unsigned>( JVET_C0024_ZERO_OUT_TH, blockHeight );
        const uint32_t      totalValues   = blkWidthNZOut * blkHeightNZOut;

        sId2NbSbb = new NbInfoSbb[ totalValues ];
        sId2NbOut = new NbInfoOut[ totalValues ];

        for( uint32_t scanId = 0; scanId < totalValues; scanId++ )
        {
          raster2id[scanId2RP[scanId].idx] = scanId;
        }

        for( unsigned scanId = 0; scanId < totalValues; scanId++ )
        {
          const int posX = scanId2RP[scanId].x;
          const int posY = scanId2RP[scanId].y;
          const int rpos = scanId2RP[scanId].idx;
          {
            //===== inside subband neighbours =====
            NbInfoSbb&     nbSbb  = sId2NbSbb[ scanId ];
            const int      begSbb = scanId - ( scanId & (groupSize-1) ); // first pos in current subblock
            int            cpos[5];

            cpos[0] = ( posX + 1 < blkWidthNZOut                              ? ( raster2id[rpos+1           ] < groupSize + begSbb ? raster2id[rpos+1           ] - begSbb : 0 ) : 0 );
            cpos[1] = ( posX + 2 < blkWidthNZOut                              ? ( raster2id[rpos+2           ] < groupSize + begSbb ? raster2id[rpos+2           ] - begSbb : 0 ) : 0 );
            cpos[2] = ( posX + 1 < blkWidthNZOut && posY + 1 < blkHeightNZOut ? ( raster2id[rpos+1+blockWidth] < groupSize + begSbb ? raster2id[rpos+1+blockWidth] - begSbb : 0 ) : 0 );
            cpos[3] = ( posY + 1 < blkHeightNZOut                             ? ( raster2id[rpos+  blockWidth] < groupSize + begSbb ? raster2id[rpos+  blockWidth] - begSbb : 0 ) : 0 );
            cpos[4] = ( posY + 2 < blkHeightNZOut                             ? ( raster2id[rpos+2*blockWidth] < groupSize + begSbb ? raster2id[rpos+2*blockWidth] - begSbb : 0 ) : 0 );

            for( nbSbb.num = 0; true; )
            {
              int nk = -1;
              for( int k = 0; k < 5; k++ )
              {
                if( cpos[k] != 0 && ( nk < 0 || cpos[k] < cpos[nk] ) )
                {
                  nk = k;
                }
              }
              if( nk < 0 )
              {
                break;
              }
              nbSbb.inPos[ nbSbb.num++ ] = uint8_t( cpos[nk] );
              cpos[nk] = 0;
            }
            for( int k = nbSbb.num; k < 5; k++ )
            {
              nbSbb.inPos[k] = 0;
            }
          }
          {
            //===== outside subband neighbours =====
            NbInfoOut&     nbOut  = sId2NbOut[ scanId ];
            const int      begSbb = scanId - ( scanId & (groupSize-1) ); // first pos in current subblock
            int            cpos[5];

            cpos[0] = ( posX + 1 < blkWidthNZOut                              ? ( raster2id[rpos+1           ] >= groupSize + begSbb ? raster2id[rpos+1           ] : 0 ) : 0 );
            cpos[1] = ( posX + 2 < blkWidthNZOut                              ? ( raster2id[rpos+2           ] >= groupSize + begSbb ? raster2id[rpos+2           ] : 0 ) : 0 );
            cpos[2] = ( posX + 1 < blkWidthNZOut && posY + 1 < blkHeightNZOut ? ( raster2id[rpos+1+blockWidth] >= groupSize + begSbb ? raster2id[rpos+1+blockWidth] : 0 ) : 0 );
            cpos[3] = ( posY + 1 < blkHeightNZOut                             ? ( raster2id[rpos+  blockWidth] >= groupSize + begSbb ? raster2id[rpos+  blockWidth] : 0 ) : 0 );
            cpos[4] = ( posY + 2 < blkHeightNZOut                             ? ( raster2id[rpos+2*blockWidth] >= groupSize + begSbb ? raster2id[rpos+2*blockWidth] : 0 ) : 0 );

            for( nbOut.num = 0; true; )
            {
              int nk = -1;
              for( int k = 0; k < 5; k++ )
              {
                if( cpos[k] != 0 && ( nk < 0 || cpos[k] < cpos[nk] ) )
                {
                  nk = k;
                }
              }
              if( nk < 0 )
              {
                break;
              }
              nbOut.outPos[ nbOut.num++ ] = uint16_t( cpos[nk] );
              cpos[nk] = 0;
            }
            for( int k = nbOut.num; k < 5; k++ )
            {
              nbOut.outPos[k] = 0;
            }
            nbOut.maxDist = ( scanId == 0 ? 0 : sId2NbOut[scanId-1].maxDist );
            for( int k = 0; k < nbOut.num; k++ )
            {
              if( nbOut.outPos[k] > nbOut.maxDist )
              {
                nbOut.maxDist = nbOut.outPos[k];
              }
            }
          }
        }

        // make it relative
        for( unsigned scanId = 0; scanId < totalValues; scanId++ )
        {
          NbInfoOut& nbOut  = sId2NbOut[scanId];
          const int  begSbb = scanId - ( scanId & (groupSize-1) ); // first pos in current subblock
          for( int k = 0; k < nbOut.num; k++ )
          {
            CHECK(begSbb > nbOut.outPos[k], "Position must be past sub block begin");
            nbOut.outPos[k] -= begSbb;
          }
          nbOut.maxDist -= scanId;
        }

        for( int chId = 0; chId < MAX_NUM_CHANNEL_TYPE; chId++ )
        {
          m_tuParameters[hd][vd][chId] = new TUParameters( *this, blockWidth, blockHeight, ChannelType(chId) );
        }
      }
    }
    m_scansInitialized = true;
  }

  void Rom::xUninitScanArrays()
  {
    if( !m_scansInitialized )
    {
      return;
    }
    for( int hd = 0; hd <= MAX_CU_DEPTH; hd++ )
    {
      for( int vd = 0; vd <= MAX_CU_DEPTH; vd++ )
      {
        NbInfoSbb*& sId2NbSbb = m_scanId2NbInfoSbbArray[hd][vd];
        NbInfoOut*& sId2NbOut = m_scanId2NbInfoOutArray[hd][vd];
        if( sId2NbSbb )
        {
          delete [] sId2NbSbb;
        }
        if( sId2NbOut )
        {
          delete [] sId2NbOut;
        }
        for( int chId = 0; chId < MAX_NUM_CHANNEL_TYPE; chId++ )
        {
          TUParameters*& tuPars = m_tuParameters[hd][vd][chId];
          if( tuPars )
          {
            delete tuPars;
          }
        }
      }
    }
    m_scansInitialized = false;
  }


  static Rom g_Rom;


  TUParameters::TUParameters( const Rom& rom, const unsigned width, const unsigned height, const ChannelType chType )
  {
    m_chType              = chType;
    m_width               = width;
    m_height              = height;
    const uint32_t nonzeroWidth  = std::min<uint32_t>(JVET_C0024_ZERO_OUT_TH, m_width);
    const uint32_t nonzeroHeight = std::min<uint32_t>(JVET_C0024_ZERO_OUT_TH, m_height);
    m_numCoeff                   = nonzeroWidth * nonzeroHeight;
    const int log2W       = floorLog2( m_width  );
    const int log2H       = floorLog2( m_height );
    m_log2SbbWidth        = g_log2SbbSize[ log2W ][ log2H ][0];
    m_log2SbbHeight       = g_log2SbbSize[ log2W ][ log2H ][1];
    m_log2SbbSize         = m_log2SbbWidth + m_log2SbbHeight;
    m_sbbSize             = ( 1 << m_log2SbbSize );
    m_sbbMask             = m_sbbSize - 1;
    m_widthInSbb  = nonzeroWidth >> m_log2SbbWidth;
    m_heightInSbb = nonzeroHeight >> m_log2SbbHeight;
    m_numSbb              = m_widthInSbb * m_heightInSbb;
    m_scanType            = SCAN_DIAG;
    SizeType        hsbb  = gp_sizeIdxInfo->idxFrom( m_widthInSbb  );
    SizeType        vsbb  = gp_sizeIdxInfo->idxFrom( m_heightInSbb );
    SizeType        hsId  = gp_sizeIdxInfo->idxFrom( m_width  );
    SizeType        vsId  = gp_sizeIdxInfo->idxFrom( m_height );
    m_scanSbbId2SbbPos    = g_scanOrder     [ SCAN_UNGROUPED   ][ m_scanType ][ hsbb ][ vsbb ];
    m_scanId2BlkPos       = g_scanOrder     [ SCAN_GROUPED_4x4 ][ m_scanType ][ hsId ][ vsId ];
    m_scanId2NbInfoSbb    = rom.getNbInfoSbb( log2W, log2H );
    m_scanId2NbInfoOut    = rom.getNbInfoOut( log2W, log2H );
    m_scanInfo            = new ScanInfo[ m_numCoeff ];
    for( int scanIdx = 0; scanIdx < m_numCoeff; scanIdx++ )
    {
      xSetScanInfo( m_scanInfo[scanIdx], scanIdx );
    }
  }


  void TUParameters::xSetScanInfo( ScanInfo& scanInfo, int scanIdx )
  {
    scanInfo.chType     = m_chType;
    scanInfo.tuWidth    = m_width;
    scanInfo.tuHeight   = m_height;
    scanInfo.sbbSize    = m_sbbSize;
    scanInfo.numSbb     = m_numSbb;
    scanInfo.scanIdx    = scanIdx;
    scanInfo.rasterPos  = m_scanId2BlkPos[scanIdx].idx;
    scanInfo.sbbPos     = m_scanSbbId2SbbPos[scanIdx >> m_log2SbbSize].idx;
    scanInfo.insidePos  = scanIdx & m_sbbMask;
    scanInfo.eosbb      = ( scanInfo.insidePos == 0 );
    scanInfo.spt        = SCAN_ISCSBB;
    if(  scanInfo.insidePos == m_sbbMask && scanIdx > scanInfo.sbbSize && scanIdx < m_numCoeff - 1 )
      scanInfo.spt      = SCAN_SOCSBB;
    else if( scanInfo.eosbb && scanIdx > 0 && scanIdx < m_numCoeff - m_sbbSize )
      scanInfo.spt      = SCAN_EOCSBB;
    scanInfo.posX = m_scanId2BlkPos[scanIdx].x;
    scanInfo.posY = m_scanId2BlkPos[scanIdx].y;
    if( scanIdx )
    {
      const int nextScanIdx = scanIdx - 1;
      const int diag        = m_scanId2BlkPos[nextScanIdx].x + m_scanId2BlkPos[nextScanIdx].y;
      if( m_chType == CHANNEL_TYPE_LUMA )
      {
        scanInfo.sigCtxOffsetNext = ( diag < 2 ? 8 : diag < 5 ?  4 : 0 );
        scanInfo.gtxCtxOffsetNext = ( diag < 1 ? 16 : diag < 3 ? 11 : diag < 10 ? 6 : 1 );
      }
      else
      {
        scanInfo.sigCtxOffsetNext = ( diag < 2 ? 4 : 0 );
        scanInfo.gtxCtxOffsetNext = ( diag < 1 ? 6 : 1 );
      }
      scanInfo.nextInsidePos      = nextScanIdx & m_sbbMask;
      scanInfo.nextNbInfoSbb      = m_scanId2NbInfoSbb[ nextScanIdx ];
      if( scanInfo.eosbb )
      {
        const int nextSbbPos  = m_scanSbbId2SbbPos[nextScanIdx >> m_log2SbbSize].idx;
        const int nextSbbPosY = nextSbbPos               / m_widthInSbb;
        const int nextSbbPosX = nextSbbPos - nextSbbPosY * m_widthInSbb;
        scanInfo.nextSbbRight = ( nextSbbPosX < m_widthInSbb  - 1 ? nextSbbPos + 1            : 0 );
        scanInfo.nextSbbBelow = ( nextSbbPosY < m_heightInSbb - 1 ? nextSbbPos + m_widthInSbb : 0 );
      }
    }
  }



  class RateEstimator
  {
  public:
    RateEstimator () {}
    ~RateEstimator() {}
    void initCtx  ( const TUParameters& tuPars, const TransformUnit& tu, const ComponentID compID, const FracBitsAccess& fracBitsAccess );

    inline const BinFracBits *sigSbbFracBits() const { return m_sigSbbFracBits; }
    inline const BinFracBits *sigFlagBits(unsigned stateId) const
    {
#if TCQ_8STATES
			return m_sigFracBits[(stateId & 1) ? 1 + ((stateId & 3) >> 1) : 0];
#else
      return m_sigFracBits[std::max(((int) stateId) - 1, 0)];
#endif
    }
    inline const CoeffFracBits *gtxFracBits(unsigned stateId) const { return m_gtxFracBits; }
    inline int32_t              lastOffset(unsigned scanIdx) const
    {
      return m_lastBitsX[m_scanId2Pos[scanIdx].x] + m_lastBitsY[m_scanId2Pos[scanIdx].y];
    }

  private:
    void  xSetLastCoeffOffset ( const FracBitsAccess& fracBitsAccess, const TUParameters& tuPars, const TransformUnit& tu, const ComponentID compID );
    void  xSetSigSbbFracBits  ( const FracBitsAccess& fracBitsAccess, ChannelType chType );
    void  xSetSigFlagBits     ( const FracBitsAccess& fracBitsAccess, ChannelType chType );
    void  xSetGtxFlagBits     ( const FracBitsAccess& fracBitsAccess, ChannelType chType );

  private:
    static const unsigned sm_numCtxSetsSig    = 3;
    static const unsigned sm_numCtxSetsGtx    = 2;
    static const unsigned sm_maxNumSigSbbCtx  = 2;
    static const unsigned sm_maxNumSigCtx     = 12;
    static const unsigned sm_maxNumGtxCtx     = 21;

  private:
    const ScanElement * m_scanId2Pos;
    int32_t             m_lastBitsX      [ MAX_TB_SIZEY ];
    int32_t             m_lastBitsY      [ MAX_TB_SIZEY ];
    BinFracBits         m_sigSbbFracBits [ sm_maxNumSigSbbCtx ];
    BinFracBits         m_sigFracBits    [ sm_numCtxSetsSig   ][ sm_maxNumSigCtx ];
    CoeffFracBits       m_gtxFracBits                          [ sm_maxNumGtxCtx ];
  };

  void RateEstimator::initCtx( const TUParameters& tuPars, const TransformUnit& tu, const ComponentID compID, const FracBitsAccess& fracBitsAccess )
  {
    m_scanId2Pos = tuPars.m_scanId2BlkPos;
    xSetSigSbbFracBits  ( fracBitsAccess, tuPars.m_chType );
    xSetSigFlagBits     ( fracBitsAccess, tuPars.m_chType );
    xSetGtxFlagBits     ( fracBitsAccess, tuPars.m_chType );
    xSetLastCoeffOffset ( fracBitsAccess, tuPars, tu, compID );
  }

  void RateEstimator::xSetLastCoeffOffset( const FracBitsAccess& fracBitsAccess, const TUParameters& tuPars, const TransformUnit& tu, const ComponentID compID )
  {
    const ChannelType chType = ( compID == COMPONENT_Y ? CHANNEL_TYPE_LUMA : CHANNEL_TYPE_CHROMA );
    int32_t cbfDeltaBits = 0;
    if( compID == COMPONENT_Y && !CU::isIntra(*tu.cu) && !tu.depth )
    {
      const BinFracBits bits  = fracBitsAccess.getFracBitsArray( Ctx::QtRootCbf() );
      cbfDeltaBits            = int32_t( bits.intBits[1] ) - int32_t( bits.intBits[0] );
    }
    else
    {
      BinFracBits bits;
      bool prevLumaCbf           = false;
      bool lastCbfIsInferred     = false;
      bool useIntraSubPartitions = tu.cu->ispMode && isLuma(chType);
      if( useIntraSubPartitions )
      {
        bool rootCbfSoFar = false;
        bool isLastSubPartition = CU::isISPLast(*tu.cu, tu.Y(), compID);
        uint32_t nTus = tu.cu->ispMode == HOR_INTRA_SUBPARTITIONS ? tu.cu->lheight() >> floorLog2(tu.lheight()) : tu.cu->lwidth() >> floorLog2(tu.lwidth());
        if( isLastSubPartition )
        {
          TransformUnit* tuPointer = tu.cu->firstTU;
          for( int tuIdx = 0; tuIdx < nTus - 1; tuIdx++ )
          {
            rootCbfSoFar |= TU::getCbfAtDepth(*tuPointer, COMPONENT_Y, tu.depth);
            tuPointer     = tuPointer->next;
          }
          if( !rootCbfSoFar )
          {
            lastCbfIsInferred = true;
          }
        }
        if( !lastCbfIsInferred )
        {
          prevLumaCbf = TU::getPrevTuCbfAtDepth(tu, compID, tu.depth);
        }
        bits = fracBitsAccess.getFracBitsArray(Ctx::QtCbf[compID](DeriveCtx::CtxQtCbf(compID, prevLumaCbf, true)));
      }
      else
      {
        bits = fracBitsAccess.getFracBitsArray(Ctx::QtCbf[compID](DeriveCtx::CtxQtCbf(compID, tu.cbf[COMPONENT_Cb])));
      }
      cbfDeltaBits = lastCbfIsInferred ? 0 : int32_t(bits.intBits[1]) - int32_t(bits.intBits[0]);
    }

    static const unsigned prefixCtx[] = { 0, 0, 0, 3, 6, 10, 15, 21 };
    uint32_t              ctxBits  [ LAST_SIGNIFICANT_GROUPS ];
    for( unsigned xy = 0; xy < 2; xy++ )
    {
      int32_t             bitOffset   = ( xy ? cbfDeltaBits : 0 );
      int32_t*            lastBits    = ( xy ? m_lastBitsY : m_lastBitsX );
      const unsigned      size        = ( xy ? tuPars.m_height : tuPars.m_width );
      const unsigned      log2Size    = ceilLog2( size );
      const bool          useYCtx     = ( xy != 0 );
      const CtxSet&       ctxSetLast  = ( useYCtx ? Ctx::LastY : Ctx::LastX )[ chType ];
      const unsigned      lastShift   = ( compID == COMPONENT_Y ? (log2Size+1)>>2 : Clip3<unsigned>(0,2,size>>3) );
      const unsigned      lastOffset  = ( compID == COMPONENT_Y ? ( prefixCtx[log2Size] ) : 0 );
      uint32_t            sumFBits    = 0;
      unsigned            maxCtxId    = g_uiGroupIdx[std::min<unsigned>(JVET_C0024_ZERO_OUT_TH, size) - 1];
      for( unsigned ctxId = 0; ctxId < maxCtxId; ctxId++ )
      {
        const BinFracBits bits  = fracBitsAccess.getFracBitsArray( ctxSetLast( lastOffset + ( ctxId >> lastShift ) ) );
        ctxBits[ ctxId ]        = sumFBits + bits.intBits[0] + ( ctxId>3 ? ((ctxId-2)>>1)<<SCALE_BITS : 0 ) + bitOffset;
        sumFBits               +=            bits.intBits[1];
      }
      ctxBits  [ maxCtxId ]     = sumFBits + ( maxCtxId>3 ? ((maxCtxId-2)>>1)<<SCALE_BITS : 0 ) + bitOffset;
      for (unsigned pos = 0; pos < std::min<unsigned>(JVET_C0024_ZERO_OUT_TH, size); pos++)
      {
        lastBits[ pos ]         = ctxBits[ g_uiGroupIdx[ pos ] ];
      }
    }
  }

  void RateEstimator::xSetSigSbbFracBits( const FracBitsAccess& fracBitsAccess, ChannelType chType )
  {
    const CtxSet& ctxSet = Ctx::SigCoeffGroup[ chType ];
    for( unsigned ctxId = 0; ctxId < sm_maxNumSigSbbCtx; ctxId++ )
    {
      m_sigSbbFracBits[ ctxId ] = fracBitsAccess.getFracBitsArray( ctxSet( ctxId ) );
    }
  }

  void RateEstimator::xSetSigFlagBits( const FracBitsAccess& fracBitsAccess, ChannelType chType )
  {
    for( unsigned ctxSetId = 0; ctxSetId < sm_numCtxSetsSig; ctxSetId++ )
    {
      BinFracBits*    bits    = m_sigFracBits [ ctxSetId ];
      const CtxSet&   ctxSet  = Ctx::SigFlag  [ chType + 2*ctxSetId ];
      const unsigned  numCtx  = ( chType == CHANNEL_TYPE_LUMA ? 12 : 8 );
      for( unsigned ctxId = 0; ctxId < numCtx; ctxId++ )
      {
        bits[ ctxId ] = fracBitsAccess.getFracBitsArray( ctxSet( ctxId ) );
      }
    }
  }

  void RateEstimator::xSetGtxFlagBits( const FracBitsAccess& fracBitsAccess, ChannelType chType )
  {
    const CtxSet&   ctxSetPar   = Ctx::ParFlag [     chType ];
    const CtxSet&   ctxSetGt1   = Ctx::GtxFlag [ 2 + chType ];
    const CtxSet&   ctxSetGt2   = Ctx::GtxFlag [     chType ];
    const unsigned  numCtx      = ( chType == CHANNEL_TYPE_LUMA ? 21 : 11 );
    for( unsigned ctxId = 0; ctxId < numCtx; ctxId++ )
    {
      BinFracBits     fbPar = fracBitsAccess.getFracBitsArray( ctxSetPar( ctxId ) );
      BinFracBits     fbGt1 = fracBitsAccess.getFracBitsArray( ctxSetGt1( ctxId ) );
      BinFracBits     fbGt2 = fracBitsAccess.getFracBitsArray( ctxSetGt2( ctxId ) );
      CoeffFracBits&  cb    = m_gtxFracBits[ ctxId ];
      int32_t         par0  = (1<<SCALE_BITS) + int32_t(fbPar.intBits[0]);
      int32_t         par1  = (1<<SCALE_BITS) + int32_t(fbPar.intBits[1]);
      cb.bits[0] = 0;
      cb.bits[1] = fbGt1.intBits[0] + (1 << SCALE_BITS);
      cb.bits[2] = fbGt1.intBits[1] + par0 + fbGt2.intBits[0];
      cb.bits[3] = fbGt1.intBits[1] + par1 + fbGt2.intBits[0];
      cb.bits[4] = fbGt1.intBits[1] + par0 + fbGt2.intBits[1];
      cb.bits[5] = fbGt1.intBits[1] + par1 + fbGt2.intBits[1];
    }
  }





  /*================================================================================*/
  /*=====                                                                      =====*/
  /*=====   D A T A   S T R U C T U R E S                                      =====*/
  /*=====                                                                      =====*/
  /*================================================================================*/


  struct PQData
  {
    TCoeff  absLevel;
    int64_t deltaDist;
  };


  struct Decision
  {
    int64_t rdCost;
    TCoeff  absLevel;
    int     prevId;
  };




  /*================================================================================*/
  /*=====                                                                      =====*/
  /*=====   P R E - Q U A N T I Z E R                                          =====*/
  /*=====                                                                      =====*/
  /*================================================================================*/

  class Quantizer
  {
  public:
    Quantizer() {}
#if TCQ_8STATES
    void  dequantBlock         ( const TransformUnit& tu, const ComponentID compID, const QpParam& cQP, CoeffBuf& recCoeff, bool enableScalingLists, int* piDequantCoef, const uint64_t stateTransTab ) const;
#else
    void  dequantBlock         ( const TransformUnit& tu, const ComponentID compID, const QpParam& cQP, CoeffBuf& recCoeff, bool enableScalingLists, int* piDequantCoef ) const;
#endif

		void  initQuantBlock       ( const TransformUnit& tu, const ComponentID compID, const QpParam& cQP, const double lambda, int gValue );
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT_VS
    inline void   preQuantCoeff( const TCoeff absCoeff, PQData *pqData, TCoeff quanCoeff ) const;
#else
    inline void   preQuantCoeff( const TCoeff absCoeff, PQData *pqData, int quanCoeff ) const;
#endif
    inline TCoeff getLastThreshold() const { return m_thresLast; }
    inline TCoeff getSSbbThreshold() const { return m_thresSSbb; }

    inline int64_t getQScale()       const { return m_QScale; }
  private:
    // quantization
    int               m_QShift;
    int64_t           m_QAdd;
    int64_t           m_QScale;
    TCoeff            m_maxQIdx;
    TCoeff            m_thresLast;
    TCoeff            m_thresSSbb;
    // distortion normalization
    int               m_DistShift;
    int64_t           m_DistAdd;
    int64_t           m_DistStepAdd;
    int64_t           m_DistOrgFact;
  };

  inline int ceil_log2(uint64_t x)
  {
    static const uint64_t t[6] = { 0xFFFFFFFF00000000ull, 0x00000000FFFF0000ull, 0x000000000000FF00ull, 0x00000000000000F0ull, 0x000000000000000Cull, 0x0000000000000002ull };
    int y = (((x & (x - 1)) == 0) ? 0 : 1);
    int j = 32;
    for( int i = 0; i < 6; i++)
    {
      int k = (((x & t[i]) == 0) ? 0 : j);
      y += k;
      x >>= k;
      j >>= 1;
    }
    return y;
  }
  void Quantizer::initQuantBlock(const TransformUnit& tu, const ComponentID compID, const QpParam& cQP, const double lambda, int gValue = -1)
  {
    CHECKD( lambda <= 0.0, "Lambda must be greater than 0" );

    const int         qpDQ                  = cQP.Qp(tu.mtsIdx[compID] == MTS_SKIP) + 1;
    const int         qpPer                 = qpDQ / 6;
    const int         qpRem                 = qpDQ - 6 * qpPer;
    const SPS&        sps                   = *tu.cs->sps;
    const CompArea&   area                  = tu.blocks[ compID ];
    const ChannelType chType                = toChannelType( compID );
    const int         channelBitDepth       = sps.getBitDepth( chType );
    const int         maxLog2TrDynamicRange = sps.getMaxLog2TrDynamicRange( chType );
    const int         nomTransformShift     = getTransformShift( channelBitDepth, area.size(), maxLog2TrDynamicRange );
    const bool        clipTransformShift    = ( tu.mtsIdx[compID] == MTS_SKIP && sps.getSpsRangeExtension().getExtendedPrecisionProcessingFlag());
    const bool    needsSqrt2ScaleAdjustment = TU::needsSqrt2Scale(tu, compID);
    const int         transformShift        = ( clipTransformShift ? std::max<int>( 0, nomTransformShift ) : nomTransformShift ) + (needsSqrt2ScaleAdjustment?-1:0);
    // quant parameters
    m_QShift                    = QUANT_SHIFT  - 1 + qpPer + transformShift;
    m_QAdd                      = -( ( 3 << m_QShift ) >> 1 );
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT_VS
    int               invShift  = IQUANT_SHIFT + 1 - qpPer - transformShift;
#else
    Intermediate_Int  invShift  = IQUANT_SHIFT + 1 - qpPer - transformShift;
#endif
    m_QScale                    = g_quantScales[needsSqrt2ScaleAdjustment?1:0][ qpRem ];
    const unsigned    qIdxBD    = std::min<unsigned>( maxLog2TrDynamicRange + 1, 8*sizeof(Intermediate_Int) + invShift - IQUANT_SHIFT - 1 );
    m_maxQIdx                   = ( 1 << (qIdxBD-1) ) - 4;
    m_thresLast                 = TCoeff((int64_t(4) << m_QShift));
    m_thresSSbb                 = TCoeff((int64_t(3) << m_QShift));
    // distortion calculation parameters
    const int64_t qScale        = (gValue==-1) ? m_QScale : gValue;
    const int nomDShift =
      SCALE_BITS - 2 * (nomTransformShift + DISTORTION_PRECISION_ADJUSTMENT(channelBitDepth)) + m_QShift + (needsSqrt2ScaleAdjustment ? 1 : 0);
    const double  qScale2       = double( qScale * qScale );
    const double  nomDistFactor = ( nomDShift < 0 ? 1.0/(double(int64_t(1)<<(-nomDShift))*qScale2*lambda) : double(int64_t(1)<<nomDShift)/(qScale2*lambda) );
    const int64_t pow2dfShift   = (int64_t)( nomDistFactor * qScale2 ) + 1;
    const int     dfShift       = ceil_log2( pow2dfShift );
    m_DistShift                 = 62 + m_QShift - 2*maxLog2TrDynamicRange - dfShift;
    m_DistAdd                   = (int64_t(1) << m_DistShift) >> 1;
    m_DistStepAdd               = (int64_t)( nomDistFactor * double(int64_t(1)<<(m_DistShift+m_QShift)) + .5 );
    m_DistOrgFact               = (int64_t)( nomDistFactor * double(int64_t(1)<<(m_DistShift+1       )) + .5 );
  }

#if TCQ_8STATES
  void Quantizer::dequantBlock( const TransformUnit& tu, const ComponentID compID, const QpParam& cQP, CoeffBuf& recCoeff, bool enableScalingLists, int* piDequantCoef, const uint64_t stateTransTab ) const
#else
  void Quantizer::dequantBlock( const TransformUnit& tu, const ComponentID compID, const QpParam& cQP, CoeffBuf& recCoeff, bool enableScalingLists, int* piDequantCoef) const
#endif
	{

    //----- set basic parameters -----
    const CompArea&     area      = tu.blocks[ compID ];
    const int           numCoeff  = area.area();
    const SizeType      hsId      = gp_sizeIdxInfo->idxFrom( area.width  );
    const SizeType      vsId      = gp_sizeIdxInfo->idxFrom( area.height );
    const CoeffScanType scanType  = SCAN_DIAG;
    const ScanElement *scan       = g_scanOrder[SCAN_GROUPED_4x4][scanType][hsId][vsId];
    const TCoeff*       qCoeff    = tu.getCoeffs( compID ).buf;
          TCoeff*       tCoeff    = recCoeff.buf;

    //----- reset coefficients and get last scan index -----
    ::memset( tCoeff, 0, numCoeff * sizeof(TCoeff) );
    int lastScanIdx = -1;
    for( int scanIdx = numCoeff - 1; scanIdx >= 0; scanIdx-- )
    {
      if (qCoeff[scan[scanIdx].idx])
      {
        lastScanIdx = scanIdx;
        break;
      }
    }
    if( lastScanIdx < 0 )
    {
      return;
    }

    //----- set dequant parameters -----
    const int         qpDQ                  = cQP.Qp(tu.mtsIdx[compID] == MTS_SKIP) + 1;
    const int         qpPer                 = qpDQ / 6;
    const int         qpRem                 = qpDQ - 6 * qpPer;
    const SPS&        sps                   = *tu.cs->sps;
    const ChannelType chType                = toChannelType( compID );
    const int         channelBitDepth       = sps.getBitDepth( chType );
    const int         maxLog2TrDynamicRange = sps.getMaxLog2TrDynamicRange( chType );
    const TCoeff      minTCoeff             = -( 1 << maxLog2TrDynamicRange );
    const TCoeff      maxTCoeff             =  ( 1 << maxLog2TrDynamicRange ) - 1;
    const int         nomTransformShift     = getTransformShift( channelBitDepth, area.size(), maxLog2TrDynamicRange );
    const bool        clipTransformShift    = ( tu.mtsIdx[compID] == MTS_SKIP && sps.getSpsRangeExtension().getExtendedPrecisionProcessingFlag());
    const bool    needsSqrt2ScaleAdjustment = TU::needsSqrt2Scale(tu, compID);
    const int         transformShift        = ( clipTransformShift ? std::max<int>( 0, nomTransformShift ) : nomTransformShift ) + (needsSqrt2ScaleAdjustment?-1:0);
    Intermediate_Int  shift                 = IQUANT_SHIFT + 1 - qpPer - transformShift + (enableScalingLists ? LOG2_SCALING_LIST_NEUTRAL_VALUE : 0);
    Intermediate_Int  invQScale             = g_invQuantScales[needsSqrt2ScaleAdjustment?1:0][ qpRem ];
    Intermediate_Int  add = (shift < 0) ? 0 : ((1 << shift) >> 1);
    //----- dequant coefficients -----
    for( int state = 0, scanIdx = lastScanIdx; scanIdx >= 0; scanIdx-- )
    {
      const unsigned  rasterPos = scan[scanIdx].idx;
      const TCoeff&   level     = qCoeff[ rasterPos ];
      if( level )
      {
        if (enableScalingLists)
          invQScale = piDequantCoef[rasterPos];//scalingfactor*levelScale
        if (shift < 0 && (enableScalingLists || scanIdx == lastScanIdx))
        {
          invQScale <<= -shift;
        }
#if TCQ_8STATES
        Intermediate_Int  qIdx      = ( level << 1 ) + ( level > 0 ? -(state&1) : (state&1) );
#else
        Intermediate_Int  qIdx      = ( level << 1 ) + ( level > 0 ? -(state>>1) : (state>>1) );
#endif
				int64_t  nomTCoeff          = ((int64_t)qIdx * (int64_t)invQScale + add) >> ((shift < 0) ? 0 : shift);
        tCoeff[rasterPos]           = (TCoeff)Clip3<int64_t>(minTCoeff, maxTCoeff, nomTCoeff);
      }
#if TCQ_8STATES
      state = int( ( stateTransTab >> ((state<<3)+((level&1)<<2)) ) & 15 );
#else
      state = ( 32040 >> ((state<<2)+((level&1)<<1)) ) & 3;   // the 16-bit value "32040" represent the state transition table
#endif
		}
  }

#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT_VS
  inline void Quantizer::preQuantCoeff(const TCoeff absCoeff, PQData *pqData, TCoeff quanCoeff) const
#else
  inline void Quantizer::preQuantCoeff(const TCoeff absCoeff, PQData *pqData, int quanCoeff) const
#endif
  {
    int64_t scaledOrg = int64_t( absCoeff ) * quanCoeff;
    TCoeff  qIdx      = std::max<TCoeff>( 1, std::min<TCoeff>( m_maxQIdx, TCoeff( ( scaledOrg + m_QAdd ) >> m_QShift ) ) );
    int64_t scaledAdd = qIdx * m_DistStepAdd - scaledOrg * m_DistOrgFact;
    PQData& pq_a      = pqData[ qIdx & 3 ];
    pq_a.deltaDist    = ( scaledAdd * qIdx + m_DistAdd ) >> m_DistShift;
    pq_a.absLevel     = ( ++qIdx ) >> 1;
    scaledAdd        += m_DistStepAdd;
    PQData& pq_b      = pqData[ qIdx & 3 ];
    pq_b.deltaDist    = ( scaledAdd * qIdx + m_DistAdd ) >> m_DistShift;
    pq_b.absLevel     = ( ++qIdx ) >> 1;
    scaledAdd        += m_DistStepAdd;
    PQData& pq_c      = pqData[ qIdx & 3 ];
    pq_c.deltaDist    = ( scaledAdd * qIdx + m_DistAdd ) >> m_DistShift;
    pq_c.absLevel     = ( ++qIdx ) >> 1;
    scaledAdd        += m_DistStepAdd;
    PQData& pq_d      = pqData[ qIdx & 3 ];
    pq_d.deltaDist    = ( scaledAdd * qIdx + m_DistAdd ) >> m_DistShift;
    pq_d.absLevel     = ( ++qIdx ) >> 1;
  }







  /*================================================================================*/
  /*=====                                                                      =====*/
  /*=====   T C Q   S T A T E                                                  =====*/
  /*=====                                                                      =====*/
  /*================================================================================*/

  class State;

  struct SbbCtx
  {
    uint8_t*  sbbFlags;
    uint8_t*  levels;
  };

  class CommonCtx
  {
  public:
#if TCQ_8STATES
    CommonCtx() : m_currSbbCtx( m_allSbbCtx ), m_prevSbbCtx( m_currSbbCtx + 8 ) {}
#else
    CommonCtx() : m_currSbbCtx( m_allSbbCtx ), m_prevSbbCtx( m_currSbbCtx + 4 ) {}
#endif

    inline void swap() { std::swap(m_currSbbCtx, m_prevSbbCtx); }

    inline void reset( const TUParameters& tuPars, const RateEstimator &rateEst)
    {
      m_nbInfo = tuPars.m_scanId2NbInfoOut;
      ::memcpy( m_sbbFlagBits, rateEst.sigSbbFracBits(), 2*sizeof(BinFracBits) );
      const int numSbb    = tuPars.m_numSbb;
      const int chunkSize = numSbb + tuPars.m_numCoeff;
      uint8_t*  nextMem   = m_memory;
#if TCQ_8STATES
      for( int k = 0; k < 16; k++, nextMem += chunkSize )
#else
      for( int k = 0; k < 8; k++, nextMem += chunkSize )
#endif
      {
        m_allSbbCtx[k].sbbFlags = nextMem;
        m_allSbbCtx[k].levels   = nextMem + numSbb;
      }
    }

    inline void update(const ScanInfo &scanInfo, const State *prevState, State &currState);

  private:
    const NbInfoOut*            m_nbInfo;
    BinFracBits                 m_sbbFlagBits[2];
#if TCQ_8STATES
    SbbCtx                      m_allSbbCtx  [16];
#else
    SbbCtx                      m_allSbbCtx  [8];
#endif
    SbbCtx*                     m_currSbbCtx;
    SbbCtx*                     m_prevSbbCtx;
#if TCQ_8STATES
    uint8_t                     m_memory[ 16 * ( MAX_TB_SIZEY * MAX_TB_SIZEY + MLS_GRP_NUM ) ];
#else
    uint8_t                     m_memory[ 8 * ( MAX_TB_SIZEY * MAX_TB_SIZEY + MLS_GRP_NUM ) ];
#endif
  };

#define RICEMAX 32
  const int32_t g_goRiceBits[4][RICEMAX] =
  {
    { 32768,  65536,  98304, 131072, 163840, 196608, 262144, 262144, 327680, 327680, 327680, 327680, 393216, 393216, 393216, 393216, 393216, 393216, 393216, 393216, 458752, 458752, 458752, 458752, 458752, 458752, 458752, 458752, 458752, 458752, 458752, 458752},
    { 65536,  65536,  98304,  98304, 131072, 131072, 163840, 163840, 196608, 196608, 229376, 229376, 294912, 294912, 294912, 294912, 360448, 360448, 360448, 360448, 360448, 360448, 360448, 360448, 425984, 425984, 425984, 425984, 425984, 425984, 425984, 425984},
    { 98304,  98304,  98304,  98304, 131072, 131072, 131072, 131072, 163840, 163840, 163840, 163840, 196608, 196608, 196608, 196608, 229376, 229376, 229376, 229376, 262144, 262144, 262144, 262144, 327680, 327680, 327680, 327680, 327680, 327680, 327680, 327680},
    {131072, 131072, 131072, 131072, 131072, 131072, 131072, 131072, 163840, 163840, 163840, 163840, 163840, 163840, 163840, 163840, 196608, 196608, 196608, 196608, 196608, 196608, 196608, 196608, 229376, 229376, 229376, 229376, 229376, 229376, 229376, 229376}
  };

  class State
  {
    friend class CommonCtx;
  public:
    State( const RateEstimator& rateEst, CommonCtx& commonCtx, const int stateId );

    template<uint8_t numIPos>
    inline void updateState(const ScanInfo &scanInfo, const State *prevStates, const Decision &decision);
#if TCQ_8STATES
    inline void updateStateEOS(const ScanInfo &scanInfo, const State *prevStates, const State *skipStates,
                               const Decision &decision, int skipOffset );
#else
    inline void updateStateEOS(const ScanInfo &scanInfo, const State *prevStates, const State *skipStates,
                               const Decision &decision);
#endif

#if TCQ_8STATES
    inline void init( int effectiveWidth, int effectiveHeight )
#else
    inline void init()
#endif
    {
      ::memset(m_absLevelsAndCtxInit, 0, sizeof(m_absLevelsAndCtxInit));
      m_rdCost        = std::numeric_limits<int64_t>::max()>>1;
      m_numSigSbb     = 0;
      m_remRegBins    = 4;  // just large enough for last scan pos
      m_refSbbCtxId   = -1;
      m_sigFracBits   = m_sigFracBitsArray[ 0 ];
      m_coeffFracBits = m_gtxFracBitsArray[ 0 ];
      m_goRicePar     = 0;
      m_goRiceZero    = 0;
#if TCQ_8STATES
      effWidth        = effectiveWidth;
      effHeight       = effectiveHeight;
#endif
    }

    void checkRdCosts( const ScanPosType spt, const PQData& pqDataA, const PQData& pqDataB, Decision& decisionA, Decision& decisionB ) const
    {
      const int32_t*  goRiceTab = g_goRiceBits[m_goRicePar];
      int64_t         rdCostA   = m_rdCost + pqDataA.deltaDist;
      int64_t         rdCostB   = m_rdCost + pqDataB.deltaDist;
      int64_t         rdCostZ   = m_rdCost;

      const TCoeff absLevelA = pqDataA.absLevel;
      const TCoeff absLevelB = pqDataB.absLevel;
      const int32_t* bits    = m_coeffFracBits.bits;

      if( m_remRegBins >= 4 )
      {
        if( absLevelA < 4 )
        {
          rdCostA += bits[absLevelA];
        }
        else
        {
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT_VS
          const TCoeff value = ( absLevelA - 4 ) >> 1;
#else
          const int value = ( absLevelA - 4 ) >> 1;
#endif
          rdCostA += bits[absLevelA - ( value << 1 )] + goRiceTab[value < RICEMAX ? value : RICEMAX - 1];
        }

        if( absLevelB < 4 )
        {
          rdCostB += bits[absLevelB];
        }
        else
        {
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT_VS
          const TCoeff value = ( absLevelB - 4 ) >> 1;
#else
          const int value = ( absLevelB - 4 ) >> 1;
#endif
          rdCostB += bits[absLevelB - ( value << 1 )] + goRiceTab[value < RICEMAX ? value : RICEMAX - 1];
        }

        const int sigBit1 = m_sigFracBits.intBits[1];

        if( spt == SCAN_ISCSBB )
        {
          rdCostA += sigBit1;
          rdCostB += sigBit1;
          rdCostZ += m_sigFracBits.intBits[0];
        }
        else if( spt == SCAN_SOCSBB )
        {
          const int sbbBit1 = m_sbbFracBits.intBits[1];

          rdCostA += sbbBit1 + sigBit1;
          rdCostB += sbbBit1 + sigBit1;
          rdCostZ += sbbBit1 + m_sigFracBits.intBits[0];
        }
        else if( m_numSigSbb )
        {
          rdCostA += sigBit1;
          rdCostB += sigBit1;
          rdCostZ += m_sigFracBits.intBits[0];
        }
        else
        {
          rdCostZ = decisionA.rdCost;
        }
      }
      else
      {
        rdCostA += ( 1 << SCALE_BITS ) + goRiceTab[absLevelA <= m_goRiceZero ? absLevelA - 1 : ( absLevelA < RICEMAX ? absLevelA : RICEMAX - 1 )];
        rdCostB += ( 1 << SCALE_BITS ) + goRiceTab[absLevelB <= m_goRiceZero ? absLevelB - 1 : ( absLevelB < RICEMAX ? absLevelB : RICEMAX - 1 )];
        rdCostZ += goRiceTab[m_goRiceZero];
      }

      if( rdCostA < decisionA.rdCost )
      {
        decisionA.rdCost = rdCostA;
        decisionA.absLevel = absLevelA;
        decisionA.prevId = m_stateId;
      }
      if( rdCostZ < decisionA.rdCost )
      {
        decisionA.rdCost = rdCostZ;
        decisionA.absLevel = 0;
        decisionA.prevId = m_stateId;
      }
      if( rdCostB < decisionB.rdCost )
      {
        decisionB.rdCost = rdCostB;
        decisionB.absLevel = absLevelB;
        decisionB.prevId = m_stateId;
      }
    }

    inline void checkRdCostStart(int32_t lastOffset, const PQData &pqData, Decision &decision) const
    {
      int64_t rdCost = pqData.deltaDist + lastOffset;
      if (pqData.absLevel < 4)
      {
        rdCost += m_coeffFracBits.bits[pqData.absLevel];
      }
      else
      {
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT_VS
        const TCoeff value = (pqData.absLevel - 4) >> 1;
#else
        const unsigned value = (pqData.absLevel - 4) >> 1;
#endif
        rdCost += m_coeffFracBits.bits[pqData.absLevel - (value << 1)] + g_goRiceBits[m_goRicePar][value < RICEMAX ? value : RICEMAX-1];
      }
      if( rdCost < decision.rdCost )
      {
        decision.rdCost   = rdCost;
        decision.absLevel = pqData.absLevel;
        decision.prevId   = -1;
      }
    }

#if TCQ_8STATES
		inline void checkRdCostSkipSbb(Decision &decision, int skipOffset) const
#else
    inline void checkRdCostSkipSbb(Decision &decision) const
#endif
    {
      int64_t rdCost = m_rdCost + m_sbbFracBits.intBits[0];
      if( rdCost < decision.rdCost )
      {
        decision.rdCost   = rdCost;
        decision.absLevel = 0;
#if TCQ_8STATES
        decision.prevId   = skipOffset + m_stateId;
#else
        decision.prevId   = 4+m_stateId;
#endif
      }
    }
#if !TU_256 || EXTENDED_LFNST
#if TCQ_8STATES
    inline void checkRdCostSkipSbbZeroOut(Decision &decision, int skipOffset) const
#else
    inline void checkRdCostSkipSbbZeroOut(Decision &decision) const
#endif
    {
      int64_t rdCost = m_rdCost + m_sbbFracBits.intBits[0];
      decision.rdCost = rdCost;
      decision.absLevel = 0;
#if TCQ_8STATES
      decision.prevId = skipOffset + m_stateId;
#else
      decision.prevId = 4 + m_stateId;
#endif
    }
#endif

  private:
    int64_t                   m_rdCost;
    uint16_t                  m_absLevelsAndCtxInit[24];  // 16x8bit for abs levels + 16x16bit for ctx init id
    int8_t                    m_numSigSbb;
    int                       m_remRegBins;
    int8_t                    m_refSbbCtxId;
    BinFracBits               m_sbbFracBits;
    BinFracBits               m_sigFracBits;
    CoeffFracBits             m_coeffFracBits;
    int8_t                    m_goRicePar;
    int8_t                    m_goRiceZero;
    const int8_t              m_stateId;
    const BinFracBits*const   m_sigFracBitsArray;
    const CoeffFracBits*const m_gtxFracBitsArray;
    CommonCtx&                m_commonCtx;
  public:
    unsigned                  effWidth;
    unsigned                  effHeight;
  };


  State::State( const RateEstimator& rateEst, CommonCtx& commonCtx, const int stateId )
    : m_sbbFracBits     { { 0, 0 } }
    , m_stateId         ( stateId )
    , m_sigFracBitsArray( rateEst.sigFlagBits(stateId) )
    , m_gtxFracBitsArray( rateEst.gtxFracBits(stateId) )
    , m_commonCtx       ( commonCtx )
  {
  }

  template<uint8_t numIPos>
  inline void State::updateState(const ScanInfo &scanInfo, const State *prevStates, const Decision &decision)
  {
    m_rdCost = decision.rdCost;
    if( decision.prevId > -2 )
    {
      if( decision.prevId >= 0 )
      {
        const State*  prvState  = prevStates            +   decision.prevId;
        m_numSigSbb             = prvState->m_numSigSbb + !!decision.absLevel;
        m_refSbbCtxId           = prvState->m_refSbbCtxId;
        m_sbbFracBits           = prvState->m_sbbFracBits;
        m_remRegBins            = prvState->m_remRegBins - 1;
        m_goRicePar             = prvState->m_goRicePar;
        if( m_remRegBins >= 4 )
        {
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT_VS
          m_remRegBins -= (decision.absLevel < 2 ? (unsigned)decision.absLevel : 3);
#else
          m_remRegBins -= (decision.absLevel < 2 ? decision.absLevel : 3);
#endif
        }
        ::memcpy( m_absLevelsAndCtxInit, prvState->m_absLevelsAndCtxInit, 48*sizeof(uint8_t) );
      }
      else
      {
        m_numSigSbb     =  1;
        m_refSbbCtxId   = -1;
        int ctxBinSampleRatio = (scanInfo.chType == CHANNEL_TYPE_LUMA) ? MAX_TU_LEVEL_CTX_CODED_BIN_CONSTRAINT_LUMA : MAX_TU_LEVEL_CTX_CODED_BIN_CONSTRAINT_CHROMA;
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT_VS
        m_remRegBins = (effWidth * effHeight *ctxBinSampleRatio) / 16 - (decision.absLevel < 2 ? (unsigned)decision.absLevel : 3);
#else
        m_remRegBins = (effWidth * effHeight *ctxBinSampleRatio) / 16 - (decision.absLevel < 2 ? decision.absLevel : 3);
#endif
        ::memset( m_absLevelsAndCtxInit, 0, 48*sizeof(uint8_t) );
      }

      uint8_t* levels               = reinterpret_cast<uint8_t*>(m_absLevelsAndCtxInit);
      levels[ scanInfo.insidePos ]  = (uint8_t)std::min<TCoeff>( 255, decision.absLevel );

      if (m_remRegBins >= 4)
      {
        TCoeff  tinit = m_absLevelsAndCtxInit[8 + scanInfo.nextInsidePos];
        TCoeff  sumAbs1 = (tinit >> 3) & 31;
        TCoeff  sumNum = tinit & 7;
#define UPDATE(k) {TCoeff t=levels[scanInfo.nextNbInfoSbb.inPos[k]]; sumAbs1+=std::min<TCoeff>(4+(t&1),t); sumNum+=!!t; }
        if (numIPos == 1)
        {
          UPDATE(0);
        }
        else if (numIPos == 2)
        {
          UPDATE(0);
          UPDATE(1);
        }
        else if (numIPos == 3)
        {
          UPDATE(0);
          UPDATE(1);
          UPDATE(2);
        }
        else if (numIPos == 4)
        {
          UPDATE(0);
          UPDATE(1);
          UPDATE(2);
          UPDATE(3);
        }
        else if (numIPos == 5)
        {
          UPDATE(0);
          UPDATE(1);
          UPDATE(2);
          UPDATE(3);
          UPDATE(4);
        }
#undef UPDATE
        TCoeff sumGt1 = sumAbs1 - sumNum;
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
        m_sigFracBits = m_sigFracBitsArray[scanInfo.sigCtxOffsetNext + std::min<TCoeff>( (sumAbs1+1)>>1, 3 )];
#else
        m_sigFracBits = m_sigFracBitsArray[scanInfo.sigCtxOffsetNext + std::min( (sumAbs1+1)>>1, 3 )];
#endif
        m_coeffFracBits = m_gtxFracBitsArray[scanInfo.gtxCtxOffsetNext + (sumGt1 < 4 ? sumGt1 : 4)];

        TCoeff  sumAbs = m_absLevelsAndCtxInit[8 + scanInfo.nextInsidePos] >> 8;
#define UPDATE(k) {TCoeff t=levels[scanInfo.nextNbInfoSbb.inPos[k]]; sumAbs+=t; }
        if (numIPos == 1)
        {
          UPDATE(0);
        }
        else if (numIPos == 2)
        {
          UPDATE(0);
          UPDATE(1);
        }
        else if (numIPos == 3)
        {
          UPDATE(0);
          UPDATE(1);
          UPDATE(2);
        }
        else if (numIPos == 4)
        {
          UPDATE(0);
          UPDATE(1);
          UPDATE(2);
          UPDATE(3);
        }
        else if (numIPos == 5)
        {
          UPDATE(0);
          UPDATE(1);
          UPDATE(2);
          UPDATE(3);
          UPDATE(4);
        }
#undef UPDATE
        int sumAll = std::max(std::min(31, (int)sumAbs - 4 * 5), 0);
        m_goRicePar = g_auiGoRiceParsCoeff[sumAll];
      }
      else
      {
        TCoeff  sumAbs = m_absLevelsAndCtxInit[8 + scanInfo.nextInsidePos] >> 8;
#define UPDATE(k) {TCoeff t=levels[scanInfo.nextNbInfoSbb.inPos[k]]; sumAbs+=t; }
        if (numIPos == 1)
        {
          UPDATE(0);
        }
        else if (numIPos == 2)
        {
          UPDATE(0);
          UPDATE(1);
        }
        else if (numIPos == 3)
        {
          UPDATE(0);
          UPDATE(1);
          UPDATE(2);
        }
        else if (numIPos == 4)
        {
          UPDATE(0);
          UPDATE(1);
          UPDATE(2);
          UPDATE(3);
        }
        else if (numIPos == 5)
        {
          UPDATE(0);
          UPDATE(1);
          UPDATE(2);
          UPDATE(3);
          UPDATE(4);
        }
#undef UPDATE
        sumAbs = std::min<TCoeff>(31, sumAbs);
        m_goRicePar = g_auiGoRiceParsCoeff[sumAbs];
        m_goRiceZero = g_auiGoRicePosCoeff0(m_stateId, m_goRicePar);
      }
    }
  }

#if TCQ_8STATES
  inline void State::updateStateEOS(const ScanInfo &scanInfo, const State *prevStates, const State *skipStates,
                                    const Decision &decision, int skipOffset )
#else
  inline void State::updateStateEOS(const ScanInfo &scanInfo, const State *prevStates, const State *skipStates,
                                    const Decision &decision)
#endif
  {
    m_rdCost = decision.rdCost;
    if( decision.prevId > -2 )
    {
      const State* prvState = 0;
#if TCQ_8STATES
      if( decision.prevId  >= skipOffset )
#else
      if( decision.prevId  >= 4 )
#endif
      {
        CHECK( decision.absLevel != 0, "cannot happen" );
#if TCQ_8STATES
        prvState    = skipStates + ( decision.prevId - skipOffset );
#else
        prvState    = skipStates + ( decision.prevId - 4 );
#endif
        m_numSigSbb = 0;
        ::memset( m_absLevelsAndCtxInit, 0, 16*sizeof(uint8_t) );
      }
      else if( decision.prevId  >= 0 )
      {
        prvState    = prevStates            +   decision.prevId;
        m_numSigSbb = prvState->m_numSigSbb + !!decision.absLevel;
        ::memcpy( m_absLevelsAndCtxInit, prvState->m_absLevelsAndCtxInit, 16*sizeof(uint8_t) );
      }
      else
      {
        m_numSigSbb = 1;
        ::memset( m_absLevelsAndCtxInit, 0, 16*sizeof(uint8_t) );
      }
      reinterpret_cast<uint8_t*>(m_absLevelsAndCtxInit)[ scanInfo.insidePos ] = (uint8_t)std::min<TCoeff>( 255, decision.absLevel );

      m_commonCtx.update( scanInfo, prvState, *this );

      TCoeff  tinit   = m_absLevelsAndCtxInit[ 8 + scanInfo.nextInsidePos ];
      TCoeff  sumNum  =   tinit        & 7;
      TCoeff  sumAbs1 = ( tinit >> 3 ) & 31;
      TCoeff  sumGt1  = sumAbs1        - sumNum;
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT
      m_sigFracBits   = m_sigFracBitsArray[ scanInfo.sigCtxOffsetNext + std::min<TCoeff>( (sumAbs1+1)>>1, 3 ) ];
#else
      m_sigFracBits   = m_sigFracBitsArray[ scanInfo.sigCtxOffsetNext + std::min( (sumAbs1+1)>>1, 3 ) ];
#endif
      m_coeffFracBits = m_gtxFracBitsArray[ scanInfo.gtxCtxOffsetNext + ( sumGt1  < 4 ? sumGt1  : 4 ) ];
    }
  }

  inline void CommonCtx::update(const ScanInfo &scanInfo, const State *prevState, State &currState)
  {
    uint8_t*    sbbFlags  = m_currSbbCtx[ currState.m_stateId ].sbbFlags;
    uint8_t*    levels    = m_currSbbCtx[ currState.m_stateId ].levels;
    std::size_t setCpSize = m_nbInfo[ scanInfo.scanIdx - 1 ].maxDist * sizeof(uint8_t);
    if( prevState && prevState->m_refSbbCtxId >= 0 )
    {
      ::memcpy( sbbFlags,                  m_prevSbbCtx[prevState->m_refSbbCtxId].sbbFlags,                  scanInfo.numSbb*sizeof(uint8_t) );
      ::memcpy( levels + scanInfo.scanIdx, m_prevSbbCtx[prevState->m_refSbbCtxId].levels + scanInfo.scanIdx, setCpSize );
    }
    else
    {
      ::memset( sbbFlags,                  0, scanInfo.numSbb*sizeof(uint8_t) );
      ::memset( levels + scanInfo.scanIdx, 0, setCpSize );
    }
    sbbFlags[ scanInfo.sbbPos ] = !!currState.m_numSigSbb;
    ::memcpy( levels + scanInfo.scanIdx, currState.m_absLevelsAndCtxInit, scanInfo.sbbSize*sizeof(uint8_t) );

    const int       sigNSbb   = ( ( scanInfo.nextSbbRight ? sbbFlags[ scanInfo.nextSbbRight ] : false ) || ( scanInfo.nextSbbBelow ? sbbFlags[ scanInfo.nextSbbBelow ] : false ) ? 1 : 0 );
    currState.m_numSigSbb     = 0;
    if (prevState)
    {
      currState.m_remRegBins = prevState->m_remRegBins;
    }
    else
    {
      int ctxBinSampleRatio = (scanInfo.chType == CHANNEL_TYPE_LUMA) ? MAX_TU_LEVEL_CTX_CODED_BIN_CONSTRAINT_LUMA : MAX_TU_LEVEL_CTX_CODED_BIN_CONSTRAINT_CHROMA;
      currState.m_remRegBins = (currState.effWidth * currState.effHeight *ctxBinSampleRatio) / 16;
    }
    currState.m_goRicePar     = 0;
    currState.m_refSbbCtxId   = currState.m_stateId;
    currState.m_sbbFracBits   = m_sbbFlagBits[ sigNSbb ];

    uint16_t          templateCtxInit[16];
    const int         scanBeg   = scanInfo.scanIdx - scanInfo.sbbSize;
    const NbInfoOut*  nbOut     = m_nbInfo + scanBeg;
    const uint8_t*    absLevels = levels   + scanBeg;
    for( int id = 0; id < scanInfo.sbbSize; id++, nbOut++ )
    {
      if( nbOut->num )
      {
        TCoeff sumAbs = 0, sumAbs1 = 0, sumNum = 0;
#define UPDATE(k) {TCoeff t=absLevels[nbOut->outPos[k]]; sumAbs+=t; sumAbs1+=std::min<TCoeff>(4+(t&1),t); sumNum+=!!t; }
        UPDATE(0);
        if( nbOut->num > 1 )
        {
          UPDATE(1);
          if( nbOut->num > 2 )
          {
            UPDATE(2);
            if( nbOut->num > 3 )
            {
              UPDATE(3);
              if( nbOut->num > 4 )
              {
                UPDATE(4);
              }
            }
          }
        }
#undef UPDATE
        templateCtxInit[id] = uint16_t(sumNum) + ( uint16_t(sumAbs1) << 3 ) + ( (uint16_t)std::min<TCoeff>( 127, sumAbs ) << 8 );
      }
      else
      {
        templateCtxInit[id] = 0;
      }
    }
    ::memset( currState.m_absLevelsAndCtxInit,     0,               16*sizeof(uint8_t) );
    ::memcpy( currState.m_absLevelsAndCtxInit + 8, templateCtxInit, scanInfo.sbbSize * sizeof( uint16_t ) );
    ::memset( currState.m_absLevelsAndCtxInit + 8 + scanInfo.sbbSize, 0, ( 16 - scanInfo.sbbSize )*sizeof(uint16_t) );
  }



  /*================================================================================*/
  /*=====                                                                      =====*/
  /*=====   T C Q                                                              =====*/
  /*=====                                                                      =====*/
  /*================================================================================*/
  class DepQuant : private RateEstimator
  {
  public:
    DepQuant();

    void    quant   ( TransformUnit& tu, const CCoeffBuf& srcCoeff, const ComponentID compID, const QpParam& cQP, const double lambda, const Ctx& ctx, TCoeff& absSum, bool enableScalingLists, int* quantCoeff );
    void    dequant ( const TransformUnit& tu, CoeffBuf& recCoeff, const ComponentID compID, const QpParam& cQP, bool enableScalingLists, int* quantCoeff );

  private:
#if TU_256 && !EXTENDED_LFNST
#if TCQ_8STATES
		template<int numStates>
#endif
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT_VS
    void    xDecideAndUpdate  ( const TCoeff absCoeff, const ScanInfo& scanInfo, TCoeff quantCoeff);
#else
    void    xDecideAndUpdate  ( const TCoeff absCoeff, const ScanInfo& scanInfo, int quantCoeff);
#endif    
#if TCQ_8STATES
		template<int numStates>
#endif
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT_VS
    void    xDecide           ( const ScanPosType spt, const TCoeff absCoeff, const int lastOffset, Decision* decisions, TCoeff quantCoeff );
#else
    void    xDecide           ( const ScanPosType spt, const TCoeff absCoeff, const int lastOffset, Decision* decisions, int quantCoeff );
#endif    
#else
#if TCQ_8STATES
		template<int numStates>
#endif
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT_VS
    void    xDecideAndUpdate  ( const TCoeff absCoeff, const ScanInfo& scanInfo, bool zeroOut, TCoeff quantCoeff);
#else
    void    xDecideAndUpdate  ( const TCoeff absCoeff, const ScanInfo& scanInfo, bool zeroOut, int quantCoeff);
#endif    
#if TCQ_8STATES
		template<int numStates>
#endif
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT_VS
    void    xDecide           ( const ScanPosType spt, const TCoeff absCoeff, const int lastOffset, Decision* decisions, bool zeroOut, TCoeff quantCoeff );
#else
		void    xDecide           ( const ScanPosType spt, const TCoeff absCoeff, const int lastOffset, Decision* decisions, bool zeroOut, int quantCoeff );
#endif
#endif

  private:
    CommonCtx   m_commonCtx;
#if TCQ_8STATES
    State       m_allStates[ 24 ];
#else
    State       m_allStates[ 12 ];
#endif
    State*      m_currStates;
    State*      m_prevStates;
    State*      m_skipStates;
    State       m_startState;
    Quantizer   m_quant;
#if TCQ_8STATES
    Decision    m_trellis[ MAX_TB_SIZEY * MAX_TB_SIZEY * 16 ];
#else
    Decision    m_trellis[ MAX_TB_SIZEY * MAX_TB_SIZEY ][ 8 ];
#endif
  };


#define TINIT(x) {*this,m_commonCtx,x}
  DepQuant::DepQuant()
    : RateEstimator ()
    , m_commonCtx   ()
#if TCQ_8STATES
    , m_allStates   {TINIT(0),TINIT(1),TINIT(2),TINIT(3),TINIT(4),TINIT(5),TINIT(6),TINIT(7),
                     TINIT(0),TINIT(1),TINIT(2),TINIT(3),TINIT(4),TINIT(5),TINIT(6),TINIT(7),
                     TINIT(0),TINIT(1),TINIT(2),TINIT(3),TINIT(4),TINIT(5),TINIT(6),TINIT(7)}
#else
    , m_allStates   {TINIT(0),TINIT(1),TINIT(2),TINIT(3),TINIT(0),TINIT(1),TINIT(2),TINIT(3),TINIT(0),TINIT(1),TINIT(2),TINIT(3)}
#endif
		, m_currStates  (  m_allStates      )
#if TCQ_8STATES
    , m_prevStates  (  m_currStates + 8 )
    , m_skipStates  (  m_prevStates + 8 )
#else
    , m_prevStates  (  m_currStates + 4 )
    , m_skipStates  (  m_prevStates + 4 )
#endif
    , m_startState  TINIT(0)
  {}
#undef TINIT

  void DepQuant::dequant( const TransformUnit& tu,  CoeffBuf& recCoeff, const ComponentID compID, const QpParam& cQP, bool enableScalingLists, int* piDequantCoef )
  {
#if TCQ_8STATES
    m_quant.dequantBlock( tu, compID, cQP, recCoeff, enableScalingLists, piDequantCoef, g_stateTransTab[ tu.cs->slice->getDepQuantEnabledIdc() ] );
#else
    m_quant.dequantBlock( tu, compID, cQP, recCoeff, enableScalingLists, piDequantCoef );
#endif
  }


#define DINIT(l,p) {std::numeric_limits<int64_t>::max()>>2,l,p}
#if TCQ_8STATES
  static const Decision startDec[2][16] = {
    {
      DINIT(-1,-2),DINIT(-1,-2),DINIT(-1,-2),DINIT(-1,-2),DINIT( 0, 4),DINIT( 0, 5),DINIT( 0, 6),DINIT( 0, 7),
      DINIT( 0, 0),DINIT( 0, 0),DINIT( 0, 0),DINIT( 0, 0),DINIT( 0, 0),DINIT( 0, 0),DINIT( 0, 0),DINIT( 0, 0)
    },{
      DINIT(-1,-2),DINIT(-1,-2),DINIT(-1,-2),DINIT(-1,-2),DINIT(-1,-2),DINIT(-1,-2),DINIT(-1,-2),DINIT(-1,-2),
      DINIT( 0, 8),DINIT( 0, 9),DINIT( 0,10),DINIT( 0,11),DINIT( 0,12),DINIT( 0,13),DINIT( 0,14),DINIT( 0,15)
    }};
#else
  static const Decision startDec[8] = {DINIT(-1,-2),DINIT(-1,-2),DINIT(-1,-2),DINIT(-1,-2),DINIT(0,4),DINIT(0,5),DINIT(0,6),DINIT(0,7)};
#endif
#undef  DINIT

#if TCQ_8STATES 
  template<int numStates>
#endif
#if TU_256 && !EXTENDED_LFNST
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT_VS
  void DepQuant::xDecide( const ScanPosType spt, const TCoeff absCoeff, const int lastOffset, Decision* decisions, TCoeff quanCoeff)
#else
  void DepQuant::xDecide( const ScanPosType spt, const TCoeff absCoeff, const int lastOffset, Decision* decisions, int quanCoeff)
#endif  
#else
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT_VS
  void DepQuant::xDecide( const ScanPosType spt, const TCoeff absCoeff, const int lastOffset, Decision* decisions, bool zeroOut, TCoeff quanCoeff)
#else
  void DepQuant::xDecide( const ScanPosType spt, const TCoeff absCoeff, const int lastOffset, Decision* decisions, bool zeroOut, int quanCoeff)
#endif  
#endif
  {
#if TCQ_8STATES
    ::memcpy( decisions, startDec[numStates>>3], (numStates<<1)*sizeof(Decision) );
#else
    ::memcpy( decisions, startDec, 8*sizeof(Decision) );
#endif
#if !TU_256 || EXTENDED_LFNST
    if( zeroOut )
    {
      if( spt==SCAN_EOCSBB )
      {
#if TCQ_8STATES
        m_skipStates[0].checkRdCostSkipSbbZeroOut( decisions[0], numStates );
        m_skipStates[1].checkRdCostSkipSbbZeroOut( decisions[1], numStates );
        m_skipStates[2].checkRdCostSkipSbbZeroOut( decisions[2], numStates );
        m_skipStates[3].checkRdCostSkipSbbZeroOut( decisions[3], numStates );
        m_skipStates[4].checkRdCostSkipSbbZeroOut( decisions[4], numStates );
        m_skipStates[5].checkRdCostSkipSbbZeroOut( decisions[5], numStates );
        m_skipStates[6].checkRdCostSkipSbbZeroOut( decisions[6], numStates );
        m_skipStates[7].checkRdCostSkipSbbZeroOut( decisions[7], numStates );
#else
        m_skipStates[0].checkRdCostSkipSbbZeroOut(decisions[0]);
        m_skipStates[1].checkRdCostSkipSbbZeroOut(decisions[1]);
        m_skipStates[2].checkRdCostSkipSbbZeroOut(decisions[2]);
        m_skipStates[3].checkRdCostSkipSbbZeroOut(decisions[3]);
#endif
      }
      return;
    }
#endif

    PQData  pqData[4];
    m_quant.preQuantCoeff( absCoeff, pqData, quanCoeff );
#if TCQ_8STATES
    if( numStates > 4 )
    {
      m_prevStates[0].checkRdCosts( spt, pqData[0], pqData[2], decisions[0], decisions[2]);
      m_prevStates[1].checkRdCosts( spt, pqData[3], pqData[1], decisions[5], decisions[7]);
      m_prevStates[2].checkRdCosts( spt, pqData[0], pqData[2], decisions[1], decisions[3]);
      m_prevStates[3].checkRdCosts( spt, pqData[3], pqData[1], decisions[6], decisions[4]);
      m_prevStates[4].checkRdCosts( spt, pqData[0], pqData[2], decisions[2], decisions[0]);
      m_prevStates[5].checkRdCosts( spt, pqData[3], pqData[1], decisions[4], decisions[6]);
      m_prevStates[6].checkRdCosts( spt, pqData[0], pqData[2], decisions[3], decisions[1]);
      m_prevStates[7].checkRdCosts( spt, pqData[3], pqData[1], decisions[7], decisions[5]);
    }
    else
    {
      m_prevStates[0].checkRdCosts( spt, pqData[0], pqData[2], decisions[0], decisions[1]);
      m_prevStates[1].checkRdCosts( spt, pqData[3], pqData[1], decisions[2], decisions[3]);
      m_prevStates[2].checkRdCosts( spt, pqData[0], pqData[2], decisions[1], decisions[0]);
      m_prevStates[3].checkRdCosts( spt, pqData[3], pqData[1], decisions[3], decisions[2]);
    }
#else
    m_prevStates[0].checkRdCosts( spt, pqData[0], pqData[2], decisions[0], decisions[2]);
    m_prevStates[1].checkRdCosts( spt, pqData[0], pqData[2], decisions[2], decisions[0]);
    m_prevStates[2].checkRdCosts( spt, pqData[3], pqData[1], decisions[1], decisions[3]);
    m_prevStates[3].checkRdCosts( spt, pqData[3], pqData[1], decisions[3], decisions[1]);
#endif
    if( spt==SCAN_EOCSBB )
    {
#if TCQ_8STATES
      m_skipStates[0].checkRdCostSkipSbb( decisions[0], numStates );
      m_skipStates[1].checkRdCostSkipSbb( decisions[1], numStates );
      m_skipStates[2].checkRdCostSkipSbb( decisions[2], numStates );
      m_skipStates[3].checkRdCostSkipSbb( decisions[3], numStates );
      m_skipStates[4].checkRdCostSkipSbb( decisions[4], numStates );
      m_skipStates[5].checkRdCostSkipSbb( decisions[5], numStates );
      m_skipStates[6].checkRdCostSkipSbb( decisions[6], numStates );
      m_skipStates[7].checkRdCostSkipSbb( decisions[7], numStates );
#else
        m_skipStates[0].checkRdCostSkipSbb( decisions[0] );
        m_skipStates[1].checkRdCostSkipSbb( decisions[1] );
        m_skipStates[2].checkRdCostSkipSbb( decisions[2] );
        m_skipStates[3].checkRdCostSkipSbb( decisions[3] );
#endif
    }
#if TCQ_8STATES
    if( numStates > 4 )
    {
      m_startState.checkRdCostStart( lastOffset, pqData[0], decisions[0] );
      m_startState.checkRdCostStart( lastOffset, pqData[2], decisions[2] );
    }
    else
    {
      m_startState.checkRdCostStart( lastOffset, pqData[0], decisions[0] );
      m_startState.checkRdCostStart( lastOffset, pqData[2], decisions[1] );
    }
#else
    m_startState.checkRdCostStart( lastOffset, pqData[0], decisions[0] );
    m_startState.checkRdCostStart( lastOffset, pqData[2], decisions[2] );
#endif
  }

#if TCQ_8STATES
	template<int numStates>
#endif
#if TU_256 && !EXTENDED_LFNST
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT_VS
  void DepQuant::xDecideAndUpdate( const TCoeff absCoeff, const ScanInfo& scanInfo, TCoeff quantCoeff )
#else
  void DepQuant::xDecideAndUpdate( const TCoeff absCoeff, const ScanInfo& scanInfo, int quantCoeff )
#endif
#else
#if JVET_R0351_HIGH_BIT_DEPTH_SUPPORT_VS
  void DepQuant::xDecideAndUpdate( const TCoeff absCoeff, const ScanInfo& scanInfo, bool zeroOut, TCoeff quantCoeff )
#else
  void DepQuant::xDecideAndUpdate( const TCoeff absCoeff, const ScanInfo& scanInfo, bool zeroOut, int quantCoeff )
#endif  
#endif
  {
#if TCQ_8STATES
    Decision* decisions = m_trellis + scanInfo.scanIdx * (numStates<<1);
#else
    Decision* decisions = m_trellis[ scanInfo.scanIdx ];
#endif

    std::swap( m_prevStates, m_currStates );
#if TU_256 && !EXTENDED_LFNST
#if TCQ_8STATES
    xDecide<numStates>( scanInfo.spt, absCoeff, lastOffset(scanInfo.scanIdx), decisions, quantCoeff );
#else
    xDecide( scanInfo.spt, absCoeff, lastOffset(scanInfo.scanIdx), decisions, quantCoeff );
#endif
#else
#if TCQ_8STATES
		xDecide<numStates>(scanInfo.spt, absCoeff, lastOffset(scanInfo.scanIdx), decisions, zeroOut, quantCoeff);
#else
    xDecide( scanInfo.spt, absCoeff, lastOffset(scanInfo.scanIdx), decisions, zeroOut, quantCoeff );
#endif
#endif
    if( scanInfo.scanIdx )
    {
      if( scanInfo.eosbb )
      {
        m_commonCtx.swap();
#if TCQ_8STATES
        m_currStates[0].updateStateEOS( scanInfo, m_prevStates, m_skipStates, decisions[0], numStates );
        m_currStates[1].updateStateEOS( scanInfo, m_prevStates, m_skipStates, decisions[1], numStates );
        m_currStates[2].updateStateEOS( scanInfo, m_prevStates, m_skipStates, decisions[2], numStates );
        m_currStates[3].updateStateEOS( scanInfo, m_prevStates, m_skipStates, decisions[3], numStates );

        if( numStates > 4 )
        {
          m_currStates[4].updateStateEOS( scanInfo, m_prevStates, m_skipStates, decisions[4], numStates );
          m_currStates[5].updateStateEOS( scanInfo, m_prevStates, m_skipStates, decisions[5], numStates );
          m_currStates[6].updateStateEOS( scanInfo, m_prevStates, m_skipStates, decisions[6], numStates );
          m_currStates[7].updateStateEOS( scanInfo, m_prevStates, m_skipStates, decisions[7], numStates );
        }

        ::memcpy(decisions + numStates, decisions, numStates * sizeof(Decision));
#else
        m_currStates[0].updateStateEOS( scanInfo, m_prevStates, m_skipStates, decisions[0] );
        m_currStates[1].updateStateEOS( scanInfo, m_prevStates, m_skipStates, decisions[1] );
        m_currStates[2].updateStateEOS( scanInfo, m_prevStates, m_skipStates, decisions[2] );
        m_currStates[3].updateStateEOS( scanInfo, m_prevStates, m_skipStates, decisions[3] );
        ::memcpy( decisions+4, decisions, 4*sizeof(Decision) );
#endif
      }
#if TU_256 && !EXTENDED_LFNST
      else
#else
      else if( !zeroOut )
#endif
      {
#if TCQ_8STATES
        switch( scanInfo.nextNbInfoSbb.num )
        {
        case 0:
          m_currStates[0].updateState<0>( scanInfo, m_prevStates, decisions[0] );
          m_currStates[1].updateState<0>( scanInfo, m_prevStates, decisions[1] );
          m_currStates[2].updateState<0>( scanInfo, m_prevStates, decisions[2] );
          m_currStates[3].updateState<0>( scanInfo, m_prevStates, decisions[3] );

          if( numStates > 4 )
          {
            m_currStates[4].updateState<0>( scanInfo, m_prevStates, decisions[4] );
            m_currStates[5].updateState<0>( scanInfo, m_prevStates, decisions[5] );
            m_currStates[6].updateState<0>( scanInfo, m_prevStates, decisions[6] );
            m_currStates[7].updateState<0>( scanInfo, m_prevStates, decisions[7] );
          }
          break;
        case 1:
          m_currStates[0].updateState<1>( scanInfo, m_prevStates, decisions[0] );
          m_currStates[1].updateState<1>( scanInfo, m_prevStates, decisions[1] );
          m_currStates[2].updateState<1>( scanInfo, m_prevStates, decisions[2] );
          m_currStates[3].updateState<1>( scanInfo, m_prevStates, decisions[3] );

          if( numStates > 4 )
          {
            m_currStates[4].updateState<1>( scanInfo, m_prevStates, decisions[4] );
            m_currStates[5].updateState<1>( scanInfo, m_prevStates, decisions[5] );
            m_currStates[6].updateState<1>( scanInfo, m_prevStates, decisions[6] );
            m_currStates[7].updateState<1>( scanInfo, m_prevStates, decisions[7] );
          }
          break;
        case 2:
          m_currStates[0].updateState<2>( scanInfo, m_prevStates, decisions[0] );
          m_currStates[1].updateState<2>( scanInfo, m_prevStates, decisions[1] );
          m_currStates[2].updateState<2>( scanInfo, m_prevStates, decisions[2] );
          m_currStates[3].updateState<2>( scanInfo, m_prevStates, decisions[3] );

          if( numStates > 4 )
          {
            m_currStates[4].updateState<2>( scanInfo, m_prevStates, decisions[4] );
            m_currStates[5].updateState<2>( scanInfo, m_prevStates, decisions[5] );
            m_currStates[6].updateState<2>( scanInfo, m_prevStates, decisions[6] );
            m_currStates[7].updateState<2>( scanInfo, m_prevStates, decisions[7] );
          }
          break;
        case 3:
          m_currStates[0].updateState<3>( scanInfo, m_prevStates, decisions[0] );
          m_currStates[1].updateState<3>( scanInfo, m_prevStates, decisions[1] );
          m_currStates[2].updateState<3>( scanInfo, m_prevStates, decisions[2] );
          m_currStates[3].updateState<3>( scanInfo, m_prevStates, decisions[3] );

          if( numStates > 4 )
          {
            m_currStates[4].updateState<3>( scanInfo, m_prevStates, decisions[4] );
            m_currStates[5].updateState<3>( scanInfo, m_prevStates, decisions[5] );
            m_currStates[6].updateState<3>( scanInfo, m_prevStates, decisions[6] );
            m_currStates[7].updateState<3>( scanInfo, m_prevStates, decisions[7] );
          }
          break;
        case 4:
          m_currStates[0].updateState<4>( scanInfo, m_prevStates, decisions[0] );
          m_currStates[1].updateState<4>( scanInfo, m_prevStates, decisions[1] );
          m_currStates[2].updateState<4>( scanInfo, m_prevStates, decisions[2] );
          m_currStates[3].updateState<4>( scanInfo, m_prevStates, decisions[3] );

          if( numStates > 4 )
          {
            m_currStates[4].updateState<4>( scanInfo, m_prevStates, decisions[4] );
            m_currStates[5].updateState<4>( scanInfo, m_prevStates, decisions[5] );
            m_currStates[6].updateState<4>( scanInfo, m_prevStates, decisions[6] );
            m_currStates[7].updateState<4>( scanInfo, m_prevStates, decisions[7] );
          }
          break;
        default:
          m_currStates[0].updateState<5>( scanInfo, m_prevStates, decisions[0] );
          m_currStates[1].updateState<5>( scanInfo, m_prevStates, decisions[1] );
          m_currStates[2].updateState<5>( scanInfo, m_prevStates, decisions[2] );
          m_currStates[3].updateState<5>( scanInfo, m_prevStates, decisions[3] );

          if( numStates > 4 )
          {
            m_currStates[4].updateState<5>( scanInfo, m_prevStates, decisions[4] );
            m_currStates[5].updateState<5>( scanInfo, m_prevStates, decisions[5] );
            m_currStates[6].updateState<5>( scanInfo, m_prevStates, decisions[6] );
            m_currStates[7].updateState<5>( scanInfo, m_prevStates, decisions[7] );
          }
        }
#else
        switch( scanInfo.nextNbInfoSbb.num )
        {
        case 0:
          m_currStates[0].updateState<0>( scanInfo, m_prevStates, decisions[0] );
          m_currStates[1].updateState<0>( scanInfo, m_prevStates, decisions[1] );
          m_currStates[2].updateState<0>( scanInfo, m_prevStates, decisions[2] );
          m_currStates[3].updateState<0>( scanInfo, m_prevStates, decisions[3] );
          break;
        case 1:
          m_currStates[0].updateState<1>( scanInfo, m_prevStates, decisions[0] );
          m_currStates[1].updateState<1>( scanInfo, m_prevStates, decisions[1] );
          m_currStates[2].updateState<1>( scanInfo, m_prevStates, decisions[2] );
          m_currStates[3].updateState<1>( scanInfo, m_prevStates, decisions[3] );
          break;
        case 2:
          m_currStates[0].updateState<2>( scanInfo, m_prevStates, decisions[0] );
          m_currStates[1].updateState<2>( scanInfo, m_prevStates, decisions[1] );
          m_currStates[2].updateState<2>( scanInfo, m_prevStates, decisions[2] );
          m_currStates[3].updateState<2>( scanInfo, m_prevStates, decisions[3] );
          break;
        case 3:
          m_currStates[0].updateState<3>( scanInfo, m_prevStates, decisions[0] );
          m_currStates[1].updateState<3>( scanInfo, m_prevStates, decisions[1] );
          m_currStates[2].updateState<3>( scanInfo, m_prevStates, decisions[2] );
          m_currStates[3].updateState<3>( scanInfo, m_prevStates, decisions[3] );
          break;
        case 4:
          m_currStates[0].updateState<4>( scanInfo, m_prevStates, decisions[0] );
          m_currStates[1].updateState<4>( scanInfo, m_prevStates, decisions[1] );
          m_currStates[2].updateState<4>( scanInfo, m_prevStates, decisions[2] );
          m_currStates[3].updateState<4>( scanInfo, m_prevStates, decisions[3] );
          break;
        default:
          m_currStates[0].updateState<5>( scanInfo, m_prevStates, decisions[0] );
          m_currStates[1].updateState<5>( scanInfo, m_prevStates, decisions[1] );
          m_currStates[2].updateState<5>( scanInfo, m_prevStates, decisions[2] );
          m_currStates[3].updateState<5>( scanInfo, m_prevStates, decisions[3] );
        }
#endif
      }

      if( scanInfo.spt == SCAN_SOCSBB )
      {
        std::swap( m_prevStates, m_skipStates );
      }
    }
  }


  void DepQuant::quant( TransformUnit& tu, const CCoeffBuf& srcCoeff, const ComponentID compID, const QpParam& cQP, const double lambda, const Ctx& ctx, TCoeff& absSum, bool enableScalingLists, int* quantCoeff )
  {
    CHECKD( tu.cs->sps->getSpsRangeExtension().getExtendedPrecisionProcessingFlag(), "ext precision is not supported" );
#if SIGN_PREDICTION
    CoeffBuf signBuff = tu.getCoeffSigns(compID);
    uint32_t uiWidth  = tu.blocks[compID].width;
    uint32_t uiHeight = tu.blocks[compID].height;

    if ( tu.cs->sps->getNumPredSigns() > 0  && uiHeight >= 4 && uiWidth >= 4)
    {
      TCoeff *coeff = signBuff.buf;
      for (uint32_t y = 0; y < SIGN_PRED_FREQ_RANGE; y++)
      {
        coeff[0] = TrQuant::SIGN_PRED_BYPASS;
        coeff[1] = TrQuant::SIGN_PRED_BYPASS;
        coeff[2] = TrQuant::SIGN_PRED_BYPASS;
        coeff[3] = TrQuant::SIGN_PRED_BYPASS;
        coeff += signBuff.stride;
      }
    }
#endif

    //===== reset / pre-init =====
#if TCQ_8STATES
    const int numStates = 2 << tu.cs->slice->getDepQuantEnabledIdc();
#endif
    const TUParameters& tuPars  = *g_Rom.getTUPars( tu.blocks[compID], compID );
    m_quant.initQuantBlock    ( tu, compID, cQP, lambda );
    TCoeff*       qCoeff      = tu.getCoeffs( compID ).buf;
    const TCoeff* tCoeff      = srcCoeff.buf;
    const int     numCoeff    = tu.blocks[compID].area();
    ::memset( tu.getCoeffs( compID ).buf, 0x00, numCoeff*sizeof(TCoeff) );
    absSum          = 0;

    const CompArea& area     = tu.blocks[ compID ];
    const uint32_t  width    = area.width;
    const uint32_t  height   = area.height;
    const uint32_t  lfnstIdx = tu.cu->lfnstIdx;

    //===== scaling matrix ====
    //const int         qpDQ = cQP.Qp + 1;
    //const int         qpPer = qpDQ / 6;
    //const int         qpRem = qpDQ - 6 * qpPer;

    //TCoeff thresTmp = thres;
    bool zeroOut = false;
#if EXTENDED_LFNST || !TU_256
    int effWidth = tuPars.m_width, effHeight = tuPars.m_height; 
#endif
#if !TU_256
    if( ( tu.mtsIdx[compID] > MTS_SKIP || (tu.cs->sps->getUseMTS() && tu.cu->sbtInfo != 0 && tuPars.m_height <= 32 && tuPars.m_width <= 32)) && compID == COMPONENT_Y )
    {
      effHeight = (tuPars.m_height == 32) ? 16 : tuPars.m_height;
      effWidth = (tuPars.m_width == 32) ? 16 : tuPars.m_width;
      zeroOut = (effHeight < tuPars.m_height || effWidth < tuPars.m_width);
    }
#endif
#if EXTENDED_LFNST
    int lfnst_threshold = 0;
    if (lfnstIdx > 0 && tu.mtsIdx[compID] != MTS_SKIP && width >= 4 && height >= 4)
    {
      const bool whge3 = width >= 8 && height >= 8;
      lfnst_threshold = whge3 ? 8 : 4;
      zeroOut = true;
      effWidth = lfnst_threshold;
      effHeight = lfnst_threshold;
    }
#endif
    bool zeroOutforThres = zeroOut;
#if !TU_256
    zeroOutforThres = zeroOut || (32 < tuPars.m_height || 32 < tuPars.m_width);
#endif

    //===== find first test position =====
    int firstTestPos = numCoeff - 1;
#if !EXTENDED_LFNST
    if (lfnstIdx > 0 && tu.mtsIdx[compID] != MTS_SKIP && width >= 4 && height >= 4)
    {
#if JVET_W0119_LFNST_EXTENSION
      firstTestPos = PU::getLFNSTMatrixDim( width, height ) - 1;
#else
      firstTestPos = ( ( width == 4 && height == 4 ) || ( width == 8 && height == 8 ) )  ? 7 : 15 ;
#endif
    }
#endif
    const TCoeff defaultQuantisationCoefficient = (TCoeff)m_quant.getQScale();
    const TCoeff thres = m_quant.getLastThreshold();
    for( ; firstTestPos >= 0; firstTestPos-- )
    {
      if (zeroOutforThres && (tuPars.m_scanId2BlkPos[firstTestPos].x >= ((tuPars.m_width == 32 && zeroOut) ? 16 : 32)
                           || tuPars.m_scanId2BlkPos[firstTestPos].y >= ((tuPars.m_height == 32 && zeroOut) ? 16 : 32)))
        continue;

#if EXTENDED_LFNST
      if( lfnst_threshold > 0 && ( tuPars.m_scanId2BlkPos[firstTestPos].x >= lfnst_threshold || tuPars.m_scanId2BlkPos[firstTestPos].y >= lfnst_threshold ) )
      {
        continue;
      }
#endif
      TCoeff thresTmp = (enableScalingLists) ? TCoeff(thres / (4 * quantCoeff[tuPars.m_scanId2BlkPos[firstTestPos].idx]))
                                             : TCoeff(thres / (4 * defaultQuantisationCoefficient));

      if (abs(tCoeff[tuPars.m_scanId2BlkPos[firstTestPos].idx]) > thresTmp)
      {
        break;
      }
    }
    if( firstTestPos < 0 )
    {
      return;
    }

    //===== real init =====
    RateEstimator::initCtx( tuPars, tu, compID, ctx.getFracBitsAcess() );
    m_commonCtx.reset( tuPars, *this );
#if !TCQ_8STATES 
    for( int k = 0; k < 12; k++ )
    {
      m_allStates[k].init();
    }
    m_startState.init();
#endif
#if TU_256
    int effectWidth  = width;
    int effectHeight = height;
#else
    int effectWidth = std::min(32, effWidth);
    int effectHeight = std::min(32, effHeight);
#endif
#if EXTENDED_LFNST
    effectWidth  = ( lfnst_threshold > 0 ) ? lfnst_threshold : effectWidth;
    effectHeight = ( lfnst_threshold > 0 ) ? lfnst_threshold : effectHeight;
#endif
#if TCQ_8STATES
    for( int k = 0; k < numStates; k++ )
    {
      m_allStates[k   ].init( effectWidth, effectHeight );
      m_allStates[k+8 ].init( effectWidth, effectHeight );
      m_allStates[k+16].init( effectWidth, effectHeight );
    }
    m_startState.init( effectWidth, effectHeight );
#else
    for (int k = 0; k < 12; k++)
    {
      m_allStates[k].effWidth = effectWidth;
      m_allStates[k].effHeight = effectHeight;
    }
    m_startState.effWidth = effectWidth;
    m_startState.effHeight = effectHeight;
#endif
    //===== populate trellis ===== // 
#if TCQ_8STATES
#if TU_256 && !EXTENDED_LFNST
		void (DepQuant::*fdecide)(const TCoeff, const ScanInfo&, int) = (numStates>4 ? &DepQuant::xDecideAndUpdate<8> : &DepQuant::xDecideAndUpdate<4>);
#else
    void (DepQuant::*fdecide)(const TCoeff,const ScanInfo&,bool,int) = ( numStates>4 ? &DepQuant::xDecideAndUpdate<8> : &DepQuant::xDecideAndUpdate<4> );
#endif
#endif
    for( int scanIdx = firstTestPos; scanIdx >= 0; scanIdx-- )
    {
      const ScanInfo& scanInfo = tuPars.m_scanInfo[ scanIdx ];

      if (enableScalingLists)
      {
        m_quant.initQuantBlock(tu, compID, cQP, lambda, quantCoeff[scanInfo.rasterPos]);
#if TU_256 && !EXTENDED_LFNST
#if TCQ_8STATES
				(this->*fdecide)(abs(tCoeff[scanInfo.rasterPos]), scanInfo, quantCoeff[scanInfo.rasterPos]);
#else
        xDecideAndUpdate( abs( tCoeff[scanInfo.rasterPos]), scanInfo, quantCoeff[scanInfo.rasterPos] );
#endif
#else
#if TCQ_8STATES
        (this->*fdecide)( abs( tCoeff[scanInfo.rasterPos]), scanInfo, (zeroOut && (scanInfo.posX >= effWidth || scanInfo.posY >= effHeight)), quantCoeff[scanInfo.rasterPos] );
#else
        xDecideAndUpdate( abs( tCoeff[scanInfo.rasterPos]), scanInfo, (zeroOut && (scanInfo.posX >= effWidth || scanInfo.posY >= effHeight)), quantCoeff[scanInfo.rasterPos] );
#endif
#endif
      }
      else
#if TU_256 && !EXTENDED_LFNST
#if TCQ_8STATES
        (this->*fdecide)( abs( tCoeff[scanInfo.rasterPos]), scanInfo, defaultQuantisationCoefficient );
#else
        xDecideAndUpdate( abs( tCoeff[scanInfo.rasterPos]), scanInfo, defaultQuantisationCoefficient );
#endif
#else
#if TCQ_8STATES
        (this->*fdecide)( abs( tCoeff[scanInfo.rasterPos]), scanInfo, (zeroOut && (scanInfo.posX >= effWidth || scanInfo.posY >= effHeight)), defaultQuantisationCoefficient );
#else
        xDecideAndUpdate( abs( tCoeff[scanInfo.rasterPos]), scanInfo, (zeroOut && (scanInfo.posX >= effWidth || scanInfo.posY >= effHeight)), defaultQuantisationCoefficient );
#endif
#endif
    }

    //===== find best path =====
    Decision  decision    = { std::numeric_limits<int64_t>::max(), -1, -2 };
    int64_t   minPathCost =  0;
#if TCQ_8STATES
    for( int8_t stateId = 0; stateId < numStates; stateId++ )
#else
    for( int8_t stateId = 0; stateId < 4; stateId++ )
#endif
    {
#if TCQ_8STATES
      int64_t pathCost = m_trellis[stateId].rdCost;
#else
      int64_t pathCost = m_trellis[0][stateId].rdCost;
#endif
      if( pathCost < minPathCost )
      {
        decision.prevId = stateId;
        minPathCost     = pathCost;
      }
    }

    //===== backward scanning =====
    int scanIdx = 0;
#if TCQ_8STATES
    for( const Decision* trellisStage = m_trellis; decision.prevId >= 0; scanIdx++, trellisStage += (numStates<<1) )
#else
    for( ; decision.prevId >= 0; scanIdx++ )
#endif
    {
#if TCQ_8STATES
      decision          = trellisStage[ decision.prevId ];
#else
      decision          = m_trellis[ scanIdx ][ decision.prevId ];
#endif
      int32_t blkpos    = tuPars.m_scanId2BlkPos[scanIdx].idx;
      qCoeff[ blkpos ]  = ( tCoeff[ blkpos ] < 0 ? -decision.absLevel : decision.absLevel );
      absSum           += decision.absLevel;
    }
  }

}; // namespace DQIntern




//===== interface class =====
DepQuant::DepQuant( const Quant* other, bool enc ) : QuantRDOQ( other )
{
  const DepQuant* dq = dynamic_cast<const DepQuant*>( other );
  CHECK( other && !dq, "The DepQuant cast must be successfull!" );
  p = new DQIntern::DepQuant();
  if( enc )
  {
    DQIntern::g_Rom.init();
  }
}

DepQuant::~DepQuant()
{
  delete static_cast<DQIntern::DepQuant*>(p);
}

void DepQuant::quant( TransformUnit &tu, const ComponentID &compID, const CCoeffBuf &pSrc, TCoeff &uiAbsSum, const QpParam &cQP, const Ctx& ctx )
{
  const bool useRegularResidualCoding = tu.cu->slice->getTSResidualCodingDisabledFlag() || tu.mtsIdx[compID] != MTS_SKIP;

#if TCQ_8STATES
  if( tu.cs->slice->getDepQuantEnabledIdc() && useRegularResidualCoding )
#else
  if( tu.cs->slice->getDepQuantEnabledFlag() && useRegularResidualCoding )
#endif
  {
    //===== scaling matrix ====
    const int         qpDQ            = cQP.Qp(tu.mtsIdx[compID] == MTS_SKIP) + 1;
    const int         qpPer           = qpDQ / 6;
    const int         qpRem           = qpDQ - 6 * qpPer;
    const CompArea    &rect           = tu.blocks[compID];
    const int         width           = rect.width;
    const int         height          = rect.height;
    uint32_t          scalingListType = getScalingListType(tu.cu->predMode, compID);
    CHECK(scalingListType >= SCALING_LIST_NUM, "Invalid scaling list");
    const uint32_t    log2TrWidth     = floorLog2(width);
    const uint32_t    log2TrHeight    = floorLog2(height);

    const bool        disableSMForLFNST = tu.cs->slice->getExplicitScalingListUsed() ? tu.cs->slice->getSPS()->getDisableScalingMatrixForLfnstBlks() : false;
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
    const bool        isLfnstApplied = tu.cu->lfnstIdx > 0 && (tu.cu->isSepTree() ? true : isLuma(compID));
#else
    const bool        isLfnstApplied = tu.cu->lfnstIdx > 0 && (CS::isDualITree(*tu.cs) ? true : isLuma(compID));
#endif
    const bool        disableSMForACT = tu.cs->slice->getSPS()->getScalingMatrixForAlternativeColourSpaceDisabledFlag() && (tu.cs->slice->getSPS()->getScalingMatrixDesignatedColourSpaceFlag() == tu.cu->colorTransform);
    const bool        enableScalingLists = getUseScalingList(width, height, (tu.mtsIdx[compID] == MTS_SKIP), isLfnstApplied, disableSMForLFNST, disableSMForACT);
    static_cast<DQIntern::DepQuant*>(p)->quant( tu, pSrc, compID, cQP, Quant::m_dLambda, ctx, uiAbsSum, enableScalingLists, Quant::getQuantCoeff(scalingListType, qpRem, log2TrWidth, log2TrHeight) );
  }
  else
  {
    QuantRDOQ::quant( tu, compID, pSrc, uiAbsSum, cQP, ctx );
  }
}

#if SIGN_PREDICTION
uint32_t DepQuant::getPredictedSigns( TransformUnit& tu, const ComponentID compID, std::vector<Position> &predSignsXY )
{
  uint32_t numPredSigns = 0;
  bool bUseSignPred = TU::getUseSignPred( tu, compID );

  if( bUseSignPred )
  {
    const int32_t  maxNumPredSigns = tu.cs->sps->getNumPredSigns();
    CoeffBuf bufQCoeffs = tu.getCoeffs( compID );
    CoeffBuf bufSigns = tu.getCoeffSigns( compID );
    TCoeff *coeff = bufQCoeffs.buf;
    TCoeff *signs = bufSigns.buf;

    for( uint32_t uiY = 0; uiY < SIGN_PRED_FREQ_RANGE; uiY++ )
    {
      for( uint32_t uiX = 0; uiX < SIGN_PRED_FREQ_RANGE; uiX++ )
      {
        uint32_t uiLevel = abs( coeff[uiX] );

        if( uiLevel )
        {
          if( signs[uiX] != TrQuant::SIGN_PRED_HIDDEN )
          {
            predSignsXY.push_back( Position( uiX, uiY ) );

            if( ++numPredSigns >= maxNumPredSigns )
            {
              break;
            }

          }
        }
      }
      if( numPredSigns == maxNumPredSigns )
      {
        break;
      }

      coeff += bufQCoeffs.stride;
      signs += bufSigns.stride;
    }
  }
  return numPredSigns;
}
#endif

void DepQuant::dequant( const TransformUnit &tu, CoeffBuf &dstCoeff, const ComponentID &compID, const QpParam &cQP )
{
  const bool useRegularResidualCoding = tu.cu->slice->getTSResidualCodingDisabledFlag() || tu.mtsIdx[compID] != MTS_SKIP;

#if TCQ_8STATES
  if( tu.cs->slice->getDepQuantEnabledIdc() && useRegularResidualCoding )
#else
  if( tu.cs->slice->getDepQuantEnabledFlag() && useRegularResidualCoding )
#endif
  {
    const int         qpDQ            = cQP.Qp(tu.mtsIdx[compID] == MTS_SKIP) + 1;
    const int         qpPer           = qpDQ / 6;
    const int         qpRem           = qpDQ - 6 * qpPer;
    const CompArea    &rect           = tu.blocks[compID];
    const int         width           = rect.width;
    const int         height          = rect.height;
    uint32_t          scalingListType = getScalingListType(tu.cu->predMode, compID);
    CHECK(scalingListType >= SCALING_LIST_NUM, "Invalid scaling list");
    const uint32_t    log2TrWidth  = floorLog2(width);
    const uint32_t    log2TrHeight = floorLog2(height);

    const bool disableSMForLFNST = tu.cs->slice->getExplicitScalingListUsed() ? tu.cs->slice->getSPS()->getDisableScalingMatrixForLfnstBlks() : false;
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
    const bool isLfnstApplied = tu.cu->lfnstIdx > 0 && (tu.cu->isSepTree() ? true : isLuma(compID));
#else
    const bool isLfnstApplied = tu.cu->lfnstIdx > 0 && (CS::isDualITree(*tu.cs) ? true : isLuma(compID));
#endif
    const bool disableSMForACT = tu.cs->slice->getSPS()->getScalingMatrixForAlternativeColourSpaceDisabledFlag() && (tu.cs->slice->getSPS()->getScalingMatrixDesignatedColourSpaceFlag() == tu.cu->colorTransform);
    const bool enableScalingLists = getUseScalingList(width, height, (tu.mtsIdx[compID] == MTS_SKIP), isLfnstApplied, disableSMForLFNST, disableSMForACT);
    static_cast<DQIntern::DepQuant*>(p)->dequant( tu, dstCoeff, compID, cQP, enableScalingLists, Quant::getDequantCoeff(scalingListType, qpRem, log2TrWidth, log2TrHeight) );
  }
  else
  {
    QuantRDOQ::dequant( tu, dstCoeff, compID, cQP );
  }
}



