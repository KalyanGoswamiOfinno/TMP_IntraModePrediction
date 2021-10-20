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

/** \file     Contexts.h
 *  \brief    Classes providing probability descriptions and contexts (header)
 */

#ifndef __CONTEXTS__
#define __CONTEXTS__

#include "CommonDef.h"
#include "Slice.h"

#include <vector>

static constexpr int     PROB_BITS   = 15;   // Nominal number of bits to represent probabilities
#if EC_HIGH_PRECISION
static constexpr int     PROB_BITS_0 = 15;   // Number of bits to represent 1st estimate
static constexpr int     PROB_BITS_1 = 15;   // Number of bits to represent 2nd estimate
#else
static constexpr int     PROB_BITS_0 = 10;   // Number of bits to represent 1st estimate
static constexpr int     PROB_BITS_1 = 14;   // Number of bits to represent 2nd estimate
#endif
static constexpr int     MASK_0      = ~(~0u << PROB_BITS_0) << (PROB_BITS - PROB_BITS_0);
static constexpr int     MASK_1      = ~(~0u << PROB_BITS_1) << (PROB_BITS - PROB_BITS_1);
static constexpr uint8_t DWS         = 8;   // 0x47 Default window sizes

struct BinFracBits
{
  uint32_t intBits[2];
};


enum BPMType
{
  BPM_Undefined = 0,
  BPM_Std,
  BPM_NUM
};

class ProbModelTables
{
protected:
#if EC_HIGH_PRECISION
	static const BinFracBits m_binFracBits[512];
#else
  static const BinFracBits m_binFracBits[256];
#endif
  static const uint8_t      m_RenormTable_32  [ 32];          // Std         MP   MPI
};

class BinProbModelBase : public ProbModelTables
{
public:
  BinProbModelBase () {}
  ~BinProbModelBase() {}
  static uint32_t estFracBitsEP ()                    { return  (       1 << SCALE_BITS ); }
  static uint32_t estFracBitsEP ( unsigned numBins )  { return  ( numBins << SCALE_BITS ); }
};

class BinProbModel_Std : public BinProbModelBase
{
public:
  BinProbModel_Std()
  {
    uint16_t half = 1 << (PROB_BITS - 1);
    m_state[0]    = half;
    m_state[1]    = half;
    m_rate        = DWS;
  }
  ~BinProbModel_Std ()                {}
public:
  void            init              ( int qp, int initId );
  void update(unsigned bin)
  {
    int rate0 = m_rate >> 4;
    int rate1 = m_rate & 15;

    m_state[0] -= (m_state[0] >> rate0) & MASK_0;
    m_state[1] -= (m_state[1] >> rate1) & MASK_1;
    if (bin)
    {
      m_state[0] += (0x7fffu >> rate0) & MASK_0;
      m_state[1] += (0x7fffu >> rate1) & MASK_1;
    }
  }
  void setLog2WindowSize(uint8_t log2WindowSize)
  {
    int rate0 = 2 + ((log2WindowSize >> 2) & 3);
    int rate1 = 3 + rate0 + (log2WindowSize & 3);
    m_rate    = 16 * rate0 + rate1;
    CHECK(rate1 > 9, "Second window size is too large!");
  }
  void estFracBitsUpdate(unsigned bin, uint64_t &b)
  {
    b += estFracBits(bin);
    update(bin);
  }
  uint32_t        estFracBits(unsigned bin) const { return getFracBitsArray().intBits[bin]; }
  static uint32_t estFracBitsTrm(unsigned bin) { return (bin ? 0x3bfbb : 0x0010c); }
#if EC_HIGH_PRECISION
	BinFracBits     getFracBitsArray() const { return m_binFracBits[state_est()]; }
#else
  BinFracBits     getFracBitsArray() const { return m_binFracBits[state()]; }
#endif
public:
#if EC_HIGH_PRECISION
	uint16_t state_est() const { return (m_state[0] + m_state[1]) >> 7; }
	uint16_t state() const { return (m_state[0] + m_state[1]) >> 1; }
	uint8_t mps() const { return state() >> 14; }
	uint8_t getLPS(unsigned range) const
	{
		uint16_t q = state();
		if (q & 0x4000)
			q = q ^ 0x7fff;
		return ((range * (q >> 6)) >> 9) + 1;
	}
#else
  uint8_t state() const { return (m_state[0] + m_state[1]) >> 8; }
  uint8_t mps() const { return state() >> 7; }
  uint8_t getLPS(unsigned range) const
  {
    uint16_t q = state();
    if (q & 0x80)
      q = q ^ 0xff;
    return ((q >> 2) * (range >> 5) >> 1) + 4;
  }
#endif
  static uint8_t  getRenormBitsLPS  ( unsigned LPS )                    { return    m_RenormTable_32  [LPS>>3]; }
  static uint8_t  getRenormBitsRange( unsigned range )                  { return    1; }
  uint16_t getState() const { return m_state[0] + m_state[1]; }
  void     setState(uint16_t pState)
  {
    m_state[0] = (pState >> 1) & MASK_0;
    m_state[1] = (pState >> 1) & MASK_1;
  }
public:
  uint64_t estFracExcessBits(const BinProbModel_Std &r) const
  {
#if EC_HIGH_PRECISION
		int n = 2 * state_est() + 1;
		return ((1024 - n) * r.estFracBits(0) + n * r.estFracBits(1) + 512) >> 10;
#else
    int n = 2 * state() + 1;
    return ((512 - n) * r.estFracBits(0) + n * r.estFracBits(1) + 256) >> 9;
#endif
  }
private:
  uint16_t m_state[2];
  uint8_t  m_rate;
};



class CtxSet
{
public:
  CtxSet( uint16_t offset, uint16_t size ) : Offset( offset ), Size( size ) {}
  CtxSet( const CtxSet& ctxSet ) : Offset( ctxSet.Offset ), Size( ctxSet.Size ) {}
  CtxSet( std::initializer_list<CtxSet> ctxSets );
public:
  uint16_t  operator()  ()  const
  {
    return Offset;
  }
  uint16_t  operator()  ( uint16_t inc )  const
  {
    CHECKD( inc >= Size, "Specified context increment (" << inc << ") exceed range of context set [0;" << Size - 1 << "]." );
    return Offset + inc;
  }
  bool operator== ( const CtxSet& ctxSet ) const
  {
    return ( Offset == ctxSet.Offset && Size == ctxSet.Size );
  }
  bool operator!= ( const CtxSet& ctxSet ) const
  {
    return ( Offset != ctxSet.Offset || Size != ctxSet.Size );
  }
public:
  uint16_t  Offset;
  uint16_t  Size;
};



class ContextSetCfg
{
public:
  // context sets: specify offset and size
  static const CtxSet   SplitFlag;
  static const CtxSet   SplitQtFlag;
  static const CtxSet   SplitHvFlag;
  static const CtxSet   Split12Flag;
#if !INTRA_RM_SMALL_BLOCK_SIZE_CONSTRAINTS
  static const CtxSet   ModeConsFlag;
#endif
  static const CtxSet   SkipFlag;
  static const CtxSet   MergeFlag;
  static const CtxSet   RegularMergeFlag;
  static const CtxSet   MergeIdx;
#if TM_MRG
  static const CtxSet   TmMergeIdx;
#endif
  static const CtxSet   PredMode;
  static const CtxSet   MultiRefLineIdx;
  static const CtxSet   IntraLumaMpmFlag;
#if SECONDARY_MPM
  static const CtxSet   IntraLumaSecondMpmFlag;
#endif
  static const CtxSet   IntraLumaPlanarFlag;
#if SECONDARY_MPM
  static const CtxSet   IntraLumaMPMIdx;
#endif
  static const CtxSet   CclmModeFlag;
  static const CtxSet   CclmModeIdx;
  static const CtxSet   IntraChromaPredMode;
  static const CtxSet   MipFlag;
#if JVET_V0130_INTRA_TMP
  static const CtxSet   TmpFlag;
#endif
#if MMLM
  static const CtxSet   MMLMFlag;
#endif
  static const CtxSet   DeltaQP;
  static const CtxSet   InterDir;
  static const CtxSet   RefPic;
  static const CtxSet   MmvdFlag;
  static const CtxSet   MmvdMergeIdx;
  static const CtxSet   MmvdStepMvpIdx;
#if JVET_W0097_GPM_MMVD_TM
  static const CtxSet   GeoMmvdFlag;
  static const CtxSet   GeoMmvdStepMvpIdx;
#endif
  static const CtxSet   SubblockMergeFlag;
  static const CtxSet   AffineFlag;
  static const CtxSet   AffineType;
  static const CtxSet   AffMergeIdx;
#if AFFINE_MMVD
  static const CtxSet   AfMmvdFlag;
  static const CtxSet   AfMmvdIdx;
  static const CtxSet   AfMmvdOffsetStep;
#endif
#if TM_MRG
  static const CtxSet   TMMergeFlag;
#endif
  static const CtxSet   Mvd;
#if MULTI_HYP_PRED
  static const CtxSet   MultiHypothesisFlag;
  static const CtxSet   MHRefPic;
  static const CtxSet   MHWeight;
#endif
  static const CtxSet   BDPCMMode;
  static const CtxSet   QtRootCbf;
  static const CtxSet   ACTFlag;
  static const CtxSet   QtCbf           [3];    // [ channel ]
  static const CtxSet   SigCoeffGroup   [2];    // [ ChannelType ]
  static const CtxSet   LastX           [2];    // [ ChannelType ]
  static const CtxSet   LastY           [2];    // [ ChannelType ]
  static const CtxSet   SigFlag         [6];    // [ ChannelType + State ]
  static const CtxSet   ParFlag         [2];    // [ ChannelType ]
  static const CtxSet   GtxFlag         [4];    // [ ChannelType + x ]
  static const CtxSet   TsSigCoeffGroup;
  static const CtxSet   TsSigFlag;
  static const CtxSet   TsParFlag;
  static const CtxSet   TsGtxFlag;
  static const CtxSet   TsLrg1Flag;
  static const CtxSet   TsResidualSign;
  static const CtxSet   MVPIdx;
  static const CtxSet   SaoMergeFlag;
  static const CtxSet   SaoTypeIdx;
#if JVET_V0094_BILATERAL_FILTER
  static const CtxSet   BifCtrlFlags;
#endif
#if JVET_W0066_CCSAO
  static const CtxSet   CcSaoControlIdc;
#endif
  static const CtxSet   TransformSkipFlag;
  static const CtxSet   MTSIdx;
  static const CtxSet   LFNSTIdx;
  static const CtxSet   PLTFlag;
  static const CtxSet   RotationFlag;
  static const CtxSet   RunTypeFlag;
  static const CtxSet   IdxRunModel;
  static const CtxSet   CopyRunModel;
  static const CtxSet   SbtFlag;
  static const CtxSet   SbtQuadFlag;
  static const CtxSet   SbtHorFlag;
  static const CtxSet   SbtPosFlag;
  static const CtxSet   ChromaQpAdjFlag;
  static const CtxSet   ChromaQpAdjIdc;
  static const CtxSet   ImvFlag;
#if ENABLE_DIMD
  static const CtxSet   DimdFlag;
#endif
#if JVET_W0123_TIMD_FUSION
  static const CtxSet   TimdFlag;
#endif
#if ENABLE_OBMC
  static const CtxSet   ObmcFlag;
#endif 
  static const CtxSet   BcwIdx;
  static const CtxSet   ctbAlfFlag;
  static const CtxSet   ctbAlfAlternative;
  static const CtxSet   AlfUseTemporalFilt;
  static const CtxSet   CcAlfFilterControlFlag;
  static const CtxSet   CiipFlag;
  static const CtxSet   SmvdFlag;
  static const CtxSet   IBCFlag;
  static const CtxSet   ISPMode;
  static const CtxSet   JointCbCrFlag;
#if INTER_LIC
  static const CtxSet   LICFlag;
#endif
#if SIGN_PREDICTION
  static const CtxSet   signPred[2];
#endif
  static const unsigned NumberOfContexts;

  // combined sets for less complex copying
  // NOTE: The contained CtxSet's should directly follow each other in the initalization list;
  //       otherwise, you will copy more elements than you want !!!
  static const CtxSet   Sao;
  static const CtxSet   Alf;
  static const CtxSet   Palette;

public:
  static const std::vector<uint8_t>&  getInitTable( unsigned initId );
private:
  static std::vector<std::vector<uint8_t> > sm_InitTables;
  static CtxSet addCtxSet( std::initializer_list<std::initializer_list<uint8_t> > initSet2d );
};



class FracBitsAccess
{
public:
  virtual BinFracBits getFracBitsArray( unsigned ctxId ) const = 0;
};



template <class BinProbModel>
class CtxStore : public FracBitsAccess
{
public:
  CtxStore();
  CtxStore( bool dummy );
  CtxStore( const CtxStore<BinProbModel>& ctxStore );
public:
  void copyFrom   ( const CtxStore<BinProbModel>& src )                        { checkInit(); ::memcpy( m_Ctx,               src.m_Ctx,               sizeof( BinProbModel ) * ContextSetCfg::NumberOfContexts ); }
  void copyFrom   ( const CtxStore<BinProbModel>& src, const CtxSet& ctxSet )  { checkInit(); ::memcpy( m_Ctx+ctxSet.Offset, src.m_Ctx+ctxSet.Offset, sizeof( BinProbModel ) * ctxSet.Size ); }
  void init       ( int qp, int initId );
  void setWinSizes( const std::vector<uint8_t>&   log2WindowSizes );
  void loadPStates( const std::vector<uint16_t>&  probStates );
  void savePStates( std::vector<uint16_t>&        probStates )  const;

  const BinProbModel& operator[]      ( unsigned  ctxId  )  const { return m_Ctx[ctxId]; }
  BinProbModel&       operator[]      ( unsigned  ctxId  )        { return m_Ctx[ctxId]; }
  uint32_t            estFracBits     ( unsigned  bin,
                                        unsigned  ctxId  )  const { return m_Ctx[ctxId].estFracBits(bin); }

  BinFracBits         getFracBitsArray( unsigned  ctxId  )  const { return m_Ctx[ctxId].getFracBitsArray(); }

private:
  inline void checkInit() { if( m_Ctx ) return; m_CtxBuffer.resize( ContextSetCfg::NumberOfContexts ); m_Ctx = m_CtxBuffer.data(); }
private:
  std::vector<BinProbModel> m_CtxBuffer;
  BinProbModel*             m_Ctx;
};



class Ctx;
class SubCtx
{
  friend class Ctx;
public:
  SubCtx( const CtxSet& ctxSet, const Ctx& ctx ) : m_CtxSet( ctxSet          ), m_Ctx( ctx          ) {}
  SubCtx( const SubCtx& subCtx )                 : m_CtxSet( subCtx.m_CtxSet ), m_Ctx( subCtx.m_Ctx ) {}
  const SubCtx& operator= ( const SubCtx& ) = delete;
private:
  const CtxSet  m_CtxSet;
  const Ctx&    m_Ctx;
};



class Ctx : public ContextSetCfg
{
public:
  Ctx();
  Ctx( const BinProbModel_Std*    dummy );
  Ctx( const Ctx&                 ctx   );

public:
  const Ctx& operator= ( const Ctx& ctx )
  {
    m_BPMType = ctx.m_BPMType;
    switch( m_BPMType )
    {
    case BPM_Std:   m_CtxStore_Std  .copyFrom( ctx.m_CtxStore_Std   );  break;
    default:        break;
    }
    ::memcpy( m_GRAdaptStats, ctx.m_GRAdaptStats, sizeof( unsigned ) * RExt__GOLOMB_RICE_ADAPTATION_STATISTICS_SETS );
    return *this;
  }

  SubCtx operator= ( SubCtx&& subCtx )
  {
    m_BPMType = subCtx.m_Ctx.m_BPMType;
    switch( m_BPMType )
    {
    case BPM_Std:   m_CtxStore_Std  .copyFrom( subCtx.m_Ctx.m_CtxStore_Std,   subCtx.m_CtxSet );  break;
    default:        break;
    }
    return std::move(subCtx);
  }

  void  init ( int qp, int initId )
  {
    switch( m_BPMType )
    {
    case BPM_Std:   m_CtxStore_Std  .init( qp, initId );  break;
    default:        break;
    }
    for( std::size_t k = 0; k < RExt__GOLOMB_RICE_ADAPTATION_STATISTICS_SETS; k++ )
    {
      m_GRAdaptStats[k] = 0;
    }
  }

  void  loadPStates( const std::vector<uint16_t>& probStates )
  {
    switch( m_BPMType )
    {
    case BPM_Std:   m_CtxStore_Std  .loadPStates( probStates );  break;
    default:        break;
    }
  }

  void  savePStates( std::vector<uint16_t>& probStates ) const
  {
    switch( m_BPMType )
    {
    case BPM_Std:   m_CtxStore_Std  .savePStates( probStates );  break;
    default:        break;
    }
  }

  void  initCtxAndWinSize( unsigned ctxId, const Ctx& ctx, const uint8_t winSize )
  {
    switch( m_BPMType )
    {
    case BPM_Std:
      m_CtxStore_Std  [ctxId] = ctx.m_CtxStore_Std  [ctxId];
      m_CtxStore_Std  [ctxId] . setLog2WindowSize   (winSize);
      break;
    default:
      break;
    }
  }

  const unsigned&     getGRAdaptStats ( unsigned      id )      const { return m_GRAdaptStats[id]; }
  unsigned&           getGRAdaptStats ( unsigned      id )            { return m_GRAdaptStats[id]; }

public:
  unsigned            getBPMType      ()                        const { return m_BPMType; }
  const Ctx&          getCtx          ()                        const { return *this; }
  Ctx&                getCtx          ()                              { return *this; }

  explicit operator   const CtxStore<BinProbModel_Std>  &()     const { return m_CtxStore_Std; }
  explicit operator         CtxStore<BinProbModel_Std>  &()           { return m_CtxStore_Std; }

  const FracBitsAccess&   getFracBitsAcess()  const
  {
    switch( m_BPMType )
    {
    case BPM_Std:   return m_CtxStore_Std;
    default:        THROW("BPMType out of range");
    }
  }

private:
  BPMType                       m_BPMType;
  CtxStore<BinProbModel_Std>    m_CtxStore_Std;
protected:
  unsigned                      m_GRAdaptStats[RExt__GOLOMB_RICE_ADAPTATION_STATISTICS_SETS];
#if ENABLE_SPLIT_PARALLELISM

public:
  int64_t cacheId;
  bool    cacheUsed;
#endif
};



typedef dynamic_cache<Ctx> CtxCache;

class TempCtx
{
  TempCtx( const TempCtx& ) = delete;
  const TempCtx& operator=( const TempCtx& ) = delete;
public:
  TempCtx ( CtxCache* cache )                     : m_ctx( *cache->get() ), m_cache( cache ) {}
  TempCtx ( CtxCache* cache, const Ctx& ctx    )  : m_ctx( *cache->get() ), m_cache( cache ) { m_ctx = ctx; }
  TempCtx ( CtxCache* cache, SubCtx&&   subCtx )  : m_ctx( *cache->get() ), m_cache( cache ) { m_ctx = std::forward<SubCtx>(subCtx); }
  ~TempCtx()                                      { m_cache->cache( &m_ctx ); }
  const Ctx& operator=( const Ctx& ctx )          { return ( m_ctx = ctx ); }
  SubCtx     operator=( SubCtx&&   subCtx )       { return m_ctx = std::forward<SubCtx>( subCtx ); }
  operator const Ctx& ()           const          { return m_ctx; }
  operator       Ctx& ()                          { return m_ctx; }
private:
  Ctx&      m_ctx;
  CtxCache* m_cache;
};



#endif
