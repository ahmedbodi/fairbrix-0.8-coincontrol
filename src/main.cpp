// Copyright (c) 2009-2010 Satoshi Nakamoto
// Copyright (c) 2009-2012 The Bitcoin developers
// Distributed under the MIT/X11 software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#include "alert.h"
#include "checkpoints.h"
#include "db.h"
#include "txdb.h"
#include "net.h"
#include "init.h"
#include "ui_interface.h"
#include "checkqueue.h"
#include <boost/algorithm/string/replace.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

using namespace std;
using namespace boost;

//
// Global state
//

CCriticalSection cs_setpwalletRegistered;
set<CWallet*> setpwalletRegistered;

CCriticalSection cs_main;

CTxMemPool mempool;
unsigned int nTransactionsUpdated = 0;

map<uint256, CBlockIndex*> mapBlockIndex;
// FBX
//uint256 hashGenesisBlock("0x12a765e31ffd4059bada1e25190f6e98c99d9714d334efa41a195a7e7e04bfe2");
//static CBigNum bnProofOfWorkLimit(~uint256(0) >> 20); // Litecoin: starting difficulty is 1 / 2^12
uint256 hashGenesisBlock("0x000000000019d6689c085ae165831e934ff763ae46a2a6c172b3f1b60a8ce26f");
static CBigNum bnProofOfWorkLimit(~uint256(0) >> 20); // fbx v0.3.x testnet default: (~uint256(0) >> 20)
CBlockIndex* pindexGenesisBlock = NULL;
int nBestHeight = -1;
uint256 nBestChainWork = 0;
uint256 nBestInvalidWork = 0;
uint256 hashBestChain = 0;
CBlockIndex* pindexBest = NULL;
set<CBlockIndex*, CBlockIndexWorkComparator> setBlockIndexValid; // may contain all CBlockIndex*'s that have validness >=BLOCK_VALID_TRANSACTIONS, and must contain those who aren't failed
int64 nTimeBestReceived = 0;
int nScriptCheckThreads = 0;
bool fImporting = false;
bool fReindex = false;
bool fBenchmark = false;
bool fTxIndex = false;
unsigned int nCoinCacheSize = 5000;

// FBX fees -- Base fee is either nMinTxFee or nMinRelayTxFee
// Minimum fee (when fees apply) 0.01. This is less than litecoins so the fee isn't
// quite so large when you send an output that is nearly a cent.
// However, there is code elsewhere that will increase this fee for very small outputs.
//
///** Fees smaller than this (in satoshi) are considered zero fee (for transaction creation) */
//int64 CTransaction::nMinTxFee = 2000000;
///** Fees smaller than this (in satoshi) are considered zero fee (for relaying) */
//int64 CTransaction::nMinRelayTxFee = 2000000;
int64 CTransaction::nMinTxFee = 1000000;
int64 CTransaction::nMinRelayTxFee = 1000000;

CMedianFilter<int> cPeerBlockCounts(8, 0); // Amount of blocks that other nodes claim to have

map<uint256, CBlock*> mapOrphanBlocks;
multimap<uint256, CBlock*> mapOrphanBlocksByPrev;

map<uint256, CTransaction> mapOrphanTransactions;
map<uint256, set<uint256> > mapOrphanTransactionsByPrev;

// Constant stuff for coinbase transactions we create:
CScript COINBASE_FLAGS;

const string strMessageMagic = "Fairbrix Signed Message:\n";

double dHashesPerSec = 0.0;
int64 nHPSTimerStart = 0;

// Settings
int64 nTransactionFee = 0;
int64 nMinimumInputValue = DUST_HARD_LIMIT;


//////////////////////////////////////////////////////////////////////////////
//
// dispatching functions
//

// These functions dispatch to one or all registered wallets


void RegisterWallet(CWallet* pwalletIn)
{
    {
        LOCK(cs_setpwalletRegistered);
        setpwalletRegistered.insert(pwalletIn);
    }
}

void UnregisterWallet(CWallet* pwalletIn)
{
    {
        LOCK(cs_setpwalletRegistered);
        setpwalletRegistered.erase(pwalletIn);
    }
}

// get the wallet transaction with the given hash (if it exists)
bool static GetTransaction(const uint256& hashTx, CWalletTx& wtx)
{
    BOOST_FOREACH(CWallet* pwallet, setpwalletRegistered)
        if (pwallet->GetTransaction(hashTx,wtx))
            return true;
    return false;
}

// erases transaction with the given hash from all wallets
void static EraseFromWallets(uint256 hash)
{
    BOOST_FOREACH(CWallet* pwallet, setpwalletRegistered)
        pwallet->EraseFromWallet(hash);
}

// make sure all wallets know about the given transaction, in the given block
void SyncWithWallets(const uint256 &hash, const CTransaction& tx, const CBlock* pblock, bool fUpdate)
{
    BOOST_FOREACH(CWallet* pwallet, setpwalletRegistered)
        pwallet->AddToWalletIfInvolvingMe(hash, tx, pblock, fUpdate);
}

// notify wallets about a new best chain
void static SetBestChain(const CBlockLocator& loc)
{
    BOOST_FOREACH(CWallet* pwallet, setpwalletRegistered)
        pwallet->SetBestChain(loc);
}

// notify wallets about an updated transaction
void static UpdatedTransaction(const uint256& hashTx)
{
    BOOST_FOREACH(CWallet* pwallet, setpwalletRegistered)
        pwallet->UpdatedTransaction(hashTx);
}

// dump all wallets
void static PrintWallets(const CBlock& block)
{
    BOOST_FOREACH(CWallet* pwallet, setpwalletRegistered)
        pwallet->PrintWallet(block);
}

// notify wallets about an incoming inventory (for request counts)
void static Inventory(const uint256& hash)
{
    BOOST_FOREACH(CWallet* pwallet, setpwalletRegistered)
        pwallet->Inventory(hash);
}

// ask wallets to resend their transactions
void static ResendWalletTransactions()
{
    BOOST_FOREACH(CWallet* pwallet, setpwalletRegistered)
        pwallet->ResendWalletTransactions();
}







//////////////////////////////////////////////////////////////////////////////
//
// CCoinsView implementations
//

bool CCoinsView::GetCoins(const uint256 &txid, CCoins &coins) { return false; }
bool CCoinsView::SetCoins(const uint256 &txid, const CCoins &coins) { return false; }
bool CCoinsView::HaveCoins(const uint256 &txid) { return false; }
CBlockIndex *CCoinsView::GetBestBlock() { return NULL; }
bool CCoinsView::SetBestBlock(CBlockIndex *pindex) { return false; }
bool CCoinsView::BatchWrite(const std::map<uint256, CCoins> &mapCoins, CBlockIndex *pindex) { return false; }
bool CCoinsView::GetStats(CCoinsStats &stats) { return false; }


CCoinsViewBacked::CCoinsViewBacked(CCoinsView &viewIn) : base(&viewIn) { }
bool CCoinsViewBacked::GetCoins(const uint256 &txid, CCoins &coins) { return base->GetCoins(txid, coins); }
bool CCoinsViewBacked::SetCoins(const uint256 &txid, const CCoins &coins) { return base->SetCoins(txid, coins); }
bool CCoinsViewBacked::HaveCoins(const uint256 &txid) { return base->HaveCoins(txid); }
CBlockIndex *CCoinsViewBacked::GetBestBlock() { return base->GetBestBlock(); }
bool CCoinsViewBacked::SetBestBlock(CBlockIndex *pindex) { return base->SetBestBlock(pindex); }
void CCoinsViewBacked::SetBackend(CCoinsView &viewIn) { base = &viewIn; }
bool CCoinsViewBacked::BatchWrite(const std::map<uint256, CCoins> &mapCoins, CBlockIndex *pindex) { return base->BatchWrite(mapCoins, pindex); }
bool CCoinsViewBacked::GetStats(CCoinsStats &stats) { return base->GetStats(stats); }

CCoinsViewCache::CCoinsViewCache(CCoinsView &baseIn, bool fDummy) : CCoinsViewBacked(baseIn), pindexTip(NULL) { }

bool CCoinsViewCache::GetCoins(const uint256 &txid, CCoins &coins) {
    if (cacheCoins.count(txid)) {
        coins = cacheCoins[txid];
        return true;
    }
    if (base->GetCoins(txid, coins)) {
        cacheCoins[txid] = coins;
        return true;
    }
    return false;
}

std::map<uint256,CCoins>::iterator CCoinsViewCache::FetchCoins(const uint256 &txid) {
    std::map<uint256,CCoins>::iterator it = cacheCoins.lower_bound(txid);
    if (it != cacheCoins.end() && it->first == txid)
        return it;
    CCoins tmp;
    if (!base->GetCoins(txid,tmp))
        return cacheCoins.end();
    std::map<uint256,CCoins>::iterator ret = cacheCoins.insert(it, std::make_pair(txid, CCoins()));
    tmp.swap(ret->second);
    return ret;
}

CCoins &CCoinsViewCache::GetCoins(const uint256 &txid) {
    std::map<uint256,CCoins>::iterator it = FetchCoins(txid);
    assert(it != cacheCoins.end());
    return it->second;
}

bool CCoinsViewCache::SetCoins(const uint256 &txid, const CCoins &coins) {
    cacheCoins[txid] = coins;
    return true;
}

bool CCoinsViewCache::HaveCoins(const uint256 &txid) {
    return FetchCoins(txid) != cacheCoins.end();
}

CBlockIndex *CCoinsViewCache::GetBestBlock() {
    if (pindexTip == NULL)
        pindexTip = base->GetBestBlock();
    return pindexTip;
}

bool CCoinsViewCache::SetBestBlock(CBlockIndex *pindex) {
    pindexTip = pindex;
    return true;
}

bool CCoinsViewCache::BatchWrite(const std::map<uint256, CCoins> &mapCoins, CBlockIndex *pindex) {
    for (std::map<uint256, CCoins>::const_iterator it = mapCoins.begin(); it != mapCoins.end(); it++)
        cacheCoins[it->first] = it->second;
    pindexTip = pindex;
    return true;
}

bool CCoinsViewCache::Flush() {
    bool fOk = base->BatchWrite(cacheCoins, pindexTip);
    if (fOk)
        cacheCoins.clear();
    return fOk;
}

unsigned int CCoinsViewCache::GetCacheSize() {
    return cacheCoins.size();
}

/** CCoinsView that brings transactions from a memorypool into view.
    It does not check for spendings by memory pool transactions. */
CCoinsViewMemPool::CCoinsViewMemPool(CCoinsView &baseIn, CTxMemPool &mempoolIn) : CCoinsViewBacked(baseIn), mempool(mempoolIn) { }

bool CCoinsViewMemPool::GetCoins(const uint256 &txid, CCoins &coins) {
    if (base->GetCoins(txid, coins))
        return true;
    if (mempool.exists(txid)) {
        const CTransaction &tx = mempool.lookup(txid);
        coins = CCoins(tx, MEMPOOL_HEIGHT);
        return true;
    }
    return false;
}

bool CCoinsViewMemPool::HaveCoins(const uint256 &txid) {
    return mempool.exists(txid) || base->HaveCoins(txid);
}

CCoinsViewCache *pcoinsTip = NULL;
CBlockTreeDB *pblocktree = NULL;

//////////////////////////////////////////////////////////////////////////////
//
// mapOrphanTransactions
//

bool AddOrphanTx(const CTransaction& tx)
{
    uint256 hash = tx.GetHash();
    if (mapOrphanTransactions.count(hash))
        return false;

    // Ignore big transactions, to avoid a
    // send-big-orphans memory exhaustion attack. If a peer has a legitimate
    // large transaction with a missing parent then we assume
    // it will rebroadcast it later, after the parent transaction(s)
    // have been mined or received.
    // 10,000 orphans, each of which is at most 5,000 bytes big is
    // at most 500 megabytes of orphans:
    unsigned int sz = tx.GetSerializeSize(SER_NETWORK, CTransaction::CURRENT_VERSION);
    if (sz > 5000)
    {
        printf("ignoring large orphan tx (size: %u, hash: %s)\n", sz, hash.ToString().c_str());
        return false;
    }

    mapOrphanTransactions[hash] = tx;
    BOOST_FOREACH(const CTxIn& txin, tx.vin)
        mapOrphanTransactionsByPrev[txin.prevout.hash].insert(hash);

    printf("stored orphan tx %s (mapsz %"PRIszu")\n", hash.ToString().c_str(),
        mapOrphanTransactions.size());
    return true;
}

void static EraseOrphanTx(uint256 hash)
{
    if (!mapOrphanTransactions.count(hash))
        return;
    const CTransaction& tx = mapOrphanTransactions[hash];
    BOOST_FOREACH(const CTxIn& txin, tx.vin)
    {
        mapOrphanTransactionsByPrev[txin.prevout.hash].erase(hash);
        if (mapOrphanTransactionsByPrev[txin.prevout.hash].empty())
            mapOrphanTransactionsByPrev.erase(txin.prevout.hash);
    }
    mapOrphanTransactions.erase(hash);
}

unsigned int LimitOrphanTxSize(unsigned int nMaxOrphans)
{
    unsigned int nEvicted = 0;
    while (mapOrphanTransactions.size() > nMaxOrphans)
    {
        // Evict a random orphan:
        uint256 randomhash = GetRandHash();
        map<uint256, CTransaction>::iterator it = mapOrphanTransactions.lower_bound(randomhash);
        if (it == mapOrphanTransactions.end())
            it = mapOrphanTransactions.begin();
        EraseOrphanTx(it->first);
        ++nEvicted;
    }
    return nEvicted;
}







//////////////////////////////////////////////////////////////////////////////
//
// CTransaction / CTxOut
//

bool CTxOut::IsDust() const
{
    // Litecoin: IsDust() detection disabled, allows any valid dust to be relayed.
    // The fees imposed on each dust txo is considered sufficient spam deterrant. 
    return false;
}

bool CTransaction::IsStandard(string& strReason) const
{
    if (nVersion > CTransaction::CURRENT_VERSION || nVersion < 1) {
        strReason = "version";
        return false;
    }

    if (!IsFinal()) {
        strReason = "not-final";
        return false;
    }

    // Extremely large transactions with lots of inputs can cost the network
    // almost as much to process as they cost the sender in fees, because
    // computing signature hashes is O(ninputs*txsize). Limiting transactions
    // to MAX_STANDARD_TX_SIZE mitigates CPU exhaustion attacks.
    unsigned int sz = this->GetSerializeSize(SER_NETWORK, CTransaction::CURRENT_VERSION);
    if (sz >= MAX_STANDARD_TX_SIZE) {
        strReason = "tx-size";
        return false;
    }

    BOOST_FOREACH(const CTxIn& txin, vin)
    {
        // Biggest 'standard' txin is a 3-signature 3-of-3 CHECKMULTISIG
        // pay-to-script-hash, which is 3 ~80-byte signatures, 3
        // ~65-byte public keys, plus a few script ops.
        if (txin.scriptSig.size() > 500) {
            strReason = "scriptsig-size";
            return false;
        }
        if (!txin.scriptSig.IsPushOnly()) {
            strReason = "scriptsig-not-pushonly";
            return false;
        }
    }
    BOOST_FOREACH(const CTxOut& txout, vout) {
        if (!::IsStandard(txout.scriptPubKey)) {
            strReason = "scriptpubkey";
            return false;
        }
        if (txout.IsDust()) {
            strReason = "dust";
            return false;
        }
    }
    return true;
}

//
// Check transaction inputs, and make sure any
// pay-to-script-hash transactions are evaluating IsStandard scripts
//
// Why bother? To avoid denial-of-service attacks; an attacker
// can submit a standard HASH... OP_EQUAL transaction,
// which will get accepted into blocks. The redemption
// script can be anything; an attacker could use a very
// expensive-to-check-upon-redemption script like:
//   DUP CHECKSIG DROP ... repeated 100 times... OP_1
//
bool CTransaction::AreInputsStandard(CCoinsViewCache& mapInputs) const
{
    if (IsCoinBase())
        return true; // Coinbases don't use vin normally

    for (unsigned int i = 0; i < vin.size(); i++)
    {
        const CTxOut& prev = GetOutputFor(vin[i], mapInputs);

        vector<vector<unsigned char> > vSolutions;
        txnouttype whichType;
        // get the scriptPubKey corresponding to this input:
        const CScript& prevScript = prev.scriptPubKey;
        if (!Solver(prevScript, whichType, vSolutions))
            return false;
        int nArgsExpected = ScriptSigArgsExpected(whichType, vSolutions);
        if (nArgsExpected < 0)
            return false;

        // Transactions with extra stuff in their scriptSigs are
        // non-standard. Note that this EvalScript() call will
        // be quick, because if there are any operations
        // beside "push data" in the scriptSig the
        // IsStandard() call returns false
        vector<vector<unsigned char> > stack;
        if (!EvalScript(stack, vin[i].scriptSig, *this, i, false, 0))
            return false;

        if (whichType == TX_SCRIPTHASH)
        {
            if (stack.empty())
                return false;
            CScript subscript(stack.back().begin(), stack.back().end());
            vector<vector<unsigned char> > vSolutions2;
            txnouttype whichType2;
            if (!Solver(subscript, whichType2, vSolutions2))
                return false;
            if (whichType2 == TX_SCRIPTHASH)
                return false;

            int tmpExpected;
            tmpExpected = ScriptSigArgsExpected(whichType2, vSolutions2);
            if (tmpExpected < 0)
                return false;
            nArgsExpected += tmpExpected;
        }

        if (stack.size() != (unsigned int)nArgsExpected)
            return false;
    }

    return true;
}

unsigned int CTransaction::GetLegacySigOpCount() const
{
    unsigned int nSigOps = 0;
    BOOST_FOREACH(const CTxIn& txin, vin)
    {
        nSigOps += txin.scriptSig.GetSigOpCount(false);
    }
    BOOST_FOREACH(const CTxOut& txout, vout)
    {
        nSigOps += txout.scriptPubKey.GetSigOpCount(false);
    }
    return nSigOps;
}


int CMerkleTx::SetMerkleBranch(const CBlock* pblock)
{
    CBlock blockTmp;

    if (pblock == NULL) {
        CCoins coins;
        if (pcoinsTip->GetCoins(GetHash(), coins)) {
            CBlockIndex *pindex = FindBlockByHeight(coins.nHeight);
            if (pindex) {
                if (!blockTmp.ReadFromDisk(pindex))
                    return 0;
                pblock = &blockTmp;
            }
        }
    }

    if (pblock) {
        // Update the tx's hashBlock
        hashBlock = pblock->GetHash();

        // Locate the transaction
        for (nIndex = 0; nIndex < (int)pblock->vtx.size(); nIndex++)
            if (pblock->vtx[nIndex] == *(CTransaction*)this)
                break;
        if (nIndex == (int)pblock->vtx.size())
        {
            vMerkleBranch.clear();
            nIndex = -1;
            printf("ERROR: SetMerkleBranch() : couldn't find tx in block\n");
            return 0;
        }

        // Fill in merkle branch
        vMerkleBranch = pblock->GetMerkleBranch(nIndex);
    }

    // Is the tx in a block that's in the main chain
    map<uint256, CBlockIndex*>::iterator mi = mapBlockIndex.find(hashBlock);
    if (mi == mapBlockIndex.end())
        return 0;
    CBlockIndex* pindex = (*mi).second;
    if (!pindex || !pindex->IsInMainChain())
        return 0;

    return pindexBest->nHeight - pindex->nHeight + 1;
}







bool CTransaction::CheckTransaction(CValidationState &state) const
{
    // Basic checks that don't depend on any context
    if (vin.empty())
        return state.DoS(10, error("CTransaction::CheckTransaction() : vin empty"));
    if (vout.empty())
        return state.DoS(10, error("CTransaction::CheckTransaction() : vout empty"));
    // Size limits
    if (::GetSerializeSize(*this, SER_NETWORK, PROTOCOL_VERSION) > MAX_BLOCK_SIZE)
        return state.DoS(100, error("CTransaction::CheckTransaction() : size limits failed"));

    // Check for negative or overflow output values
    int64 nValueOut = 0;
    BOOST_FOREACH(const CTxOut& txout, vout)
    {
        if (txout.nValue < 0)
            return state.DoS(100, error("CTransaction::CheckTransaction() : txout.nValue negative"));
        if (txout.nValue > MAX_MONEY)
            return state.DoS(100, error("CTransaction::CheckTransaction() : txout.nValue too high"));
        nValueOut += txout.nValue;
        if (!MoneyRange(nValueOut))
            return state.DoS(100, error("CTransaction::CheckTransaction() : txout total out of range"));
    }

    // Check for duplicate inputs
    set<COutPoint> vInOutPoints;
    BOOST_FOREACH(const CTxIn& txin, vin)
    {
        if (vInOutPoints.count(txin.prevout))
            return state.DoS(100, error("CTransaction::CheckTransaction() : duplicate inputs"));
        vInOutPoints.insert(txin.prevout);
    }

    if (IsCoinBase())
    {
        if (vin[0].scriptSig.size() < 2 || vin[0].scriptSig.size() > 100)
            return state.DoS(100, error("CTransaction::CheckTransaction() : coinbase script size"));
    }
    else
    {
        BOOST_FOREACH(const CTxIn& txin, vin)
            if (txin.prevout.IsNull())
                return state.DoS(10, error("CTransaction::CheckTransaction() : prevout is null"));
    }

    return true;
}

int64 CTransaction::GetMinFee(unsigned int nBlockSize, bool fAllowFree,
                              enum GetMinFee_mode mode) const
{
    // Base fee is either nMinTxFee or nMinRelayTxFee
    int64 nBaseFee = (mode == GMF_RELAY) ? nMinRelayTxFee : nMinTxFee;

    unsigned int nBytes = ::GetSerializeSize(*this, SER_NETWORK, PROTOCOL_VERSION);
    unsigned int nNewBlockSize = nBlockSize + nBytes;

// FBX fees
    int smallTxOutCount = 0;
// This makes large sized transactions cost more than before.
//    int64 nMinFee = (1 + (int64)nBytes / 1000) * nBaseFee;
    int64 nMinFee = (1 + (int64)nBytes / 500) * nBaseFee;

    if (fAllowFree)
    {
        if (nBlockSize == 1)
        {
            // Transactions under 10K are free
            // (about 4500 BTC if made of 50 BTC inputs)
            if (nBytes < 10000)
                nMinFee = 0;
        }
        else
        {
            // Free transaction area
// FBX fees
// since blocks are faster than in bitcoin, reserve less space for free transactions.
// (see also 'nBlockPrioritySize')
//            if (nNewBlockSize < 27000)
            if (nNewBlockSize < 12000)
                nMinFee = 0;
        }
    }

    // Litecoin
    // To limit dust spam, add nBaseFee for each output less than DUST_SOFT_LIMIT
// FBX fees
//    BOOST_FOREACH(const CTxOut& txout, vout)
//        if (txout.nValue < DUST_SOFT_LIMIT)
//            nMinFee += nBaseFee;
    BOOST_FOREACH(const CTxOut& txout, vout)  {
        if (txout.nValue < CENT/100) { // outputs smaller than 0.0001
            nMinFee += nBaseFee * 100;  // fee of 1
            smallTxOutCount++;
        }
    else if ((txout.nValue < CENT)) {
            nMinFee += nBaseFee;
            smallTxOutCount++;
        }
    }

    // Raise the price as the block approaches full
    if (nBlockSize != 1 && nNewBlockSize >= MAX_BLOCK_SIZE_GEN/2)
    {
        if (nNewBlockSize >= MAX_BLOCK_SIZE_GEN)
            return MAX_MONEY;
        nMinFee *= MAX_BLOCK_SIZE_GEN / (MAX_BLOCK_SIZE_GEN - nNewBlockSize);
    }

    if (!MoneyRange(nMinFee))
        nMinFee = MAX_MONEY;

// FBX fees
    // This is the core change to limit dust spam.  Instead of a flat fee for small outputs charge a fee
    // for each small output.  If there are more than 15 small outputs than don't allow the transaction at all.
            if(smallTxOutCount > 15)
                nMinFee = MAX_MONEY;

    return nMinFee;
}

void CTxMemPool::pruneSpent(const uint256 &hashTx, CCoins &coins)
{
    LOCK(cs);

    std::map<COutPoint, CInPoint>::iterator it = mapNextTx.lower_bound(COutPoint(hashTx, 0));

    // iterate over all COutPoints in mapNextTx whose hash equals the provided hashTx
    while (it != mapNextTx.end() && it->first.hash == hashTx) {
        coins.Spend(it->first.n); // and remove those outputs from coins
        it++;
    }
}

bool CTxMemPool::accept(CValidationState &state, CTransaction &tx, bool fCheckInputs, bool fLimitFree,
                        bool* pfMissingInputs)
{
    if (pfMissingInputs)
        *pfMissingInputs = false;

    if (!tx.CheckTransaction(state))
        return error("CTxMemPool::accept() : CheckTransaction failed");

    // Coinbase is only valid in a block, not as a loose transaction
    if (tx.IsCoinBase())
        return state.DoS(100, error("CTxMemPool::accept() : coinbase as individual tx"));

    // To help v0.1.5 clients who would see it as a negative number
    if ((int64)tx.nLockTime > std::numeric_limits<int>::max())
        return error("CTxMemPool::accept() : not accepting nLockTime beyond 2038 yet");

    // Rather not work on nonstandard transactions (unless -testnet)
    string strNonStd;
    if (!fTestNet && !tx.IsStandard(strNonStd))
        return error("CTxMemPool::accept() : nonstandard transaction (%s)",
                     strNonStd.c_str());

    // is it already in the memory pool?
    uint256 hash = tx.GetHash();
    {
        LOCK(cs);
        if (mapTx.count(hash))
            return false;
    }

    // Check for conflicts with in-memory transactions
    CTransaction* ptxOld = NULL;
    for (unsigned int i = 0; i < tx.vin.size(); i++)
    {
        COutPoint outpoint = tx.vin[i].prevout;
        if (mapNextTx.count(outpoint))
        {
            // Disable replacement feature for now
            return false;

            // Allow replacing with a newer version of the same transaction
            if (i != 0)
                return false;
            ptxOld = mapNextTx[outpoint].ptx;
            if (ptxOld->IsFinal())
                return false;
            if (!tx.IsNewerThan(*ptxOld))
                return false;
            for (unsigned int i = 0; i < tx.vin.size(); i++)
            {
                COutPoint outpoint = tx.vin[i].prevout;
                if (!mapNextTx.count(outpoint) || mapNextTx[outpoint].ptx != ptxOld)
                    return false;
            }
            break;
        }
    }

    if (fCheckInputs)
    {
        CCoinsView dummy;
        CCoinsViewCache view(dummy);

        {
        LOCK(cs);
        CCoinsViewMemPool viewMemPool(*pcoinsTip, *this);
        view.SetBackend(viewMemPool);

        // do we already have it?
        if (view.HaveCoins(hash))
            return false;

        // do all inputs exist?
        // Note that this does not check for the presence of actual outputs (see the next check for that),
        // only helps filling in pfMissingInputs (to determine missing vs spent).
        BOOST_FOREACH(const CTxIn txin, tx.vin) {
            if (!view.HaveCoins(txin.prevout.hash)) {
                if (pfMissingInputs)
                    *pfMissingInputs = true;
                return false;
            }
        }

        // are the actual inputs available?
        if (!tx.HaveInputs(view))
            return state.Invalid(error("CTxMemPool::accept() : inputs already spent"));

        // Bring the best block into scope
        view.GetBestBlock();

        // we have all inputs cached now, so switch back to dummy, so we don't need to keep lock on mempool
        view.SetBackend(dummy);
        }

        // Check for non-standard pay-to-script-hash in inputs
        if (!tx.AreInputsStandard(view) && !fTestNet)
            return error("CTxMemPool::accept() : nonstandard transaction input");

        // Note: if you modify this code to accept non-standard transactions, then
        // you should add code here to check that the transaction does a
        // reasonable number of ECDSA signature verifications.

        int64 nFees = tx.GetValueIn(view)-tx.GetValueOut();
        unsigned int nSize = ::GetSerializeSize(tx, SER_NETWORK, PROTOCOL_VERSION);

        // Don't accept it if it can't get into a block
        int64 txMinFee = tx.GetMinFee(1000, true, GMF_RELAY);
        if (fLimitFree && nFees < txMinFee)
            return error("CTxMemPool::accept() : not enough fees %s, %"PRI64d" < %"PRI64d,
                         hash.ToString().c_str(),
                         nFees, txMinFee);

        // Continuously rate-limit free transactions
        // This mitigates 'penny-flooding' -- sending thousands of free transactions just to
        // be annoying or make others' transactions take longer to confirm.
        if (fLimitFree && nFees < CTransaction::nMinRelayTxFee)
        {
// FBX fees
            static double dFreeRelay;
            static double dPartialRelay;
            static double dNewFreeCount;

            static double dFreeCount;
            static int64 nLastTime;
            int64 nNow = GetTime();

            LOCK(cs);

            // Use an exponentially decaying ~10-minute window:
            dFreeCount *= pow(1.0 - 1.0/600.0, (double)(nNow - nLastTime));
            nLastTime = nNow;
            // -limitfreerelay unit is thousand-bytes-per-minute

// FBX fees
// Why: another spam attack mitigation. Don't relay more than (on average) 5000 bytes
// of free transactions a minute.  This is roughly equivilant to 20 normal sized transactions per minute.
// At default rate it would take several months to fill 1GB
//
//            // At default rate it would take over a month to fill 1GB
//            if (dFreeCount >= GetArg("-limitfreerelay", 15)*10*1000)
//            return error("CTxMemPool::accept() : free transaction rejected by rate limiter");
            dFreeRelay = GetArg("-limitfreerelay", 5)*10*1000;
            dPartialRelay = dFreeRelay * 0.75;
            dNewFreeCount = dFreeCount + nSize;
            if( !( dNewFreeCount <= dFreeRelay
                || dFreeCount < dPartialRelay
// fixme                || IsFromMe(tx)
              )
            )
                return error("CTxMemPool::accept() : free transaction rejected by rate limiter");

            if (fDebug)
                printf("Rate limit dFreeCount: %g => %g\n", dFreeCount, dFreeCount+nSize);
            dFreeCount += nSize;
        }

        // Check against previous transactions
        // This is done last to help prevent CPU exhaustion denial-of-service attacks.
        if (!tx.CheckInputs(state, view, true, SCRIPT_VERIFY_P2SH | SCRIPT_VERIFY_STRICTENC))
        {
            return error("CTxMemPool::accept() : ConnectInputs failed %s", hash.ToString().c_str());
        }
    }

    // Store transaction in memory
    {
        LOCK(cs);
        if (ptxOld)
        {
            printf("CTxMemPool::accept() : replacing tx %s with new version\n", ptxOld->GetHash().ToString().c_str());
            remove(*ptxOld);
        }
        addUnchecked(hash, tx);
    }

    ///// are we sure this is ok when loading transactions or restoring block txes
    // If updated, erase old tx from wallet
    if (ptxOld)
        EraseFromWallets(ptxOld->GetHash());
    SyncWithWallets(hash, tx, NULL, true);

    printf("CTxMemPool::accept() : accepted %s (poolsz %"PRIszu")\n",
           hash.ToString().c_str(),
           mapTx.size());
    return true;
}

bool CTransaction::AcceptToMemoryPool(CValidationState &state, bool fCheckInputs, bool fLimitFree, bool* pfMissingInputs)
{
    try {
        return mempool.accept(state, *this, fCheckInputs, fLimitFree, pfMissingInputs);
    } catch(std::runtime_error &e) {
        return state.Abort(_("System error: ") + e.what());
    }
}

bool CTxMemPool::addUnchecked(const uint256& hash, const CTransaction &tx)
{
    // Add to memory pool without checking anything.  Don't call this directly,
    // call CTxMemPool::accept to properly check the transaction first.
    {
        mapTx[hash] = tx;
        for (unsigned int i = 0; i < tx.vin.size(); i++)
            mapNextTx[tx.vin[i].prevout] = CInPoint(&mapTx[hash], i);
        nTransactionsUpdated++;
    }
    return true;
}


bool CTxMemPool::remove(const CTransaction &tx, bool fRecursive)
{
    // Remove transaction from memory pool
    {
        LOCK(cs);
        uint256 hash = tx.GetHash();
        if (fRecursive) {
            for (unsigned int i = 0; i < tx.vout.size(); i++) {
                std::map<COutPoint, CInPoint>::iterator it = mapNextTx.find(COutPoint(hash, i));
                if (it != mapNextTx.end())
                    remove(*it->second.ptx, true);
            }
        }
        if (mapTx.count(hash))
        {
            BOOST_FOREACH(const CTxIn& txin, tx.vin)
                mapNextTx.erase(txin.prevout);
            mapTx.erase(hash);
            nTransactionsUpdated++;
        }
    }
    return true;
}

bool CTxMemPool::removeConflicts(const CTransaction &tx)
{
    // Remove transactions which depend on inputs of tx, recursively
    LOCK(cs);
    BOOST_FOREACH(const CTxIn &txin, tx.vin) {
        std::map<COutPoint, CInPoint>::iterator it = mapNextTx.find(txin.prevout);
        if (it != mapNextTx.end()) {
            const CTransaction &txConflict = *it->second.ptx;
            if (txConflict != tx)
                remove(txConflict, true);
        }
    }
    return true;
}

void CTxMemPool::clear()
{
    LOCK(cs);
    mapTx.clear();
    mapNextTx.clear();
    ++nTransactionsUpdated;
}

void CTxMemPool::queryHashes(std::vector<uint256>& vtxid)
{
    vtxid.clear();

    LOCK(cs);
    vtxid.reserve(mapTx.size());
    for (map<uint256, CTransaction>::iterator mi = mapTx.begin(); mi != mapTx.end(); ++mi)
        vtxid.push_back((*mi).first);
}




int CMerkleTx::GetDepthInMainChain(CBlockIndex* &pindexRet) const
{
    if (hashBlock == 0 || nIndex == -1)
        return 0;

    // Find the block it claims to be in
    map<uint256, CBlockIndex*>::iterator mi = mapBlockIndex.find(hashBlock);
    if (mi == mapBlockIndex.end())
        return 0;
    CBlockIndex* pindex = (*mi).second;
    if (!pindex || !pindex->IsInMainChain())
        return 0;

    // Make sure the merkle branch connects to this block
    if (!fMerkleVerified)
    {
        if (CBlock::CheckMerkleBranch(GetHash(), vMerkleBranch, nIndex) != pindex->hashMerkleRoot)
            return 0;
        fMerkleVerified = true;
    }

    pindexRet = pindex;
    return pindexBest->nHeight - pindex->nHeight + 1;
}


int CMerkleTx::GetBlocksToMaturity() const
{
    if (!IsCoinBase())
        return 0;
// FBX -- mined balance available after 151 confirmations
//    return max(0, (COINBASE_MATURITY+20) - GetDepthInMainChain());
    return max(0, (COINBASE_MATURITY+20) - GetDepthInMainChain());
}


bool CMerkleTx::AcceptToMemoryPool(bool fCheckInputs, bool fLimitFree)
{
    CValidationState state;
    return CTransaction::AcceptToMemoryPool(state, fCheckInputs, fLimitFree);
}



bool CWalletTx::AcceptWalletTransaction(bool fCheckInputs)
{
    {
        LOCK(mempool.cs);
        // Add previous supporting transactions first
        BOOST_FOREACH(CMerkleTx& tx, vtxPrev)
        {
            if (!tx.IsCoinBase())
            {
                uint256 hash = tx.GetHash();
                if (!mempool.exists(hash) && pcoinsTip->HaveCoins(hash))
                    tx.AcceptToMemoryPool(fCheckInputs, false);
            }
        }
        return AcceptToMemoryPool(fCheckInputs, false);
    }
    return false;
}


// Return transaction in tx, and if it was found inside a block, its hash is placed in hashBlock
bool GetTransaction(const uint256 &hash, CTransaction &txOut, uint256 &hashBlock, bool fAllowSlow)
{
    CBlockIndex *pindexSlow = NULL;
    {
        LOCK(cs_main);
        {
            LOCK(mempool.cs);
            if (mempool.exists(hash))
            {
                txOut = mempool.lookup(hash);
                return true;
            }
        }

        if (fTxIndex) {
            CDiskTxPos postx;
            if (pblocktree->ReadTxIndex(hash, postx)) {
                CAutoFile file(OpenBlockFile(postx, true), SER_DISK, CLIENT_VERSION);
                CBlockHeader header;
                try {
                    file >> header;
                    fseek(file, postx.nTxOffset, SEEK_CUR);
                    file >> txOut;
                } catch (std::exception &e) {
                    return error("%s() : deserialize or I/O error", __PRETTY_FUNCTION__);
                }
                hashBlock = header.GetHash();
                if (txOut.GetHash() != hash)
                    return error("%s() : txid mismatch", __PRETTY_FUNCTION__);
                return true;
            }
        }

        if (fAllowSlow) { // use coin database to locate block that contains transaction, and scan it
            int nHeight = -1;
            {
                CCoinsViewCache &view = *pcoinsTip;
                CCoins coins;
                if (view.GetCoins(hash, coins))
                    nHeight = coins.nHeight;
            }
            if (nHeight > 0)
                pindexSlow = FindBlockByHeight(nHeight);
        }
    }

    if (pindexSlow) {
        CBlock block;
        if (block.ReadFromDisk(pindexSlow)) {
            BOOST_FOREACH(const CTransaction &tx, block.vtx) {
                if (tx.GetHash() == hash) {
                    txOut = tx;
                    hashBlock = pindexSlow->GetBlockHash();
                    return true;
                }
            }
        }
    }

    return false;
}






//////////////////////////////////////////////////////////////////////////////
//
// CBlock and CBlockIndex
//

static CBlockIndex* pblockindexFBBHLast;
CBlockIndex* FindBlockByHeight(int nHeight)
{
    CBlockIndex *pblockindex;
    if (nHeight < nBestHeight / 2)
        pblockindex = pindexGenesisBlock;
    else
        pblockindex = pindexBest;
    if (pblockindexFBBHLast && abs(nHeight - pblockindex->nHeight) > abs(nHeight - pblockindexFBBHLast->nHeight))
        pblockindex = pblockindexFBBHLast;
    while (pblockindex->nHeight > nHeight)
        pblockindex = pblockindex->pprev;
    while (pblockindex->nHeight < nHeight)
        pblockindex = pblockindex->pnext;
    pblockindexFBBHLast = pblockindex;
    return pblockindex;
}

bool CBlock::ReadFromDisk(const CBlockIndex* pindex)
{
    if (!ReadFromDisk(pindex->GetBlockPos()))
        return false;
    if (GetHash() != pindex->GetBlockHash())
        return error("CBlock::ReadFromDisk() : GetHash() doesn't match index");
    return true;
}

uint256 static GetOrphanRoot(const CBlockHeader* pblock)
{
    // Work back to the first block in the orphan chain
    while (mapOrphanBlocks.count(pblock->hashPrevBlock))
        pblock = mapOrphanBlocks[pblock->hashPrevBlock];
    return pblock->GetHash();
}

int64 static GetBlockValue(int nHeight, int64 nFees)
{
// FBX
//    int64 nSubsidy = 50 * COIN;
//
//    // Subsidy is cut in half every 840000 blocks, which will occur approximately every 4 years
//    nSubsidy >>= (nHeight / 840000); // Litecoin: 840k blocks in ~4 years
    int64 nSubsidy = 25 * COIN;

    return nSubsidy + nFees;
}

// FBX
//static const int64 nTargetTimespan = 3.5 * 24 * 60 * 60; // Litecoin: 3.5 days
//static const int64 nTargetSpacing = 2.5 * 60; // Litecoin: 2.5 minutes
static const int64 nTargetTimespan = 7 * 24 * 60 * 60; // one week
static const int64 nTargetSpacing = 5 * 60;
static const int64 nInterval = nTargetTimespan / nTargetSpacing;

//
// minimum amount of work that could possibly be required nTime after
// minimum work required was nBase
//
unsigned int ComputeMinWork(unsigned int nBase, int64 nTime)
{
    // Testnet has min-difficulty blocks
    // after nTargetSpacing*2 time between blocks:
    if (fTestNet && nTime > nTargetSpacing*2)
        return bnProofOfWorkLimit.GetCompact();

    CBigNum bnResult;
    bnResult.SetCompact(nBase);
    while (nTime > 0 && bnResult < bnProofOfWorkLimit)
    {
        // Maximum 400% adjustment...
        bnResult *= 4;
        // ... in best-case exactly 4-times-normal target time
        nTime -= nTargetTimespan*4;
    }
    if (bnResult > bnProofOfWorkLimit)
        bnResult = bnProofOfWorkLimit;
    return bnResult.GetCompact();
}

unsigned int static GetNextWorkRequired(const CBlockIndex* pindexLast, const CBlockHeader *pblock)
{
    unsigned int nProofOfWorkLimit = bnProofOfWorkLimit.GetCompact();

    // Genesis block
    if (pindexLast == NULL)
        return nProofOfWorkLimit;

    // Only change once per interval
    if ((pindexLast->nHeight+1) % nInterval != 0)
    {
        // Special difficulty rule for testnet:
        if (fTestNet)
        {
            // If the new block's timestamp is more than 2* 10 minutes
            // then allow mining of a min-difficulty block.
            if (pblock->nTime > pindexLast->nTime + nTargetSpacing*2)
                return nProofOfWorkLimit;
            else
            {
                // Return the last non-special-min-difficulty-rules-block
                const CBlockIndex* pindex = pindexLast;
                while (pindex->pprev && pindex->nHeight % nInterval != 0 && pindex->nBits == nProofOfWorkLimit)
                    pindex = pindex->pprev;
                return pindex->nBits;
            }
        }

        return pindexLast->nBits;
    }

    // Litecoin: This fixes an issue where a 51% attack can change difficulty at will.
    // Go back the full period unless it's the first retarget after genesis. Code courtesy of Art Forz
    int blockstogoback = nInterval-1;
    if ((pindexLast->nHeight+1) != nInterval)
        blockstogoback = nInterval;

    // Go back by what we want to be 14 days worth of blocks
    const CBlockIndex* pindexFirst = pindexLast;
    for (int i = 0; pindexFirst && i < blockstogoback; i++)
        pindexFirst = pindexFirst->pprev;
    assert(pindexFirst);

    // Limit adjustment step
    int64 nActualTimespan = pindexLast->GetBlockTime() - pindexFirst->GetBlockTime();
    printf("  nActualTimespan = %"PRI64d"  before bounds\n", nActualTimespan);
    if (nActualTimespan < nTargetTimespan/4)
        nActualTimespan = nTargetTimespan/4;
    if (nActualTimespan > nTargetTimespan*4)
        nActualTimespan = nTargetTimespan*4;

    // Retarget
    CBigNum bnNew;
    bnNew.SetCompact(pindexLast->nBits);
    bnNew *= nActualTimespan;
    bnNew /= nTargetTimespan;

    if (bnNew > bnProofOfWorkLimit)
        bnNew = bnProofOfWorkLimit;

    /// debug print
    printf("GetNextWorkRequired RETARGET\n");
    printf("nTargetTimespan = %"PRI64d"    nActualTimespan = %"PRI64d"\n", nTargetTimespan, nActualTimespan);
    printf("Before: %08x  %s\n", pindexLast->nBits, CBigNum().SetCompact(pindexLast->nBits).getuint256().ToString().c_str());
    printf("After:  %08x  %s\n", bnNew.GetCompact(), bnNew.getuint256().ToString().c_str());

    return bnNew.GetCompact();
}

bool CheckProofOfWork(uint256 hash, unsigned int nBits)
{
    CBigNum bnTarget;
    bnTarget.SetCompact(nBits);

    // Check range
    if (bnTarget <= 0 || bnTarget > bnProofOfWorkLimit)
        return error("CheckProofOfWork() : nBits below minimum work");

    // Check proof of work matches claimed amount
    if (hash > bnTarget.getuint256())
        return error("CheckProofOfWork() : hash doesn't match nBits");

    return true;
}

// Return maximum amount of blocks that other nodes claim to have
int GetNumBlocksOfPeers()
{
    return std::max(cPeerBlockCounts.median(), Checkpoints::GetTotalBlocksEstimate());
}

bool IsInitialBlockDownload()
{
    if (pindexBest == NULL || fImporting || fReindex || nBestHeight < Checkpoints::GetTotalBlocksEstimate())
        return true;
    static int64 nLastUpdate;
    static CBlockIndex* pindexLastBest;
    if (pindexBest != pindexLastBest)
    {
        pindexLastBest = pindexBest;
        nLastUpdate = GetTime();
    }
    return (GetTime() - nLastUpdate < 10 &&
            pindexBest->GetBlockTime() < GetTime() - 24 * 60 * 60);
}

void static InvalidChainFound(CBlockIndex* pindexNew)
{
    if (pindexNew->nChainWork > nBestInvalidWork)
    {
        nBestInvalidWork = pindexNew->nChainWork;
        pblocktree->WriteBestInvalidWork(CBigNum(nBestInvalidWork));
        uiInterface.NotifyBlocksChanged();
    }
    printf("InvalidChainFound: invalid block=%s  height=%d  log2_work=%.8g  date=%s\n",
      pindexNew->GetBlockHash().ToString().c_str(), pindexNew->nHeight,
      log(pindexNew->nChainWork.getdouble())/log(2.0), DateTimeStrFormat("%Y-%m-%d %H:%M:%S",
      pindexNew->GetBlockTime()).c_str());
    printf("InvalidChainFound:  current best=%s  height=%d  log2_work=%.8g  date=%s\n",
      hashBestChain.ToString().c_str(), nBestHeight, log(nBestChainWork.getdouble())/log(2.0),
      DateTimeStrFormat("%Y-%m-%d %H:%M:%S", pindexBest->GetBlockTime()).c_str());
    if (pindexBest && nBestInvalidWork > nBestChainWork + (pindexBest->GetBlockWork() * 6).getuint256())
        printf("InvalidChainFound: Warning: Displayed transactions may not be correct! You may need to upgrade, or other nodes may need to upgrade.\n");
}

void static InvalidBlockFound(CBlockIndex *pindex) {
    pindex->nStatus |= BLOCK_FAILED_VALID;
    pblocktree->WriteBlockIndex(CDiskBlockIndex(pindex));
    setBlockIndexValid.erase(pindex);
    InvalidChainFound(pindex);
    if (pindex->pnext) {
        CValidationState stateDummy;
        ConnectBestBlock(stateDummy); // reorganise away from the failed block
    }
}

bool ConnectBestBlock(CValidationState &state) {
    do {
        CBlockIndex *pindexNewBest;

        {
            std::set<CBlockIndex*,CBlockIndexWorkComparator>::reverse_iterator it = setBlockIndexValid.rbegin();
            if (it == setBlockIndexValid.rend())
                return true;
            pindexNewBest = *it;
        }

        if (pindexNewBest == pindexBest || (pindexBest && pindexNewBest->nChainWork == pindexBest->nChainWork))
            return true; // nothing to do

        // check ancestry
        CBlockIndex *pindexTest = pindexNewBest;
        std::vector<CBlockIndex*> vAttach;
        do {
            if (pindexTest->nStatus & BLOCK_FAILED_MASK) {
                // mark descendants failed
                CBlockIndex *pindexFailed = pindexNewBest;
                while (pindexTest != pindexFailed) {
                    pindexFailed->nStatus |= BLOCK_FAILED_CHILD;
                    setBlockIndexValid.erase(pindexFailed);
                    pblocktree->WriteBlockIndex(CDiskBlockIndex(pindexFailed));
                    pindexFailed = pindexFailed->pprev;
                }
                InvalidChainFound(pindexNewBest);
                break;
            }

            if (pindexBest == NULL || pindexTest->nChainWork > pindexBest->nChainWork)
                vAttach.push_back(pindexTest);

            if (pindexTest->pprev == NULL || pindexTest->pnext != NULL) {
                reverse(vAttach.begin(), vAttach.end());
                BOOST_FOREACH(CBlockIndex *pindexSwitch, vAttach) {
                    boost::this_thread::interruption_point();
                    try {
                        if (!SetBestChain(state, pindexSwitch))
                            return false;
                    } catch(std::runtime_error &e) {
                        return state.Abort(_("System error: ") + e.what());
                    }
                }
                return true;
            }
            pindexTest = pindexTest->pprev;
        } while(true);
    } while(true);
}

void CBlockHeader::UpdateTime(const CBlockIndex* pindexPrev)
{
    nTime = max(pindexPrev->GetMedianTimePast()+1, GetAdjustedTime());

    // Updating time can change work required on testnet:
    if (fTestNet)
        nBits = GetNextWorkRequired(pindexPrev, this);
}











const CTxOut &CTransaction::GetOutputFor(const CTxIn& input, CCoinsViewCache& view)
{
    const CCoins &coins = view.GetCoins(input.prevout.hash);
    assert(coins.IsAvailable(input.prevout.n));
    return coins.vout[input.prevout.n];
}

int64 CTransaction::GetValueIn(CCoinsViewCache& inputs) const
{
    if (IsCoinBase())
        return 0;

    int64 nResult = 0;
    for (unsigned int i = 0; i < vin.size(); i++)
        nResult += GetOutputFor(vin[i], inputs).nValue;

    return nResult;
}

unsigned int CTransaction::GetP2SHSigOpCount(CCoinsViewCache& inputs) const
{
    if (IsCoinBase())
        return 0;

    unsigned int nSigOps = 0;
    for (unsigned int i = 0; i < vin.size(); i++)
    {
        const CTxOut &prevout = GetOutputFor(vin[i], inputs);
        if (prevout.scriptPubKey.IsPayToScriptHash())
            nSigOps += prevout.scriptPubKey.GetSigOpCount(vin[i].scriptSig);
    }
    return nSigOps;
}

void CTransaction::UpdateCoins(CValidationState &state, CCoinsViewCache &inputs, CTxUndo &txundo, int nHeight, const uint256 &txhash) const
{
    // mark inputs spent
    if (!IsCoinBase()) {
        BOOST_FOREACH(const CTxIn &txin, vin) {
            CCoins &coins = inputs.GetCoins(txin.prevout.hash);
            CTxInUndo undo;
            assert(coins.Spend(txin.prevout, undo));
            txundo.vprevout.push_back(undo);
        }
    }

    // add outputs
    assert(inputs.SetCoins(txhash, CCoins(*this, nHeight)));
}

bool CTransaction::HaveInputs(CCoinsViewCache &inputs) const
{
    if (!IsCoinBase()) {
        // first check whether information about the prevout hash is available
        for (unsigned int i = 0; i < vin.size(); i++) {
            const COutPoint &prevout = vin[i].prevout;
            if (!inputs.HaveCoins(prevout.hash))
                return false;
        }

        // then check whether the actual outputs are available
        for (unsigned int i = 0; i < vin.size(); i++) {
            const COutPoint &prevout = vin[i].prevout;
            const CCoins &coins = inputs.GetCoins(prevout.hash);
            if (!coins.IsAvailable(prevout.n))
                return false;
        }
    }
    return true;
}

bool CScriptCheck::operator()() const {
    const CScript &scriptSig = ptxTo->vin[nIn].scriptSig;
    if (!VerifyScript(scriptSig, scriptPubKey, *ptxTo, nIn, nFlags, nHashType))
        return error("CScriptCheck() : %s VerifySignature failed", ptxTo->GetHash().ToString().c_str());
    return true;
}

bool VerifySignature(const CCoins& txFrom, const CTransaction& txTo, unsigned int nIn, unsigned int flags, int nHashType)
{
    return CScriptCheck(txFrom, txTo, nIn, flags, nHashType)();
}

bool CTransaction::CheckInputs(CValidationState &state, CCoinsViewCache &inputs, bool fScriptChecks, unsigned int flags, std::vector<CScriptCheck> *pvChecks) const
{
    if (!IsCoinBase())
    {
        if (pvChecks)
            pvChecks->reserve(vin.size());

        // This doesn't trigger the DoS code on purpose; if it did, it would make it easier
        // for an attacker to attempt to split the network.
        if (!HaveInputs(inputs))
            return state.Invalid(error("CheckInputs() : %s inputs unavailable", GetHash().ToString().c_str()));

        // While checking, GetBestBlock() refers to the parent block.
        // This is also true for mempool checks.
        int nSpendHeight = inputs.GetBestBlock()->nHeight + 1;
        int64 nValueIn = 0;
        int64 nFees = 0;
        for (unsigned int i = 0; i < vin.size(); i++)
        {
            const COutPoint &prevout = vin[i].prevout;
            const CCoins &coins = inputs.GetCoins(prevout.hash);

            // If prev is coinbase, check that it's matured
            if (coins.IsCoinBase()) {
                if (nSpendHeight - coins.nHeight < COINBASE_MATURITY)
                    return state.Invalid(error("CheckInputs() : tried to spend coinbase at depth %d", nSpendHeight - coins.nHeight));
            }

            // Check for negative or overflow input values
            nValueIn += coins.vout[prevout.n].nValue;
            if (!MoneyRange(coins.vout[prevout.n].nValue) || !MoneyRange(nValueIn))
                return state.DoS(100, error("CheckInputs() : txin values out of range"));

        }

        if (nValueIn < GetValueOut())
            return state.DoS(100, error("CheckInputs() : %s value in < value out", GetHash().ToString().c_str()));

        // Tally transaction fees
        int64 nTxFee = nValueIn - GetValueOut();
        if (nTxFee < 0)
            return state.DoS(100, error("CheckInputs() : %s nTxFee < 0", GetHash().ToString().c_str()));
        nFees += nTxFee;
        if (!MoneyRange(nFees))
            return state.DoS(100, error("CheckInputs() : nFees out of range"));

        // The first loop above does all the inexpensive checks.
        // Only if ALL inputs pass do we perform expensive ECDSA signature checks.
        // Helps prevent CPU exhaustion attacks.

        // Skip ECDSA signature verification when connecting blocks
        // before the last block chain checkpoint. This is safe because block merkle hashes are
        // still computed and checked, and any change will be caught at the next checkpoint.
        if (fScriptChecks) {
            for (unsigned int i = 0; i < vin.size(); i++) {
                const COutPoint &prevout = vin[i].prevout;
                const CCoins &coins = inputs.GetCoins(prevout.hash);

                // Verify signature
                CScriptCheck check(coins, *this, i, flags, 0);
                if (pvChecks) {
                    pvChecks->push_back(CScriptCheck());
                    check.swap(pvChecks->back());
                } else if (!check()) {
                    if (flags & SCRIPT_VERIFY_STRICTENC) {
                        // For now, check whether the failure was caused by non-canonical
                        // encodings or not; if so, don't trigger DoS protection.
                        CScriptCheck check(coins, *this, i, flags & (~SCRIPT_VERIFY_STRICTENC), 0);
                        if (check())
                            return state.Invalid();
                    }
                    return state.DoS(100,false);
                }
            }
        }
    }

    return true;
}




bool CBlock::DisconnectBlock(CValidationState &state, CBlockIndex *pindex, CCoinsViewCache &view, bool *pfClean)
{
    assert(pindex == view.GetBestBlock());

    if (pfClean)
        *pfClean = false;

    bool fClean = true;

    CBlockUndo blockUndo;
    CDiskBlockPos pos = pindex->GetUndoPos();
    if (pos.IsNull())
        return error("DisconnectBlock() : no undo data available");
    if (!blockUndo.ReadFromDisk(pos, pindex->pprev->GetBlockHash()))
        return error("DisconnectBlock() : failure reading undo data");

    if (blockUndo.vtxundo.size() + 1 != vtx.size())
        return error("DisconnectBlock() : block and undo data inconsistent");

    // undo transactions in reverse order
    for (int i = vtx.size() - 1; i >= 0; i--) {
        const CTransaction &tx = vtx[i];
        uint256 hash = tx.GetHash();

        // check that all outputs are available
        if (!view.HaveCoins(hash)) {
            fClean = fClean && error("DisconnectBlock() : outputs still spent? database corrupted");
            view.SetCoins(hash, CCoins());
        }
        CCoins &outs = view.GetCoins(hash);

        CCoins outsBlock = CCoins(tx, pindex->nHeight);
        // The CCoins serialization does not serialize negative numbers.
        // No network rules currently depend on the version here, so an inconsistency is harmless
        // but it must be corrected before txout nversion ever influences a network rule.
        if (outsBlock.nVersion < 0)
            outs.nVersion = outsBlock.nVersion;
        if (outs != outsBlock)
            fClean = fClean && error("DisconnectBlock() : added transaction mismatch? database corrupted");

        // remove outputs
        outs = CCoins();

        // restore inputs
        if (i > 0) { // not coinbases
            const CTxUndo &txundo = blockUndo.vtxundo[i-1];
            if (txundo.vprevout.size() != tx.vin.size())
                return error("DisconnectBlock() : transaction and undo data inconsistent");
            for (unsigned int j = tx.vin.size(); j-- > 0;) {
                const COutPoint &out = tx.vin[j].prevout;
                const CTxInUndo &undo = txundo.vprevout[j];
                CCoins coins;
                view.GetCoins(out.hash, coins); // this can fail if the prevout was already entirely spent
                if (undo.nHeight != 0) {
                    // undo data contains height: this is the last output of the prevout tx being spent
                    if (!coins.IsPruned())
                        fClean = fClean && error("DisconnectBlock() : undo data overwriting existing transaction");
                    coins = CCoins();
                    coins.fCoinBase = undo.fCoinBase;
                    coins.nHeight = undo.nHeight;
                    coins.nVersion = undo.nVersion;
                } else {
                    if (coins.IsPruned())
                        fClean = fClean && error("DisconnectBlock() : undo data adding output to missing transaction");
                }
                if (coins.IsAvailable(out.n))
                    fClean = fClean && error("DisconnectBlock() : undo data overwriting existing output");
                if (coins.vout.size() < out.n+1)
                    coins.vout.resize(out.n+1);
                coins.vout[out.n] = undo.txout;
                if (!view.SetCoins(out.hash, coins))
                    return error("DisconnectBlock() : cannot restore coin inputs");
            }
        }
    }

    // move best block pointer to prevout block
    view.SetBestBlock(pindex->pprev);

    if (pfClean) {
        *pfClean = fClean;
        return true;
    } else {
        return fClean;
    }
}

void static FlushBlockFile(bool fFinalize = false)
{
    LOCK(cs_LastBlockFile);

    CDiskBlockPos posOld(nLastBlockFile, 0);

    FILE *fileOld = OpenBlockFile(posOld);
    if (fileOld) {
        if (fFinalize)
            TruncateFile(fileOld, infoLastBlockFile.nSize);
        FileCommit(fileOld);
        fclose(fileOld);
    }

    fileOld = OpenUndoFile(posOld);
    if (fileOld) {
        if (fFinalize)
            TruncateFile(fileOld, infoLastBlockFile.nUndoSize);
        FileCommit(fileOld);
        fclose(fileOld);
    }
}

bool FindUndoPos(CValidationState &state, int nFile, CDiskBlockPos &pos, unsigned int nAddSize);

static CCheckQueue<CScriptCheck> scriptcheckqueue(128);

void ThreadScriptCheck() {
    RenameThread("bitcoin-scriptch");
    scriptcheckqueue.Thread();
}

bool CBlock::ConnectBlock(CValidationState &state, CBlockIndex* pindex, CCoinsViewCache &view, bool fJustCheck)
{
    // Check it again in case a previous version let a bad block in
    if (!CheckBlock(state, !fJustCheck, !fJustCheck))
        return false;

    // verify that the view's current state corresponds to the previous block
    assert(pindex->pprev == view.GetBestBlock());

    // Special case for the genesis block, skipping connection of its transactions
    // (its coinbase is unspendable)
    if (GetHash() == hashGenesisBlock) {
        view.SetBestBlock(pindex);
        pindexGenesisBlock = pindex;
        return true;
    }

    bool fScriptChecks = pindex->nHeight >= Checkpoints::GetTotalBlocksEstimate();

    // Do not allow blocks that contain transactions which 'overwrite' older transactions,
    // unless those are already completely spent.
    // If such overwrites are allowed, coinbases and transactions depending upon those
    // can be duplicated to remove the ability to spend the first instance -- even after
    // being sent to another address.
    // See BIP30 and http://r6.ca/blog/20120206T005236Z.html for more information.
    // This logic is not necessary for memory pool transactions, as AcceptToMemoryPool
    // already refuses previously-known transaction ids entirely.
    // This rule was originally applied all blocks whose timestamp was after October 1, 2012, 0:00 UTC.
    // Now that the whole chain is irreversibly beyond that time it is applied to all blocks,
    // this prevents exploiting the issue against nodes in their initial block download.
    bool fEnforceBIP30 = true;

    if (fEnforceBIP30) {
        for (unsigned int i=0; i<vtx.size(); i++) {
            uint256 hash = GetTxHash(i);
            if (view.HaveCoins(hash) && !view.GetCoins(hash).IsPruned())
// FBX
//                return state.DoS(100, error("ConnectBlock() : tried to overwrite transaction"));
            {
                if (pindex->nHeight <= 3075)
                    printf("ConnectBlock() : tried to overwrite transaction (allowed for early blocks)\n");
                else
                    return state.DoS(100, error("ConnectBlock() : tried to overwrite transaction"));
            }
        }
    }

    // BIP16 didn't become active until Oct 1 2012
    int64 nBIP16SwitchTime = 1349049600;
    bool fStrictPayToScriptHash = (pindex->nTime >= nBIP16SwitchTime);

    unsigned int flags = SCRIPT_VERIFY_NOCACHE |
                         (fStrictPayToScriptHash ? SCRIPT_VERIFY_P2SH : SCRIPT_VERIFY_NONE);

    CBlockUndo blockundo;

    CCheckQueueControl<CScriptCheck> control(fScriptChecks && nScriptCheckThreads ? &scriptcheckqueue : NULL);

    int64 nStart = GetTimeMicros();
    int64 nFees = 0;
    int nInputs = 0;
    unsigned int nSigOps = 0;
    CDiskTxPos pos(pindex->GetBlockPos(), GetSizeOfCompactSize(vtx.size()));
    std::vector<std::pair<uint256, CDiskTxPos> > vPos;
    vPos.reserve(vtx.size());
    for (unsigned int i=0; i<vtx.size(); i++)
    {
        const CTransaction &tx = vtx[i];

        nInputs += tx.vin.size();
        nSigOps += tx.GetLegacySigOpCount();
        if (nSigOps > MAX_BLOCK_SIGOPS)
            return state.DoS(100, error("ConnectBlock() : too many sigops"));

        if (!tx.IsCoinBase())
        {
            if (!tx.HaveInputs(view))
                return state.DoS(100, error("ConnectBlock() : inputs missing/spent"));

            if (fStrictPayToScriptHash)
            {
                // Add in sigops done by pay-to-script-hash inputs;
                // this is to prevent a "rogue miner" from creating
                // an incredibly-expensive-to-validate block.
                nSigOps += tx.GetP2SHSigOpCount(view);
                if (nSigOps > MAX_BLOCK_SIGOPS)
                     return state.DoS(100, error("ConnectBlock() : too many sigops"));
            }

            nFees += tx.GetValueIn(view)-tx.GetValueOut();

            std::vector<CScriptCheck> vChecks;
            if (!tx.CheckInputs(state, view, fScriptChecks, flags, nScriptCheckThreads ? &vChecks : NULL))
                return false;
            control.Add(vChecks);
        }

        CTxUndo txundo;
        tx.UpdateCoins(state, view, txundo, pindex->nHeight, GetTxHash(i));
        if (!tx.IsCoinBase())
            blockundo.vtxundo.push_back(txundo);

        vPos.push_back(std::make_pair(GetTxHash(i), pos));
        pos.nTxOffset += ::GetSerializeSize(tx, SER_DISK, CLIENT_VERSION);
    }
    int64 nTime = GetTimeMicros() - nStart;
    if (fBenchmark)
        printf("- Connect %u transactions: %.2fms (%.3fms/tx, %.3fms/txin)\n", (unsigned)vtx.size(), 0.001 * nTime, 0.001 * nTime / vtx.size(), nInputs <= 1 ? 0 : 0.001 * nTime / (nInputs-1));

    if (vtx[0].GetValueOut() > GetBlockValue(pindex->nHeight, nFees))
        return state.DoS(100, error("ConnectBlock() : coinbase pays too much (actual=%"PRI64d" vs limit=%"PRI64d")", vtx[0].GetValueOut(), GetBlockValue(pindex->nHeight, nFees)));

    if (!control.Wait())
        return state.DoS(100, false);
    int64 nTime2 = GetTimeMicros() - nStart;
    if (fBenchmark)
        printf("- Verify %u txins: %.2fms (%.3fms/txin)\n", nInputs - 1, 0.001 * nTime2, nInputs <= 1 ? 0 : 0.001 * nTime2 / (nInputs-1));

    if (fJustCheck)
        return true;

    // Write undo information to disk
    if (pindex->GetUndoPos().IsNull() || (pindex->nStatus & BLOCK_VALID_MASK) < BLOCK_VALID_SCRIPTS)
    {
        if (pindex->GetUndoPos().IsNull()) {
            CDiskBlockPos pos;
            if (!FindUndoPos(state, pindex->nFile, pos, ::GetSerializeSize(blockundo, SER_DISK, CLIENT_VERSION) + 40))
                return error("ConnectBlock() : FindUndoPos failed");
            if (!blockundo.WriteToDisk(pos, pindex->pprev->GetBlockHash()))
                return state.Abort(_("Failed to write undo data"));

            // update nUndoPos in block index
            pindex->nUndoPos = pos.nPos;
            pindex->nStatus |= BLOCK_HAVE_UNDO;
        }

        pindex->nStatus = (pindex->nStatus & ~BLOCK_VALID_MASK) | BLOCK_VALID_SCRIPTS;

        CDiskBlockIndex blockindex(pindex);
        if (!pblocktree->WriteBlockIndex(blockindex))
            return state.Abort(_("Failed to write block index"));
    }

    if (fTxIndex)
        if (!pblocktree->WriteTxIndex(vPos))
            return state.Abort(_("Failed to write transaction index"));

    // add this block to the view's block chain
    assert(view.SetBestBlock(pindex));

    // Watch for transactions paying to me
    for (unsigned int i=0; i<vtx.size(); i++)
        SyncWithWallets(GetTxHash(i), vtx[i], this, true);

    return true;
}

bool SetBestChain(CValidationState &state, CBlockIndex* pindexNew)
{
    // All modifications to the coin state will be done in this cache.
    // Only when all have succeeded, we push it to pcoinsTip.
    CCoinsViewCache view(*pcoinsTip, true);

    // Find the fork (typically, there is none)
    CBlockIndex* pfork = view.GetBestBlock();
    CBlockIndex* plonger = pindexNew;
    while (pfork && pfork != plonger)
    {
        while (plonger->nHeight > pfork->nHeight) {
            plonger = plonger->pprev;
            assert(plonger != NULL);
        }
        if (pfork == plonger)
            break;
        pfork = pfork->pprev;
        assert(pfork != NULL);
    }

    // List of what to disconnect (typically nothing)
    vector<CBlockIndex*> vDisconnect;
    for (CBlockIndex* pindex = view.GetBestBlock(); pindex != pfork; pindex = pindex->pprev)
        vDisconnect.push_back(pindex);

    // List of what to connect (typically only pindexNew)
    vector<CBlockIndex*> vConnect;
    for (CBlockIndex* pindex = pindexNew; pindex != pfork; pindex = pindex->pprev)
        vConnect.push_back(pindex);
    reverse(vConnect.begin(), vConnect.end());

    if (vDisconnect.size() > 0) {
        printf("REORGANIZE: Disconnect %"PRIszu" blocks; %s..\n", vDisconnect.size(), pfork->GetBlockHash().ToString().c_str());
        printf("REORGANIZE: Connect %"PRIszu" blocks; ..%s\n", vConnect.size(), pindexNew->GetBlockHash().ToString().c_str());
    }

    // Disconnect shorter branch
    vector<CTransaction> vResurrect;
    BOOST_FOREACH(CBlockIndex* pindex, vDisconnect) {
        CBlock block;
        if (!block.ReadFromDisk(pindex))
            return state.Abort(_("Failed to read block"));
        int64 nStart = GetTimeMicros();
        if (!block.DisconnectBlock(state, pindex, view))
            return error("SetBestBlock() : DisconnectBlock %s failed", pindex->GetBlockHash().ToString().c_str());
        if (fBenchmark)
            printf("- Disconnect: %.2fms\n", (GetTimeMicros() - nStart) * 0.001);

        // Queue memory transactions to resurrect.
        // We only do this for blocks after the last checkpoint (reorganisation before that
        // point should only happen with -reindex/-loadblock, or a misbehaving peer.
        BOOST_FOREACH(const CTransaction& tx, block.vtx)
            if (!tx.IsCoinBase() && pindex->nHeight > Checkpoints::GetTotalBlocksEstimate())
                vResurrect.push_back(tx);
    }

    // Connect longer branch
    vector<CTransaction> vDelete;
    BOOST_FOREACH(CBlockIndex *pindex, vConnect) {
        CBlock block;
        if (!block.ReadFromDisk(pindex))
            return state.Abort(_("Failed to read block"));
        int64 nStart = GetTimeMicros();
        if (!block.ConnectBlock(state, pindex, view)) {
            if (state.IsInvalid()) {
                InvalidChainFound(pindexNew);
                InvalidBlockFound(pindex);
            }
            return error("SetBestBlock() : ConnectBlock %s failed", pindex->GetBlockHash().ToString().c_str());
        }
        if (fBenchmark)
            printf("- Connect: %.2fms\n", (GetTimeMicros() - nStart) * 0.001);

        // Queue memory transactions to delete
        BOOST_FOREACH(const CTransaction& tx, block.vtx)
            vDelete.push_back(tx);
    }

    // Flush changes to global coin state
    int64 nStart = GetTimeMicros();
    int nModified = view.GetCacheSize();
    assert(view.Flush());
    int64 nTime = GetTimeMicros() - nStart;
    if (fBenchmark)
        printf("- Flush %i transactions: %.2fms (%.4fms/tx)\n", nModified, 0.001 * nTime, 0.001 * nTime / nModified);

    // Make sure it's successfully written to disk before changing memory structure
    bool fIsInitialDownload = IsInitialBlockDownload();
    if (!fIsInitialDownload || pcoinsTip->GetCacheSize() > nCoinCacheSize) {
        // Typical CCoins structures on disk are around 100 bytes in size.
        // Pushing a new one to the database can cause it to be written
        // twice (once in the log, and once in the tables). This is already
        // an overestimation, as most will delete an existing entry or
        // overwrite one. Still, use a conservative safety factor of 2.
        if (!CheckDiskSpace(100 * 2 * 2 * pcoinsTip->GetCacheSize()))
            return state.Error();
        FlushBlockFile();
        pblocktree->Sync();
        if (!pcoinsTip->Flush())
            return state.Abort(_("Failed to write to coin database"));
    }

    // At this point, all changes have been done to the database.
    // Proceed by updating the memory structures.

    // Disconnect shorter branch
    BOOST_FOREACH(CBlockIndex* pindex, vDisconnect)
        if (pindex->pprev)
            pindex->pprev->pnext = NULL;

    // Connect longer branch
    BOOST_FOREACH(CBlockIndex* pindex, vConnect)
        if (pindex->pprev)
            pindex->pprev->pnext = pindex;

    // Resurrect memory transactions that were in the disconnected branch
    BOOST_FOREACH(CTransaction& tx, vResurrect) {
        // ignore validation errors in resurrected transactions
        CValidationState stateDummy;
        if (!tx.AcceptToMemoryPool(stateDummy, true, false))
            mempool.remove(tx, true);
    }

    // Delete redundant memory transactions that are in the connected branch
    BOOST_FOREACH(CTransaction& tx, vDelete) {
        mempool.remove(tx);
        mempool.removeConflicts(tx);
    }

    // Update best block in wallet (so we can detect restored wallets)
    if ((pindexNew->nHeight % 20160) == 0 || (!fIsInitialDownload && (pindexNew->nHeight % 144) == 0))
    {
        const CBlockLocator locator(pindexNew);
        ::SetBestChain(locator);
    }

    // New best block
    hashBestChain = pindexNew->GetBlockHash();
    pindexBest = pindexNew;
    pblockindexFBBHLast = NULL;
    nBestHeight = pindexBest->nHeight;
    nBestChainWork = pindexNew->nChainWork;
    nTimeBestReceived = GetTime();
    nTransactionsUpdated++;
    printf("SetBestChain: new best=%s  height=%d  log2_work=%.8g  tx=%lu  date=%s progress=%f\n",
      hashBestChain.ToString().c_str(), nBestHeight, log(nBestChainWork.getdouble())/log(2.0), (unsigned long)pindexNew->nChainTx,
      DateTimeStrFormat("%Y-%m-%d %H:%M:%S", pindexBest->GetBlockTime()).c_str(),
      Checkpoints::GuessVerificationProgress(pindexBest));

    // Check the version of the last 100 blocks to see if we need to upgrade:
    if (!fIsInitialDownload)
    {
        int nUpgraded = 0;
        const CBlockIndex* pindex = pindexBest;
        for (int i = 0; i < 100 && pindex != NULL; i++)
        {
            if (pindex->nVersion > CBlock::CURRENT_VERSION)
                ++nUpgraded;
            pindex = pindex->pprev;
        }
        if (nUpgraded > 0)
            printf("SetBestChain: %d of last 100 blocks above version %d\n", nUpgraded, CBlock::CURRENT_VERSION);
        if (nUpgraded > 100/2)
            // strMiscWarning is read by GetWarnings(), called by Qt and the JSON-RPC code to warn the user:
            strMiscWarning = _("Warning: This version is obsolete, upgrade required!");
    }

    std::string strCmd = GetArg("-blocknotify", "");

    if (!fIsInitialDownload && !strCmd.empty())
    {
        boost::replace_all(strCmd, "%s", hashBestChain.GetHex());
        boost::thread t(runCommand, strCmd); // thread runs free
    }

    return true;
}


bool CBlock::AddToBlockIndex(CValidationState &state, const CDiskBlockPos &pos)
{
    // Check for duplicate
    uint256 hash = GetHash();
    if (mapBlockIndex.count(hash))
        return state.Invalid(error("AddToBlockIndex() : %s already exists", hash.ToString().c_str()));

    // Construct new block index object
    CBlockIndex* pindexNew = new CBlockIndex(*this);
    assert(pindexNew);
    map<uint256, CBlockIndex*>::iterator mi = mapBlockIndex.insert(make_pair(hash, pindexNew)).first;
    pindexNew->phashBlock = &((*mi).first);
    map<uint256, CBlockIndex*>::iterator miPrev = mapBlockIndex.find(hashPrevBlock);
    if (miPrev != mapBlockIndex.end())
    {
        pindexNew->pprev = (*miPrev).second;
        pindexNew->nHeight = pindexNew->pprev->nHeight + 1;
    }
    pindexNew->nTx = vtx.size();
    pindexNew->nChainWork = (pindexNew->pprev ? pindexNew->pprev->nChainWork : 0) + pindexNew->GetBlockWork().getuint256();
    pindexNew->nChainTx = (pindexNew->pprev ? pindexNew->pprev->nChainTx : 0) + pindexNew->nTx;
    pindexNew->nFile = pos.nFile;
    pindexNew->nDataPos = pos.nPos;
    pindexNew->nUndoPos = 0;
    pindexNew->nStatus = BLOCK_VALID_TRANSACTIONS | BLOCK_HAVE_DATA;
    setBlockIndexValid.insert(pindexNew);

    if (!pblocktree->WriteBlockIndex(CDiskBlockIndex(pindexNew)))
        return state.Abort(_("Failed to write block index"));

    // New best?
    if (!ConnectBestBlock(state))
        return false;

    if (pindexNew == pindexBest)
    {
        // Notify UI to display prev block's coinbase if it was ours
        static uint256 hashPrevBestCoinBase;
        UpdatedTransaction(hashPrevBestCoinBase);
        hashPrevBestCoinBase = GetTxHash(0);
    }

    if (!pblocktree->Flush())
        return state.Abort(_("Failed to sync block index"));

    uiInterface.NotifyBlocksChanged();
    return true;
}


bool FindBlockPos(CValidationState &state, CDiskBlockPos &pos, unsigned int nAddSize, unsigned int nHeight, uint64 nTime, bool fKnown = false)
{
    bool fUpdatedLast = false;

    LOCK(cs_LastBlockFile);

    if (fKnown) {
        if (nLastBlockFile != pos.nFile) {
            nLastBlockFile = pos.nFile;
            infoLastBlockFile.SetNull();
            pblocktree->ReadBlockFileInfo(nLastBlockFile, infoLastBlockFile);
            fUpdatedLast = true;
        }
    } else {
        while (infoLastBlockFile.nSize + nAddSize >= MAX_BLOCKFILE_SIZE) {
            printf("Leaving block file %i: %s\n", nLastBlockFile, infoLastBlockFile.ToString().c_str());
            FlushBlockFile(true);
            nLastBlockFile++;
            infoLastBlockFile.SetNull();
            pblocktree->ReadBlockFileInfo(nLastBlockFile, infoLastBlockFile); // check whether data for the new file somehow already exist; can fail just fine
            fUpdatedLast = true;
        }
        pos.nFile = nLastBlockFile;
        pos.nPos = infoLastBlockFile.nSize;
    }

    infoLastBlockFile.nSize += nAddSize;
    infoLastBlockFile.AddBlock(nHeight, nTime);

    if (!fKnown) {
        unsigned int nOldChunks = (pos.nPos + BLOCKFILE_CHUNK_SIZE - 1) / BLOCKFILE_CHUNK_SIZE;
        unsigned int nNewChunks = (infoLastBlockFile.nSize + BLOCKFILE_CHUNK_SIZE - 1) / BLOCKFILE_CHUNK_SIZE;
        if (nNewChunks > nOldChunks) {
            if (CheckDiskSpace(nNewChunks * BLOCKFILE_CHUNK_SIZE - pos.nPos)) {
                FILE *file = OpenBlockFile(pos);
                if (file) {
                    printf("Pre-allocating up to position 0x%x in blk%05u.dat\n", nNewChunks * BLOCKFILE_CHUNK_SIZE, pos.nFile);
                    AllocateFileRange(file, pos.nPos, nNewChunks * BLOCKFILE_CHUNK_SIZE - pos.nPos);
                    fclose(file);
                }
            }
            else
                return state.Error();
        }
    }

    if (!pblocktree->WriteBlockFileInfo(nLastBlockFile, infoLastBlockFile))
        return state.Abort(_("Failed to write file info"));
    if (fUpdatedLast)
        pblocktree->WriteLastBlockFile(nLastBlockFile);

    return true;
}

bool FindUndoPos(CValidationState &state, int nFile, CDiskBlockPos &pos, unsigned int nAddSize)
{
    pos.nFile = nFile;

    LOCK(cs_LastBlockFile);

    unsigned int nNewSize;
    if (nFile == nLastBlockFile) {
        pos.nPos = infoLastBlockFile.nUndoSize;
        nNewSize = (infoLastBlockFile.nUndoSize += nAddSize);
        if (!pblocktree->WriteBlockFileInfo(nLastBlockFile, infoLastBlockFile))
            return state.Abort(_("Failed to write block info"));
    } else {
        CBlockFileInfo info;
        if (!pblocktree->ReadBlockFileInfo(nFile, info))
            return state.Abort(_("Failed to read block info"));
        pos.nPos = info.nUndoSize;
        nNewSize = (info.nUndoSize += nAddSize);
        if (!pblocktree->WriteBlockFileInfo(nFile, info))
            return state.Abort(_("Failed to write block info"));
    }

    unsigned int nOldChunks = (pos.nPos + UNDOFILE_CHUNK_SIZE - 1) / UNDOFILE_CHUNK_SIZE;
    unsigned int nNewChunks = (nNewSize + UNDOFILE_CHUNK_SIZE - 1) / UNDOFILE_CHUNK_SIZE;
    if (nNewChunks > nOldChunks) {
        if (CheckDiskSpace(nNewChunks * UNDOFILE_CHUNK_SIZE - pos.nPos)) {
            FILE *file = OpenUndoFile(pos);
            if (file) {
                printf("Pre-allocating up to position 0x%x in rev%05u.dat\n", nNewChunks * UNDOFILE_CHUNK_SIZE, pos.nFile);
                AllocateFileRange(file, pos.nPos, nNewChunks * UNDOFILE_CHUNK_SIZE - pos.nPos);
                fclose(file);
            }
        }
        else
            return state.Error();
    }

    return true;
}


// FBX proof of stake voting test
std::string strPosxExactMatch;
int64 posxValueDiffIn;    // aggregate money inflow for all addresses matching the test string
int64 posxValueDiffOut;   // aggregate money outflow for all addresses matching the test string
int64 posxTxValue;      // output of one tx
int posxErrCount;
int posxCritCount;
int posxSanCount;
int posxSig;            // inflow or outflow
static const int posx_ReverseBase58[128+1] = { 0,
                                            0, 0, 0, 0, 0, 0, 0, 0,
                                            0, 0, 0, 0, 0, 0, 0, 0,
                                            0, 0, 0, 0, 0, 0, 0, 0,
                                            0, 0, 0, 0, 0, 0, 0, 0,
                                            0, 0, 0, 0, 0, 0, 0, 0,
                                            0, 0, 0, 0, 0, 0, 0, 0,
                                            1, 2, 3, 4, 5, 6, 7, 8,
                                            9, 0, 0, 0, 0, 0, 0, 0,          // 57..64
                                            10, 11, 12, 13, 14, 15, 16, 17,  // 65..72, 'A'..'H'
                                            0, 18, 19, 20, 21, 22, 0, 23,    // 73..80, 'I'..'P' ('I' and 'O' skipped)
                                            24, 25, 26, 27, 28, 29, 30, 31,  // 81..88, 'Q'..'X'
                                            32, 33, 0, 0, 0, 0, 0, 0,        // 89..96
                                            34, 35, 36, 37, 38, 39, 40, 41,  // 97..104, 'a'..'h'
                                            42, 43, 44, 45, 0, 46, 47, 48,   // 105..112, 'i'..'p' ('l' skipped)
                                            49, 50, 51, 52, 53, 54, 55, 56,
                                            57, 58, 0, 0, 0, 0, 0, 0 };
#define POSX_ROUND_MAX 2    // allow 2 trading rounds at the same time for rollover
#define POSX_PAIR_MAX 4     // can be extended to about 12 without much code change
#define POSX_ANSWER_MAX 7   // 8 possible answers of the oracle -- 7 choices to bet on
#define POSX_ASK_MAX 2      // bid + ask
#define POSX_PROB_MAX 19    // 5% minimum price tick, can be extended to 49 (2% minimum price tick) without much code change
std::string strPosxAddresses[POSX_ROUND_MAX][POSX_PAIR_MAX][POSX_ANSWER_MAX][POSX_ASK_MAX][POSX_PROB_MAX] = {
    {
// even trading round, pair 1
{{{"fXB111a99YrhZCn9yi2jM8Rv6dyQXfh297", "fXB112nHRwbtpEH38z3JfvqjSQ767U4Ev6", "fXB113AmpAwimBZkvfpQ4zgRBTPfcN8QvE", "fXB114TvBJhai8oiBjb8txZ7oJU7Mch4dv", "fXB115DU5fD39nXSmSztg3GxeeGVTxknPe", "fXB116UNKk9GtVam2bymm4k7yw8vVHumXf", "fXB117MMTYvZyRX3pmm3xhRFjpPBbZyWpy", "fXB1187vMoyQrngR6P2s6skLf4a29wsSYp", "fXB119zus6hqSeGfyhhtS8YxgsSRT5f1uS", "fXB11ALky91FeKYaoGFbjoLq7ecd6eYbU3", "fXB11BgrEwTs8iEMoheYQALgkSgUjyyaDG", "fXB11CqoUgbH6r3RiFvGLouNtLj5oKsFeJ", "fXB11DKEYkifpMApSCBVkZiAj3VqrKcFum", "fXB11EJSMBEpvSvjLBkK8CRo7GHpZ73WXN", "fXB11Fdb8RAeCveF9G9ebyyhBrW8g3S8FD", "fXB11G9d6oDBw3dzHPCQ6tgz42zEy1qBJe", "fXB11HUJZpb9DzJzc9SCPoP8NmQJVUaxrd", "fXB11JFLYo9uukHu2FrE68Ds7PbWvoTk4r", "fXB11KSeoVARKAzM9GGBasyRyvrsSP17cD"},
  {"fXA111tNcnnCpBzEyLRrfRTE4KCJHKx4pa", "fXA112TUpCDUXTU8XJtj2tzJ4imyYdXT6s", "fXA113yJAJqfX4DNBZbwRAeYqDbKLrJ8Gq", "fXA114Nm7usp6hRcfCq1S6zM3tyk2s7bSU", "fXA1158kbrJpSYjJniKjapi8nXW7omeYws", "fXA116bVPs2dxi6DvfmBAu3yvMdZ1p72yh", "fXA117SYdTrhB68Jh1mnaK99FNFjWz8wHm", "fXA118KBcWAgTeKpMzYHG4DFbdMh9kUAzm", "fXA119bvkGNyevJ5HNNyrqnRDsDkc9kuKQ", "fXA11AK8L8zeDER6Vh917btTqirQCqQR5Q", "fXA11BLcLv9EmHDYMRGn87P5mu1PxPyw8i", "fXA11Cxtc51Xd3JSV1wNmW5YtU3wK112Ke", "fXA11Dn63ARCEquMep3fCpALHx1eqYrAsz", "fXA11E2MSxZQavNnktMksjYuR7LriEzK5V", "fXA11FumJskaEZCvJ5QfHqKot71f52meub", "fXA11GPidcqu7nTCs77TZcsvF36mGr4x2d", "fXA11HoEoPAgQqcTrP9vj2fe4LN5aQkuAa", "fXA11JY8fKFhYDk4EKcT5QHaAfWMoZG6Mx", "fXA11KdtRfkz5r3wXHT1HUxXow1n9Tsa6T"}},
 {{"fXB121s7Tbb2e4Lj4PrBBnW6ZxWkZ279fn", "fXB122fyVSTvqmyn9FMrykMYdWsatUKeRZ", "fXB123oVAa1oPUtwkMMRP7eG9yiUjykAHN", "fXB124mbUWVns5MuxSz2FaLtD6m1p3iHp9", "fXB125swbSSMU9vN18kz1M2GrPpYbAxqmW", "fXB1263en2B1CoipLWsASNyxts4rR2Xq4t", "fXB127n1jKWStAUuGZ8Kh7hCb5GvT3bUxc", "fXB1281q3reusXWppUTACmh6ocKgeJMydU", "fXB129qcYdT6Azb3JPcoPsjnGYcJRESpKt", "fXB12A6HNzPymxMctW9LoLZZDFm6fiUnXN", "fXB12BMsRsduV81nyFRMsZnsnAf3nUSGWz", "fXB12CJXqwwo6pLVTeDM8po88LPHxCfqQo", "fXB12DmcBGRn4ydBdu8EaHd1vqCuxMDs3e", "fXB12EnoZNiEPfv8nsobrUe7Z9YBezXbMc", "fXB12FtxVX9R1HBxJf2P5BhGws7yweypYj", "fXB12GR7Kstt2iTv4oix3BBKJQGrbDbJPe", "fXB12HZscbiDJF2qyt5tJjF5qSxG55kLrV", "fXB12JQPftd5o2zteRA1Pq9aM3UJWmSC81", "fXB12KVsr5T9z6T6aMH2LEHwAyfN4mshNj"},
  {"fXA121bxPwKDmFfvCv8tDgAQpjGAZPnqqA", "fXA122XhUPzqNmYq9pDzniWochKfr4n7K6", "fXA1232xh3Eb46ALqbP7WbsD2yq4gJBsKN", "fXA124k9KwkpeQqhaydybvX6A5Z5WH63nb", "fXA125w1Ta2wq3voFSvyff7ux2kYKiYyT3", "fXA126NEJXnfod5VfZERVTnMTxo7jXknh6", "fXA127kd95mnpSWwf3nYCmoaKWrMAW41Rw", "fXA128fRiyoX9Mkwag4rQXekak4L5wMGBt", "fXA129wVtE7BaC57qWUj1gfrh8F7hb85xW", "fXA12APXDPTGiByZQvYtQMfPe3Mk1fWNSM", "fXA12Bg95WkcLve8LxC57B5euKBYtKYt6M", "fXA12CBiPUMxQ9kvTkbntrkz3PheQiPgNU", "fXA12Dn1S8tkwj3DDE3a1uZJMU2TG4bXRP", "fXA12E4t33wi46mT6GUqToASbSFesjbRkM", "fXA12FQQekHvkFrBph529EgyeHYm4HPauF", "fXA12GEkd8LM9CgQFZohbt47rzEWgtL7HV", "fXA12HoGXWqwRUyTiUwXKfzE77apeQD5fF", "fXA12JYbnAXakyGU7tLcQAJseupLrdNCwz", "fXA12KLKN5xRLRXnTudaTpg2tGkM371nUJ"}},
 {{"fXB131BfcgewKmbLVuU48Q83zQrgbAMX3S", "fXB132rT6hZE2d2jgwd5kYV8GXg4d2Jk9s", "fXB133tiRoyg1pi7XMtmpxxaSa7p6xYq7a", "fXB134yxP6yVKC8r2i71fjSaNfY87NFUm7", "fXB1351W4UjwDKChHrPqLqedASRaf1qWNC", "fXB13655QjekftEXx3v1yCWEQMV2PmDbCx", "fXB137TAgXmPmLJvenrHHQR5Wd6oq9KiMV", "fXB138HhmUhTSQ6hQupVvQwTzGdHu6Kt7e", "fXB139WyrgJbX3NDmbWxpbQmBBPKtPuUt6", "fXB13ABH16VTZAcMttNDH4VjsDNvYTyVCk", "fXB13BtvrSZmAZDgiusEASUAHorqrRf2cA", "fXB13CvBxxK24JEHNxXq7T3dqn98sHdigF", "fXB13DGA9JZbJ7vTAoJZQAcQAzZZ4mcWSj", "fXB13E2UofFW1beSeT2am2HAEWJxrhgJBs", "fXB13F1amjizyxFpSUxTDa6dnJo1qjyHKU", "fXB13Gb4gJDUmk2Mscf42ZgWMhqbDc3LnA", "fXB13HwGLFwT4CtV1Uh9buVJ7oevoNuLfS", "fXB13Jf6iv3RBeBUjKaHR2KUqJ2RYJQCJB", "fXB13K8En9bn7bQrUKYDHUqrxr1RjRm1fM"},
  {"fXA131qp7WtGZTH3mtpLW14vT3ect7qxSV", "fXA132qS9qaBHBG43Utv5vHsFzY7FCTw1y", "fXA133AReTTKytmwtc3p4tivt49xvu61MM", "fXA134j2SK881oBQhTGjxuw2rCmfHzoP9g", "fXA1358XXsZnTFyUSX2tr2tBNAYpzmWTS9", "fXA136hvzJvZS1LmArXUE2VGjQkvys7xj8", "fXA137kcCWx5SwDzcnZBvkX7Kh1Sv7DRJH", "fXA138qpxWRQs9cDEjMejXD2MxP4fwZkGb", "fXA139w9AerpXzmD16XjHiWqPxJnWcKVKC", "fXA13A7qvq3EdXX2n4FsdxUkcvM2bA5hPT", "fXA13ByvcNJ3kkemnywt9DCaXghrECyL7Z", "fXA13CmyCBHAdSVR8QYCC7o89S4C7VymLC", "fXA13DuSHRu79YaoNQxG3cL5eiurnj6pFB", "fXA13EHENMNfSBhUFRBSvBEg8pxvMYEQ7G", "fXA13Fh32Wc8zq1fxpfga6Cdk56RaPzGqs", "fXA13Gnwb15xL7DAm9w2weH7Y6aV1Nn444", "fXA13HVeY1myAd7KMS2AW8gbZnPPsT9mhd", "fXA13Jrhxb59U5jhGdBitX6oPDoF6jwjvH", "fXA13KthzWYRi1E5fBYtpde5oVDUw6PGGz"}},
 {{"fXB141wP7La1jcD2AeXphTzLxfKkQ1LGPs", "fXB142wTuUhcJ5UbW7oD6e4LtqvBjpSgLg", "fXB1438Z1hvwLemiobjDw49uoqjW2Z5fdv", "fXB144hN8vC9L93yqmTdFjU1bwpGtE4sHg", "fXB145H3JtBd1kYG2FFv2MySYC4bKmKBA7", "fXB146s12hQjT7qn1p6wXuCupP5x2SShhX", "fXB147EjNj8VkbHV3CRq7VCTk85Cy4MV3t", "fXB148rW4cPpwtyfdgrm9ubcddLgxYk9sk", "fXB149rEh6pLC35DrWq5x5vsSAowCmCikC", "fXB14AFvEmEt2KxzpjWrPs5WF5k2VKV18c", "fXB14BDWBcn86dy5sMcFNQUiiT7fezS8Q3", "fXB14CX59JN98fRkzzEHxJr8ytQsEshsL2", "fXB14Dopy6F7qaZc8pMgj6rgt79kZ21YZc", "fXB14EyibofSL8gGfEEmnkLg53Dt3FaZda", "fXB14FNRSq8xDshZGTYwx7cvbzwHDrD4vb", "fXB14GrU1GLUisZFwBdfHfFW65kTS3kbB9", "fXB14HoGTyRzPrra9HRnyzeX7DB8n33UrB", "fXB14Jm2dMzhMAXAfVXRSd8n4U7DZpGDb3", "fXB14K1zPCG2bm7VCPn3Sjb4M33EdsWc73"},
  {"fXA141CV4Y7UgraNbNn8RbRZsF3p1aydwz", "fXA142Uh72u61SpuG6KYnC1FTDPNPi1fJV", "fXA143K7WsZ1y7i9kZ3337wtf6aP7dN4sC", "fXA1442yvQTFFe2EPkawb3fAhcDFHTxQ6u", "fXA145t9zpWfwHaojSiARY7Dvtsa6g7k9y", "fXA146AX9ThfcsGCZgPfHi7N7cbVHpxqRs", "fXA147o1NG3EyZvgTWiMFY4fACs51EPSgG", "fXA148pNp9DhQbbcCJdDwBGhbHXpYaCQYm", "fXA149dKLd5yLceAnt1cMYiSxAVxbXvYzX", "fXA14A6XSba6Lg8yKTKgjXodZSQgJsDpgr", "fXA14BtSpxpHMBwRtTQKyoo1bUu5Y553e1", "fXA14Cw9fke4NCvgM7Wn3MTyD9nxRf664b", "fXA14Dk7e5o3mQNbPovtzDqfqMvZt86wHi", "fXA14EdTmYXmUDJcz9qxKJvFUnSdSh6gHX", "fXA14FwvavyGoKLfhJTWNi4xynMWp6N2Sg", "fXA14Gh48QimKoMRSEzKQ6teWffj1VJbdV", "fXA14HhD7k8MheYbdR3Poek33CbRVQnrVP", "fXA14JpK3gHRWme6ae8VFroxCzJ2XteaUJ", "fXA14K71s85jQQb2JPvPCTb3iojeAVUPjY"}},
 {{"fXB151J2QamKxPBQdzAspVoBB8BjwiEBQb", "fXB152ytQ5saxW53udDwrmPWeGguLzwXdA", "fXB153xbtWG5JqrwEqbcBLc25MWc4Ls3Ma", "fXB154G6kQhWDzoZuUdS5oEZgt4wTcF4MF", "fXB155AXe5vd8FqJkLZk6XTjBDibS5uBFN", "fXB156XLvoLWwsvK4WaArKwpxCz9cQy5tZ", "fXB157kkE8E2uFk1e7dUj9GkujUfy5LEsr", "fXB1584jwTB7S4eo5BX2kUXTpi5mHRGoov", "fXB159JancH4SkrxrraUHepVbjndWcC48Y", "fXB15ARrkTLTbjfLVLz8CAoguxSJo9QFT1", "fXB15BTyD9wj9vFUEsfHiPmJ3Hi7CwCBPM", "fXB15CR7ooG2nEiJkCtzpibAUpqvB3FNNH", "fXB15D7cQ1Xv9f9FHg7zzZ3eHea25kL5H4", "fXB15EhkRmPzBDDAHh7XEJy2TWxtQjqXfX", "fXB15FWT5Fhf7dS4cL3Vrd6pAHpySnyYaX", "fXB15GjdpWM6q35A8wkfe2JKSrqivzfPvv", "fXB15H1jq5yra2e4rjHuoNmeQ6M9rcjkfo", "fXB15Jm8u7rQhGZBkcCB88A1WvBeK67sfZ", "fXB15Kp5MoTntHAHfpUA8a8ttuWBn9KtFh"},
  {"fXA151o2pMTeE2Hqhvp64B31LdhgXgBaJm", "fXA152KR39u36DjP3R9K6qLQfPiVkXWNKZ", "fXA153eed1iXRMsCrSWLfNFfb43AVu5HT8", "fXA154nW1eYpWyCtrJ1ecYL3y2nYwciwg3", "fXA155p9AWrqKt7hNTdyoKcw43S6td6fGs", "fXA1566Fe3tNsqk4ubSEzmuicqEA63L4Xz", "fXA157wjScucsBTfhxg9xqM4hci1rMwPrF", "fXA158pXD9Hh3z39cVvLgMa1tYaco2PjFk", "fXA159gkRG9eeSUz5hUktGK5dVU6f4ThDa", "fXA15ASdFESbKTS1uuVYy6MnceoPWvxF7k", "fXA15BYCpqrmLdcbwrBSEj2xdAGR5SfEvc", "fXA15CGugosaaymo2mncasYYJiJ7TJ8WfY", "fXA15D5mMXTUswSbECxuU6jADAu9s5wBdC", "fXA15EedjKPM91zjuNbfG3u2SQRZbUq1LV", "fXA15Fa1W3KrWgNeHWBZ3ZCKeLd4nMbLfW", "fXA15GWrSpe6mQ3Fmb49jbPkU89Jxi9SXx", "fXA15H3QMk76KrGzmJJMUrZJGrYf2txy6X", "fXA15JEcjfrGmac65SGYimGQXqsGX2d69g", "fXA15KYdwEFU6hBQa2aZMqguTk6o85inKM"}},
 {{"fXB1614dqZaf8cCXKsVAioq85ZMT9ZToQo", "fXB162FFxeMcLXFbFVyskxpjMQiwMsEM1t", "fXB163SwvdHcn9Mh7PcmMWZtf13YD4kSMw", "fXB1647phfGGFn665kF5oK4PRZyhvqhssB", "fXB165ojeeku1NdA8P13ivmAREsFMPSKhY", "fXB166CZeL23t2erb12LgRnFvxh2kpW4x5", "fXB167hBeBrbGPXJdawedANZQKqjZpVDBq", "fXB168M1Y2MkwZndaSGeUuyqK3TKVph6Af", "fXB169a9gofLZhgysz9ThAGezBz4yMd5au", "fXB16A6SyAg6pxxEwstHBMMkcb6UQj4xor", "fXB16BGnHrUugp74CoKeR3iaBZvdgB5g6P", "fXB16C6QH8rdQyMpB6XKKnUQoCWmvnX8Eo", "fXB16DRkrzgweKjAjW3kV4oBUvaxYHAWA6", "fXB16Emj8tE1YWKwXXqtLVCxjAajBzBcXJ", "fXB16FGYdsap88fPNEgraTJaq2xHwzmkAF", "fXB16GzmdtG3FoHAuYftdsjb6eFA6b9Xim", "fXB16Hd5pgdCFhbEDNe91jPmH9XALEp9Xi", "fXB16JnSP9baUaYmV4rW1Cu7at9ZHTRJjF", "fXB16KB6u9vf1syLhWi68jckGF6CVMCQLB"},
  {"fXA161y7buhnR6Dq1XxK1M52mUgXjwH6Si", "fXA1627AVzfZ3L7QF4fE7ktJKLkJ1QSTYu", "fXA163ngGo81xwjKmyZdSebieFnEHfWiEH", "fXA164zdA8t5WfWX8T637Q4ucdW1MSqSmJ", "fXA165CCiTpdQncESfM23FFkEwbDpSCGas", "fXA16689RU3yUjpHQeqfzc5FCA6edbDG7v", "fXA167YuwDp5kfbMDXD9AA8ws2QLZCsCxF", "fXA168gnsCAbNmuTkgv4Ye8SHgeV37wRUs", "fXA169z1YizXadv5tVSrXSzyBp4QZfGGZk", "fXA16AefLHcYKJm1797v2EkkZQ1eshXprs", "fXA16BiW5EfGwdQ59Hfn5Y1kV6o59bPZmp", "fXA16Ccr68CMEqT5jsRUjazQk35LgRu9aS", "fXA16DXbg3ksLXEjVAySsqSrZZfTR8sa7s", "fXA16EYxaejt3iqkuJmXTQSX4QCnoDH7A2", "fXA16FkD6ozoESvbqex7xi2aVDxfgK6SLZ", "fXA16Ge7hzxcd6RZXgRgmH5WVz47AaYhWz", "fXA16HLFZNGzdRZeN5LyJHHMSKmmhkoxHF", "fXA16JdcsTm6RUxwQtQa4o88knn2Haai9p", "fXA16KwujdrrMhkTWkNFdUFyzcrAc3wk8B"}},
 {{"fXB171Ev496vHHo3UasMHfssxkZaDWv431", "fXB172qjGWJfRVTM8LhvdxmmcJdRw2SL8x", "fXB1733mnzCFyN8firtSL3GtBUrTJc8GSs", "fXB174c6WnoZvXZpRSfjEJFeSRrLSVxhtj", "fXB175dxN4BRKKd53GTMFoHyt9U23Bn46S", "fXB176bwQBwBhgfo9fsR2dMMTAk1B3e7r4", "fXB177GLBTN9TKqqeXktftpz5TcK2Y3ni2", "fXB178X9vzBvCJSXGZnTMYXK3H5gYBrZgt", "fXB17984N4aG1qyjLV17hJckVRZJRR17kV", "fXB17AXMAUTLRPNnUDGZ51BvxSmyTbGnJB", "fXB17BsD37qtgsfqURDUX82bZbJcAvQJFt", "fXB17CoMLNLCwFSLtpZmot2dmAs19SWQo2", "fXB17DmnAgfS6JW2JbYmW1SSKXQ7641xTm", "fXB17E3X7M4HjYmoYxKRoD3fvyYhN9MTsV", "fXB17Fx4XMgB4Rfzg4mGJJdHp2FKWM4xqJ", "fXB17GstDtwZxBCkU4uRF8gHbJNQYcU7UP", "fXB17HK2SYCSZEUqPSQyuegbQZjaFLupCi", "fXB17JfgwxkbJBnWH2yCMdFp3ZshBU9R5Q", "fXB17KkFQgMm2voLR1m2mK6FRNccq2Jzjd"},
  {"fXA171kM99K5dEHzjWLGNqDaKmWKip1M7r", "fXA172MXm3iqgsSy9LuUCXVBa328oBicW4", "fXA173CLp2RoEJJy2Z1oYeHxEBpbeYQV96", "fXA174vf7i9brQpZXa2Y6bCgpe42ihNXQa", "fXA175CpQ4RZRYqWfiSR9Z3FZma1u7s1NJ", "fXA176BtWWLBCBH3L8Ckk56DUTRVgQiCVn", "fXA177b2pyuDDQBLVPVLNoyg3Y8EecHHsP", "fXA178mKgSwhnfh7WZGvEfNfXDqBNRHuLX", "fXA179tMVvytVhXQTSsksKXSjDHdSfeKHK", "fXA17AwWguY9E4Z3C9dg7DgA8nMAsC6rJ1", "fXA17BVsZdsqTShCN6hyXiEXbmxkF3UHjy", "fXA17CFtoYtNmcp7MeqHCoL3vneksLz3vt", "fXA17DSbMUWASszkSNHUL1C1NCUA6v7Ncw", "fXA17E856Qj7woqGn9QatUMs1VdQReBEPN", "fXA17FGTNbm82NmBPGxMHcB2qCwmAmjZkx", "fXA17GqEqDmBm4cdcsC2caRfkhoV53XMGt", "fXA17HZj5SphFXX8HN7Rekf749AtJRJgPK", "fXA17JMUUkWJ7Km8K1CFmytbbfKAwnn5Fi", "fXA17KgCgXf8rPtqCE3MDkKh9D4kxKq9YK"}}},
// even trading round, pair 2
{{{"fXB211ZTx8jy1ZtSoRaeTF7Wyrj1znff5i", "fXB212WfdNgvJxRK9oXM5FDqQw64m1eP9i", "fXB213yNrjB5pvcXL4UsANuts51eJ7UHTq", "fXB214yPcXVpoHahD2BgN8CF7yvCXW5jzv", "fXB215okcmEqfggWZfYvtkYkJ9uo4Sbzi1", "fXB216v95M8mNm11zhgduCLDBWgpLfQhpK", "fXB217y58S8e8gqE5NvehgV2ACq1SdMMU2", "fXB218qghxREcg2Aui6UVHJFaB4ZwTwC62", "fXB219jLnZXvetsqaPcaBVtBJpPr9Nvyfe", "fXB21Ate54uouiJqkzDk9ZJEqycFXnyov6", "fXB21Bnj41C9oVNPvubp7cLeYfkPVeHTzc", "fXB21CsZ2AzC3xB95yBLX5dHwMLT3SggP5", "fXB21Dyavh9bvRgpvRNJXGvwBW9bwK2Bdh", "fXB21EcjoVScyNHFuVfyM6VHJ945TMemdu", "fXB21F2XYtr7nFqfzZhkuRbkZaafVEy6tW", "fXB21GbYYRAybs1ydm8Eb3NyeukvGnt8cg", "fXB21HZa6kDJJdC8uieoKBSF1ByCPdQSMh", "fXB21JYbSWbbYeFB6yWjpfcGajnEcdCpro", "fXB21K789UEf238d7Rk9vSREpXB7KNVRBr"},
  {"fXA211t1FWZHYGAR7i1Eqb5E9x16SEhp9U", "fXA212vNHutVXx6sqV77q6u99LWJuiDEpu", "fXA213JqbijLGjCbUXLuQWzzaLWU3ypNk7", "fXA214VUt1U6JjE8PfaQyiNdoLd2bSvFBs", "fXA215CEMvxQXDdWVrF8xYLPt3UyViAJa8", "fXA216iERGfAd8YECezkTXRTMPeoWdXqgh", "fXA217QX7AdYtVhyB5bN69p4ofb5sfi9Q2", "fXA218YGaRbx1EQdkzoHwDpbQVHpazJpxh", "fXA219gQcJKUe8WjJMBszmqLeKgBi9G2Uw", "fXA21AM1KMnR8VPYLMWcHz8Sk3F6RUC9Gy", "fXA21BLupch3RUUWCuuiBrpgumcascPHX1", "fXA21CP9PdnGjho9KDe7shE1BLV945ahJB", "fXA21DHBSFPTop7yCkij9Ng8d8SvD6rRJK", "fXA21ERfWGmi84baRP9iLAaiCFehbQCs9M", "fXA21FQ5bXxzXSNxy3dRt2qRKNT3b5rC5Q", "fXA21G6MH6VvtzAewNTrPFsLraQhD5gD4B", "fXA21HWWmimXC4czfthQCoQjpMZYvX8S6j", "fXA21JrdPAA1p54tnxRJD96UAbBrmz59vK", "fXA21KgTirUg7PTAmucRHHSj8gh3CzcrFZ"}},
 {{"fXB221BWmCiATq6xM696DcJ6BChGx62b5u", "fXB222FjJMkwtXzoy7BCfJH4eUJnYU4r9x", "fXB223oybQCM6HQs1MYy8RA3aapbQYJc1m", "fXB224oScQdDmRRvhXVJFXozmfuzQDoZ6a", "fXB225ovtCYXEfHFsjhHvX71kmdXmn63F4", "fXB226L61JjKg7sXwCXVktSFv5SNv8JfLg", "fXB227P7gj9MT6qzjrx9sA1YtbtCAQmsir", "fXB228uY7WeZkTtojFZof7aLfwqXvBif9A", "fXB229cxZozNU4YhyzZmsSJWNkZtuuH4Lo", "fXB22A7YRuDFrtzYofhSka9wdhFHEPK56i", "fXB22BxkpN322yYv37pM9xGBEMoaRaRqDt", "fXB22CmkwHR5CGYgfM13KfBT5Li1v3te8q", "fXB22DhJsKsA3XmTHZBH7tyZjEEgZeRZUW", "fXB22ENhJKZbaJToR4u7AUweaeNpXaS8EK", "fXB22FYPuTLJa1zGAwhxXvweXL6TJyHxQY", "fXB22GhBhUgeCqgYmb2PEW72nZq8cXgS7v", "fXB22H681QMcGsv5bLNCCuBduGH3hJEPvA", "fXB22JRwo9aPvqoFMS9FQWdxnJUbZs2Ktq", "fXB22K5Eo3fSpNsBKjsqg79MneDaQ3KPZY"},
  {"fXA2215vxVmuFsh7qmpVsEooi4kzzWiiQy", "fXA222DFYojBgqPU1jHarnDGF2dnWG7KRq", "fXA223b1tF6W2SNttEtmx6kFys69Jb8zty", "fXA2242vwTgKEmqgjGXPAuqTKLeBDwFcCX", "fXA225Y4PQyqKMiu3rDT3UCXf5eifXuVxh", "fXA226hsAy7ui72SmNzpfdc99KzGznNRdB", "fXA227hZWRgWVSgfkJfT3Uqon1fWHWcz5g", "fXA228Es3uu14f6fmWDPZXKv95yELzcxgW", "fXA2294Ga4MRdoKDogDHDYvCzeDxjedPxX", "fXA22AUjgKcnnophVY19769YRgEadft3Hp", "fXA22B3SMVXwufvSVSaJKqMnq29TFUH95u", "fXA22Cawv2uS9B3W7MBmsBg9qASj9gNcur", "fXA22Dj85SajDgBfJ4DnpA2EqXJxtnrEqA", "fXA22EzCaphr9HWSSKTmrPvgS5aGFtJYRq", "fXA22FGE6CCjzYAm5CQf3mDQLoDhZxBjow", "fXA22Gp78ctwKBPJj7ys1sKr4yATpfY5K5", "fXA22HkE95UrjMbre1WDRkJzhtaBzPha9w", "fXA22J7P1eauZuaas6E4LgsespiDod8Csg", "fXA22KdKSHhSgtDRnvzsxH15xM4c3nZYQq"}},
 {{"fXB2318ojbeQeYjtjfZDGhfLNAtk5fuktr", "fXB232EFsKxUnxfnwxsmrY5qM2id235mjr", "fXB233SWnr1LsV7kTNEBXmt4EEHMYETyZU", "fXB234fzLUZYMr6SUG5qZXWRcNns12SSCP", "fXB235x12kntQFaqBP3MB2RxN3jEbLF5u6", "fXB236QsD8U1kisu2fEpMz4fQZYi6i3ye7", "fXB237frby9fcsNTfHxybUQML8qSeLsRkr", "fXB238QjZt8yds5AXfy6JyQJUc2W76GTPU", "fXB239Po5E3TpsCppFXbNa4fLPFGAj3ZWD", "fXB23AMiJJVAStWTkme1Y8FexW9njPz6FU", "fXB23BSTvSK9BQH3czZwhGJf1xi9gUiQEu", "fXB23CHAQd46dub7vJQnPwizY1g2rLSnNs", "fXB23DfKZtvBUWDwTEbznF2hawX25RqSLA", "fXB23ESJjn878swqKjh5JVCnhCE3etgPps", "fXB23Fu72pKUi3m2i3Ey7vNnM4uuRCKvf9", "fXB23GjgtgYD1Giv3czgRvJJMmdKDh4fN1", "fXB23HDpymEBhN1kNgRVDnMyVjyKTYSWRG", "fXB23JaYPSFgtv5Wt9CAqaA5rVn9wRQX6b", "fXB23Kj7guuPSzzWHxNEvZv7pDtmmpDJfJ"},
  {"fXA231oSqR6Ckhvdt7LjSvzf78ZdXCZt3U", "fXA232eWKCmNYuWbfAki9s1NqGu3hdzYSf", "fXA233tgsVgSY4xmm4BLaU7qhz8s93qLsh", "fXA234idtYxNH9XTw8hYNeJnuTsfNWb6Ni", "fXA235KqnWTxqMtxtXfWniWUnn5iPYx2nt", "fXA236vFG63gtCMTscqBaxqk7aM9qRTS42", "fXA237NV7ACnBsk4ijZ9kB7DQcaHysMvjT", "fXA238Q5xu2wtXTPCBFFcpT9Bxv4LLEgph", "fXA239fZTWC1QoKnpYKWuYdJ9UZ8QrLrh7", "fXA23APuyM5XzFRpkSJdf2pTRbbUPSf6FU", "fXA23Bgp1UhmxbuSKM2mbKbo6us2L6qFNU", "fXA23CNWbNQogoEYKSWZMb2BNHRQbt62Xu", "fXA23D2c7U2xWwwJPrmjWtgsUUk8jks62A", "fXA23ET24pcZ19maMuaRYjAmRP6aVrJffN", "fXA23F7rQgFht2LrcTcEzCcEYiZnoqk7kC", "fXA23GsEBGGbpGN6GKAXR4NDkP2HZr6p7A", "fXA23H96KpJR32co25hHNrMdTHTSiYZ3At", "fXA23JuUNPqHY2HpG7aseeViVPns6u7hRt", "fXA23KbqJQ7vxkEtnSWis87b1ShchSE53B"}},
 {{"fXB2412JCjcp6hN7nYi3GNJDAxzbLoJJHo", "fXB2429KywXjwVdn5WPbRe4Vwiqv1xtpCP", "fXB243Lf5xtieF5NM7LCy9U2QSon78Gdo1", "fXB244aCZoyDYHWxxYKu1HZrNxRXHJQ8pg", "fXB2457CEpWoieetCRwXGZfdMuqFQ68wEq", "fXB246uYEL9y6koV5HPuMP63YPYL3iCPsR", "fXB247js8qWj1LgCHa8gPVrH9f3rEdwVdD", "fXB248iuCG2KC8VFL6n6HgfeRauDY5Jmaw", "fXB249CWd2GWqHGwzFRPdhUhecC9VknQ1k", "fXB24A41vP4EU7t39JJbio6nDpn8Uv73RK", "fXB24B6CDZ9g6LDUPLJ3mFM2CGwrLH85zL", "fXB24CYbcVb5fsPppD6PqgWEP53UxqLWVZ", "fXB24DD1UvVhZLGSePMmSCaPuZXMxi9Xf1", "fXB24E7uDYPuuhaaeHu9csQsLTEyckFWt4", "fXB24FxQyeMJhZiWfi8uHMYLHcEVTaRn76", "fXB24GwLbcM4cYdhAGtXeMDHMmLJNtFZbo", "fXB24HxoiPsLwiAba9P1nJwDxgguvDZZDc", "fXB24JtNy8ifcvx1fQTbWez3F82ZpmZmBS", "fXB24KYxWLQ4fJ6QxAZKmXR6ZD3ujgLFoT"},
  {"fXA2416KeQWUVVAxc7Et3krKFhx6btAnti", "fXA242eF3qNyn5fjtpkdGF9SGJdyPcmmQa", "fXA243XR4H84LEHdoP4gEbz3Po9NrMWSLg", "fXA244fW9VX5STY6AhZamBgeMJyrZToZJK", "fXA245JifeDv3UuVxZTDpuWbUQ3WjLnJSR", "fXA246ahbgQFwyjvPMh4gVwHdxnGHnrfL6", "fXA2475enKVXN1x8FAvQwNi974pGkkcg9v", "fXA248qU3DJvZCXaRj83ohAzn8qAMzDvTf", "fXA2499nciVTRw6ZnjSzK1Yf9EGs2CukwS", "fXA24AW9gjvzP5Zz2SDqAhaHJ2uAkqQ1Fz", "fXA24BTt5QXYxTHNmNzYwixS5Jf3kmQGz6", "fXA24CsUK8eEf4WdHoHr4A6RzLGQivnhN6", "fXA24DTVnuyMyMvtTwhAPcJhXJJ5dYcRyM", "fXA24EP8XofhPfkNFiRNBk9xusEmAe9XEW", "fXA24FhVmCHrLWnET72kVWPMFYB5h72caU", "fXA24G59YCrkzBi54j3uztkuhcHhTWF6Ls", "fXA24H2XyXs31DRX1ktS6PAEpsoHnT1bjn", "fXA24JbPAFoPwhZpPpWs6F76oWQN69jfMm", "fXA24KhM1tt33MXrLEoXCcqKfdYwTVRqyH"}},
 {{"fXB251zc7qjpzmeBvTPHSAwy9vwWERxMBU", "fXB252bHdcsnVUFL1zkrJQDMAb8NptWT2t", "fXB253Nksytfy5k7Xkbf4eSbQZ2SSvJYb8", "fXB254nP83SadeNQaPXjeam5YtbarkM5Xc", "fXB255nvRcucns6FDYMnfmkG2peBuiY4QT", "fXB2562v2yxvKfpgL1n6PWKKuiotmrVS9s", "fXB257yuvqD1uBYeeuxzShsHusYmJZ578o", "fXB258i39DuPUWrPpVkcJ28faKX9dKMpJP", "fXB259kS4932nC41wgg9Wq9gvyYMCbZ2bk", "fXB25AbPvjxqzuGpfxUqcNCN3k1Hm6aE2s", "fXB25BtMVC9oCNbGMkyQ1TzV923Qe84a4P", "fXB25CPD2PwjA2QsWxLfYgnGooAmsD4F6D", "fXB25DdvWpLy8DutWwv8Tpn8KLMBX75qri", "fXB25E849dC9tz8rPX5xZ66EtcLocPAboW", "fXB25FVdB4xHvjxJX7CQZnUros9xGDedVM", "fXB25G6G7uD4o6feurc8727Z4eC6ivsrjo", "fXB25HC1aNgLbsnm56ifG8T5aXb3b9gmQb", "fXB25JfJ3cR4shQkUUgpkR6QSMfXynPoj7", "fXB25K2Bs8XRFSv8G6gcDTqJyvEVoEZnid"},
  {"fXA251Ndn6Jfg98PPzrWJsAFaYJPV4G5Jz", "fXA252YWZaDJ2gweduv5ZVgyMQgYnRTsC3", "fXA2533cJHcoMrJsLmPC6QUEHjA1caDzt9", "fXA254U26bVDybLnet7dmM48AWZRKGBdTE", "fXA2558KaYtA9YWKzQDALqhGtgbJhWmHZk", "fXA256CPwnJrSwSpRWx85MhuVe1ZkaSktY", "fXA257Rp1CRBMvPvMjUKJaxDPgjnQtmEdJ", "fXA258CSu792Bibc4AnA786sdeoVwGkhH1", "fXA259TW6yAHJVgurmJtpsmuR87ecH4Dep", "fXA25A4kw9iHDV5RdT3Pdn6ZgsgmePLotq", "fXA25B4oSvxvW12x1dKAXvpKMCUcAoAYKG", "fXA25Cvdf2LvJej2CPocDV46VuacszY4Xb", "fXA25DCgvEjo6AoBQTFfNDjUgpjGmVBTTE", "fXA25ETo2JoURQzqRf7aCrpm6eWvvS5KyG", "fXA25F4GzkWrrmBtyYBLnFQUqqivN4FwwJ", "fXA25GXBFKM4GPnQdGz9wfmoJ6dHQXkhpj", "fXA25H8C3VcgrwPTV5eT9rDmmCbVCu5Qxa", "fXA25JaYLm2hT67wo83VfT3hr1GRNc8ECp", "fXA25KXYXP9jdqVMb8KKqU5H12sqFEcwWb"}},
 {{"fXB261YyEUh4tvoajQbeNnZoBEYZmtw3Dr", "fXB262B93WYrwMYS13tRvmnJxeYgAfaGY9", "fXB263qz1sPqeW4dt7wbXbUPpoBc3vW1FW", "fXB264x266QoXv9pet37huMJbeYECYYzHs", "fXB265mCjPwRkknmKSv3ViwZn6qjW7V9mm", "fXB266Cjh24Rzp5vDVWtSD9SdwT4VaA4qb", "fXB267LwM8ojif4mKqQP4zZo3GQ6CFmiag", "fXB268cmrfm2worNc16kmcuzDLU3YVUeTs", "fXB269BEZcwXdhJ2FPESvVaG7iM7XRqwg5", "fXB26AwPihSNM22SuUNc9HgD1dyRigYhPk", "fXB26Bd6vLvDYv7ggFxERXX9gGC5C8QH6S", "fXB26CtqVZ5dUXgtKFAJZbxKqicExZCg2Z", "fXB26DMJGGE9XDQ3wtuMde2SuCfXD5ai3W", "fXB26EkaD9iPXdo9c2SoeGu7n6hQVTqaND", "fXB26FsELJ1KRXjWHyQvubWYvo6bkcCapd", "fXB26G13DNPQ5Qhy19WwapBhQKgajmkiaF", "fXB26HP593D5gjL85hAudgsxmEAzWpLhzd", "fXB26JpqFBtLY3axf6KN5CxxzJ7z3uMRZB", "fXB26KDR9fFfN7BfPR8JRDu8JELSP2PG8m"},
  {"fXA2613Rs5nsuoesGjchAy71LYevJgv51s", "fXA262pdLbQXeHXjZQVfxsPHNJU5uv6iD3", "fXA263y88U1rTKpm3hBchUJmMwPajDGccf", "fXA264FxCmaFNBogrxgpz1mjNS9y2pT3CU", "fXA265jrVZR7hzfCBtZUd4nh2xWofsxRif", "fXA266whT5u8cg6dZ5vuetS6YXiPvuYAUq", "fXA2672y4WPga1Lmkr1XJkoWhmPccuPe7K", "fXA2688dVzhmaQGCxq4wSD7wQserEgdHsz", "fXA269bYww7fEhYojPyamqrn19DEQaeXA9", "fXA26AY2j9G6C5mhFZEh17wnspAGueXn8c", "fXA26BiwvPbpxYzfkbNP9bX2bkczokzwGt", "fXA26C6162bQMQSeSZEg1xY63ra98T3bC1", "fXA26Dy6YEdbd2BpYZJ8BEu1StZA92S83P", "fXA26EPbbzCX6Z5VFDk1LqBaMyLhsRLcb9", "fXA26FdjzjLKD7mDz6WRDJZ4ZzZehskQRc", "fXA26Gucs85iG1X6uVRNnw2Ejpr1Hj7Vid", "fXA26HoKNBPBbXR2ooZ4MsybCQYcXxMvZ5", "fXA26JPydqLRA41aYRJcCY2Kf9rT94HJqo", "fXA26KyzXy8YvdqnabEnBgTfzoMbKiP5ZZ"}},
 {{"fXB271dnQcMyovvdQU7Jjpzdxo8L9b4onC", "fXB272LJ3SyXhNaTAxfZ431e1WqLNk2fku", "fXB273HW67jTguj55mzNBVnD5ES98uoLh7", "fXB274yNqoHfTcGRTUhRdLqRqbpKj8X3cq", "fXB275kXkD7XVr1JxdhCMiBX1cB6efhNmk", "fXB276GTXpy3x6J1yXbCKvo3Von5tSWswm", "fXB277JYJ54ATfmKeS367gSXxjhRqdEJms", "fXB278UVSxCjmoJNPFkb4AvMmfphrrm6SB", "fXB279LDNcapzMGJAHpozYabMe781SqByh", "fXB27AyVLtNisaiBxV1PjcE9EU5utGALvF", "fXB27B2BPe7vkjVZVgCzDN3R1vLwERy5xT", "fXB27CGkTp8JBwdySRD59qsefxwVyjkF6o", "fXB27DwiWmKte48CJQBR6YG8Ew5ZFNf7NF", "fXB27E1EqE8PjXvhssjixgpaEsTY1fK65Y", "fXB27FbFPDAaMU314V3F6e3HYgF1QhGAua", "fXB27Gmk3ozSzCnQ9mPjFNchECpZvYVoc7", "fXB27HMwyDKfRG1hfbzQzDFEr36Lf6jZaE", "fXB27Jq5dYgz5Wt4bpoVEHYKcNKR2EHVwR", "fXB27KzDTsgRBBkQgQ1UigqmFKFCeNyYDB"},
  {"fXA271AU6Hz8tQzuyQ1Qzt6nWE6XDs83Ja", "fXA27247rsujkCyjfm2sWErBpEFHmKYbn3", "fXA273dYc399mKHrxLYk5FeU63uA4bZywR", "fXA274uQ2Hqu8tVmJBRn2aefJjKWYiCgEb", "fXA275foKkct2KMZiG3qERezcPczdnrnTb", "fXA276SjBccaeGkDEU7mb7EXWEbBJQDUTo", "fXA277yjGfjYLkg4kK8v5H6qprv2iGNM9i", "fXA278YTdrswhVnAWmxzZ3X7jRGxPPo2hj", "fXA279v1xWimKPqfuxCh79by2YwjRV1WbD", "fXA27A1a5pnr4rNiVRGAr4chr65XS6BRrq", "fXA27BmYDSZDMctAXYfQF7snjQNQWT2krY", "fXA27CWFD18MmAgp9xPpRDgrTRHTVwCyaH", "fXA27DttV8NqyuEKBba4JyQmiGxRtBaMNu", "fXA27EgJK9BiEaCCDGZ69PDQJ2HPQwrybq", "fXA27FQP8yy8znJCChPUdfjMd6Dmv9MV9f", "fXA27G4yMVpUo8hv5dfAPfZsR6AJFVpYb2", "fXA27HjaiXn9WqKS6ye3StpyuL6ShqP14m", "fXA27JfMJNH9QYGTwwECtNRaJjgks15tMk", "fXA27KU2FuzhhagnWwWRrvPxW3HKQ2q13C"}}},
// even trading round, pair 3
{{{"fXB311WqWepwmzRTfhGoTemPoF6aAufyn8", "fXB312BttWuwwqN5CPdey3nJdLtMFNtaYt", "fXB3136y6RF5Ny62ceQr9wUrVwZvzHNAd9", "fXB314a7djzjkCsQH4qMEyUd6GFUoRdGs7", "fXB315qZVwKx4i84bFwgp3RUW1BK1yo87f", "fXB316oGtp855Wk5QRqdSVG4mJKiktWSEn", "fXB317csG71DJLoFtTZCm1ZVYyts5zWt79", "fXB31862pc82jzUDGokQp1jh481kLfiPrm", "fXB3199dSLwKUBuNWMzjheiPxtX9ybU81d", "fXB31ASKr1LRGxLRtEr95VekRC4bKF27iL", "fXB31B17QftSQpyS1kZR8QCkqU7cihx8RJ", "fXB31Co5gASHHKgWEHHsVsDhTzAhkXAv14", "fXB31DhyUbp6e4Yjp8qfDHAywx1K1KJdTX", "fXB31EDNHfk7GDX78gxckzHSbmm9ApBU7r", "fXB31FvB4FMEhGcbwELb8HJ5fy1f6nTTkX", "fXB31GzeJb7sK9TUMmJ9G8UyARns6BY3i7", "fXB31Hr2aqxtB9LRPNw2uqobXMLdgbhzM5", "fXB31JXnigtSqiFZCRm8Lbzvp3BvsoCXBj", "fXB31KuoXUdMtVxaqGonRft6hkYL1ptgQ9"},
  {"fXA311yw4KzGYdgtCEVgazz5RhzCGRYhXj", "fXA312L4uwDxHHwxzN28toBVEQXW81GVKr", "fXA313XvWiki9Cx93Vav2EzmqfW5tmLY4Q", "fXA314MnaZNTakbjy836thA7shxMiPyr9y", "fXA3158HtMgrgjZwXGAKRkCKowtixT9pjU", "fXA316k9X2a73edYAuo7w4BRNbEMn32NkE", "fXA317x8BhsAd7y4dohxhJUUoL9VaK5eMB", "fXA318qZ7d7JKRMrviscoueyEbnPQGAuLS", "fXA319kGVBrzp2V2bWwT7rdh4se5p7CtYV", "fXA31AB6uhSDy62f5cZiTB2H3y2DVvy3fA", "fXA31ByYTWgCq875quPEoaRyie882XLYY9", "fXA31CYoruwBxFgFTZaigzuYDKsp5ZADqD", "fXA31DUvMU6J9yCiaER5wAX4uQqXFikNNN", "fXA31EtJHTV9tSAaRh1UwbKKCLZJrPw22k", "fXA31FkcqovGaKLzeaUVctQGZGhjEnKaxb", "fXA31GUoZxgXTAK6qwjFjzuAHqbdYys1Nn", "fXA31H7JvymRGLjqssJsP9s8kwostCbYS7", "fXA31Jkw8FjxBZAyoJRgguVfaUDwUW8x8X", "fXA31Kp4f8DC1S7j5RA5aqJJe6rx4hyT5q"}},
 {{"fXB321Bsh1bpbVan3bxLFbdF3ZAk5gkqMo", "fXB322JoHCc3E8jFWvBSuHhmSzC76K1Eqx", "fXB323JTzMZbj3ToZTnb9FUvGBxL5jNqvm", "fXB324nW2vbSdvjsqkCoPG5J9KALYgRmPv", "fXB325MZmfY8mVXExy3LtegC41A8QaPg8B", "fXB326KS7xt7x8sHpraLutjpLWdbm3P1Vr", "fXB327Lai9kiJHmHnuaPiVLxSjFuYfyAQY", "fXB328Gt7zVVUdvxC4sB73yosQwkRpPRsN", "fXB329gdpWjPstUpXK273GdhmTYXJrLvqp", "fXB32A7PX33U6537EdgHcBSp6N4b7VpGNX", "fXB32BENUGA1kjRMVX8L4EDPABjzMw4h2R", "fXB32Cb1kd2A7cPSZyFEQWxoQgTBmgFXSm", "fXB32DGE9dhj7yDGCW8esrhJPQUs9CJ6j6", "fXB32E2Zfh1a2etHRwW9fD8ixXnduBpvEs", "fXB32FA8YD1fyrNUhxCd5boH2LgA9c5yAE", "fXB32GSEfGJizkib8Rqj3QeFtUnKDWRKck", "fXB32H65LdPrkUcEBbcX3RKwuHjfdgp2NR", "fXB32JWq12rcjnapdhYiT4a1qnJxbgUnKz", "fXB32KgVAkYEkyShFnNDTKxE5dSU3fLpGo"},
  {"fXA3219CjNGLbJRVFzFxyA8mSrFDZZkFy3", "fXA322swufZuLAr5BHKEeCespbB1Jx4KRE", "fXA323LffDHUmj55rKitCPPz9vAyAB1Nft", "fXA324PiAywE9Ybt7cP94i99ZQtw1TRao2", "fXA325wkoGzAgoMPNoFb2FsSbKasaxjwb1", "fXA326eWNwXLEHGHULzxxrkcpwb4Ev51yK", "fXA327wDzFbzk6d1anL8uFoaeGTPzYU7sB", "fXA328qju65uvAFAy237MTwrdkYrPPZn7T", "fXA329pVWUjz7DHe53E7wtphr59MmbqHNR", "fXA32ASL4hu7ZYTYgaKkYZQB8f45u2VrKB", "fXA32BNwnamDBeL3Vijc4R6gqJtTEo5kGn", "fXA32CL7GWoaBicP3DXb7qZEen6GfS9dQb", "fXA32DVHtym3tcNBVw8Q6ix1XZLxeUbbkD", "fXA32EXA85n44eV76UsxXwJLkc1sRezLSA", "fXA32FZd3XGYEzH6Eg2RmQSWpg2CevoehJ", "fXA32GaWrpoUKuHzE6opY3M3f4LoKjA2oX", "fXA32HxekQFeyjkGmo9bhdUt3drSuv5869", "fXA32JWYgdD2HGmHWSrqznf7EGWEpsy9qH", "fXA32KACbCBVEE3UmybRSY93ps9jKgGhGZ"}},
 {{"fXB331bWD1CNM7grxwPzzoQWuZy28vPeEz", "fXB332GG4RJzmKHj5uzYsjqmJPCCeT6k5Z", "fXB333wHyqd4j1eEGRS6d4oHJJ5JjWFEvU", "fXB334DGpduo2Ebb5NC8q8KKHpq5dsbYMn", "fXB335Bq4Niofx42978WCk9Mnagze9YRta", "fXB336Hkxcq88dVYAa8kmmGtcKS9QG2tQy", "fXB337tiU3zMLdxk15y269pbXit8d5mRwZ", "fXB3387TKoKZaSHnJ8ueRSJZwPSQLFzsu7", "fXB3392AwS7hXrr5A8iv8Rowg89QbnxkW5", "fXB33AkhL5mFbW68YSMu6ZJSoGAg9y3qVJ", "fXB33BhfFCNcpm1PUdWMLUXMno5R3fZjxY", "fXB33CL5M4MhToyigMKnckmpqf11X7NQNW", "fXB33DgEuhi7kFshrvsKRhceLZfscb24DP", "fXB33EiErUD9GjB5QBgyR8s3N9KzvENMGo", "fXB33FXGjS3tKxtPAaCLfGGVo2Y5LgtYPr", "fXB33GU7wzfYhry4XVWXcH623aJWs35eNZ", "fXB33HjWpTB7MccrTDf819sdTnG3iUYJt9", "fXB33JkYJMRkcjUCWVRiyNVE5gpS8KCQv9", "fXB33KcbpYgQR7QgqgDzGpLxYrZvazn7n4"},
  {"fXA331EVtoc5j8BYvQXrVFsDgePaUBMvXo", "fXA332JuSjMsvWLo7XUCwKTdUciSKgZJ2C", "fXA3334gbaPPCNwRCFJBL7jaXAEyXyETMh", "fXA3341uyEX7eHh6YutLT5JTNXNBVohe1F", "fXA335CCvVEiYgZAMcW14PDzjYFAUcRaXn", "fXA336jRYydR1h4ZmrZedfaFCp8AnG8zP6", "fXA3374umTMY5yCg1Ep3EL3MnMXJuhV7Pk", "fXA338LVj8uK5tyoFxwoHSBXLKSUagVMPA", "fXA339E1SdfNoaQMxEtbBMmMAJqXpiwgpd", "fXA33ACuStMEBiVYHU5SmQr4bqKeWii1EM", "fXA33B9TngLtgAZQFC94UsEog8wuQgRwPH", "fXA33CkopfdbCQNJ4DZAoRjr342kdKtS4C", "fXA33DDUP4YA67K3JQhSsJTQaWXHvc5tuH", "fXA33ECJuFYjYjzj46KjZCaNMuSKP3FJub", "fXA33FHX7MM6VLVHcZpecBn3XviCcu5ywX", "fXA33GZ4NZimn9DF39PFu7sRKmrfr2shUh", "fXA33HxBNfg1ukswFanDCFRdfoLqW326GU", "fXA33JHELhFgzvKUDTprM4gGNirjeDXKG9", "fXA33Kd911ij512wdCbJFnYmGrcw13b4SW"}},
 {{"fXB341RdaJ6FqaVsvhAHGM9DWaKVbCPQwQ", "fXB342Ea3gJQ7HiuwQRBtEXoWS5qDVmz7W", "fXB343VXuZnTeWigHVewE6Fo5wFjksreVz", "fXB344s5j8a8asSUREwJwsSJ3X2mjWVhzR", "fXB345GJ8oaxvfx5b2HcLDdd8y5Ln2idzc", "fXB3469sYCy6NvthgMtSb6dz1iia7qkoDw", "fXB347fwvNVWv2WUm4g7HQsx8HysRL9MEY", "fXB348Uv34CgZD2U5395RsnSWKrC55wRN8", "fXB349gysEnXvgAjR1rnqj4ZF5hpp957nn", "fXB34A3jZTGuHvwFMWrxxDTqhQU9XGk7kx", "fXB34B73n3RytotkKB9XDfk4YxfegmMDBh", "fXB34CaD3KLnjSbMbfxkVUtV24aPSsZScZ", "fXB34DCyihJro4inYti4Zqwm2KAVNE4vEm", "fXB34EUEJtmWNCU7fuL883BAKSbQWdt7rK", "fXB34FDfTFr94taxsVrSHHTuXD3FENBfVs", "fXB34GVM5KgaGwVVEXScqL4fqLXQkfHuaH", "fXB34HVnr4msQEQbuk3oRx8EToKsKzuKCv", "fXB34J4eEsoHBUtrQ41HA1szChQeuGxGGH", "fXB34Kmvgks6rpX5sAW8f1eZHVjNKG7Hib"},
  {"fXA341sXquCC6wypq3wGBm6C9GHDZtf9ak", "fXA342tP6ztdnokurf3cMTm6Vaa6KryJef", "fXA3433eNpovsXUkQAPrjp1hYAWDvHLMgo", "fXA344qWgh6e6xSL2xHjAdbgMEC4hxAgLB", "fXA345D6sP2MikorxQe47rLsY1FePN1Dcy", "fXA3462v5tkQTc5ZqbwzLpJJ6ZWv1WZGMc", "fXA34715rGTTYUuV9FvyeUw2u6yWBaHcEK", "fXA348MHuzvTwUoWRBHmrnc7s1YCQ2jmYS", "fXA349feJhDwfBiuQTe653tPeyzzMZHoor", "fXA34AmZxzen9GJp3PhS69G3TKJuBsFkom", "fXA34BmBXHrfvQ24W3BPFiZpgGizTUjuVN", "fXA34Cpg3sN6CzsRH3Xhx977Esp42LdzxX", "fXA34DnWuofG7MnmZdJSaHf8ECxbnnpsbU", "fXA34EoN6gY2zadMRARvof3rC81g47Rfvx", "fXA34FQPqcw5wJTkqJnrBgmgHkjwBNUyJM", "fXA34G3UdQhHdojWccQ2x2iWDPemz7ZPZS", "fXA34HYacJ5oZfR5yMSzRjtLao62vqtjSr", "fXA34JFPPimPS7wdbj3nfNGNH9v6niTvch", "fXA34KhwGejisEhuz1eyTzA3BPUhf82eUW"}},
 {{"fXB351Mu1RVXfvAGGegTcfNeFPh1B8rbSF", "fXB352rf4CKqFy3sEHdToVvVkgryNGHdQt", "fXB3536DA4Ab7dNHqekZ8qaDwGGpmoenTu", "fXB3548eKhNzmBn8AKTc4m5HV72gqSWs2D", "fXB355UPQHBQqvwhgFbVwtPSZRqb9XJxPy", "fXB356YAmQDbSuTuSQgnHeFH61jHTQCpSu", "fXB357v43gmZPKNGc1yG3WocAfAsPt7t1a", "fXB358vfG1GrCdjGzG57A57dZc7D6Gzvsn", "fXB359vPH11NDpRsx8QnT1fwY3VVtiXiWh", "fXB35ARteQ6xYXrmxqBbWBTgdCR2uZxrAX", "fXB35BTiZk5kmRLadfEyQHz2nmchS1wSuL", "fXB35CkCEDcKp619AK8ux67vvSLmzSW8c7", "fXB35DnpzX9Zd9Vf7VdnyTLDFe9tqAsE2g", "fXB35EVnGwuZushhB5RsukEE9NciwMvVTv", "fXB35FGsP8Bu5BtkMEhAeZWSCs3owv1Df6", "fXB35GzX3PQTt1izYSPZE67cfXavWLPdpm", "fXB35Hn6rcgRYZywf1JxP6XwuRaRrEnsHE", "fXB35JULnpiVjn7sqMpVasuAZD6et8Lxjc", "fXB35Kbi1JT2fsZvaLHJDGCfsmjrG8YuPS"},
  {"fXA351DqotLf3aJrFT6upKcivZbcU3mbzb", "fXA352Cjwrr3WCWoFMsSs2hgP1Ryhzgvdp", "fXA3538b2WnSNehxS6KYD5qULDbKb9WZPG", "fXA3543hF7591ZPkBMsnnABBGVfiUazKEV", "fXA355skXKxrxkkaTzpHedm8bfn9XSNrF8", "fXA3563zcwoEDPc7Tz7Ky8KmJD2Aq1kXnN", "fXA357aTXxidBrnGGDXcMaN9uzEZQgTGrN", "fXA358WXFX3J9dr1b9XTQpUJ8c4uQVFQTo", "fXA359c2Mf3BGiVz6Euqca4o4Y9VowKzRp", "fXA35AHi91qnvPt1n3HBaNyMJLqUCBkJE5", "fXA35B1zNArpsbt588qCvJh6oXgUC223C5", "fXA35CNe3rgjTw4jAXGxxZg45qodF5nnpA", "fXA35DuMsNW1EK6d5YJDCCCZDajeCFak6f", "fXA35ExvNLd1USnMp4soNjttU2LnuE8Sq2", "fXA35FH3rn3MMrmojfM8hs5U7kK4DAMWPR", "fXA35G9kfQecRxKPamFMbczVMDrZhyJUHE", "fXA35Hi1JBDyz1qmJf97vJrA5yAtZr7AZq", "fXA35J7nQzS3sQGMGPi2tVRUmo49gEmDmg", "fXA35KroDBcR69KxGJ41PbHkx2DWyjdGQK"}},
 {{"fXB3611XnV2h31uLNHyBf3NTTFE8Z73AoL", "fXB362B4P7ec6JmuBYTgkHuVfT3Zwh1EQA", "fXB363KuYpnxEbUmwD7BJs2xEFQFeW97fG", "fXB364zGGeD6Qnu7u2ARDNqZnsTZAVn8n4", "fXB365Q4aaVAieNxdwXnRe5hTmf199BjxD", "fXB366D49XY5jSJvbDSi6URxPJqQ1KteNf", "fXB3676TnGqS4gPqvBU7viHjzRUwrYB9K5", "fXB3684XEMkUuxAbCsjJqE8HHj6Z32sYLB", "fXB369KXuKmWzaDRjyviVRtmJipPgfzHxM", "fXB36AAhBoJuhaDP66zEGNmPhXU5b828kb", "fXB36BsnQ1Cd4VUbiHoS5GXiYGFkkBBHNe", "fXB36C3NdoUPhqbWuUX9EYRYbJPiqg6ThE", "fXB36DZWdJwN2erD7xunyaLoDPkpkuh3Qp", "fXB36EZwUCcUBGtrE65FUCSh8JQ7orrBJs", "fXB36FG7VxcdwtrVbvSYqKpUXSqnkvhxde", "fXB36G6YQhrFQHLUuWiEZqHVh8afZtnqgK", "fXB36H5mQmCJSD2KPzTuUXQWG1UVZDTgpv", "fXB36JkuMaXvbtuciq7VxA35oFc4EFZfMM", "fXB36K4EJFRGN4EtoM23QeMNfmyw964Pdj"},
  {"fXA361akv4qWCwyFBqtwXAQCCRxGqZxiTw", "fXA362FXzjxnitX4MRPk81sbCTmddgK46z", "fXA363K1nSZqSvFFaaiSbN3RJaEaRdMBiD", "fXA364CFnQdAnMesayoWK84sroKGndUpjC", "fXA365BALYrVtPA7mDy8HqWXprV3Q3TAuZ", "fXA366cP2PNsRYEYHvea5URZZg4Wz14Xqf", "fXA367caaVT4yWznfnNjfGmnXbRKEHrG43", "fXA368Dp9NE9aqPaawW79WvUvo4e426E7o", "fXA369jWwRZ5CA3Psp4Kzz1CqHj1cs1Enq", "fXA36A3LssUvCMNVVafRvvS4NDq7fa2goT", "fXA36BDETFKPZ77nWrTyoiNBrJS2sgmmJ7", "fXA36CwkgcBSikFyUuYr2tSV1ELJSuT6AT", "fXA36DqhhhSPSm1yRHDXoAkvdTnUf5Gqzs", "fXA36EHyZ33xNgLvWxc3CdXF12sWZ3eY61", "fXA36FUuehxoE5NAXhr4hdcF88buyWXJCb", "fXA36Gjh22xegoa3XyqdEafTaUbyDSycFn", "fXA36HLmvaqC7Rcgh6PpBRhYkM3QCaXxQk", "fXA36JuChX4eGwwWb8bceA5f9uB1DqtCp2", "fXA36Kgqf2udVmpptgDaCVfeW5DhKSL3Yc"}},
 {{"fXB3713fgZoe5B7gj7AQ813k9QKfjXq1ME", "fXB372z1uKLExgH4rx4PKddFS19ZXqo7qc", "fXB3735LDmBjqVQverejmKxSRoFnAZ3Qfq", "fXB374RLoXvbavsNKdpzRSHn7dZ1mpubbw", "fXB375NEboFHXZHLbNfHhxuMYDBTo3QYzF", "fXB376iFV5VLffCdxruafDNX8KLFoUgmKJ", "fXB377KsXG82rWtPAtdu1HguFzYkvYcRLU", "fXB378ERYPQQxrZEzL7mYMSBYeKwcpxyZj", "fXB379RRRNhZhUyqFRdCxumfrqzmq9duef", "fXB37AAENUGMpZzFXXtCzJ3Td3dnYkRAXK", "fXB37BT45iEZ3mN3zjG4EupT275UbSAqbD", "fXB37CE4Cjcc1riE9hogEC1YnsL58EoGRR", "fXB37D7BACmiLDxtvj8EqxGuHoxnPZ69Xv", "fXB37EN3JWB1Ns15QBVeEZrcur4C6yZccZ", "fXB37F2UvfDWgFf3X4ns2yaQXwzSHYNFPo", "fXB37GaVrMrwKS5zCfeKQjXp4PvXKiG478", "fXB37HqSKbzkVkWV6EMyNgFFuYnEo3ZWuu", "fXB37JhRiio9L6uYQQF4Pe43W1nMyZxMxp", "fXB37KPsHcYJnya8owpb9bZC2phvM47sRt"},
  {"fXA371VAjnPLaBFuvs24zbfLQyNpVrnafT", "fXA3728WV1x2JeLTmuTbDo9mCdafwva9hg", "fXA373aBWqRa1cFMUHs6sJp7GEvXwrHGPT", "fXA374coa8z3NvDZx27qRMGVeFRbxKWGoN", "fXA375riGZgy8AuAYoaC2hrW7Db6Yy4jee", "fXA376ssJcQtuMLEbtVZ1zPkTww71QPAEw", "fXA377vzmJ7RTjWUP9f1qNFvcBErUJUj6A", "fXA378fV81Uncd6Kc1v6TcMcWrHTDVpDnF", "fXA379mWfggkL7SVhEo19aTujQz9tdyGdU", "fXA37AiB5zTV7HmAeTj7P7nnRB3r7ye3HH", "fXA37B42xRPTbgz8nGi75ZfJWNpLAJEX2t", "fXA37CQCRmdhG7gMYE66sQ6uQzNo4zAkMF", "fXA37DorRDugoF3pS7r8qC2G5ktfUGjnKm", "fXA37EkLKXEQ18oKmxaLCSvwTXddtjBVNB", "fXA37FwRfzHFba1dptPDa6qD6qsf49z2RY", "fXA37GjQC9MEriYZQ3TFtB6QgSPFwSppkZ", "fXA37HzMxPs3ZmyccDyX84qDigA5sTuZuC", "fXA37JxtwFiZ1EnQBkBgdRamYKifEpkoqi", "fXA37KZdmyBNSbKXS71Ky4NySWBbxvt7Sa"}}},
// even trading round, pair 4
{{{"fXB411nM3LFSvErLMiN4WRUcstUpfTZiGt", "fXB4121jwfwkimFtfZgWU5pvB8bA1Awd8Z", "fXB413tx4Mrps5zFLjqoZStVKtMUVVrvCz", "fXB414JAcWeNAHxX1TjNDJ3w847vMQ9giH", "fXB415aoRuftdws27wTTVNFsvCGKazYbaL", "fXB416XqoNgvGS3eJhLbLiDenVM2Utz2Ve", "fXB417mYVZTy4HkTtVxcwgF23c5dksNNjK", "fXB418Cjyi4SF9J6QsnpCskUvscPsPEy18", "fXB419GX1HbbMnbmYGpqetfWskMMvqCh22", "fXB41AuWakGg3brSRKwESejdc8hhhzkM3r", "fXB41By7FY4oonSHN6cqsGfLwYmjTgZzsf", "fXB41CfwvrKYCUpupMbsLG1tCuTY2FgpEh", "fXB41DkWNeaA4HrM1N4QNPjQ3Q6hnk3xYY", "fXB41EVuPeXGvegQ7GHLYpX1Hw37SfTCMR", "fXB41F5euu7vtwrYHP29Fq7cJYPxKvCTNC", "fXB41GmkHkXq9nw95Z4pQL9NTFTBjATTrT", "fXB41HNvLz3axAyQBHMFxLRbRomPfnf6W4", "fXB41JnSJFDBKj1WngKPR3cagvMz4oRWWB", "fXB41Ku83owQW6o4ijNcmhTkhb4bjsDC1o"},
  {"fXA411YGZev41mv8vN5Wwnyux8BSwGXAeT", "fXA412tKMLYQFH3hqoS2wyCwVsuibJkPi8", "fXA413tD3USmQiwGjFmMTWoXNCWEfMNoWE", "fXA4144u6H4sRyoSaqYHdTv4DxKx8Ffan3", "fXA415jnfKsvB5RM1raYYeTawDhfg1vVFY", "fXA4163Xoi8wvGb74MWp5ytwDc1KNHKnCA", "fXA4172rA8mJcA2cq52L5LaZ3eCzqYFn1G", "fXA418Z1yDonWeGdE62g6F1UskUNzPY8Ns", "fXA419tdiYezkhuRpB6b2qHoKXZ16WaXc4", "fXA41ArLAf6FMXvGKbAF2RAK3yg3fssy3o", "fXA41BWXDnJPGBKT3vwnZkkpgYK3mwf3Z3", "fXA41C3PaEKDak8gV9U3QtDjFetGayFcS8", "fXA41DqGFaDwA634raYRz3apiQ3zhDWnhN", "fXA41E9VSEtvn7M3ZXKQwHHeAdiRmp4zn5", "fXA41FxrtnerigHyR8jMQs9gDRxg4vCjC5", "fXA41GCBbAbhfPjS6P73YnCrJAJJwK7kPs", "fXA41HAZusN9N4Hdwia8LBWsd6pz8VBj2H", "fXA41JAd3L81P7p4SYUbDcpigkEVQzPZ2D", "fXA41K7LdFriFQH39GhBQ6nQRENXQ3vUBR"}},
 {{"fXB4219BfMfRo3WVnxpuZz6hhiB4ztGrBk", "fXB422JD6LSVCozP1UkvQfCPbKffs89xRX", "fXB423JMstQT2fsvxVWpEserXp9ZdsHjbf", "fXB424RGctFjT5AGrj7YmVLECz6ABHNy47", "fXB425RrXxKzxfgKAyGhzkTQVZ3acndoby", "fXB426NAJh7NCD4s9uGs4p1x1KXtPUtZfT", "fXB427XzjLT73NKZvNaVSJaQRSZdfnH9f1", "fXB428HR7iYbzYGDLjPxCXM8kday9ZVJQs", "fXB429w1TnMotrDhksHhZyq7LBYiVAMYuJ", "fXB42ATCWA42yWjWQP3sLeJJP7WbT5cdiN", "fXB42B1b7DYuJQSWdu46sqcdJkXvk3SSj5", "fXB42CCHPE6uXMFhiWrNsV6LkWSe2JUNFz", "fXB42DDvMcDb2Gtk1ZvJsE9MCCR8ZsWkbJ", "fXB42EbWoRPgWXUxCQ4TVZnxrt3UcZfPqM", "fXB42FaoUGQyszohPGTVh8u7BNURimaBzt", "fXB42Gbfsm4fpoRCH1hYQ364kqhzESTs8F", "fXB42H1k61WpRxxuuaZY8zMyZPcA5P2V2V", "fXB42JrVzhirudrYDbtykFeJ8yJdKtVgHG", "fXB42KtF7ARwjeqjwX92vL9A6d8XVfna42"},
  {"fXA421TEnzukoPyxTAm4DApT4GYVk7tUkp", "fXA422t3pS7Ct9Pv6aT8RvXrYBGr1Lqyjd", "fXA423qP8Gbi3Z59E759dVnqkUoarS7KZC", "fXA424VsmHxtUSpNpk2JT6czZRhvxhNXAi", "fXA4251xfo7s6TxT8ysdkXNz1SKoUaagtP", "fXA426b2tCBmFv4uh2trLNkrtzGWetwAer", "fXA427uKHjvKyRdKjffJt29G7vqhrNPQMW", "fXA4286dbmaRTLvw61nq9ESB5HKtJLkJgm", "fXA4295hoR3oBkiufrQM26nEhJEn46UBob", "fXA42A6yspE2FGiuPtu8Lv2NMC9hyQvwSi", "fXA42B4PXvDQYYgScgySbZxDUSE7sAwTdY", "fXA42CyDZXn62iCZL37j7ivH56Ypqnumny", "fXA42DhW7tJJtBPZn5be4ggEfJbcY7rkxg", "fXA42EkUcHQ1gqqH6aNeVRtSykybRDTmNC", "fXA42FyazmscHWav2AsdyCZNugfvR3jcyq", "fXA42G961MB41tfm3uFvQdppoQmb3qHe8X", "fXA42HtBQeRX6TUia4Lo6fS8k7kcVFNNeB", "fXA42Jj6WdDRDpCSMXxRUTWUArus8UKN2C", "fXA42KeCeGad2t3P3FwWEyhB4cVqDbBXfU"}},
 {{"fXB4313yEQckLeET4NU4gv5ceVHpHu8V6s", "fXB432k9Fs2eqGGkhjaW6AGgu2MNtY6TsK", "fXB4335jYuNP6fTcP6RqAxLxKog3VWK1Zv", "fXB434j2Af84oqyoLUyS3vqvRNv363RkRf", "fXB435YTBB2XfT5jHeYpGPBciTPkdyfq3y", "fXB436enr5cHEL6e3kyZ3BH1Ju9yw7DWGX", "fXB437nhJSHS7efjiqR7RDE3gEnZH1SxCi", "fXB438E5SJk6na8xJZLWNZ8CFVi4J9Nkni", "fXB439jEsVv9p5SPtwXMDNAKKjD8WxEMeh", "fXB43AWBq5GaDMxSm7JcuVMNpL52h8fCCv", "fXB43BprNPyNjwSuP8B3JHj1fGcbVFVqta", "fXB43CqWyBim4tioZf7hbtSfYjJKVBAkCy", "fXB43DNmVbg65D7eBTXXNFfkhsozvB1FSo", "fXB43ExrDQY4Ryjb1dRi2ZidEaaA9NuBth", "fXB43FnYaBhzh4sHKMW37rgqNKqMSjQCo5", "fXB43G8NEdyMpBZQXfQVrwkhbiHzq5254M", "fXB43HMVFeMbQjf3vUZzrA718vDDtdvuUr", "fXB43JWr9LWZ1fwqJELZFm1fUvdaYm4kYx", "fXB43KcXMCgKqQY88pV7ppGMAPEqRb5XMr"},
  {"fXA431WrZhagW7btQLGRrRPm27ToZaHhHG", "fXA4323PSBj1rxd8Ne8bE4F1CW9PVYCx3Y", "fXA433Bw7x2ZfxhjjLMpXKuXFtGq3teCJ2", "fXA434HvDYGT3poNoUzjcPbDNRxEtQLYuh", "fXA435PPK2gPywSoakD9VSwWADHmgVvqKX", "fXA436XNGY4tnYtDUndywv9ZAhsYA5w9wJ", "fXA437nfuidRRZRcbyXQQzwUH9LmBcEe56", "fXA438WykpRidHTSrht1fLtkdugr6TRWnH", "fXA439ZB1mXUPyMfyweQxZzpnPWu2CYmS4", "fXA43AY3JvZQzeY9uMDuX8zNi3AmuNsQ4T", "fXA43B3AtmKW3TEukWcedFguA9RXRtdR67", "fXA43CHHnDScoTirEbgd4bKiZggzxBjLwt", "fXA43Dao4gjgbdB7S4tW29VRNWyZB85htM", "fXA43ErtG4nNFWUobUw7e3P2VVpmKeNLNd", "fXA43FTTLD7HUFPpfvXZJarRYoFVjY3VCf", "fXA43GLPbAGsWb5UD4stm6e8rJjRkcfzmd", "fXA43HzeCr5749tXkDHeEotWp1bTQ51MBP", "fXA43JcjMdcuTWYopdR2UEJ4SvD375GNjX", "fXA43K2nJPUyvy3yMG2rVmys792MeMuzea"}},
 {{"fXB441mDd41pB9Ww57T8m2hfgjQu6KZUFN", "fXB442JhuYtGPWrAyXJEde2NhuBLnnG1Kf", "fXB443MEqDAr1JLYfWhmGH8TSyVtBCzmUZ", "fXB444oMuXuZFR49M9LNEPdMPvCYdza7eK", "fXB445FpKUHLWiGh5tQLw5FFFkhAkJxEHs", "fXB446GUEELkxNS6tpQFEzHgwg2BmmMnHe", "fXB447HmCwEJTk4LBcMcL3Je3w1V2LzyRv", "fXB448h9p5v1HxK6VueUCfW25ynysg7vpe", "fXB449yZU8oeXjDNbmHZvaVFsg9oq8caGq", "fXB44A6yy1UGLSHprGa9KQvtjoBEJ3wmvM", "fXB44BiZESWNqxM5swfCsBDywrk7Lh1QNe", "fXB44CDTNCydH1tuGhvLvpSf9tB7ZvJk3v", "fXB44DYYqon2uFSiRmVZZ22HTADD2X5TyV", "fXB44EDWBSQWoKrGE8cHEXRDWibwnomtJW", "fXB44FWJKaEDZ6MBm28ghqCsjcZKUnSkPm", "fXB44G1LbMyb5Hxetoi4Nmsu1ZnL6LSYBY", "fXB44HnHuQP5q3hajNAoaegua7jEgcgZpb", "fXB44JwS72mYScj1PJFJzM2ECqYiWPH6xn", "fXB44Kt64ASPYrorzSzGTo91kVwhYakWH6"},
  {"fXA441YPSxsG8WejkQ97ssNoz2sWD39chv", "fXA442n3hjH1mWA9ymyiNV6izgxdcKdNg8", "fXA443SoLEDoJhGNv9FDQQhU8U5WEUtEbt", "fXA444C9134fLDxk3xi1wFEpARvPm42b8G", "fXA44542pdFvz9ojrek6wWf5anoNhRGn5K", "fXA446JAS1J9xs4ZZf5b71HheVPHS8vqA8", "fXA447J66ZGyWr3zX4XfNL37TbDVdKJYWo", "fXA4485tKi6Upn6dJ18Ag1Snvmrfotv715", "fXA4497ufXBMrMMKdTPhkYNFysUwvKJeJZ", "fXA44AjV8xJD3H2bfDxNxPaAEi2yxs69RN", "fXA44BAVA8jFyHX5UFcMbeLEqeBJjAYyhc", "fXA44C1iKfoKiqugnfe4wxbgP1TVNPiwR3", "fXA44DBn5bczK6tqiiJCYgvwaNkUTgL2oT", "fXA44EFFZUSzpzWzB2r2Pnj2pKuVeSLCfX", "fXA44F9sYTEsV79u3hDaoyX9KKb7YDR48p", "fXA44Gdte8bvzjKjsRzrjTrX51sLp3AT2p", "fXA44HxRHcUwaQPmtzkNmkqqLQnyggxdg5", "fXA44JtJDKJYGTHfwqzqmVwP52KVGwmDge", "fXA44KsXEu1RJUXx5Qg9htjDUvayfhc4F3"}},
 {{"fXB451ZzckZ2sk2FVsKZQ4ZMX7wcMJDRKZ", "fXB452rmEXe8wn6ATzt6q9qzE2g461vb9Y", "fXB453AwRnUWFmBXyhgQaHAZRLZcPjowZS", "fXB454i7g2bKWwmGJERrRBWBbgDBgo8zNx", "fXB455kmX9d7gimtU44KVywa91i7UMqrZN", "fXB4566NSKfKcCet5CEVsK9TwRqAgqE3ED", "fXB457gyptWJFJoLEeGvqrEkNbMpue2owf", "fXB458a4GipjoUEdw88xDLUcoWCie8frSA", "fXB459YYPiCowFa41pRjc6ZBNuo5vypgxV", "fXB45ATuACzsVUatbXZEyDd7FeTwAwM5ft", "fXB45BDdCjYAK1FGcmpVJqW7M4uR2wrnGv", "fXB45CTq9KUmzS2tRZbTND3KBRHJ6gpQm1", "fXB45DJY763sAGhSitZNLoKNybUyn7uSR8", "fXB45EUnPhNiasr3M9CMvu8NSyKKqB4JUi", "fXB45FGzQosPGYzfrkmBgH8qQZ7m5EEoL2", "fXB45GQpUCaqUakXVt7PJPPFgUrQDFZYH2", "fXB45H45cVhrVeTfHwvAStroQscNZYL5nP", "fXB45JtoS5GQhu6PiGtc1FnHpiCSq2VAtS", "fXB45Krn9eaqL3hPgrGwJ4ZaS4dNUa2kHG"},
  {"fXA451st4rsv6mhe4nK7mvWYtahdSdvVNz", "fXA452FUbMH9okmwgtKSUk1tbfwnzgFYxv", "fXA453Bcu2d372FwWK5Vc2MMKrRzgUqu2o", "fXA4547D5YUmztS9M47tRSroiXn1CnMwhM", "fXA455qNoB9MwojHwf7mSYxkzopwceftnZ", "fXA456Zq7JWFRjQfu64UGTdRbGVatyR85M", "fXA457uB5uoX5xXZcUoWzu17ehLkS8NWyY", "fXA45886ckEL313j4Yz4G3XHEhig3ETVy5", "fXA459aEer6hwNqqbyX1tWLmbBVDisVYmG", "fXA45ANu33u3zcmLDJTRjegG1VYRvDwUPo", "fXA45BrQ78bCsjpPCwJFQyvm2qNLN963xu", "fXA45CGgeRzQ4opZA6HjQcdWBCXd2duAPx", "fXA45DAx1VQnhw7cJ9TW8Bqf9bWhUdgpB8", "fXA45EhKTMozq5DSuBoEArCWDF4b8MApXe", "fXA45F2aExGhHhxYvwoy5QXBymudFjfiwd", "fXA45GA9h6JGpJQBoCopV3SwZigJVcDLZW", "fXA45HFVVcZQRZ49k13LUkRJbNMczfQ9sV", "fXA45Jr3Ta2UiptULPyKQMFpJkd26gDCYG", "fXA45KEAkbJWMcUBL9HRa7N8DQ2u3hNUum"}},
 {{"fXB461zLMQXNg5ufvLxbvA239uHtb3hBu9", "fXB462rqGf7jFT3ssPV1KzN5ScJxk8TciF", "fXB4633zZAJ7fvo4tufwLZTQQQB6MyMLDx", "fXB4648sHtdNfL7Eu6TNSYDJ3dYS6UezkZ", "fXB465XTzF1Wx6JnwcxmNHu3a6Yhvsa1Rg", "fXB46614AH7NViXKgdWUkQLMLGH6NnV6Co", "fXB467Ai93oAqvt5MZ7AfBDGaY3Nj2Cp7V", "fXB468z35k3dbpTeaeWp1bgkGTsTEKj9Jx", "fXB469k26imJw4RsFG3MzXurV7NWh19VtE", "fXB46AbKi4RkYipvDgjXfBZMG5YcEYsGoq", "fXB46BzsDY2wqGZuCFyX5Gwun8ZZYGMfZy", "fXB46C44EYrpNhv7ByyAzbpMZt2YKwtW8U", "fXB46DzGH8hudkcwDWXrzRepXBBa7j5Z57", "fXB46EuqPnW8nfrNYAnAK9mNL4Y6M25zGM", "fXB46FqrZHjyDeBPuhrC1ewtvKHEEgD3te", "fXB46GTr2rv2vPWRSRszhRu1UdtbwnD2aR", "fXB46HnvDQY6YuZJPGryg4odtdshE92Afj", "fXB46JmzLw43BM9A5iBYSedC1vTKSpaAwt", "fXB46KiSdzDezin5Hu421FYRUmjQKmyXC6"},
  {"fXA461jGMbP66SLiAjuXEcUzZKXvSNVjWr", "fXA4622mxtZsFfGzQm1iicHMZksbK1NUoX", "fXA463gSbJDupHpm72ndb8PLgKWY5PjsAR", "fXA4648W7V9EegBkjbRmwHPF6UW6ZdHAT5", "fXA465SRjcf6HYzQ8MgbnqYPdZ7FouZcjY", "fXA466uTFotpp4EEu1QsfKus2GaBwpGZ92", "fXA467gM9RW4y3cMg3Uu7sEoyynpcDcnss", "fXA468ixrJYpxvfUUMKgF3guJPnpvVJv7k", "fXA469MTRdRFDyiJ2TrE7hCUEAUA9Tnqat", "fXA46ABLCodxQjfiV5Gg4EDFWaAWvRHxwP", "fXA46B3CGMd29BgQPNpFLYqmKjJfxM3Gdt", "fXA46CLdnaoXgUk3XVjEfBo2T487UyDAxz", "fXA46DgjgrZ78S7QUAMpr799egU8drqRDE", "fXA46EqjoZEAUzy5jbf1ThPRQYTVS6hoWf", "fXA46FzVRCckTWhKp1yvyYoPVuJukN2i4h", "fXA46GfTn8GB9FmvayhUEeeXXR9XBizqRG", "fXA46HfX3uqLfAcrRbqtxoeUHJPA6y453X", "fXA46Ja2J1BLfUKHEX7oQWxN8TpVmS3MNY", "fXA46KFmwBR3W58QazWtFmcQug1LDbTQs3"}},
 {{"fXB471ZsAMiHKJ3VM2zCBeDR4SZAz2UaV7", "fXB472bdZXB4psChQnGWgU1ic2Cij87HUR", "fXB473C2ynHTpogPnea4mNEvyBw7nVy6He", "fXB474A933fz4a5bm2pyXwQC4Den4QBovH", "fXB475DVEyLyXQNv29GYhtWYBf8awHKvmB", "fXB476eC13ixBB4cLE5tiQvMeYZ143tV2M", "fXB477TEWxu9nZ61mbwposn6FgJhBxi2Fi", "fXB478ZehFmeD8you4vizRgqFWMHXcwgr6", "fXB479cHiq9ZaXuEiaDkdkByRSbxPAVnRv", "fXB47AG37an1PU4Ln5nM2j919eie79NUrA", "fXB47Bm9X8FK6323ReopLpY69JCNPj1zAr", "fXB47C6hf9BemMMtZdDxH66GJxWY1pCcGY", "fXB47DHsbWLfHRxJMd1WnnqZxaHGUpdXCw", "fXB47EwhqXhgPWQQ8DhSm1mgFLWhvQjgfv", "fXB47Ffj45N6Nc4YpVc2Zej6ByQUeqYyDf", "fXB47Gqz2iNXVNGUQHUAR3PptW3ZnwgAt4", "fXB47HouHG5ZQLkg37kjSqCgEiUcGpH1b2", "fXB47Jk4p9JRqZxkSG9XZYiqL8Mew9zX8G", "fXB47KKzXTQWoXf5kWAmM2CwieYZRMxNoo"},
  {"fXA471GK6D1G1QksUsZSU6hf1qdLuCF595", "fXA472rXuBYtJh71d2z8xTNE8q9VSsK4MA", "fXA473BsNvu8vqntNbHE2dgEfvwG6Z8LoM", "fXA474Z1kCkZNcCVJdkaAbHMFdZyyEgk1k", "fXA475tbkEcNhaFBQ4kL7gz7iuLqCBudSP", "fXA476yWyhFLusJZhhgvpQ5MAsP4S2BUm5", "fXA477LNz1RoWkBq41WPdkyD895MqDiX6o", "fXA478h8pMWDXbxH6cLQCzd1TKJJSjC7rd", "fXA479rr5m2yRatSR69MPCYREAGXExbaow", "fXA47AcSMxLjqs3SEN6nt3ftHbJeTTUZ8H", "fXA47BZPw7NrFRTSAfg7kEyuNndDyxHzkS", "fXA47CpyCDiCZ6yJYsUoP4n2yagTva6PHL", "fXA47DvvfNwZYbQZDesMRmx53goANAMtpw", "fXA47ETGqJjJmGpC5szCUQvRD9Z2X1tcsW", "fXA47FTA3W2rgJzAzyKefqzzBwmUWD18GM", "fXA47GeipixNK226R5HimgVGHki4tQae1o", "fXA47Hu4EwmzQ79pGa2kQs6Nja7kdYyrmn", "fXA47JdNfptQj2YSesyVJVB5U8jPcN1qTJ", "fXA47KpsW59qfD7zezNhH9xiJRcsGQiwR4"}}}
    },
    {
// odd trading round, pair 1
{{{"fXb111WjeHjLwiiFWC1VPJr2fYJTU9aHDm", "fXb112cZNh38R7JsoTwuTGLc4cLE2pZX39", "fXb113gaiw3D8DBR2k9V7zVawTFzcwPHaX", "fXb114R8VRj1ugKzF9qxcfcEnZ65kHF2GJ", "fXb115vamRGpttGaLho3Hup2rmSfzRKKAy", "fXb116wtPJhpWrjGoMtRuTtUZBMvQzzpct", "fXb117kZ6dCanJ62Ep9g8JNCkkw8bSAzhe", "fXb118hsUWnoYcGzXVJBPXyksgSXZ3WyXW", "fXb119x78zCAfa9GXxQD1to9kYeV21Zc8n", "fXb11A42NGA1yJ7sT2ak8bK6Q2qv2z81Hb", "fXb11BwvgupUQxMENP5tkCgmyRJov4hMSC", "fXb11CNS97kTPdnRSyTusBiuxt2eCWpNra", "fXb11DiXQH7DVLjiwrWxCueMGuawQUTXn5", "fXb11EyaDrTxJm3viGHr1kpXFQGrFabp6X", "fXb11FK3hy7MAMZU6UEqs6GK5cM8Dittqi", "fXb11GZZh9pMnamgrtcSpQr2wwjE4FCrUb", "fXb11Hm7gh9rFUtVbWnk72NfoGg38zN34e", "fXb11Ja4YkjemSMceVzAQdFJ1JKsSreJNL", "fXb11KeCtV5HxA5VpxoDcgxuspymngS1Ye"},
  {"fXa111PuZiLRMHEeyAAHUD9DeT8Lfrei6Q", "fXa1125HnPCrVAGPVDys1SWgpXq473CPVP", "fXa113qhEtQ8nDCHqEYFddHmBSkPDRxrxm", "fXa1141ZF8ve3JtWuVhPFK7ceX8kvENpTA", "fXa115qo656JFjWdrsV2RDQ1sq8rkpZwDz", "fXa1163mP7nF92VUhinvSeheqmswECvMaS", "fXa117Wo63g2f4oB1mPQ2718cqRcNkfd3W", "fXa118L5qM7jvDXcPifZ1GzcPqe4KQmAea", "fXa119MeyDEzFjNtCVaB3rNwuMkgo8vjq1", "fXa11A5pU1KHLQ4aoVptJXeDkC5uMDeWho", "fXa11B83aL2WCSJCvFrqXcgTbq2ZZDFfaw", "fXa11CqYj2YHYTakbexfMCGaXPRYNcpheH", "fXa11DegHbY38N1Bx6ucVUxjUWAgMS2j5j", "fXa11ES4THDFoVJcRyrqcuUTx4figh418H", "fXa11FeMnB2UGoAxS85v2dUijPCRXuMFFh", "fXa11GwyqftSPuMdZLjJU9v6jhZ62Bj4KC", "fXa11HAqwpQ9AdQxB4NdHWLvhTFQg1U6Rd", "fXa11J2rgkbH8V9HEauJaviA1ELDMcpGKU", "fXa11KdNWw3Pw1pTiTXtH6JpQXThEd6Yfz"}},
 {{"fXb121FLgL432K6hUekbr3qVMGEHxuiXYM", "fXb122J1JsAhyhrLVqBWYdfitz1fXvUnZ9", "fXb123zKdXkvvJQXBahZZMg6ugUeNGptu5", "fXb124dwcVoZCNAQVopv73amAmVVY5VDAz", "fXb125DhfqhS6rZYdvLhRCvv2HyctvVPpK", "fXb126HU56PGSqus4E6aVw2H2ouFcb48ez", "fXb127i3psTpPWPAJJef2yJoZYNxVyiDmf", "fXb128R473vnC1mFkLxhSkZYCvtYmLnDoV", "fXb129PspBJFyz3ss3TWhbsShbqd1az71k", "fXb12APyUrJg4XahQQmg2WnLdv4YQDe8EM", "fXb12BMLA42jngJW2JF5uhanFLXXVUzoge", "fXb12CxnrpPbrDiQYLgPZHFSUq66YFpRUg", "fXb12DY5hijBE7PgnceXY5xWuEAsD52nrD", "fXb12EwNw1fdubtZCoDaK8MRBFLywt3bQP", "fXb12FCTUeSKeEmkCt7Q87V7gL9U351tHr", "fXb12GrbHmwECFsxaNY5u1QqCssz1fCddT", "fXb12HR8CNXxvwWEX9C4AW2U5U4NXGxAZp", "fXb12JYCmy5EXMQg4wNdahzYBSexRANdnT", "fXb12KSv1q1EeeTGHcNYQcL9Dc7Bnm5ykH"},
  {"fXa121Zgv1E8Ju8CHnusvZaYpvnLpDjZEJ", "fXa122NzSk18vDyYYJvywyTYAkYDM5FPZL", "fXa123bSxGRaSWqhf4TEDW8Kg6aPmXFK5S", "fXa124fu4WyJ3zN5kYpkqLnb2K6fC658bG", "fXa125GGJ8FjPFrT1pRyGC5JWMk8CVZGnv", "fXa126yxCSffWmgSfRDX9tvix6gTVyX2Eb", "fXa127KJYDRdBzyxA4Q8iWz9S67mWoQUW1", "fXa128YSG336qCGK2JGqrQAPRWQ98sAxkr", "fXa129kG2kZz8Hz7E4MePByZ1Po1ZVXHLQ", "fXa12AeZnPo3nrp6qczCM85paDxWcTTnf2", "fXa12BKjDDaXKRwzpSmvwi8GTCePTqM9Pw", "fXa12Cs1tupt8M1gfWpqi8qqdAQi9V8QWk", "fXa12DyJM5Qf1oPELhgk1WudgbsbEYSCWU", "fXa12EKTxkzYi7enAQcRFWiepzCfGfCj1A", "fXa12FKXhT5aFjUEHUQRVBapCRSLFxVnuo", "fXa12GqxkV5SmFVcYKk59xLva9aaPhdfNi", "fXa12H7peHZqFPad2UxKRk4hcaeYjouWwJ", "fXa12J6vuJ95vDwkqWuZG3EcXcwL7XNNGS", "fXa12Kks9te1J9pLQcKfNbZoFLZojE5Gst"}},
 {{"fXb131x6eExG3WxSK58suJHE7DLMSb3pwW", "fXb132TTXsu2kxEkJZEfFryiXy1EaijYyw", "fXb133a5wDdvZ23ZtZ2caiwDTKw5HqWiME", "fXb134x4FqRmo4Kr3dxN77dMj6VGAbDzf6", "fXb135xXqgiJPVToQTtFG8Wj2rHtKmyyhj", "fXb136VY89c1E5jCwKSS8Hrx7W2eGVF8L1", "fXb137z8G6jXieCTpbpYxFwgy8PWgWYDq7", "fXb1384mMZpDBJ7a9gd3dwp9mitF3SPQeB", "fXb139cc6LJvrom5fXS35nh3UWApomMmdv", "fXb13AMgHcMhiQfg8KT8nRicJPmwRYFK1d", "fXb13BCrC5md91ebLAmhcFWxdMtEaXMTkr", "fXb13C4fBmzdH6AjEbyvDy8jz7Zmi8hk8Q", "fXb13DWwFBc7pJQaMBgtUm9foFa1qhNork", "fXb13EBSrfPRBePkmiKUEM43W3XTqE6rh7", "fXb13FeNSfP17nhQQ3k13am4hqiCTAzFMD", "fXb13GLdDvDdXwhxhR5WRxVjtqbtCy5EMR", "fXb13HDE6JJsRKoZgkdtHU32YVcvDbSasN", "fXb13J1QP9cz5SEsGvNsQWSRy8ERvWXzLd", "fXb13KU5xh47QpJZXHBC6kfumoZyDuPrjj"},
  {"fXa131ZRBHQiZBqizgYNjk5yJuA3cGXE9V", "fXa1321NvsNqU32YT4kof9GP7MB7qbgAdL", "fXa133rG7EsP2SURZbdXhCES6ozUCBAXfk", "fXa134q7Q8tjAyF3Tr4yqPknS5NKR5jDKG", "fXa13515MsJpXT1QQBm13ohpeZdWPZEf6n", "fXa136rMdA1SwhmQfPuy1THZUYNm4TS8yG", "fXa137DANqhSqUNaXxkqxEYfs5CgUvTY8b", "fXa13813gGtZYzXXnugvypnX2ZCNEFZcKV", "fXa13928zx24ke5YmgRvtYqFjM6Kiv9xuY", "fXa13AfXSnxeso5wCr3Wi9QKKc2ftz51cT", "fXa13BbfnpS2gL9wkz5dv6rLBfHbSPic42", "fXa13CuzKqUvf7UWgC75YNA3nf2x1FwN4R", "fXa13D1m5dzNUW1fPTzN5gQjmYVgz2YTKg", "fXa13EQHdtJtEqz2Pb75pTNpBDwXKyDENz", "fXa13FZapmxT8EV3thQJsKihVzbSsBFQos", "fXa13Gkcr63mqRuPpRMtEWy2wBk8bUBttL", "fXa13H2NsGWZMbWDCESYXEk9hs6gnfnqok", "fXa13JSYqtpub88VzPPFE6UPWv2pk4B7TJ", "fXa13KU1ch4ycytLqGLPSNgJGLTiKCrNPs"}},
 {{"fXb141tbJKCDZLMExVLgBJRTWpRUNfxdWA", "fXb142Pm5Tvm971GuSk4SrBUBkqWVS9TCv", "fXb143CJbMjGLWCugjMFVAAdNUNGUDkT3d", "fXb144XNcTP1UNjRbq5WBbT7dBdCS6cbjE", "fXb14541UAUvxHYzJEZzu7D4LzEi9pNtBm", "fXb146XxsaUAZFZyk3RjUxbRaDb7o3ek5L", "fXb147HMBa6DHwrEWdjCKF67TkWgZDkUAH", "fXb148mfrS9bdzGRzJoE5kNpqbA6RQb9wP", "fXb149jaWTcvpWwpsNNowPUuVKqg1VVnJi", "fXb14AD7MGCNR6c1xhetzn1T93yv5c5YM7", "fXb14B1WVc6fNpHK8hgG2BbiiEz3AusQWg", "fXb14CwqVsVCSAVaysxVDigZaCBF8jdDUG", "fXb14DBSdCnmvM3xF6TR9z6LgeGUa1fTup", "fXb14EQ9YtzXmiEdZMjBK6vYr1wYGDf5Ku", "fXb14Fgkd9bAoZsFBauWWXJohJGky2J9yV", "fXb14GAwPJVmPXCGo1eC7bXdMvCxejsyFa", "fXb14HDu3pyKRSGGuk5UnS74pG6k5eoVpt", "fXb14JfbHMb7YJTLJz7dV7eVfvg9PHxYYL", "fXb14KegN9QnheGq3uQ9xbmLNCFMycus3r"},
  {"fXa141SeArvhP4FEFFit2i4XxdXGVZsQRi", "fXa142QReah62Hkt1fYQqvtXLNrpU2XSNu", "fXa143npHgE24A2PutDrdFh2ji4difNmzm", "fXa1449wdAM4JvYNFWxwbnHxjWF8Q3HJa4", "fXa145HaiWYgieowVLtrTn2gSkMdvk5hkg", "fXa146EJyVjC1o1cPScy99w64Z6aCdWquu", "fXa147WbJknhhSmB3bfGfAhN15HMfxJcZb", "fXa148iFEvANiM2kyPDreQ9FUSdHoo9PTq", "fXa1496yfW8ftmhYu9fxPTxGkyx4Me1oQ1", "fXa14APiui6kJeDdEushH3TaHFoRhUyVoP", "fXa14BM1RmQBgwVDCJ3LWjny5jfQzXLqcU", "fXa14CdL2FrktPP7XosKqJ7pJChu3TdbYP", "fXa14DQQpkyuPoGbp9QcD2oyXAgdEUoZeQ", "fXa14EG7NtSdKn55xAC4fNbY7zYadwVTGZ", "fXa14FzTas5faSL72xFS3HC79xtzYdiLta", "fXa14Gi8dETAcQo9Hu1bR2DjNDkR512APa", "fXa14HYozMchAEtC9cw4dBoutjXjbZ825f", "fXa14JzL6Lmfq4H6H3UoPqB3LSTDR2sTRe", "fXa14KGuDJ9WYoVxF4NUsexUQwQfEGh8aM"}},
 {{"fXb151DaWs1Gbb5P5R1Q6akC3Ac389kahu", "fXb152pujCyydEHVehpiuvZuA4SVtmLrEc", "fXb153BuJE5okFga3pRCp74CqZ3gbFMfc7", "fXb1546QjWHaiXSuK4S4UBtq7jNYJymYCu", "fXb155ffq8UvMqQV898ebXbHxb7PYvzjCd", "fXb156wLYZgAwGMMbqdcAtf1ZyHekw8Wh4", "fXb157xhjB7WDLDPToJcpUBHgBGgay8ZMB", "fXb158CMfXB8ASpRpsn96YaRpi18CqJjAp", "fXb159uWNwXgfq5G84kn9NNkdwAELYZg8r", "fXb15AHbHf833geJFqnFv511PfVXjiyuoi", "fXb15BQJQVx64fRpzS8goHPKUt57k9Knmq", "fXb15Cn18oo4Lv8HHNgW4rjo5if37Qu5A3", "fXb15Dfy7HXCRKE7ecqGBiWoLHzGNMEQuP", "fXb15E8p3sbTHUR8o7RcCLZD86qqg3NM8h", "fXb15FcP9NGzWYt6wzybzjhsp1jkTMhHLq", "fXb15G3xG1EEQZLGFUfEYo1VcKRSeHLcoF", "fXb15HxYpiqVvEYDYQ1KCXk6fYLFQTr8Ws", "fXb15JaShG5mCTFKHREHYNx3QHY7rjwyeQ", "fXb15Km5KgzoCAUGLLkQAKukgX5YvXNLn5"},
  {"fXa1514U65kSfkXD5KSEFrds8KykMpoS6D", "fXa152qzvfhfaH1fxEbxmmcqRabGy5uTgQ", "fXa153kbfdvuaMy28cM58yJ69XQWAQjYnm", "fXa1546cBbxAqBRcVJqDzDqsTaNwoyG7M2", "fXa155owvgtEzdvzcimWfAnoRzk6mK1uhp", "fXa156Fyz6yDxRfvPYtReNz2CvWKRjMjx6", "fXa157tGBRYCmSyRxNkRLgtbtQjwrCzWye", "fXa158tyNANfJhFTpuR9UX1dctVdVWLATv", "fXa159U3bnHou5eTJ6iNGDC2zL4g6M6WTJ", "fXa15A2LATojpNVKAooZgaa2zjZG2z6GPS", "fXa15BC3pip7aLBQdHvwg29UyjHRbJcazD", "fXa15CfbpcCmKDne4KEKe7BVMQFBaeLQgQ", "fXa15D4TSqSoaRs9AEDSuhZLdR5SY1mT8J", "fXa15Ef2VbhCkvzV3N6iQPdCuVhYRPoBbj", "fXa15FXnMqHqSq3xtjBZziwCoLsFBKwLB4", "fXa15GnS5KptGFRkUTVTDhzPubV51GFTo2", "fXa15HbrJnJPxDn8h3k1GBvgkFQZmHCseK", "fXa15JA9sebXyUWmd7F38FJ37H9hjhiiqn", "fXa15KsSyPAUdpPCRVfKd4AVGfgQNZGvt8"}},
 {{"fXb161CGUeeiiGaieGCttvzdE8aUavQZbG", "fXb162fPzhLr7nYXtnEVwQCGa6MXVPMVuM", "fXb163NdGLYLkfZ3hYCr5kZNavUVBzhB7S", "fXb164Rd2Dv5Zaw2WvdcSUfaFfzDgRFsuR", "fXb1655a1Ydi4b4BbPJtkCVdYZg2rokRxX", "fXb166mqTfqan8pQs6mTMjSTPiiCfUnnPs", "fXb167zCMZeGDVVYBvR3yd2zprJpGudEK4", "fXb16897ncyRpgBGhmmGqZViZEoSHu87NL", "fXb169hF3iqHJvVPNyeP87Wo539KdFDexx", "fXb16ARCEed9vT4S4Fxstt6SgVVzjtGqeM", "fXb16BUZqmK5WH3CFxaUZYVpRVGsVUZD1c", "fXb16CmdUz1uTj8FJk3qT3tuM5AZqKk1Sk", "fXb16DXqMeRFQh7Y3a7CmUEijst4ooPNzi", "fXb16Ed5kYZ1hiixG5bFBohha2KVUmxu1j", "fXb16F8ZRNc9NhmQ64iPkMVWYsfRGq4174", "fXb16GbtCezwH3YycHHzmHhesvwFiS22eZ", "fXb16HHg4yfhbisMxBP9mLydh6qt3MPPnm", "fXb16JoHrN5L7t5zu51MCLGnGs6afMGyiu", "fXb16K2Pmafu5RY54orHdSSuQQx75WZNni"},
  {"fXa161PSxxLQd1wXDY5C2RprqRjAcCKE4F", "fXa162zDZeSdj7awuzpk68Fh6HsTzNFe2t", "fXa163PHuTTw34nZ6mx6USCySSL9kq524A", "fXa164kgPtLGvMkTGrt7UPpqXd3pzB5bx3", "fXa165p5eyJ4WuVvmtdPnG8MuJgzaWfyAh", "fXa166pMdVum79Roy4kBEZd5K917BNXo1f", "fXa167TSQdCFv1XqmfG7PkX4goCoqQ3voZ", "fXa168GjHRRduTcqRD8V84TuTQbdF9vZLZ", "fXa169LDApksBe3LDyFRAENdbnf3XvsvpL", "fXa16Ayn6teZMsku6cEa15keHPzVBJZ5ig", "fXa16BLY4B8ZbLhQ7bq7ye9yNF2tzL2hCL", "fXa16CQb2FXCpBiy1kPMVv95dHnPMQJkxa", "fXa16DG2DVWx7MTukbxdRpdttp9TMNiPaT", "fXa16EXSnSCqJZ7k2tuRrryGYT4ufYL3Y6", "fXa16FZHZwhsSGx5GMmGdqCgNwErNaLKMz", "fXa16GE9gQ8r4sKUuei5UXdpVVgJ9EvfSg", "fXa16HWnQXJQdLkY78bxBpHZ2nUrM2WovC", "fXa16JCSrSmNG2Z9jNoXypzeTskjaMnTgT", "fXa16KHDEYX84ZcHfQr1SNhoZ4tJSsWp1j"}},
 {{"fXb171gKJxawBcoqSyyujbtvZNC6QAKPr5", "fXb172rSFvPZbsHWaoPGTSm33ALcH4PHsA", "fXb173XsUA5DA4UcLHepsdcxYfSWrPzHHq", "fXb174xkybXUmYzud16QK79tXSgdFk581L", "fXb175AbsTUNV2zbMxJokTAxfpXx6svugq", "fXb176fzBhJzYbKzLUSertiX9SyfLwGTti", "fXb177Bm4LA4RuUxRoo74PipGegjcrHX5E", "fXb178WUYgXrt9h9m2wJhDvNqUTiY1eCYA", "fXb179tqXbY6UkS1UGHBVDae2MEyMFAvRG", "fXb17Aak5GZRS9SA7FoBZQwUwEv3se2Wxq", "fXb17BDkGWcZ2xrTnfG8tCouC5wootfwP1", "fXb17CTb1H2upwEhDL56YYDyo1QCqcPabo", "fXb17DR9Kj2DdxDyKoFT2bjeGzcQFQDcb8", "fXb17EBwMCpBKPUwXnBCCCQiJH9PfXwreD", "fXb17FGxkP7jd9jmtaFbUYWXgmy2WFjQTv", "fXb17GCnzP2PSx2MNJfL296a8LtphpssRu", "fXb17HUXFWrkCJs7fsNN92hzLMo3aSMsWL", "fXb17JqURcWKMwD7yt6dF6q9bh3FMgqq1D", "fXb17KdTw2viH7w73X4pffoujF8vW6RBhk"},
  {"fXa1716rJmahbAFD9VaS3QarAk3fBk4tjX", "fXa172D8JwxnUmH2m4jMQxmDtS3pRy1Ktg", "fXa1737UrPU3JEeHbDpZxNbXPCzaHy8ftG", "fXa174UpJZJMbYAa4a8nKaZgkZ5H8bJ7TL", "fXa1753Ay1WTVhWZqvnpjJp8K4LS3QNyAp", "fXa176FLpZSWx199qRwRxzeuXxKDJFTCqA", "fXa177UtKjndDuy5E1AarHeaitHD85RjGC", "fXa178uDwtCmgaPvEdZtgZSqjxEKgaUiBC", "fXa179tZ47iZE6MnbXM9FhBNvo6CvV8QaT", "fXa17ArLwhbSBRx86ndXh2q53XSL67cqXJ", "fXa17BKvJMBJKrrLXGpv7fUZiZ38U1DPWK", "fXa17C3meJdNxiZtMRHU8VQ4Ui7rX97rd2", "fXa17DUKoxztUWKn1h8UFnoqpjKjysmjgv", "fXa17EZ33p7AqkNg3iCupKWeiRWX2wZi37", "fXa17FrYVo9XJT2i5KNJEBNDKy4xFiZG2W", "fXa17G85t9t48JerNwATn6QdHBDDqujLAr", "fXa17HSzwfJJf8AY3xmNqfonmHXm9GNY1Z", "fXa17Jq8Lih6WGJAGosTaLrMt4fqZUa3L4", "fXa17KJz5Q1GkKfkRnBaf3mZB3TKF5PvVg"}}},
// odd trading round, pair 2
{{{"fXb211z6Fg8pbFjXCr7Ge5A47zcVkzRzRS", "fXb2123rQFuXiBmBrmZQFZtJCnLSvYLN9B", "fXb213Qg1bdQY2SFMd9tVFq122ag62avhN", "fXb2141ass4CM1rokmhHjjQPpyU38F9g7y", "fXb215U3Y9r5oacbB1STH6KB268mJYtqbJ", "fXb216KwyL1dA22fCkrshUmdbjPxz4WFB3", "fXb21793TR9SWHtL77BXRaDAv5kJPjQyEs", "fXb218Rm6egUgx2GRtW9RBrRBryaLu46eB", "fXb219jtq3GZcX6dsQ2DLo7xZx6VpyAszD", "fXb21Aa5GBHZKg8ZmmEUsRyWwC76FXnF8Y", "fXb21BhU269E7FK8BWrfVv2twhaUBVQ7Dd", "fXb21CNUBqSkb9qW78wBaYhinCZ4dKNNAp", "fXb21DLhUr8JHaELsV8iNGcRwVKw73a6cz", "fXb21Eg7ZwJV9ad3n9bLNJiumktEoSjU4F", "fXb21FYdzrCKSSbymym8cES3fYNhLMvhTK", "fXb21GyCedmKikDSv6uFrEEr1GX5YY387j", "fXb21HP2g4LGdS9mWa9xFWhu28daWZRF1V", "fXb21JNnTofAKbeyhPUe2NTTba7T2gU3ew", "fXb21KrcYsQzNhVYemzUAcCZihQYXzn6TB"},
  {"fXa2112aC34Fh9D6JDNLwXRcJ7oSEXXFhH", "fXa212mg5skQzmgYMQVpHyzkoD2UmTPQoT", "fXa213BVFmdtiZQwMp8QSgATAmqVzp6m1J", "fXa214xX4iWg9fpT79YfKVogo8Zj62rpj1", "fXa215UyDSn41EAjHbfQoiQLSxWeN7sKvM", "fXa216mHJrxwUB5L2LpSPzdnrMWQLcErKA", "fXa2171U3fGYpi17xqV5aWbH1AEzEwNhap", "fXa218noCgSGoaEUcpjUTFSeH7C9gjZXg7", "fXa219BcB2xqKBNm64aV72Pazy7GBPYH2R", "fXa21AimtELz427bPzyY3YdHjkq7p71H2o", "fXa21BM1oMQaVryVo9kuBTF7xYTVv6hrm6", "fXa21CeE694msLWHSR1VnYKWGsKzpMSH8X", "fXa21DqjunHS5XpBCNoSvyYMTqJiNokWhr", "fXa21Egn74PYxpMZi3MK2vFh8Cx887H29S", "fXa21FZHEqNaGpm9DD4VCv1KBGFsNLjZ7g", "fXa21Gn3426KAqYsiYAhmQKcNHADdeubBJ", "fXa21H5CRfCikn6dvwbzhcb9nZJTYb3rL7", "fXa21Jq5Y8w9Jt2Ns1CAftJHMF8Q8x26vW", "fXa21KZ8cU7g6voVzd79CzrVJvmA6eUeZy"}},
 {{"fXb221QZaowpw9uDsAee7uXpUs1tT3A6zh", "fXb222JfBrWKApSxr6RVq3LPeeLLAF5TmD", "fXb223sxj7uhs3mYuvsKG3GSJmXr1V8PVE", "fXb224V3pzfDW9NheYWiZ2oN6sbugmnD3v", "fXb2255GwPSivEpREuovXx6Eg1DrCRskd2", "fXb2267PWkuPaCnWzds1ciMLGAiQ7kvZFu", "fXb227tobLn6xGMH2ybHs9q3ug1zbVt9bV", "fXb228KVd4cGB4zJV7UQW2tzYJg27f9xF2", "fXb229WmY9BUqckodLgYssZAvfGWESBQRK", "fXb22ARjbiVrNFDDosbxcXxA85B3B3y3wX", "fXb22BAsxVyWjU2Ne1EfnBEWKhQ8ctt2di", "fXb22CmtvkdsfMVFscB559mS4HpUz3Xxqa", "fXb22DweMmzSU7VFJc1pRvRnbFi9uqZQkM", "fXb22EfeHDzU2L6KzRfLt6hAr4EhYhRv11", "fXb22FZbfhoZC3VM1YR1M8ook6zoTmYzKZ", "fXb22GbGvN6ARWT1c5MdUE8cCFQ6ZFooxk", "fXb22HbZTnqdiftgosesSNuz82ga1Gp2zY", "fXb22JrDETGnQp912asXNw7cdM4GiBL8iQ", "fXb22K7vGDvAaqy1TG1twPb2tnCNmmpB8m"},
  {"fXa221qqBZGjsFCKqgzRbwR5e54qPnMY6S", "fXa2224y6UfGnuFpEcafA1eVMZd5QXc34f", "fXa223eJKT66AEsi2949Vw38zroyD5Azev", "fXa224BMebd6kX4rLvQqh4pAQMe5ZhckNE", "fXa225uh4haWBEdo8XsoR69H5Ts4AoAe2a", "fXa226m2igzmG487TY6xareTDaMsymdG2X", "fXa227XJ4QYZdAmxDYgYzpqjrhPCxfrjac", "fXa228hTb8kUBTciy4FJ2VWLo5kiXpGAvC", "fXa229He2f8XSW2ukkFdGLs5WCxGtt35bs", "fXa22AXFiTBZmkUiftRiRyRsMibp3x1shi", "fXa22BuTmGPTykwmMrujwjx4KuVRS8kmqw", "fXa22Cp4BkaHNeM7V19WX4xK6NhhhGXF98", "fXa22DG7voFS4QoUyKQ3XLQLWEVgiYbYVT", "fXa22EiKrGGShwNtiJbvAYf7gZbje4JCEW", "fXa22FGe2ZfhpJk5jEWw7ADBT251peU6Vj", "fXa22GfhyoF2kWgT7yvGHQ84vviL1Y1s1K", "fXa22HbHG3jbLLQyrq93kPEYdKND2rpUtq", "fXa22Jd1gLKptHMS7eu3eVAYTMxfitMgnC", "fXa22KdWxF4yAzmMNhAzj57fVdHV8u6Brh"}},
 {{"fXb2312C6UPVwHZQdNizx8tTs72GTKFatn", "fXb232jP5RKXPWsGW2253iFhkDSAt31FLS", "fXb233JqCMPPbkZZSU19ENcoVpBPd5yHok", "fXb234MDSDdBjv6spkpv5KtnPi1oqEH81b", "fXb235tXbDVE2EgUzSTWU1aRXMMWe7iC79", "fXb236PxFnLyNncGsMwFSaRZJJqDj6Pp1q", "fXb2372qyBWDnC4qzvEa4pzLYapJGorZq3", "fXb238HZLaSqKsum8NTvrF1Jk8XvcMZ3HQ", "fXb239qJLrBWtLcLEW6vDN7vkTxnNPJCpJ", "fXb23AAKiKP9a7yhm47UVygUuJKGg5YuwP", "fXb23BWWvQRbZbFgxL1LxYd8LCtX9RBDTw", "fXb23CrQ9ToTz8dcukDtcjSjRosBvjT6Kc", "fXb23DeDxP6hrGPL7C4W5qpAA5dnxZCFjL", "fXb23ERjuW8XmLg2CyiSCDGkisKc8fvzZF", "fXb23F4fewbzQDECuoLrddWj7ka4u7TYyB", "fXb23GCo6a52v2UaJixGYocJovyfYiyTAp", "fXb23HenCyncepEzpLM5WfR9wEwk8XNk4Z", "fXb23JSsvsqe7RfewMV6qKMyFYarD8Zgtb", "fXb23K7CD4YiB82VBRNMx9LCmVRbQjPtP1"},
  {"fXa231MUtg2BUvVzpTw57m1M29GDxfv7r7", "fXa232a7PTkpwuxmSdxhUWWwue6A6yfoJz", "fXa233HhNQPccMDcX1VUQmJdyWvDubjKzZ", "fXa234999YfiwcGTxJTag3r3FAyjBSHjyt", "fXa235uVnPQWt58wHH14cxtt6EhcTyqpqx", "fXa2366PfJqmUbBhpHiA9oM4Sdm4D58pxK", "fXa237uRtgAnzha19u5biiqNYZPCMuzQty", "fXa238kDzZLHnhLC4QwZqWPxj4Pm9Pam2E", "fXa239JWwuUWq7Fo3ZkguQEpGUx2yVqQEd", "fXa23Ag7yL1ZTm5dT88r7biLyLWfbsCDK2", "fXa23B7CCFzvDqbRkftT7F1fT2KZCvJuUc", "fXa23CS17LDqdCtwL3a7paDpXuSTsYRbb6", "fXa23DXZFF21ujygLxj6eSY2WaGkeEkybP", "fXa23EjSBbcQs8dGFsutdQCGnGyXeZhuzd", "fXa23F2MBYpy5WFZckNbtBMAsu1jjsWtwD", "fXa23GudZwRptQrr7H5rsD3LwGKAaRhTBL", "fXa23HHngQJvwJhMVCyB3YkUSKcMVtb54X", "fXa23Jdh5HxLP5fYsrnA7eL7LQfRrUzb55", "fXa23KRzS7R5UxgRwjJDXHCPtu3rRJFKj9"}},
 {{"fXb241mLMGiC443CsC7t7h8CmzNXt6V22c", "fXb242HvoJmwkcKp32XKxKXHjzAVaZMuQA", "fXb243MeJ8TB4N1XGQHuJkdgSTomAYhKNX", "fXb244Bk9egxENt5B9uHXjiJceRTQwsadS", "fXb245AfKkdDfFDhEuKuZ8PfhHmrMusXZA", "fXb2462FkFNYxMvNXxHBQUYEnEJg1ziZcK", "fXb247ZW7z2rt64mn52RCrTLuyXoGaFpy9", "fXb248U1PhTt4hJ1GBsLwnCtLdxrakMXQp", "fXb249Y36PNXSPSpgvYMQL1PxCnLiiHGJ8", "fXb24Acz6XT9oGyP1YTjmFVsFTiB6megWP", "fXb24BaZyRKMo1Dwfys53BnMcsH3SEP83T", "fXb24CWLrqiNyhH9JSkNBjg2kK4QFRiE2g", "fXb24DGkb7HKV9gZYWvdn8Eceuyxy37avC", "fXb24EqV9jnnfWAWKPzbCj4w8hExMGEM9j", "fXb24F3TnZKAQdCJtGPzJFdxrWVRvyraMh", "fXb24GS66fVkiDFCSccZrBaNkoeRMpda3r", "fXb24HeoDCyhcfD3SJZL91u5saUX48bgw1", "fXb24JkPjshRow6AE3LtuAA69Vfi62dtHH", "fXb24KoMe8Q23FCWHnLc2T7pZ5yGyTEA9z"},
  {"fXa241ZoASfwFeSjXhKXiMsWCdTqgWMEG8", "fXa242PjTHFKUACKGc1sA5BoiZtPaLjig1", "fXa243tosGdGoebn1SLRfAzZg9hgss5gNX", "fXa244DB5TAgsLanQfD38yAWykjFCCsK44", "fXa245tAr1vhtTefR7sSACEWVFmudbiDwp", "fXa246YLLxucpAWAo1zPMWBbxshtKivw66", "fXa247CmXQR5wF7ai99KnwNL8ErVaQmeDN", "fXa248LgT9krM5MASGL5LVTPARW6pds8LM", "fXa249uMbjMT4tGVWP9cxRKpA9MgMDjoWL", "fXa24ApQCmqSaxJBMRa1YeBxPTc2m5FFqf", "fXa24BUztiLZN1Dr8YmGc4Zm3rzSRbqsTp", "fXa24CE32p6YoqVVUJg2wAXHiuGNq2gLmS", "fXa24D43D2gDLocwPjWUSSKMHehiR1J8P7", "fXa24Efjx4FzJUHYtuMTGa3RfBZo8VXVju", "fXa24F4WeWku6yCx5rEF7u3EVthdLzCc47", "fXa24GDznfdoSGzsrqXJQKwTC51vboyHTj", "fXa24HtcipM17kbCZf81fkQdToAfct8Amf", "fXa24JAqbUeWfSNsV17cGH3PDMiqDDhhWE", "fXa24KB2yxAqzUUUTQpjmuXnQJ27AUPPej"}},
 {{"fXb251N4vsUQMXbAuEeGdzLQTJ459ZnfsQ", "fXb2524YMExCAuJkaxvwECChNtSZxvXVfY", "fXb253Wb9oC4hjgy8ikUHjam1eurTmYecY", "fXb254Rp2KVcMFVV2ibfiC7b9Ttsvfq8j3", "fXb255f5UMjcpPV8uKGmZreTR21htF9ASu", "fXb256QgSE9ERSxSo41tpm2bm4hrF1bATL", "fXb257aVcTR2Hkp7tSRESrrreUNkKa2ais", "fXb258svaznEkapkSzSRssVqEoxX7kUQrk", "fXb259mrgiHBNPLXeFZsApBwwEbNcH8N4M", "fXb25AJGCgd1xdDx7Pgfcro7ENuFeggYTT", "fXb25BSwZKxfpYmGuSg1gk58FEohmfY1vc", "fXb25CZUyagtwyCQrTrC2Qw8NXVSZLrURH", "fXb25DxtWRMC9czXo16Nk6boqFKwwp5h5w", "fXb25EVMKgBMNRoNfmEM8D5s83YtEaBmaT", "fXb25Fot5UbBWXj3SgrmkJ7aEMCQfzyeDK", "fXb25GA2JxSmE6sjcQFiUHoMFhhbVWCm8q", "fXb25HYRnpVzowpXRnsLoMc1HxAUVVh2Uk", "fXb25JvVdvmZ3BX9Q88WxwAi4Kzoqn8j6A", "fXb25K2Re5vnS8PosmNX3dE6SwKhf4s6uE"},
  {"fXa251dCBvHcYKnnEZzEaeALHrZDAMaN9C", "fXa252iWUZkeKSgMugTToER6VE3aKwxxPJ", "fXa253JmKG8LUHyqLJhCxoRKduJ98byvvm", "fXa254aN9pMPemGyhvVLX1hoxKPJxm6WEV", "fXa255TAwqrnPxTHoCHr1Kr7vRmR3ExJCP", "fXa2565hRhX7KYQcjNrRhiejfT3NAayxEW", "fXa257zcfKsxCXUk2c3wp1yogxMX2G4Xe3", "fXa2589NqzWncHRiZGn2Zb1u1TPUHHvcdR", "fXa259KPiLn8dCkCBnu1WUVsui28P4N8Fx", "fXa25ASjpxm3oz31UzdE2YcYi1AT88qmbg", "fXa25B2ZCYjxWkYcwg47oKtkm4Jpz6q3aS", "fXa25CSkzwvhGomg1BUcfNLy7F5EAQ4VR1", "fXa25DKzq2m86r3uzuGy54JADtKUPVg6jd", "fXa25EY7U96QZpWSiVzotmy6oREtK6ejem", "fXa25FL9CyZXu3REVBUpnXoGiB94Kn4QAn", "fXa25GfPtmiwd5QGRiz82v7FYFithgUNee", "fXa25HBZKPo7DrGTfApmvewSy9tFYnJQ2m", "fXa25JSZiQU5eZsoEgkTb5Fmy39s2TC2k9", "fXa25KpfzmKyoMAyRKutXtbkZcoohP4vmd"}},
 {{"fXb2617pNEtgsXyRJT18ppiR8a5WEg2vPW", "fXb262JpecybfBAiBmmjevn2AQKXFJvCB9", "fXb263U389yAM3K59SAqzHr8tVVPR6xpCL", "fXb2642eXrfkwBXew3kRFDgKiebTVTHDV1", "fXb265ugzVY85QKo3oF1kv3T3M2wJSuBKT", "fXb266TkoMCVceuv9NoBZmg3ijqge6yNrz", "fXb267V8tQ6LiNtUfFik8MsXUXEwTG1YpB", "fXb268N6Kvz4hQzGxNFvTcKmwoYNRvvcGL", "fXb269JKsRDmiHBeNyhRAwQeaGjQYHKszy", "fXb26AdSQswhUm6B7TbeT6QneEmfRPHLep", "fXb26BvRVtioMQ2d3XngfyMkzhXp6Bjscp", "fXb26CCq4mhiz2qGtdJT6LiYuKwfH1Miv6", "fXb26DpB4GiYxKvsFAeYZtJp4NZdZYJAmS", "fXb26EGYvrogmVppFB4qFGk73JYWw8cV1f", "fXb26FSRyUxJHDX8DeuRtLC23VW7ei9uhM", "fXb26G11uiZstPsf1dQ37tGrkaXKDwAxGe", "fXb26H6maqa5iVC4yyEbzxuzMLmCJZ5Drg", "fXb26JAL5hsGcKpW8sGwCCKkouKuWAR2Kj", "fXb26KY6cPD9zxxekGf75rthFnUF3dAAKo"},
  {"fXa261nAUYfENQc9sanfgYJv1BpcR9pKsC", "fXa262h4YTzTJuTC2Vyyf4P2DpjntUY22m", "fXa26316TUQRYGg8Bd8192UYdqM1qw3QC6", "fXa264SAmLVkJQA6qX4cXRSm4ue5W8NhVF", "fXa265PumoEepjyunYH5Ti9hyZr1F6h2R3", "fXa266A3atRTjTgoy8jAhZ7g9BsRypeTzM", "fXa267dZkmWUGgoURHyDv8gW86dTh6Cdnd", "fXa268saEZCenhm9ZmH3QTthat98pB5ZSZ", "fXa269U2W41Vj2nnDfVF98Ewq6Yg7fyefx", "fXa26AriHF18748urM51F748jiYehgMv93", "fXa26BFFjUPaGPPt2JCqmXHXFEi9pxm24Z", "fXa26CGkMaMmPLgksTYG8PABkUzD48c3Jk", "fXa26Dow2v3QQn2BZPoBEr1RnvUqNAYJLm", "fXa26E4tKMzp7LY4ghyUcdVGZBUz4s8Sxo", "fXa26F9z4DyrHUuBuY3JvsSXaGMYx1Tkjg", "fXa26GTc9qbehMJtstYXBQw6kCKLNtbGZA", "fXa26Hj5GgrUPnnq3TEz5kuG4z6gtoM3kA", "fXa26JKrcRGDLaCbBTSWuDoe9vsVvqMqRh", "fXa26KZn6WPJc2vYoGEZ8MAZCoj2h7fNty"}},
 {{"fXb2714qRaaSpvcfiRkLgMcYqs51KYgsGz", "fXb272Jx7kdrLdhf8LkoRxgW2s3RZEfwZP", "fXb273zd62JawG7mrWqFdwNx3jwp5Y9jEs", "fXb274csZwpyYkcbnagK2ES3zW8sPyuvZH", "fXb275ZoMmFLqqEnzgckjf2HtctosjhnvB", "fXb276THAsUWMx4ySxTAwTa2DH3nRbQkt7", "fXb2772FxoC88XCRm9ju3pTW629xrBt9Fh", "fXb278J62pnayhVnLvGHTvoK8tU1ELD7UA", "fXb2794WQwUacibyTgZSqQz5FjvZjfDWga", "fXb27A7Hz7vtEzhB2dd8nM3PYHwR7eTdsd", "fXb27BGLgawAVQ5f59zHPnkeTKzM37UfxC", "fXb27CmkxCzBdJLdipHgBuVi5AQyKe4Guv", "fXb27DqaiPDroCL3gruP3MijU3d3vQ6i2c", "fXb27EqQcvhdoWTLgmazdgApL6KksrjuVu", "fXb27Fn3zp1y8Ya4siykbM7GxRrGGfqdkW", "fXb27GgcUkYw6WaLqWp5TVM6U3fiBTiuTb", "fXb27HRiLiUq2iVTVhQzULn2tfnhEztaP5", "fXb27J3YrmHsu6ZmtMzvshcTLxMzBMJ5Ts", "fXb27K4uVByHTpdF1M5pUCSWvkb7t4chP5"},
  {"fXa271R5acUzQ4bHYQx4JYcAXJeQ73Kyv6", "fXa272N6bvHMvY66MNM5BfiRL9kutrVbh6", "fXa273hBKvkFz5sXiLoFWbZNZqwsbqD6fX", "fXa274hyR3GZyw7K8ym27A1E2wmkNWM4Sx", "fXa275NygwxY33VNuDsfwMpchzVJ4BAM9d", "fXa276qXSeogB5roXbiGaTFZ3W1SfMFJhr", "fXa277snRcoNNHTaBLo7US87VoxoT8fen1", "fXa278dmbvDVhihuiG7qJrU9Pbv5zEjUEy", "fXa279bGCosa1uN38mQZiJzCF9TRsAfZ8c", "fXa27AmbDXQiFz3ubuKT5NQsU6kLYq7kNY", "fXa27B6i9REqBaJJ4YtmxuShnBnkUhY2vj", "fXa27C9RusVCMinrt52wjgfPgBKBM3byVi", "fXa27D7eMRa92mvV9cK29DA3BeKXnSQPKr", "fXa27EFq5ZgvDVJETRKks5wHJUVe9dYBB5", "fXa27FE3aRboYJmsWdf31dRAVfE6i4VR2n", "fXa27G7HJaHA4eiYGw5473YujoGraMfu6N", "fXa27HRg6VytFcScMscZT3dyhVdgHpRwPN", "fXa27JfPqEmTGzeFE7m8xNFJQ1v49LYy5Q", "fXa27KrTZ3HQ27TMgGikoVDg5RpboqwJLC"}}},
// odd trading round, pair 3
{{{"fXb311cQBBrzGyTBgeSmUNvksexN9b9g87", "fXb312TnzvzCtt1RK93a3iUvVjTJ2PH7dn", "fXb313qYwAchMUCvbaCMnGFM8sWfkB9HAq", "fXb314fNharhE1v6DP7jFApWNn2GX294Tz", "fXb3158aZAKtR1qj3jfhFGjiWbdr8CUMSR", "fXb316txD7ZXgfY57HkNFqQtFKyy9e2E6j", "fXb317pNdpZyS571FGMz4NkZHs6P1hhQTq", "fXb318kma9rQdYZWsboSHk8BxXth7We6TY", "fXb319qWK46ozDkd3HqXRUwfYsChzN8XQU", "fXb31A2aWHoNADpto2LHTdNJnnkwAmrcD7", "fXb31BwVBFEhmqSUt7noGR3Qr9m9dHEUo7", "fXb31CTE5aafru1JHawHCNDLddngHyUGVX", "fXb31DdTgYLXqbL9wnKxnZxAq1tKXzVb1S", "fXb31EmhQ6YDDMA4bxECHb5VJsb4g7LTKG", "fXb31FEmuoc7Q3QavtFfmUiS6jrvoNP78f", "fXb31GUDBPqj6pCLYA1bsGz8He5vse32Rd", "fXb31HtT8RRAYn8Yo4T5rCmhn6T8y4u9Dv", "fXb31Jwni7qXVuL2NZBgCZ765Hqryqo3Dr", "fXb31KmZXMma6Jd28nPZn8S14KUj8sx4Yn"},
  {"fXa311nXpHgLXct7kMcBnnzuQGEswNYjeT", "fXa312Vr4DxVXckJg29giPZGDLjPeWjf5T", "fXa313g2cpxhs9V3nZmrykijEYDYAdgFWd", "fXa3148k1BcNNwBDmQURkbsCtaaUEdKVHD", "fXa3151yq9sRVNDkKhFPXrfkSuVg9dcvqp", "fXa316nQi1PRc9uZxHwdGwjCpJiBuxi3LH", "fXa317jyoXZx62vwfWf8EonN1vnh6FdLA5", "fXa318VfKKC2yGJ9FXKZ5UimHwa1M7maYC", "fXa319SUV5ZUXghB89eLK5fxYTPyFqKuUS", "fXa31AzKLtiCTi2NAk8fWVXxevUdaH1D7h", "fXa31BEGNRwsqMh6hgLioaz7XkAGP8ju8Z", "fXa31Cg4R56LWzPcoZymd4KEnMDZmCcw6M", "fXa31DsoukTA4FZErsHrgAmRMyQ1eTXYA2", "fXa31ENh5CPsddrL4RdpTfNFjFT57xtA2w", "fXa31FG2mvJkVeVLwg1jf87NVKJnh4rSAK", "fXa31GPf1a8qXkVm9Rw8NsBQpFXU61KD9L", "fXa31H1fKEo25YtHyjbB3NEUf7Hfdeg6bB", "fXa31J6dih3w16d1DziLSFWVMjMZ349P16", "fXa31KeoMuXDWEFxKdwkHsCJL94VpxLZ68"}},
 {{"fXb321YfVW3J6216ZASAqbYd9fWNRzQTL6", "fXb322vopUGSXUZw3tXiQDsuv26J1piqJ8", "fXb323DE1EWB1KZMom5bd9YVcWiJ49x6NR", "fXb324A91QyFf5KVFxNYZ695xd2hKWthPZ", "fXb3254LwKQLzyhdpJtqXzkhExLLjoLVYV", "fXb326TtXUjqPvRscnTEUU9Vie4KkvTLBi", "fXb327ZuMAYrrDGAC9TQjryHLWabHqtJJP", "fXb328m5sywKqQBRxWRd2rbsNuN7GfuEqK", "fXb329HEMUuESH9VDt1MUQTNS9zgVDz5sC", "fXb32Aw7eJDEmuxxJQG92fpo61WoxFyLjc", "fXb32BbtoTXQ44XjoSoAgNqD3HQ115vnBQ", "fXb32CwYXGxVS6rnTrN2NrW6QLNqxP3Upo", "fXb32DmvpErJdo5Ro2SGqkyb4hMrr7DkpS", "fXb32EZ8vDF4uaa7bSfK449LgHmf2dikHQ", "fXb32FFLivbztHP7bdadtP4caQaFja7rqX", "fXb32GHKzb8sLXpGcHN95tYtqMwZfZscoX", "fXb32HkU62vJtYbSzmjHeuSqJG5CT5rJ53", "fXb32JBDNPz1bbuanMyMzVPye2wiJqBqTy", "fXb32KJeubMwbQEXweVPDoyExuLbLcKoKJ"},
  {"fXa321YCQt9cdVfv2QbMDDmPZCZRjHWsUf", "fXa322YANQhi54ij3t5U9a53MhWnxCjNVw", "fXa323mPc5dYc21GwiT8ghJgfsuhef2Sic", "fXa324zdCGdyG65jQsi2yUz8u4JDpjLjxe", "fXa325BGNpBL96cPJb2NrZG2YHaYPWFs1H", "fXa326R8F8Tgf5tjmoB3knwaagaCAePKqx", "fXa327QN2mn1b5i5QoZXocp3ZDkbrBrp2n", "fXa328PTatxTaoavLawQpWpQzLin6BNWG7", "fXa329bkYcCGJTsvEP5JVJ49okKnDVJ3ED", "fXa32AX6TTL1hT1JTfKMpsFtQND93JQGEg", "fXa32BRN9GHrxqGAPbMiFS11T5aRfkRPct", "fXa32CZvZ989kz2nwgshhnTozY2QugZPxR", "fXa32DoATTsjz4C2xHLbqUwVYM9764kdFf", "fXa32EuvSsFiLzFuB1Q53ntXYgW9awcLLr", "fXa32FkKBUL1eeEHyvqM69CtnuFfkXKVuo", "fXa32G9iJRWivzRTAzPseiCab8Y8DnX5Yx", "fXa32Hz5Dhk8knH8oEkwHNeg3EqkpVzCnu", "fXa32JtE5P9TJ5uQ7DFWgMAjujLDswBGdp", "fXa32KvRUwWef5KZhovmocMognhbubzGmW"}},
 {{"fXb3317FpAjC8m8GZUJRWThG3epeQDvxrP", "fXb332qgN7bR3ubZhKA7NGzDbh1Mgp5dL9", "fXb3334YvL4Rvc1PKntWW91LYdYyvAeqHv", "fXb334S4MSX6TSsxTnB6dVCUh3oxJ25m6Q", "fXb335bd3R47kQdVfvV7ukvCk6X21zMeKz", "fXb336CNBTVV4tXHB2rVRcUKjmiNuQziKN", "fXb337c5VccDcZUyTvGbkC4UdW7uVRVYEP", "fXb3389W9Mzw3as2okWJkyiB98BknGc5ty", "fXb339YuSEFq2tA561rENANbQD5SQLDi7t", "fXb33AL1Zj9hPKSMNrJdDUkc2QGW8nhQV6", "fXb33BTqUYzrUkjJ1FqmSvM8XdxV7oiKY7", "fXb33CPmSbp9vRhbEDfQWv7WDujhGPydEs", "fXb33DW7BZNC1eAzjfuwb3tg3uUxA1R2ei", "fXb33EUJ5zGNxgV8sR58s6awNK3MJEKYpS", "fXb33FhRNxQL1oVzoeZmeuTmHsj8MDWpEr", "fXb33GsxxLXXf53GZ4WjukZfLgH4X8CToh", "fXb33H5JGTexuFaQv11dyuBsKVpmtEJqaT", "fXb33JMJiFCbk8qgazXG76o2ryNfKaTRbd", "fXb33K8AQq77L7u5GEUFwCAdxkunC8BfFT"},
  {"fXa331HeQjwNBAXGzr5hjaKGgUeJingzeD", "fXa332jTabyc9jur1sFDmqRvsZkNcjtSx2", "fXa333wEVhQgvZY19dxzZynKGUozjbwdgC", "fXa334qSTa5uVLsYadHNiQpHSMJzqQNE86", "fXa3351whuNqx4GKiR8EVub9Dw2b4ebfYH", "fXa336qfh541eQ4WwYCBzrFtFrDFgyeyKn", "fXa337GHQCHBckJS2BK3xNSiDShmqRYJDR", "fXa338ESpadPQbt4KQs1cjUwfUZfeYEVGN", "fXa339dzQ85TUmuHMma5sUXvu26SsXkLrq", "fXa33A7UQDgrZAPJqTcckBUweCDHVzUeth", "fXa33BPTM7ge6uJ4AXcsteuepFXBbm5rAc", "fXa33CYCe5anTcP5T5Yoy3xiZ1s3fqa4v2", "fXa33Dt9K6Gp525otkKGQ6XWFVoWgmadGC", "fXa33Evj5A1sj1d432RCpPVD6AiombKkEP", "fXa33F4wNCdA9bgbgZPyXY8mUn1FbT2xsX", "fXa33Go5T6KMrpcYDDbAQDPZGnq4oHp5GZ", "fXa33HFf3A3JCezDCFmyzWAcE5BU1iTpzn", "fXa33JwCiJTojeidovVR1GdP2myC6NVCHj", "fXa33Kb9U9KdYPHW7gC8Ysxbdzy9QXJGwq"}},
 {{"fXb341yA6xn7Yccyn2iq2fi6kYaShDUgVw", "fXb342VkAdSJKYFQEwaNsLmHQPKL18oTeS", "fXb343c9eVdixQpHcLpinB6QzMZTrXk3Wb", "fXb344XMA8MYchqrJeQ4TSwpzHAYsYxHrB", "fXb345Dyb5uuvZhF1bEFyAhaZkvkEprqbX", "fXb346DHaHJiysySfo4TLu3NK8YAFytphH", "fXb347ZeQWkyyAWgj7Fys3kAgDtuLQ81Do", "fXb348EU3922Q3BzgzpBUYeaxg9ZQQB4bG", "fXb349Rf3TTzADVT7xyRdQN3ogZ4T4fgAa", "fXb34Ah5inzJsEcRCr3GVSdFwGmQfmcP1b", "fXb34BzR9HQctGWBeypadSjoDiXVu1TwBd", "fXb34CTpYEyFZP58HwSrA8SyMg7Crv54pW", "fXb34DDZuD7Dgqv1qqyYgpRXE9P5JmMLdh", "fXb34EvcFjuMryAtYK72thzzP3juGnRLyz", "fXb34FHEfyz2c97CgsQhHzVkSdrR4TgD8a", "fXb34Gjnx1sEpynwGXJR2aB8JyhCPG8uaH", "fXb34HTHxzkPbDB5ocjJeMdEodQF89PfGf", "fXb34JMSQgcpsbPxB1c4AqtsVtpGsgmpbd", "fXb34K8PKxoJt9MQJYJ7KKMnCfzM7KMEjq"},
  {"fXa341aw4agdgFyCsy46zRyDbbumWzymFP", "fXa342kzRJb5MPNATRerHJNge48kYsqdty", "fXa343xciWXU7YVbmkZuH48Lyay7HAE8ji", "fXa3441Au6TpVv5xDQyacGNLL4Cp4YUw1K", "fXa3454ohqRjQo8K3zpMWoyezXYUT7o5i2", "fXa346wnWZajhCa1gXsqUhGVAimRPMc3tD", "fXa347Wbs8udpg8xGNTUeMx5WvHRVwNC6H", "fXa348un9zD8TxNqpN7KfSsiexLDAVfmyk", "fXa349kRZSkZXZrz9kjfRjTiJfDbjYfStP", "fXa34A53TzPYCdshdTdrHPWJRTGZJUydd5", "fXa34BA9FMSNBxN1BeBFxjUacg7UwNFoz3", "fXa34CJfdaSPmawcKUQhc63575oRGv8brM", "fXa34DPLQEuQJbeFjuaMtPasqrwtyq3NY8", "fXa34ELrVfivoWWyJqJyZmjXhWSSTxPR2J", "fXa34Fyb3VzZBR1bjYaPHrZwPhVanaEkL6", "fXa34G7fYYkwaVxS42mwemeHugtYRXx4dr", "fXa34HvA4f5TaSb5d1ynachwznZLWj5yKv", "fXa34J8n1QUcfg8btBbGmHSt1edikZmZ6C", "fXa34KS5CjH6iEBPCcuCa6xPky57tLyVkq"}},
 {{"fXb351ggtjPwDdNAKLkjkNsgXpWfGahiHP", "fXb352rFBvQPjB1wz52ZE1pHeH916eEJZ1", "fXb353zQnpGCz7ccXJt4apgBzKRoHM8FBa", "fXb354Bv5xzGJHD8Cp6py7Q5jGSmRf9Qs3", "fXb355Qw1dZZThwvWCGY7M7qk7DGuFMpf1", "fXb356eHH2JuhW9nkkPdcS1gbx3JLyUs7Y", "fXb357gZhgg6N79SY69WzhoZmpkW7k7L2o", "fXb358xSDMcJ9iGy5Pxo66S5uRq2h1wjUP", "fXb359z8MetF9wKumskNQBEmF3G1ehtKfj", "fXb35Aa9iXMG4wCk7Udw6A6NTCPEcgLrW9", "fXb35BkgXVCCVrcRtuahpFWBKitGHx6Lno", "fXb35Ccd7oKmmDrSEqUgYS1mvH2oAqGufs", "fXb35DvifRKCDZemrN7hHd967dPZZMzGjn", "fXb35ECQ788ZzJQ1DvKgbThpZJVSDqVYao", "fXb35FuFsbmGwu9nAdG7bPrgC7opc3quoo", "fXb35GfzUG6dGSZftouwtKhupXsukS9zzg", "fXb35HYPoYx3XFdY4u5dWEZCQDYdSvYUd9", "fXb35JsYnqe6jAqkd51YxSzCpGAZHdUde9", "fXb35KNB8TJaPU9YyEUXJ1jovsJJZ24ZVU"},
  {"fXa3516UktXUnEFacuFVLnue1P2kpy6QsJ", "fXa352ujyZbVYsoTjzrJRBnNGHh931K8Ae", "fXa353Pw5SGuU16widYSaWUZKa7jzVYn6h", "fXa3542ae29V2GRZ8HKXgQttpmbpP2dC8E", "fXa355WXU7KEP7GWqe1FR2bkrFByzF2bGC", "fXa356Urh46zfqDCBQM2PvugEgn2x7QShN", "fXa357gTpxBPmk25BF76zNLVdv2YiBTHm4", "fXa358aYEFgLQwk7yLuqvcDB9hvT4gHRUt", "fXa359CWN3EgSvg7zRQ6LvkUh7XQ6rSbGf", "fXa35AZnFtdEFhX5DxeCb59seRXZKafU8u", "fXa35BvaTxY8Y84fvzKKGZgCfL5yp24vbN", "fXa35CmjynfcEGLoNxix6e5RyKUXyE5x1K", "fXa35DKqf7qDABbrJpo7GbYjDrPmpKCJGj", "fXa35EFZUKfGPwuHcBsTSmhF2uYyPMuLqh", "fXa35FAGS5EN6w3QRrVrcc3P5NJTaREmmP", "fXa35GGw8vk2A9aKpZzYYxxdTirtMA5oiD", "fXa35HHQHvGUrCXYARndeDfeeAxPht8619", "fXa35Jq3bGBGDDbA7ZreKBjiKMUYxA8Xpv", "fXa35KE7atKuPiDkEF7zt6rm5TmFRRvDqk"}},
 {{"fXb361Fxyw1TpFq6r2E4y38YaHvBTczYrZ", "fXb362dhQk1AN1SLkXpUhBSRd3opP9BF7q", "fXb363qoB99q3eybTX3Vbyc99VntfrMkVS", "fXb364Cd6hDWV9nt3EZb1hViXPVEsKcS3q", "fXb365BjdwFRziNaX9SNKCGurQ1udaEMop", "fXb366sQZyYxTfj1ZTxMpp2jCXrzD9xiSB", "fXb367sDBsVKVqcDfG87Ye7s4n5q8FhK25", "fXb368JZWZEA7RTsf1fJ2WNxpExcMtHxTc", "fXb369PEJq76Av7kFo9GzqipeBCRixRmRg", "fXb36AneqXZ6hv141BQ7Hvj75bMbMZAfzK", "fXb36B4bKEUfnjJRE5fiuueakrpRV79sa7", "fXb36CL6feZFtMxSW6SGZMDoSFAyPfE2Cn", "fXb36DpDLsuj2beSwYvqeYftNLqRKvF6E3", "fXb36EDwSTiLix6xaiK3zJMD2mYqU9SegN", "fXb36FDVpfztF6BuWtzBh2VXaXoJYdWvpA", "fXb36G3DMvGGrPx6NB2LgQmW5DthdxYMMS", "fXb36H9snSh8ogkN4zKgYNdzzMioQGJ7CB", "fXb36JryKdTKMKeL15Z3rsZPuxvagfNtyM", "fXb36KtnTMDPjvNqdQHew6oERM4cy1SBoS"},
  {"fXa36135edqLt6mceSjULSZcM3ij4n2Upw", "fXa362Hneyyt1YLCc18Bn88s7BTCQWRwLM", "fXa363yWJRMmUu6wMjAfMFf41YvARxAzN8", "fXa364MrXDdcWRf1FJfv7TmyXq4bu4XEVN", "fXa365b6WZCML9bVUAmzEF8iTuiD6M59fF", "fXa366FcCE39Vx9KdfEdFWSkgsevfu9UP9", "fXa3672iSB6Ri9BTDCs1vay2N6S5XdwaTh", "fXa3683HB9iizkqNYJAs2r6wy4eXNfqgKF", "fXa369DAvDwnJgsPokhEBxZDUT4uDsUtht", "fXa36AVrifGbcsiJRVzCPV1b7EFa3BbbWM", "fXa36BdwGzTMxPNc1PWXQfVs7utfedPEJB", "fXa36CdfwRW1opYR5mQKvJLXfGvwi5ZvaR", "fXa36D9Yb3VZ4w55GeLNYaDEGLfbWSWaya", "fXa36Em9rEDGqVd5HfiwV1cUWr5ipYTCo2", "fXa36F4zy6HjkrePbWyxZXKhrPYbERXbxz", "fXa36GpT6BemUjf3TAkAMtC6CDk5a5EM3e", "fXa36HXcJMdFSfdpMQ9RdtDUsBy2bC3w6m", "fXa36JvBdACr2hRhS3rTKzJWf9tZmPDFnz", "fXa36KFms6AY67sqHHzLEt7Q7v8ZvvYytJ"}},
 {{"fXb371V16caGyWZek5W58aGFhAUgVAcPer", "fXb372uSKYuT9DSqpYLG6Q9GXTNWQf3zV7", "fXb373MxiR5KspCtfwsxqtXzjqqw2N3gBK", "fXb3745tN7gBLaoEYZaNVk2HaGVBE4C8zE", "fXb375MiJdiJY8ZU8Z9Tur11L1s3xpmdea", "fXb376wuaScKZAXo5qbPQfTNuiZAkhGfhT", "fXb377viYtpjKbSLNggFudELrZg63r7KxQ", "fXb378iqYKwAH7r2BvYedhpKZoC9FEaLaR", "fXb379BhdhHsZHgqcyp14Dj6uu1BLmvhfc", "fXb37AfqDcJxLgS4EWN529qFbujAqLVLGT", "fXb37BkbDQaaPfWkhEujWobfsNX9vWbyzV", "fXb37CEAHqL78dRR4ia3nyZ3CekWM3oTNj", "fXb37DFtSNHSP3DiFtxVoLjidDdqpLhiQA", "fXb37Ea1woLr93uHRjprXRfJp7rgUntYUX", "fXb37FKQ91LN46ivqMMzTQcfQx17k4ENov", "fXb37G7svjMZERQ6CCeX3hzf2WRjN2rLvS", "fXb37HFEyAhpvn9UyHFGTXRJSThRVBhtWm", "fXb37JmmZvLBvsCBoEM2Gt612nu8gFvCwv", "fXb37K6pmxow8H5bdEYLntXqjXsStW7gW1"},
  {"fXa371NRkPEm82xvTcvCSN9Hasrsd1La1Z", "fXa3722D2WyNayguPap9JCevU94S9gvF9c", "fXa373S3JmAsKPPLp3Ac6XHkUDE51XnFUw", "fXa374S7EPB3giEzWLYg1PSFp6Uy9EP9xs", "fXa37549m1hTfKBZma3Cj7U68ohD9UWEQD", "fXa376tTPyz1vcdAkmcHAP8KtaifTwGJZf", "fXa3779X5MdzoGhza6AyHYw3pW76yMKwP3", "fXa378Ko6cbWdfDBAetqUE3aGFG6BU76CA", "fXa3792dw3ZCQavAgrM2UBDzWtz8MpmU4H", "fXa37AjdcgRuNiPgjdk79AgjpXXJUuCAGT", "fXa37BHwnXBGPZv1b1baDnZpjEHNrRBLxV", "fXa37C7iMsxtg56SuGJQcsN6NbTGmAKE4X", "fXa37DCyzW9SCV7D1BVTzqApi9DD9VC2Zi", "fXa37Erj8HJjgaYW86w2ew1QaznWdZyPf9", "fXa37FF5tqSBL9pTSV7TCfVWVqejS7Wrao", "fXa37GcGwuCfTbq8jYweRrfCWXFFxkBbiC", "fXa37H2KJDsBBAzvciJHLX8rnSjtBbPsVL", "fXa37JqdMCXrA4LNZu7zzRpptFyZz1f9nX", "fXa37KcsCTgzzr2JeHpiyedYmu3U4RLGKW"}}},
// odd trading round, pair 4
{{{"fXb411fY41Lh94JGVxSMMTcE3ss1tRkdfz", "fXb412pdtA3NAJFXXzLozJYue2LD7qUaQ7", "fXb413SwJomudeFg1Q1oVzMnYA3FjvLbbV", "fXb414HTNSfBejz1K7QVQhTNh68kq3ov4d", "fXb4157RnGyaHMEbViKBzBEq3Y2LRs6e1j", "fXb416636hSdxcyaGt7yDn9snPqsGFRK6S", "fXb417vgxxSYG1tTFEUgU7Q6VCTGHQggkr", "fXb418Veex8GQFj4mGW8t1wgNsbvXKU471", "fXb419A5eDb9AMvRjoGs6LeWzyTz8sCKNM", "fXb41AUEz7krwTW1TgenpGPNMT7Usr4Xqc", "fXb41B6bz8ymMq6CkwBZobdA6SzwhPANN7", "fXb41C2cqrfT92XFBBpC2yjyj5RgyfndQe", "fXb41DczZ8uzcW8cZPXmcKay9dNhvMWCEN", "fXb41EcCTWhkPg1uBzsSQevBizKmpAHx6V", "fXb41FjH8SzZPxgLj946mPWNWNbEvRoXWn", "fXb41Gy7joQhfs9jL66c8pGrPFChcaEd3A", "fXb41HdW1YxZLR2NVZ8qSwAb8KznWUXpev", "fXb41Jggz37uh6xYVtvW9dhX5SMjmmJpod", "fXb41Kp2hDWKQbPwpiwVScbFGNHZaCpJda"},
  {"fXa411ubUsgvKbK91yhRteTX8UmZUqteXL", "fXa412VuRZBNoqXHqLHNzH6jRJcS5Se4wv", "fXa413KDVofgKDrypixTDq5rMYxnfwBBYt", "fXa414GkdDDoq7aDsyyj45psGdxnWLL1dE", "fXa415wCEJoaVcinnRUpkQen9XsL1ccEQi", "fXa416DpWdYamsWnf9NsjfpHr7wy7P3dk2", "fXa4178iF5dGjjhN8jrAVD7ihFLsr3UXzw", "fXa4186ShRornsWfyjteBTty8xsZdBVt6V", "fXa419Z8AmpF5kBG529BDMJSK2p5r5gJBq", "fXa41AkGiBkNQPCrr8VWcMDaa3SAaGGF72", "fXa41Bp5VGG1o8BtLphmwTFEasNR5xbqpo", "fXa41CF1SqF2j15UYe8QZBiFYjo4SfbFQ1", "fXa41D4CPQyPXita8QEC3nvEUmxEu4eNNN", "fXa41E522Kigt5iiM6HCEs1Jp2Vq2UG2Yn", "fXa41FzchVZd1EufSTSPXMroAwGNSGo9po", "fXa41G7YAyiJ3NMGAwUsKfGSppjv2Dntwx", "fXa41HsvBNs2VyBjo8N7howYVrKZP6RLLx", "fXa41Jm5qh4iHqq6zwKgDFTgteMju5vReP", "fXa41Kn51PkJwwbgU5iiJpSxCxsKnWyg5i"}},
 {{"fXb421Mwfm6ShjuM7oxFrwgN2ZUgcZdE57", "fXb422tmqVTmgkF6ZoPxuXYRjHXxPc9bEw", "fXb42326MuPQgUHXdvjb99WwYc9fZz7jj7", "fXb424mtC2Jecg8kt9RtYnWuJY4qAt5hrz", "fXb4257LfToiPE3TMRvyK8PbHWU4TCDSQM", "fXb426TC9YR17aoi9nXrjoxSpESiZz3eS4", "fXb427bve6jHgPfwMG34fVTBAu1cwiDFYv", "fXb428GZFcdNwBuSUojryZdYn5RTJzXmSc", "fXb429hGayPQbSzfLvS2nPHdBVcYsZ9w25", "fXb42ASfDEknhFR8ZhzECPtFedNB7mSboe", "fXb42BQYMCGqCKZNCg46jHdvxsp55imUxm", "fXb42C66ydv1jxjLA4asC5mJUwd7xbTbAP", "fXb42DTqBAsiRFupVnaJFU5vkWTXBqjfMh", "fXb42EXzBiwdEc9t5TJNK7rqfPA9jwKEQY", "fXb42Fh1f2SCtCNkHV3YMTtc3BXYdg7K89", "fXb42GuqiQBVZBG8eRGqLDmdnpF7oh6S3N", "fXb42H89qu6iffuD1xjiziLQTxfpTWb5jp", "fXb42Jggni4bQC8rR2AC6F9m6bL7LKRcpi", "fXb42K6WmppDg8VZiZMaBxWvqMjw3XWLD2"},
  {"fXa421dA7tekcVHKX6BaBpBePWsqpVf4uh", "fXa422XMgARjiEGwxkPrf1ABKDPohgLu83", "fXa423EcnDKNyoarYTEYqjFctQSbRpdWd4", "fXa424mqbLq7MHiS7ZovtPyL2XvYDDjeAB", "fXa425yZTz8cEbDhDiJ7uqHe25VjyrjEGd", "fXa426VokhuuiNdn9uMd5ZPNoWQMFm4zHD", "fXa427ngq7EzjWmmeRsNr2DXGkDcsjNaDg", "fXa428P8H2crrsWrncUq1CuU8bmU6k5NmF", "fXa429qWYM53Ahit893PhrawKdWsRvzMu7", "fXa42AF6rk5tVsKEGyhd8dT79kh5L7K3BA", "fXa42BZiSi73YuRdDLjNQsfiscx7SQ9n1C", "fXa42CnB5LMXPoSCMx1fXKwA7eot5DZeo6", "fXa42DHYBNhizXV1VihAn9zyuezW2md1nJ", "fXa42EZexFNWW27AA7xce3tgPDSKdD5xQy", "fXa42F5BAxTtdo1w5VJwkuBrECde8EEyMF", "fXa42GEsEiSzyB71biSvkgGzFxYUSsRuzE", "fXa42HyMswHM6PFUexAiwWnECRFk7D4Gmn", "fXa42JzGKLdsJomB8WuiNozxPCYFJEJZEB", "fXa42KELms9up7GTjaryU14W6ANSiRBFsU"}},
 {{"fXb4317cNJuoEJjrvh8arwSUMuxQz3v1kE", "fXb432ss9SRwyoSkxL13vHonAcz4NwSThR", "fXb433bmanzNFuJD25XBYUnXSwGrS1XTUA", "fXb434jMx5PHHhpAWVPro7dLQqrYsPC7dd", "fXb435gc9Jqer3sMm2b7jCRyYR9AwRAcsq", "fXb436Nw36x5HPbx8HzBAxryESBXtGXky3", "fXb4371nDxMavDkyMmAVQKUhjLeWhNKwiZ", "fXb4388bzjg2jctMZ1jKKivFrkpbZDTqZs", "fXb439zsKVqKYbVWWACCpw1cFFp8nsdJjY", "fXb43AkxjZSAojaEWYpresMd5CJbUMju5d", "fXb43BcTYVNUSnx2TWZJ676ufFQMGwyHpj", "fXb43CkGmBL3CRQbGQLWZcTM2FbiwMRrN5", "fXb43D5YMVapHgPJKC6CrozbtQdPpy9H8e", "fXb43EXuC3xLfLmAer9tCx8EmUqinjAKB6", "fXb43Fzh7V3BwvZ8SWPV9MU8GiqwNQFVNn", "fXb43GGgvWdyoZRTfeGtBdoFJpS2SKNnNH", "fXb43HtiXYn9s8eSmy4nTureD4b8q9sQ52", "fXb43JbxFVacx6ides5y9VpzY1ZDt2skR2", "fXb43KRaRKgrqb8LGR55qZJm65PSS3YTK1"},
  {"fXa431Jsn8rBP9CGqzngBYn9Hd9TpWWQtk", "fXa432Lb4zN1pJzr6G2U9o87iCjNuzyma7", "fXa433MqvuK2cyxAqVFxrUbtfM4wZKSfGQ", "fXa434agCWxCJftEgM9S56KTFDWJbiNK1z", "fXa435cg8yJ7ZytGSWHxwn91tGXWNuGgFg", "fXa436dvZLDxrVL82pAVzCNuzxTm3j9hrL", "fXa437Ex7os2Hin3adUgXeXAp4tY9HZwsv", "fXa438mGMqzMjpWMyn7YvXkTX1BW79bME5", "fXa439mK4ooZhGMGacYSSGAbPHyR6wQzms", "fXa43ATzeRbkhNTs1bhikfm2hun8MGztXx", "fXa43BcsbRCoJVRBeqQx9fxVYztTiB8QNe", "fXa43CurowDbrbEXS7rBw2YXEtcNEu84Ls", "fXa43DJu1gB7oVJU9HRVDL5Sb7WHgcAN3D", "fXa43E4n3WMV7qTMF8B8RHM5SsT8hYYbQu", "fXa43FPC52bP7ZnDxK3MqguiYeZ48vY2jB", "fXa43Ggq4PXpSj8yF9st6G34LhQZETaysY", "fXa43HC2a4ZUweWTbSCiBfKK8MPSBabSrS", "fXa43JiLpkEdy5gumhF5WQ1DSJkZrdePdi", "fXa43KoWGbV9VyyuhyLKUcmDAg1vDQCNji"}},
 {{"fXb441dQjPLBJvUsUnbbYX46QwscQr5nDo", "fXb442LoYftjHta1p6X3kd1Ss6RgJfq5Bk", "fXb443fDMFM58XNzcSjiwZxU5AVTDgxibv", "fXb444S6pvCdwXtXwT2iMjzNkL4P8CH1pQ", "fXb445U34NUFfzJ1HXN9bFXBwUrWid1imo", "fXb446gASrtNJeAqar5rBs7ScsbCCN6kL8", "fXb447qo74yn3fHjZ6ESqBXB6azERoLsoN", "fXb448oDcMSm5TCtasZ7hCCy7YzwTn7RZP", "fXb44967qm8GhTL4hPkZYrgPaBVMtPHUic", "fXb44A2mNSXdvp516MCTpMfgdEWQcmYK5v", "fXb44BqTTaSwaG7ZNRptYFiksSF3yV1V5Y", "fXb44CQqrBh4TDQffjhgcKEooQtcBcYMnR", "fXb44D7FdGZp8z2ruQPYA9t98zVtjPXa34", "fXb44EuJ55wwcXiNjSwpzn1TzfMQhKawW4", "fXb44FTWExfRoU1wnVer1RNsin3ApYcPW3", "fXb44GmiExXfr3JwASbgqC13Z4eBtKLKRJ", "fXb44HRyiF2nAz5r2FpSJ7KgXWejF9AUAx", "fXb44JfPKQZx4dg2Tv25mWS71Gqmj6w57J", "fXb44KJ6eKgVPHfPNQvvcCdmU6BpBQbvQN"},
  {"fXa441qkF6Lsw1Ec61wR96iivKgP2UfppU", "fXa4421yTdeyPSTupDYyVoFmswqnLAxxUs", "fXa443aJW4bBvv6NQJUr9KGnYduedRiTs6", "fXa444LjDXi1RQ8HQ5d42GKLEWVHKZHouH", "fXa445PUTNfMcFD1CJB5acWjNXXsULX6Lt", "fXa446VY9g8ygNMANaDntJHr7Uu66JyA7x", "fXa447b6k6aszQx1yE6JtQTU384vsmVUgn", "fXa448h9DfFQTFdnpZMU3jLc5Tkwxzii94", "fXa449R1NEruCb3dNP16BTGdWU2RXfLGVa", "fXa44AdF8fEoATJhiZ4Sf5aDSMkfCCiQzk", "fXa44BgdJYQxLUkL5H1tx5Zvi1iMVR1ahm", "fXa44CDpANEC6BsH79p5eWbZXgey5J5GZX", "fXa44DEFqpApFvnrbnenTZ51Njbtu5bHja", "fXa44EiKWVi85wcxie9RrQ2NnXvnyR69wW", "fXa44FMnvwqJUYaoMCC9gUym8xGmDxujUX", "fXa44GaW9Vw5xhieTCea2gFmACLS9XfnyZ", "fXa44Hc9rnqBaUJVkkMzkhpyeHeb8DmdBB", "fXa44JU1sZ53pb8qxgUytQyE2RB1NzuSQw", "fXa44KDzuavRQw7FH1Cyin9YeNkqWNx2jJ"}},
 {{"fXb4511j5rdEMxb8EvWvgzUcz9YSM3Eixw", "fXb452ay7SLEtkDJAwDgvshdaXSs6na5gE", "fXb453DbA8GTcCAcRhGJGg5ZVUUnS2XmZi", "fXb454L5j8WDpp9zi8VcU9ZDU64Be2q6Dq", "fXb4551WUGTt8anUWDfJkorED9J6sEX2A9", "fXb456vcAXdfPtXzUHTwUAeCvmTzAQqg3c", "fXb457NAZ6qx3JhuSA7Yb9vnCSf3ivb21L", "fXb458yU9HcMnUhD1QmoKSdWDSSfnHcpV3", "fXb459BLjnJqZv7hmjpr48WX9isX5U1Pb2", "fXb45AZRLyoLH3UZwGa2yYRV5nN8pHc3Nb", "fXb45B5EGqnPxCpoy8aNs71UTKXD3F17PY", "fXb45CnEbnespYeGUte6y6Su4hWnfAgLf6", "fXb45D7uPQh2BMPU824S7CDT2hcbxo7bHN", "fXb45E6WPdKG73FiG7dHT8m2fUmqDRZNjB", "fXb45FytgwDQqqi1SNE4u5udKjMTPtP5DY", "fXb45GwbQkN3J6K82CY53p4vfWRsqtgaD6", "fXb45HnaqEJH1EaiBd7bKmoyC8T42vFRbG", "fXb45JwAmQ7ucGqeNKq9WLoNCf6U3cebF6", "fXb45KUaZSqsMVcxnH68V8reesZBnXxXWe"},
  {"fXa4516QSS5g5EJNS18gJNEjGR9N5BgvaZ", "fXa452hxVgVn1xh4xyqUE7kx5ixXPS6xJw", "fXa453R3B7K5aNRa6cegJd9FZch4CXqbmL", "fXa4544xGccG7sp1YbtT315QbzaeLn4TWY", "fXa455RDzTwj5zpqr4KhcRWeQLjzKNRZqP", "fXa456buWL6FxUL1Lth14JUJoPdWGB4dB9", "fXa457rVeUv7dQJgj2AdgQZcw2WfZv3mKE", "fXa458fecsfixY5wiqNw8expCnXnLepaK4", "fXa459iE4GnRNgiy8fsMj7sTrPnxbFtXwS", "fXa45AQZdPvAvQpVeMSF1iRKJWgzgwAzA8", "fXa45B1Z1KaXvX5TAF8jpjze8YAr3J3et4", "fXa45CvqzGzjx6puiAgVYMWYVU4JWZMRzL", "fXa45DUEQ28Gpvzm5FUzra8WBiRsZT9jmb", "fXa45EMjfmULBPkZdfkSFYNsjxZopo2KR9", "fXa45FFfcnvATvrQd2DG9pzxDHdRdcT7EQ", "fXa45GuJU4iaN7GMEjtQ3DuJGCfRezDtx5", "fXa45HxKxbDv7g49wrNvzRT39F7SemTFj5", "fXa45JwHno5EXBboDmEanmwBqt6uBBkQwd", "fXa45K9Ccaap7cCSWFbe93LguxxqeNTjpv"}},
 {{"fXb4618NQKdThX8Tfs9d5oUcFoyyVBcQgJ", "fXb462w8iBQdFrix4NyhS3amANd1BgLw9x", "fXb463rdKXgdJq9hbB4jad6AAB3Dmj8mgn", "fXb4649DTdtTZGxBtuJgzRobb9qSPv6Prn", "fXb465dxcoBBjUdqLpBX76ZJkKjbWHaoHM", "fXb466o19HpGr2PEdF5h4UUGLeVuBBpkjj", "fXb467ptWpsk8iJkxc44Vj1F2gJDxWq9ie", "fXb4681g7e2HwYVgDHxojdJ2r7cvyuqMcc", "fXb469djq9AYT1xKpQXP2Vp2X1z3ZpfHdw", "fXb46AtbZbW33t92wPaYxZP3LuKVYuVhQs", "fXb46B5tXEJwRNaePV9sowUfyf9wi4S3pT", "fXb46CeXax6i8dR1g8csvjQZwAFXeiG6tQ", "fXb46Db7rtdtikMTTV9AaF52vioBoiCN2H", "fXb46EVBKAEKFvPEFBHg9jojiTFYdWBdgA", "fXb46FUaojj6cawuJz7LQkyFRyqNjxwxgk", "fXb46GabPGwdnEBgm8RKMdvwGCvSU5bgXN", "fXb46HzVGxQBgF4cPYM4GYu77WLPk7K6sx", "fXb46JBQpeU7B85Z8RgDHszP4TdyL3gqH3", "fXb46KpWyGDuieVtwHdf1CqTk3svKcnzh9"},
  {"fXa4611ikYP8np9KKCSPX2rYXLVXM4uJuz", "fXa462jz9HfQwmY3onJyydTuoMECpeC76L", "fXa463R6PRtnYfAH2CFzUr6FjPXCz8qp3E", "fXa464QA9TVYwAzRkoNYB1g5yo8zM1etYX", "fXa465mWgcrvXGuLCKWU8qrK17RqyRkiYg", "fXa466tjLZLeZBf2d7sts5qw1qe2QZRkcQ", "fXa467yfo3RjJ6YuhsF9kqv2uG7btpY4dV", "fXa468Diza4wAFNR3VP7QFRgd6SCNhkiwY", "fXa469WTYEAghjo7h66uccGB4UKErj2GdG", "fXa46Audo9u3ikFYkVi9ofCiPAnzQA9NNC", "fXa46BBTwZ6bH9xzRSHWha4NzMj4YRYd1a", "fXa46C28mcetVTdAQMW79puoaTLz61g525", "fXa46DHazxTXg2NEtmFBh8pcVT7yxnDAL4", "fXa46EPeUGX8Cg2FP2ryeSpmH2kj9b7eYg", "fXa46Fud3c19UxitJ4YxtXNYKKRH1svkDr", "fXa46GJWDnyrUpJB9KNECNUXjeZ9c84fz2", "fXa46HEeoXSyMPbWXhbqdTjvVXijuLJ522", "fXa46J8Dz1zDX5ZRmwWKuz7hgVLizSVeo6", "fXa46KqmTZHMrffE95UenoZ5FaXgCBEbtR"}},
 {{"fXb471s1YoS5C6JfXJtzVhYeQzr6fUAyU8", "fXb472PcKvpW9tu7XTYZcR5gLMfLj9D2TV", "fXb473kessceZmQJm6wfUGwysrRmCGpkRe", "fXb474dzHNpUCaTRvvmnfy5XkAQ65TVrS1", "fXb475i7Ma2nMz5zZgAirgKDy6vc6zQj9h", "fXb476kV5qYpmc3rMKUMPjv5BGu39wJujZ", "fXb477v4eyCLDZgo2tnC2JL6xRmL2ekeWb", "fXb478R5FHumVEjYcaj9Sim59sJo9Baf9K", "fXb4793XavvAMpvDXJkqjkhJAGWyXJXiB7", "fXb47A1wNsbqZMpMZPrm8x8hTDqgaZELjd", "fXb47BhmSaM7L7fT3xxyLmL7RyXktGF682", "fXb47CmXnnL6pSRnk9EUVTVpHhezpCT1B9", "fXb47DU7NRdR7KkdBx3JigX7Zvyky3jnYQ", "fXb47E6CYBxxxR5XZSUPkXLrfyWve19Sxd", "fXb47FtruZUX5WptrLhYhywR1RGyj5z14J", "fXb47GNzxtNEMUzoC1sSEpovNsu5TfbxW4", "fXb47HzoxHYzaq951wJYCeM3mCF9RMSA25", "fXb47JnzBnFShVncj2rCtn43mfXvM1a8MB", "fXb47KV2dHoMi3fojv4B6mZC1rrBP2YyiG"},
  {"fXa471bxVxfaNbscLxwcj3nezrwvttHhtb", "fXa472YX6t4j1iaQAy4fkuY3vfJNngPxzn", "fXa473YiFaWhckkauJcBSwqo8k1n2TcG6V", "fXa474HkDvpbg6q78YRrXCvk1Bsuu8CuqG", "fXa475HHuxZPDx8RdtJzro7pvFUxYgevL1", "fXa476K2jSvAGXPVRpbK9ByPkNuFTRoQks", "fXa477KRxgcg7XgGjXgt1oFtGW5s8eMTve", "fXa478SzdiyHg7Xv3cqGjRpPGyFdheFALa", "fXa4798YjWWKysCcQDCRRfoB8W7EPVGvjD", "fXa47AQQvUJu7BkCQSt9nXjDNeLbTMW89x", "fXa47B5kFbEW1Fsw2Kq7GykaJJE4DGGsAV", "fXa47CuSCeddyYstKBs6wUfRsdoqmTDXub", "fXa47DFq314LHStDSWe7uHKZZz1Wj9Ldyw", "fXa47EXcZmXNAiGnW3t5QiNkE988N4k1K3", "fXa47FSXBgHPUtMNYGeVkpDA3qgCCsqjYV", "fXa47GVF2wdSgKJvEJTPmRHEr2FTe7ht9F", "fXa47HHLsPxU9jvQLsve34HHsy8zLTJ2Vs", "fXa47JvsfHS9Xic613kPmvefW7M6LtakSk", "fXa47KT1sAjzJyHPQSDRdMB6PXH3B8g92z"}}}
    }
};

static bool posxTestScriptPubKey(const CScript& scriptPubKey)
{
    txnouttype type;
    vector<CTxDestination> addresses;
    int nRequired;

    if (!ExtractDestinations(scriptPubKey, type, addresses, nRequired))
    {
        posxCritCount++;
        return false;
    }

    bool fmatch = false;
    int n = 0;
    BOOST_FOREACH(const CTxDestination& addr, addresses)
    {
        string s = CBitcoinAddress(addr).ToString();

        if ((s[1] == 'X') &&
            ((s[2] == 'A') || (s[2] == 'a') || (s[2] == 'B') || (s[2] == 'b')))  // bid+ask, and 2 trading rounds at the same time
        {
            fmatch = true;
            int iround = (s[2] == 'a' || s[2] == 'b') ? 1 : 0;
            int iask = (s[2] == 'a' || s[2] == 'A') ? 1 : 0;

            // (s[3]<0 || s[3]>=128) 'is always false due to limited range of data type [-Wtype-limits]'
            int ipair = posx_ReverseBase58[(int)s[3]] - 1;
            int ianswer = s[4] - '1';                           // always <10
            int iprob = posx_ReverseBase58[(int)s[5]] - 1;
// printf("ProcessBlock() : xmode : testing %s iround %d ipair %d ianswer %d iask %d iprob %d\n", s.c_str(), iround, ipair, ianswer, iask, iprob);
            if (iround >= 0 && iround < POSX_ROUND_MAX &&
                iask >= 0 && iask < POSX_ASK_MAX &&
                ipair >= 0 && ipair < POSX_PAIR_MAX &&
                ianswer >= 0 && ianswer < POSX_ANSWER_MAX &&
                iprob >= 0 && iprob < POSX_PROB_MAX)
            {
                if (s.compare(strPosxAddresses[iround][ipair][ianswer][iask][iprob]) == 0)
                    strPosxExactMatch = s;
            }
        }
        else
            fmatch = false;

        n++;

        if (fmatch)
        {
            if (posxSig > 0)
                posxValueDiffIn += posxTxValue;
            else
                posxValueDiffOut += posxTxValue;
        }
    }
    // notice outputs with many addresses
    if (n > 1)
    {
        if (fmatch) posxErrCount++;
        return false;
    }
    return true;
}

static bool posxTestBlock(CBlock block)
{
    CMerkleTx txGen(block.vtx[0]);
    txGen.SetMerkleBranch(&block);

    BOOST_FOREACH(const CTransaction&tx, block.vtx)
    {
        // code from TxToJSON
        BOOST_FOREACH(const CTxIn& txin, tx.vin)
        {
            if (tx.IsCoinBase())
                // always do this (error 'No information available about transaction' if not)
                HexStr(txin.scriptSig.begin(), txin.scriptSig.end());
            else
            {
                int64 vout2 = (boost::int64_t)txin.prevout.n;

                // backtrace transactions if not mined (need to know where the coins coming from)
                CTransaction tx2;
                uint256 hashBlock2 = 0;
                // error 'No information available about transaction'
                if (!GetTransaction(txin.prevout.hash, tx2, hashBlock2, true))
                    return false;

                // check outputs of each 'previous' transaction (to count coins sent *from* special addresses)
                posxTxValue = 0;
                posxSig = -1;
                for (unsigned int i = 0; i < tx2.vout.size(); i++)
                {
                    const CTxOut& txout4 = tx2.vout[i];
                    posxTxValue = txout4.nValue;

                    if (i != vout2) continue;

                    posxTestScriptPubKey(txout4.scriptPubKey);
                }
                // end processing 'previous' tx

            }
        }

        // check outputs of each transaction (to count coins sent *to* special addresses)
        posxTxValue = 0;
        posxSig = 1;
        for (unsigned int i = 0; i < tx.vout.size(); i++)
        {
            const CTxOut& txout = tx.vout[i];
            posxTxValue = txout.nValue;
            posxTestScriptPubKey(txout.scriptPubKey);
        }
    }
    return true;
}

bool CBlock::CheckBlock(CValidationState &state, bool fCheckPOW, bool fCheckMerkleRoot) const
{
    // These are checks that are independent of context
    // that can be verified before saving an orphan block.

    // Size limits
    if (vtx.empty() || vtx.size() > MAX_BLOCK_SIZE || ::GetSerializeSize(*this, SER_NETWORK, PROTOCOL_VERSION) > MAX_BLOCK_SIZE)
        return state.DoS(100, error("CheckBlock() : size limits failed"));

    // Litecoin: Special short-term limits to avoid 10,000 BDB lock limit:
    if (GetBlockTime() < 1376568000)  // stop enforcing 15 August 2013 00:00:00
    {
        // Rule is: #unique txids referenced <= 4,500
        // ... to prevent 10,000 BDB lock exhaustion on old clients
        set<uint256> setTxIn;
        for (size_t i = 0; i < vtx.size(); i++)
        {
            setTxIn.insert(vtx[i].GetHash());
            if (i == 0) continue; // skip coinbase txin
            BOOST_FOREACH(const CTxIn& txin, vtx[i].vin)
                setTxIn.insert(txin.prevout.hash);
        }
        size_t nTxids = setTxIn.size();
        if (nTxids > 4500)
            return error("CheckBlock() : 15 August maxlocks violation");
    }

    // Check proof of work matches claimed amount
    if (fCheckPOW && !CheckProofOfWork(GetPoWHash(), nBits))
        return state.DoS(50, error("CheckBlock() : proof of work failed"));

    // Check timestamp
    if (GetBlockTime() > GetAdjustedTime() + 2 * 60 * 60)
        return state.Invalid(error("CheckBlock() : block timestamp too far in the future"));

    // First transaction must be coinbase, the rest must not be
    if (vtx.empty() || !vtx[0].IsCoinBase())
        return state.DoS(100, error("CheckBlock() : first tx is not coinbase"));
    for (unsigned int i = 1; i < vtx.size(); i++)
        if (vtx[i].IsCoinBase())
            return state.DoS(100, error("CheckBlock() : more than one coinbase"));

    // Check transactions
    BOOST_FOREACH(const CTransaction& tx, vtx)
        if (!tx.CheckTransaction(state))
            return error("CheckBlock() : CheckTransaction failed");

    // Build the merkle tree already. We need it anyway later, and it makes the
    // block cache the transaction hashes, which means they don't need to be
    // recalculated many times during this block's validation.
    BuildMerkleTree();

    // Check for duplicate txids. This is caught by ConnectInputs(),
    // but catching it earlier avoids a potential DoS attack:
    set<uint256> uniqueTx;
    for (unsigned int i=0; i<vtx.size(); i++) {
        uniqueTx.insert(GetTxHash(i));
    }
    if (uniqueTx.size() != vtx.size())
        return state.DoS(100, error("CheckBlock() : duplicate transaction"));

    unsigned int nSigOps = 0;
    BOOST_FOREACH(const CTransaction& tx, vtx)
    {
        nSigOps += tx.GetLegacySigOpCount();
    }
    if (nSigOps > MAX_BLOCK_SIGOPS)
        return state.DoS(100, error("CheckBlock() : out-of-bounds SigOpCount"));

    // Check merkle root
    if (fCheckMerkleRoot && hashMerkleRoot != BuildMerkleTree())
        return state.DoS(100, error("CheckBlock() : hashMerkleRoot mismatch"));

    return true;
}

bool CBlock::AcceptBlock(CValidationState &state, CDiskBlockPos *dbp)
{
    // Check for duplicate
    uint256 hash = GetHash();
    if (mapBlockIndex.count(hash))
        return state.Invalid(error("AcceptBlock() : block already in mapBlockIndex"));

    // Get prev block index
    CBlockIndex* pindexPrev = NULL;
    int nHeight = 0;
    if (hash != hashGenesisBlock) {
        map<uint256, CBlockIndex*>::iterator mi = mapBlockIndex.find(hashPrevBlock);
        if (mi == mapBlockIndex.end())
            return state.DoS(10, error("AcceptBlock() : prev block not found"));
        pindexPrev = (*mi).second;
        nHeight = pindexPrev->nHeight+1;

        // Check proof of work
        if (nBits != GetNextWorkRequired(pindexPrev, this))
            return state.DoS(100, error("AcceptBlock() : incorrect proof of work"));

        // Check timestamp against prev
        if (GetBlockTime() <= pindexPrev->GetMedianTimePast())
            return state.Invalid(error("AcceptBlock() : block's timestamp is too early"));

        // Check that all transactions are finalized
        BOOST_FOREACH(const CTransaction& tx, vtx)
            if (!tx.IsFinal(nHeight, GetBlockTime()))
                return state.DoS(10, error("AcceptBlock() : contains a non-final transaction"));

        // Check that the block chain matches the known block chain up to a checkpoint
        if (!Checkpoints::CheckBlock(nHeight, hash))
            return state.DoS(100, error("AcceptBlock() : rejected by checkpoint lock-in at %d", nHeight));

        // Reject block.nVersion=1 blocks when 95% (75% on testnet) of the network has upgraded:
        if (nVersion < 2)
        {
            if ((!fTestNet && CBlockIndex::IsSuperMajority(2, pindexPrev, 950, 1000)) ||
                (fTestNet && CBlockIndex::IsSuperMajority(2, pindexPrev, 75, 100)))
            {
                return state.Invalid(error("AcceptBlock() : rejected nVersion=1 block"));
            }
        }
        // Enforce block.nVersion=2 rule that the coinbase starts with serialized block height
        if (nVersion >= 2)
        {
            // if 750 of the last 1,000 blocks are version 2 or greater (51/100 if testnet):
            if ((!fTestNet && CBlockIndex::IsSuperMajority(2, pindexPrev, 750, 1000)) ||
                (fTestNet && CBlockIndex::IsSuperMajority(2, pindexPrev, 51, 100)))
            {
                CScript expect = CScript() << nHeight;
                if (vtx[0].vin[0].scriptSig.size() < expect.size() ||
                    !std::equal(expect.begin(), expect.end(), vtx[0].vin[0].scriptSig.begin()))
                    return state.DoS(100, error("AcceptBlock() : block height mismatch in coinbase"));
            }
        }
    }

    // Write block to history file
    try {
        unsigned int nBlockSize = ::GetSerializeSize(*this, SER_DISK, CLIENT_VERSION);
        CDiskBlockPos blockPos;
        if (dbp != NULL)
            blockPos = *dbp;
        if (!FindBlockPos(state, blockPos, nBlockSize+8, nHeight, nTime, dbp != NULL))
            return error("AcceptBlock() : FindBlockPos failed");
        if (dbp == NULL)
            if (!WriteToDisk(blockPos))
                return state.Abort(_("Failed to write block"));
        if (!AddToBlockIndex(state, blockPos))
            return error("AcceptBlock() : AddToBlockIndex failed");
    } catch(std::runtime_error &e) {
        return state.Abort(_("System error: ") + e.what());
    }

    // Relay inventory, but don't relay old inventory during initial block download
    int nBlockEstimate = Checkpoints::GetTotalBlocksEstimate();
    if (hashBestChain == hash)
    {
        LOCK(cs_vNodes);
        BOOST_FOREACH(CNode* pnode, vNodes)
            if (nBestHeight > (pnode->nStartingHeight != -1 ? pnode->nStartingHeight - 2000 : nBlockEstimate))
                pnode->PushInventory(CInv(MSG_BLOCK, hash));
    }

    return true;
}

bool CBlockIndex::IsSuperMajority(int minVersion, const CBlockIndex* pstart, unsigned int nRequired, unsigned int nToCheck)
{
    // Litecoin: temporarily disable v2 block lockin until we are ready for v2 transition
    return false;
    unsigned int nFound = 0;
    for (unsigned int i = 0; i < nToCheck && nFound < nRequired && pstart != NULL; i++)
    {
        if (pstart->nVersion >= minVersion)
            ++nFound;
        pstart = pstart->pprev;
    }
    return (nFound >= nRequired);
}

bool ProcessBlock(CValidationState &state, CNode* pfrom, CBlock* pblock, CDiskBlockPos *dbp)
{
    // Check for duplicate
    uint256 hash = pblock->GetHash();
    if (mapBlockIndex.count(hash))
        return state.Invalid(error("ProcessBlock() : already have block %d %s", mapBlockIndex[hash]->nHeight, hash.ToString().c_str()));
    if (mapOrphanBlocks.count(hash))
        return state.Invalid(error("ProcessBlock() : already have block (orphan) %s", hash.ToString().c_str()));

    // Preliminary checks
    if (!pblock->CheckBlock(state))
        return error("ProcessBlock() : CheckBlock FAILED");

    CBlockIndex* pcheckpoint = Checkpoints::GetLastCheckpoint(mapBlockIndex);
    if (pcheckpoint && pblock->hashPrevBlock != hashBestChain)
    {
        // Extra checks to prevent "fill up memory by spamming with bogus blocks"
        int64 deltaTime = pblock->GetBlockTime() - pcheckpoint->nTime;
        if (deltaTime < 0)
        {
            return state.DoS(100, error("ProcessBlock() : block with timestamp before last checkpoint"));
        }
        CBigNum bnNewBlock;
        bnNewBlock.SetCompact(pblock->nBits);
        CBigNum bnRequired;
        bnRequired.SetCompact(ComputeMinWork(pcheckpoint->nBits, deltaTime));
        if (bnNewBlock > bnRequired)
        {
            return state.DoS(100, error("ProcessBlock() : block with too little proof-of-work"));
        }
    }


    // If we don't already have its previous block, shunt it off to holding area until we get it
    if (pblock->hashPrevBlock != 0 && !mapBlockIndex.count(pblock->hashPrevBlock))
    {
        printf("ProcessBlock: ORPHAN BLOCK, prev=%s\n", pblock->hashPrevBlock.ToString().c_str());

        // Accept orphans as long as there is a node to request its parents from
        if (pfrom) {
            CBlock* pblock2 = new CBlock(*pblock);
            mapOrphanBlocks.insert(make_pair(hash, pblock2));
            mapOrphanBlocksByPrev.insert(make_pair(pblock2->hashPrevBlock, pblock2));

            // Ask this guy to fill in what we're missing
            pfrom->PushGetBlocks(pindexBest, GetOrphanRoot(pblock2));
        }
        return true;
    }

// FBX proof of stake voting test
    if (nBestHeight>180000 && GetBoolArg("-xmode", false) && GetBoolArg("-txindex", false))
    {
        strPosxExactMatch = "";
        posxValueDiffIn = posxValueDiffOut = 0;
        posxErrCount = posxCritCount = posxSanCount = 0;
        CBlock block2 = * pblock; // use copy of the block as to not break something
        if (!posxTestBlock(block2) || posxCritCount || posxErrCount || posxSanCount)
            printf("ProcessBlock() : xmode : error, block is invalid, hash %s\n", block2.GetHash().ToString().c_str());
        else if (posxValueDiffIn || posxValueDiffOut)
            printf("ProcessBlock() : xmode : access of potential order book addr, coins in %d, out %d, hash %s\n",
               (int)(posxValueDiffIn/COIN), (int)(posxValueDiffOut/COIN), block2.GetHash().ToString().c_str());
        else
            printf("ProcessBlock() : xmode : test passed, hash %s\n", block2.GetHash().ToString().c_str());

        if (strPosxExactMatch.length() > 0)
            printf("ProcessBlock() : xmode : exact match, addr %s\n", strPosxExactMatch.c_str());
    }

    // Store to disk
    if (!pblock->AcceptBlock(state, dbp))
        return error("ProcessBlock() : AcceptBlock FAILED");

    // Recursively process any orphan blocks that depended on this one
    vector<uint256> vWorkQueue;
    vWorkQueue.push_back(hash);
    for (unsigned int i = 0; i < vWorkQueue.size(); i++)
    {
        uint256 hashPrev = vWorkQueue[i];
        for (multimap<uint256, CBlock*>::iterator mi = mapOrphanBlocksByPrev.lower_bound(hashPrev);
             mi != mapOrphanBlocksByPrev.upper_bound(hashPrev);
             ++mi)
        {
            CBlock* pblockOrphan = (*mi).second;
            // Use a dummy CValidationState so someone can't setup nodes to counter-DoS based on orphan resolution (that is, feeding people an invalid block based on LegitBlockX in order to get anyone relaying LegitBlockX banned)
            CValidationState stateDummy;
            if (pblockOrphan->AcceptBlock(stateDummy))
                vWorkQueue.push_back(pblockOrphan->GetHash());
            mapOrphanBlocks.erase(pblockOrphan->GetHash());
            delete pblockOrphan;
        }
        mapOrphanBlocksByPrev.erase(hashPrev);
    }

    printf("ProcessBlock: ACCEPTED\n");
    return true;
}








CMerkleBlock::CMerkleBlock(const CBlock& block, CBloomFilter& filter)
{
    header = block.GetBlockHeader();

    vector<bool> vMatch;
    vector<uint256> vHashes;

    vMatch.reserve(block.vtx.size());
    vHashes.reserve(block.vtx.size());

    for (unsigned int i = 0; i < block.vtx.size(); i++)
    {
        uint256 hash = block.vtx[i].GetHash();
        if (filter.IsRelevantAndUpdate(block.vtx[i], hash))
        {
            vMatch.push_back(true);
            vMatchedTxn.push_back(make_pair(i, hash));
        }
        else
            vMatch.push_back(false);
        vHashes.push_back(hash);
    }

    txn = CPartialMerkleTree(vHashes, vMatch);
}








uint256 CPartialMerkleTree::CalcHash(int height, unsigned int pos, const std::vector<uint256> &vTxid) {
    if (height == 0) {
        // hash at height 0 is the txids themself
        return vTxid[pos];
    } else {
        // calculate left hash
        uint256 left = CalcHash(height-1, pos*2, vTxid), right;
        // calculate right hash if not beyong the end of the array - copy left hash otherwise1
        if (pos*2+1 < CalcTreeWidth(height-1))
            right = CalcHash(height-1, pos*2+1, vTxid);
        else
            right = left;
        // combine subhashes
        return Hash(BEGIN(left), END(left), BEGIN(right), END(right));
    }
}

void CPartialMerkleTree::TraverseAndBuild(int height, unsigned int pos, const std::vector<uint256> &vTxid, const std::vector<bool> &vMatch) {
    // determine whether this node is the parent of at least one matched txid
    bool fParentOfMatch = false;
    for (unsigned int p = pos << height; p < (pos+1) << height && p < nTransactions; p++)
        fParentOfMatch |= vMatch[p];
    // store as flag bit
    vBits.push_back(fParentOfMatch);
    if (height==0 || !fParentOfMatch) {
        // if at height 0, or nothing interesting below, store hash and stop
        vHash.push_back(CalcHash(height, pos, vTxid));
    } else {
        // otherwise, don't store any hash, but descend into the subtrees
        TraverseAndBuild(height-1, pos*2, vTxid, vMatch);
        if (pos*2+1 < CalcTreeWidth(height-1))
            TraverseAndBuild(height-1, pos*2+1, vTxid, vMatch);
    }
}

uint256 CPartialMerkleTree::TraverseAndExtract(int height, unsigned int pos, unsigned int &nBitsUsed, unsigned int &nHashUsed, std::vector<uint256> &vMatch) {
    if (nBitsUsed >= vBits.size()) {
        // overflowed the bits array - failure
        fBad = true;
        return 0;
    }
    bool fParentOfMatch = vBits[nBitsUsed++];
    if (height==0 || !fParentOfMatch) {
        // if at height 0, or nothing interesting below, use stored hash and do not descend
        if (nHashUsed >= vHash.size()) {
            // overflowed the hash array - failure
            fBad = true;
            return 0;
        }
        const uint256 &hash = vHash[nHashUsed++];
        if (height==0 && fParentOfMatch) // in case of height 0, we have a matched txid
            vMatch.push_back(hash);
        return hash;
    } else {
        // otherwise, descend into the subtrees to extract matched txids and hashes
        uint256 left = TraverseAndExtract(height-1, pos*2, nBitsUsed, nHashUsed, vMatch), right;
        if (pos*2+1 < CalcTreeWidth(height-1))
            right = TraverseAndExtract(height-1, pos*2+1, nBitsUsed, nHashUsed, vMatch);
        else
            right = left;
        // and combine them before returning
        return Hash(BEGIN(left), END(left), BEGIN(right), END(right));
    }
}

CPartialMerkleTree::CPartialMerkleTree(const std::vector<uint256> &vTxid, const std::vector<bool> &vMatch) : nTransactions(vTxid.size()), fBad(false) {
    // reset state
    vBits.clear();
    vHash.clear();

    // calculate height of tree
    int nHeight = 0;
    while (CalcTreeWidth(nHeight) > 1)
        nHeight++;

    // traverse the partial tree
    TraverseAndBuild(nHeight, 0, vTxid, vMatch);
}

CPartialMerkleTree::CPartialMerkleTree() : nTransactions(0), fBad(true) {}

uint256 CPartialMerkleTree::ExtractMatches(std::vector<uint256> &vMatch) {
    vMatch.clear();
    // An empty set will not work
    if (nTransactions == 0)
        return 0;
    // check for excessively high numbers of transactions
    if (nTransactions > MAX_BLOCK_SIZE / 60) // 60 is the lower bound for the size of a serialized CTransaction
        return 0;
    // there can never be more hashes provided than one for every txid
    if (vHash.size() > nTransactions)
        return 0;
    // there must be at least one bit per node in the partial tree, and at least one node per hash
    if (vBits.size() < vHash.size())
        return 0;
    // calculate height of tree
    int nHeight = 0;
    while (CalcTreeWidth(nHeight) > 1)
        nHeight++;
    // traverse the partial tree
    unsigned int nBitsUsed = 0, nHashUsed = 0;
    uint256 hashMerkleRoot = TraverseAndExtract(nHeight, 0, nBitsUsed, nHashUsed, vMatch);
    // verify that no problems occured during the tree traversal
    if (fBad)
        return 0;
    // verify that all bits were consumed (except for the padding caused by serializing it as a byte sequence)
    if ((nBitsUsed+7)/8 != (vBits.size()+7)/8)
        return 0;
    // verify that all hashes were consumed
    if (nHashUsed != vHash.size())
        return 0;
    return hashMerkleRoot;
}







bool AbortNode(const std::string &strMessage) {
    strMiscWarning = strMessage;
    printf("*** %s\n", strMessage.c_str());
    uiInterface.ThreadSafeMessageBox(strMessage, "", CClientUIInterface::MSG_ERROR);
    StartShutdown();
    return false;
}

bool CheckDiskSpace(uint64 nAdditionalBytes)
{
    uint64 nFreeBytesAvailable = filesystem::space(GetDataDir()).available;

    // Check for nMinDiskSpace bytes (currently 50MB)
    if (nFreeBytesAvailable < nMinDiskSpace + nAdditionalBytes)
        return AbortNode(_("Error: Disk space is low!"));

    return true;
}

CCriticalSection cs_LastBlockFile;
CBlockFileInfo infoLastBlockFile;
int nLastBlockFile = 0;

FILE* OpenDiskFile(const CDiskBlockPos &pos, const char *prefix, bool fReadOnly)
{
    if (pos.IsNull())
        return NULL;
    boost::filesystem::path path = GetDataDir() / "blocks" / strprintf("%s%05u.dat", prefix, pos.nFile);
    boost::filesystem::create_directories(path.parent_path());
    FILE* file = fopen(path.string().c_str(), "rb+");
    if (!file && !fReadOnly)
        file = fopen(path.string().c_str(), "wb+");
    if (!file) {
        printf("Unable to open file %s\n", path.string().c_str());
        return NULL;
    }
    if (pos.nPos) {
        if (fseek(file, pos.nPos, SEEK_SET)) {
            printf("Unable to seek to position %u of %s\n", pos.nPos, path.string().c_str());
            fclose(file);
            return NULL;
        }
    }
    return file;
}

FILE* OpenBlockFile(const CDiskBlockPos &pos, bool fReadOnly) {
    return OpenDiskFile(pos, "blk", fReadOnly);
}

FILE* OpenUndoFile(const CDiskBlockPos &pos, bool fReadOnly) {
    return OpenDiskFile(pos, "rev", fReadOnly);
}

CBlockIndex * InsertBlockIndex(uint256 hash)
{
    if (hash == 0)
        return NULL;

    // Return existing
    map<uint256, CBlockIndex*>::iterator mi = mapBlockIndex.find(hash);
    if (mi != mapBlockIndex.end())
        return (*mi).second;

    // Create new
    CBlockIndex* pindexNew = new CBlockIndex();
    if (!pindexNew)
        throw runtime_error("LoadBlockIndex() : new CBlockIndex failed");
    mi = mapBlockIndex.insert(make_pair(hash, pindexNew)).first;
    pindexNew->phashBlock = &((*mi).first);

    return pindexNew;
}

bool static LoadBlockIndexDB()
{
    if (!pblocktree->LoadBlockIndexGuts())
        return false;

    boost::this_thread::interruption_point();

    // Calculate nChainWork
    vector<pair<int, CBlockIndex*> > vSortedByHeight;
    vSortedByHeight.reserve(mapBlockIndex.size());
    BOOST_FOREACH(const PAIRTYPE(uint256, CBlockIndex*)& item, mapBlockIndex)
    {
        CBlockIndex* pindex = item.second;
        vSortedByHeight.push_back(make_pair(pindex->nHeight, pindex));
    }
    sort(vSortedByHeight.begin(), vSortedByHeight.end());
    BOOST_FOREACH(const PAIRTYPE(int, CBlockIndex*)& item, vSortedByHeight)
    {
        CBlockIndex* pindex = item.second;
        pindex->nChainWork = (pindex->pprev ? pindex->pprev->nChainWork : 0) + pindex->GetBlockWork().getuint256();
        pindex->nChainTx = (pindex->pprev ? pindex->pprev->nChainTx : 0) + pindex->nTx;
        if ((pindex->nStatus & BLOCK_VALID_MASK) >= BLOCK_VALID_TRANSACTIONS && !(pindex->nStatus & BLOCK_FAILED_MASK))
            setBlockIndexValid.insert(pindex);
    }

    // Load block file info
    pblocktree->ReadLastBlockFile(nLastBlockFile);
    printf("LoadBlockIndexDB(): last block file = %i\n", nLastBlockFile);
    if (pblocktree->ReadBlockFileInfo(nLastBlockFile, infoLastBlockFile))
        printf("LoadBlockIndexDB(): last block file info: %s\n", infoLastBlockFile.ToString().c_str());

    // Load nBestInvalidWork, OK if it doesn't exist
    CBigNum bnBestInvalidWork;
    pblocktree->ReadBestInvalidWork(bnBestInvalidWork);
    nBestInvalidWork = bnBestInvalidWork.getuint256();

    // Check whether we need to continue reindexing
    bool fReindexing = false;
    pblocktree->ReadReindexing(fReindexing);
    fReindex |= fReindexing;

    // Check whether we have a transaction index
    pblocktree->ReadFlag("txindex", fTxIndex);
    printf("LoadBlockIndexDB(): transaction index %s\n", fTxIndex ? "enabled" : "disabled");

    // Load hashBestChain pointer to end of best chain
    pindexBest = pcoinsTip->GetBestBlock();
    if (pindexBest == NULL)
        return true;
    hashBestChain = pindexBest->GetBlockHash();
    nBestHeight = pindexBest->nHeight;
    nBestChainWork = pindexBest->nChainWork;

    // set 'next' pointers in best chain
    CBlockIndex *pindex = pindexBest;
    while(pindex != NULL && pindex->pprev != NULL) {
         CBlockIndex *pindexPrev = pindex->pprev;
         pindexPrev->pnext = pindex;
         pindex = pindexPrev;
    }
    printf("LoadBlockIndexDB(): hashBestChain=%s  height=%d date=%s\n",
        hashBestChain.ToString().c_str(), nBestHeight,
        DateTimeStrFormat("%Y-%m-%d %H:%M:%S", pindexBest->GetBlockTime()).c_str());

    return true;
}

bool VerifyDB() {
    if (pindexBest == NULL || pindexBest->pprev == NULL)
        return true;

    // Verify blocks in the best chain
    int nCheckLevel = GetArg("-checklevel", 3);
    int nCheckDepth = GetArg( "-checkblocks", 288);
    if (nCheckDepth == 0)
        nCheckDepth = 1000000000; // suffices until the year 19000
    if (nCheckDepth > nBestHeight)
        nCheckDepth = nBestHeight;
    nCheckLevel = std::max(0, std::min(4, nCheckLevel));
    printf("Verifying last %i blocks at level %i\n", nCheckDepth, nCheckLevel);
    CCoinsViewCache coins(*pcoinsTip, true);
    CBlockIndex* pindexState = pindexBest;
    CBlockIndex* pindexFailure = NULL;
    int nGoodTransactions = 0;
    CValidationState state;
    for (CBlockIndex* pindex = pindexBest; pindex && pindex->pprev; pindex = pindex->pprev)
    {
        boost::this_thread::interruption_point();
        if (pindex->nHeight < nBestHeight-nCheckDepth)
            break;
        CBlock block;
        // check level 0: read from disk
        if (!block.ReadFromDisk(pindex))
            return error("VerifyDB() : *** block.ReadFromDisk failed at %d, hash=%s", pindex->nHeight, pindex->GetBlockHash().ToString().c_str());
        // check level 1: verify block validity
        if (nCheckLevel >= 1 && !block.CheckBlock(state))
            return error("VerifyDB() : *** found bad block at %d, hash=%s\n", pindex->nHeight, pindex->GetBlockHash().ToString().c_str());
        // check level 2: verify undo validity
        if (nCheckLevel >= 2 && pindex) {
            CBlockUndo undo;
            CDiskBlockPos pos = pindex->GetUndoPos();
            if (!pos.IsNull()) {
                if (!undo.ReadFromDisk(pos, pindex->pprev->GetBlockHash()))
                    return error("VerifyDB() : *** found bad undo data at %d, hash=%s\n", pindex->nHeight, pindex->GetBlockHash().ToString().c_str());
            }
        }
        // check level 3: check for inconsistencies during memory-only disconnect of tip blocks
        if (nCheckLevel >= 3 && pindex == pindexState && (coins.GetCacheSize() + pcoinsTip->GetCacheSize()) <= 2*nCoinCacheSize + 32000) {
            bool fClean = true;
            if (!block.DisconnectBlock(state, pindex, coins, &fClean))
                return error("VerifyDB() : *** irrecoverable inconsistency in block data at %d, hash=%s", pindex->nHeight, pindex->GetBlockHash().ToString().c_str());
            pindexState = pindex->pprev;
            if (!fClean) {
                nGoodTransactions = 0;
                pindexFailure = pindex;
            } else
                nGoodTransactions += block.vtx.size();
        }
    }
    if (pindexFailure)
        return error("VerifyDB() : *** coin database inconsistencies found (last %i blocks, %i good transactions before that)\n", pindexBest->nHeight - pindexFailure->nHeight + 1, nGoodTransactions);

    // check level 4: try reconnecting blocks
    if (nCheckLevel >= 4) {
        CBlockIndex *pindex = pindexState;
        while (pindex != pindexBest) {
            boost::this_thread::interruption_point();
            pindex = pindex->pnext;
            CBlock block;
            if (!block.ReadFromDisk(pindex))
                return error("VerifyDB() : *** block.ReadFromDisk failed at %d, hash=%s", pindex->nHeight, pindex->GetBlockHash().ToString().c_str());
            if (!block.ConnectBlock(state, pindex, coins))
                return error("VerifyDB() : *** found unconnectable block at %d, hash=%s", pindex->nHeight, pindex->GetBlockHash().ToString().c_str());
        }
    }

    printf("No coin database inconsistencies in last %i blocks (%i transactions)\n", pindexBest->nHeight - pindexState->nHeight, nGoodTransactions);

    return true;
}

void UnloadBlockIndex()
{
    mapBlockIndex.clear();
    setBlockIndexValid.clear();
    pindexGenesisBlock = NULL;
    nBestHeight = 0;
    nBestChainWork = 0;
    nBestInvalidWork = 0;
    hashBestChain = 0;
    pindexBest = NULL;
}

bool LoadBlockIndex()
{
// FBX
//    if (fTestNet)
//    {
//        pchMessageStart[0] = 0xfc;
//        pchMessageStart[1] = 0xc1;
//        pchMessageStart[2] = 0xb7;
//        pchMessageStart[3] = 0xdc;
//        hashGenesisBlock = uint256("0xf5ae71e26c74beacc88382716aced69cddf3dffff24f384e1808905e0188f68f");
//    }
    pchMessageStart[0] = 0xf9;
    pchMessageStart[1] = 0xdb;
    pchMessageStart[2] = 0xf9;
    pchMessageStart[3] = 0xdb;
    hashGenesisBlock = uint256("0x002a91713910bc96eb0edf237fcd2799d7a01186e1e96023e860bc70b3916200");

    //
    // Load block index from databases
    //
    if (!fReindex && !LoadBlockIndexDB())
        return false;

    return true;
}


bool InitBlockIndex() {
    // Check whether we're already initialized
    if (pindexGenesisBlock != NULL)
        return true;

    // Use the provided setting for -txindex in the new database
    fTxIndex = GetBoolArg("-txindex", false);
    pblocktree->WriteFlag("txindex", fTxIndex);
    printf("Initializing databases...\n");

    // Only add the genesis block if not reindexing (in which case we reuse the one already on disk)
    if (!fReindex) {
        // Genesis Block:
        // CBlock(hash=12a765e31ffd4059bada, PoW=0000050c34a64b415b6b, ver=1, hashPrevBlock=00000000000000000000, hashMerkleRoot=97ddfbbae6, nTime=1317972665, nBits=1e0ffff0, nNonce=2084524493, vtx=1)
        //   CTransaction(hash=97ddfbbae6, ver=1, vin.size=1, vout.size=1, nLockTime=0)
        //     CTxIn(COutPoint(0000000000, -1), coinbase 04ffff001d0104404e592054696d65732030352f4f63742f32303131205374657665204a6f62732c204170706c65e280997320566973696f6e6172792c2044696573206174203536)
        //     CTxOut(nValue=50.00000000, scriptPubKey=040184710fa689ad5023690c80f3a4)
        //   vMerkleTree: 97ddfbbae6

// FBX
          // Genesis block
//        const char* pszTimestamp = "NY Times 05/Oct/2011 Steve Jobs, Apple’s Visionary, Dies at 56";
//        CTransaction txNew;
//        txNew.vin.resize(1);
//        txNew.vout.resize(1);
//        txNew.vin[0].scriptSig = CScript() << 486604799 << CBigNum(4) << vector<unsigned char>((const unsigned char*)pszTimestamp, (const unsigned char*)pszTimestamp + strlen(pszTimestamp));
//        txNew.vout[0].nValue = 50 * COIN;
//        txNew.vout[0].scriptPubKey = CScript() << ParseHex("040184710fa689ad5023690c80f3a49c8f13f8d45b8c857fbcbc8bc4a8e4d3eb4b10f4d4604fa08dce601aaf0f470216fe1b51850b4acf21b179c45070ac7b03a9") << OP_CHECKSIG;
//        CBlock block;
//        block.vtx.push_back(txNew);
//        block.hashPrevBlock = 0;
//        block.hashMerkleRoot = block.BuildMerkleTree();
//        block.nVersion = 1;
//        block.nTime    = 1317972665;
//        block.nBits    = 0x1e0ffff0;
//        block.nNonce   = 2084524493;
//
//        if (fTestNet)
//        {
//            block.nTime    = 1317798646;
//            block.nNonce   = 385270584;
//        }
        // Genesis block
        const char* pszTimestamp = "\"nytimes.com 10/1/2011 - Police Arrest Over 700 Protesters on Brooklyn Bridge\"";
        CTransaction txNew;
        txNew.vin.resize(1);
        txNew.vout.resize(1);
        int nBits = 486604799;
        int extra = 3;
        txNew.vin[0].scriptSig = CScript() << nBits << CBigNum(++extra) << vector<unsigned char>((const unsigned char*)pszTimestamp, (const unsigned char*)pszTimestamp + strlen(pszTimestamp));
        txNew.vout[0].nValue = 50 * COIN;
        txNew.vout[0].scriptPubKey = CScript() << ParseHex("04678afdb0fe5548271967f1a67130b7105cd6a828e03909a67962e0ea1f61deb649f6bc3f4cef38c4f35504e51ec112de5c384df7ba0b8d578a4c702b6bf11d5f") << OP_CHECKSIG;
        CBlock block;
        block.vtx.push_back(txNew);
        block.hashPrevBlock = 0;
        block.hashMerkleRoot = block.BuildMerkleTree();
        block.nVersion = 1;
        block.nTime    = 1317529878;
        block.nBits    = 0x1e0ffff0;
        block.nNonce   = 385610221;

        //// debug print
        uint256 hash = block.GetHash();
        printf("%s\n", hash.ToString().c_str());
        printf("%s\n", hashGenesisBlock.ToString().c_str());
        printf("%s\n", block.hashMerkleRoot.ToString().c_str());
// FBX
//        assert(block.hashMerkleRoot == uint256("0x97ddfbbae6be97fd6cdf3e7ca13232a3afff2353e29badfab7f73011edd4ced9"));
        assert(block.hashMerkleRoot == uint256("0x40a627262ed716f0f3d5104315fe0b600bf8e32a021929299163f74151fa52b1"));

        block.print();
        assert(hash == hashGenesisBlock);

        // Start new block file
        try {
            unsigned int nBlockSize = ::GetSerializeSize(block, SER_DISK, CLIENT_VERSION);
            CDiskBlockPos blockPos;
            CValidationState state;
            if (!FindBlockPos(state, blockPos, nBlockSize+8, 0, block.nTime))
                return error("LoadBlockIndex() : FindBlockPos failed");
            if (!block.WriteToDisk(blockPos))
                return error("LoadBlockIndex() : writing genesis block to disk failed");
            if (!block.AddToBlockIndex(state, blockPos))
                return error("LoadBlockIndex() : genesis block not accepted");
        } catch(std::runtime_error &e) {
            return error("LoadBlockIndex() : failed to initialize block database: %s", e.what());
        }
    }

    return true;
}



void PrintBlockTree()
{
    // pre-compute tree structure
    map<CBlockIndex*, vector<CBlockIndex*> > mapNext;
    for (map<uint256, CBlockIndex*>::iterator mi = mapBlockIndex.begin(); mi != mapBlockIndex.end(); ++mi)
    {
        CBlockIndex* pindex = (*mi).second;
        mapNext[pindex->pprev].push_back(pindex);
        // test
        //while (rand() % 3 == 0)
        //    mapNext[pindex->pprev].push_back(pindex);
    }

    vector<pair<int, CBlockIndex*> > vStack;
    vStack.push_back(make_pair(0, pindexGenesisBlock));

    int nPrevCol = 0;
    while (!vStack.empty())
    {
        int nCol = vStack.back().first;
        CBlockIndex* pindex = vStack.back().second;
        vStack.pop_back();

        // print split or gap
        if (nCol > nPrevCol)
        {
            for (int i = 0; i < nCol-1; i++)
                printf("| ");
            printf("|\\\n");
        }
        else if (nCol < nPrevCol)
        {
            for (int i = 0; i < nCol; i++)
                printf("| ");
            printf("|\n");
       }
        nPrevCol = nCol;

        // print columns
        for (int i = 0; i < nCol; i++)
            printf("| ");

        // print item
        CBlock block;
        block.ReadFromDisk(pindex);
        printf("%d (blk%05u.dat:0x%x)  %s  tx %"PRIszu"",
            pindex->nHeight,
            pindex->GetBlockPos().nFile, pindex->GetBlockPos().nPos,
            DateTimeStrFormat("%Y-%m-%d %H:%M:%S", block.GetBlockTime()).c_str(),
            block.vtx.size());

        PrintWallets(block);

        // put the main time-chain first
        vector<CBlockIndex*>& vNext = mapNext[pindex];
        for (unsigned int i = 0; i < vNext.size(); i++)
        {
            if (vNext[i]->pnext)
            {
                swap(vNext[0], vNext[i]);
                break;
            }
        }

        // iterate children
        for (unsigned int i = 0; i < vNext.size(); i++)
            vStack.push_back(make_pair(nCol+i, vNext[i]));
    }
}

bool LoadExternalBlockFile(FILE* fileIn, CDiskBlockPos *dbp)
{
    int64 nStart = GetTimeMillis();

    int nLoaded = 0;
    try {
        CBufferedFile blkdat(fileIn, 2*MAX_BLOCK_SIZE, MAX_BLOCK_SIZE+8, SER_DISK, CLIENT_VERSION);
        uint64 nStartByte = 0;
        if (dbp) {
            // (try to) skip already indexed part
            CBlockFileInfo info;
            if (pblocktree->ReadBlockFileInfo(dbp->nFile, info)) {
                nStartByte = info.nSize;
                blkdat.Seek(info.nSize);
            }
        }
        uint64 nRewind = blkdat.GetPos();
        while (blkdat.good() && !blkdat.eof()) {
            boost::this_thread::interruption_point();

            blkdat.SetPos(nRewind);
            nRewind++; // start one byte further next time, in case of failure
            blkdat.SetLimit(); // remove former limit
            unsigned int nSize = 0;
            try {
                // locate a header
                unsigned char buf[4];
                blkdat.FindByte(pchMessageStart[0]);
                nRewind = blkdat.GetPos()+1;
                blkdat >> FLATDATA(buf);
                if (memcmp(buf, pchMessageStart, 4))
                    continue;
                // read size
                blkdat >> nSize;
                if (nSize < 80 || nSize > MAX_BLOCK_SIZE)
                    continue;
            } catch (std::exception &e) {
                // no valid block header found; don't complain
                break;
            }
            try {
                // read block
                uint64 nBlockPos = blkdat.GetPos();
                blkdat.SetLimit(nBlockPos + nSize);
                CBlock block;
                blkdat >> block;
                nRewind = blkdat.GetPos();

                // process block
                if (nBlockPos >= nStartByte) {
                    LOCK(cs_main);
                    if (dbp)
                        dbp->nPos = nBlockPos;
                    CValidationState state;
                    if (ProcessBlock(state, NULL, &block, dbp))
                        nLoaded++;
                    if (state.IsError())
                        break;
                }
            } catch (std::exception &e) {
                printf("%s() : Deserialize or I/O error caught during load\n", __PRETTY_FUNCTION__);
            }
        }
        fclose(fileIn);
    } catch(std::runtime_error &e) {
        AbortNode(_("Error: system error: ") + e.what());
    }
    if (nLoaded > 0)
        printf("Loaded %i blocks from external file in %"PRI64d"ms\n", nLoaded, GetTimeMillis() - nStart);
    return nLoaded > 0;
}










//////////////////////////////////////////////////////////////////////////////
//
// CAlert
//

extern map<uint256, CAlert> mapAlerts;
extern CCriticalSection cs_mapAlerts;

string GetWarnings(string strFor)
{
    int nPriority = 0;
    string strStatusBar;
    string strRPC;

    if (GetBoolArg("-testsafemode"))
        strRPC = "test";

    if (!CLIENT_VERSION_IS_RELEASE)
        strStatusBar = _("This is a pre-release test build - use at your own risk - do not use for mining or merchant applications");

    // Misc warnings like out of disk space and clock is wrong
    if (strMiscWarning != "")
    {
        nPriority = 1000;
        strStatusBar = strMiscWarning;
    }

    // Longer invalid proof-of-work chain
    if (pindexBest && nBestInvalidWork > nBestChainWork + (pindexBest->GetBlockWork() * 6).getuint256())
    {
        nPriority = 2000;
        strStatusBar = strRPC = _("Warning: Displayed transactions may not be correct! You may need to upgrade, or other nodes may need to upgrade.");
    }

    // Alerts
    {
        LOCK(cs_mapAlerts);
        BOOST_FOREACH(PAIRTYPE(const uint256, CAlert)& item, mapAlerts)
        {
            const CAlert& alert = item.second;
            if (alert.AppliesToMe() && alert.nPriority > nPriority)
            {
                nPriority = alert.nPriority;
                strStatusBar = alert.strStatusBar;
            }
        }
    }

    if (strFor == "statusbar")
        return strStatusBar;
    else if (strFor == "rpc")
        return strRPC;
    assert(!"GetWarnings() : invalid parameter");
    return "error";
}








//////////////////////////////////////////////////////////////////////////////
//
// Messages
//


bool static AlreadyHave(const CInv& inv)
{
    switch (inv.type)
    {
    case MSG_TX:
        {
            bool txInMap = false;
            {
                LOCK(mempool.cs);
                txInMap = mempool.exists(inv.hash);
            }
            return txInMap || mapOrphanTransactions.count(inv.hash) ||
                pcoinsTip->HaveCoins(inv.hash);
        }
    case MSG_BLOCK:
        return mapBlockIndex.count(inv.hash) ||
               mapOrphanBlocks.count(inv.hash);
    }
    // Don't know what it is, just say we already got one
    return true;
}




// The message start string is designed to be unlikely to occur in normal data.
// The characters are rarely used upper ASCII, not valid as UTF-8, and produce
// a large 4-byte int at any alignment.
// FBX
//unsigned char pchMessageStart[4] = { 0xfb, 0xc0, 0xb6, 0xdb }; // Litecoin: increase each by adding 2 to bitcoin's value.
unsigned char pchMessageStart[4] = { 0xf9, 0xbe, 0xb4, 0xd9 };

void static ProcessGetData(CNode* pfrom)
{
    std::deque<CInv>::iterator it = pfrom->vRecvGetData.begin();

    vector<CInv> vNotFound;

    while (it != pfrom->vRecvGetData.end()) {
        // Don't bother if send buffer is too full to respond anyway
        if (pfrom->nSendSize >= SendBufferSize())
            break;

        // Don't waste work on slow peers until they catch up on the blocks we
        // give them. 80 bytes is just the size of a block header - obviously
        // the minimum we might return.
        if (pfrom->nBlocksRequested * 80 > pfrom->nSendBytes)
            break;

        const CInv &inv = *it;
        {
            boost::this_thread::interruption_point();
            it++;

            if (inv.type == MSG_BLOCK || inv.type == MSG_FILTERED_BLOCK)
            {
                // Send block from disk
                map<uint256, CBlockIndex*>::iterator mi = mapBlockIndex.find(inv.hash);
                pfrom->nBlocksRequested++;
                if (mi != mapBlockIndex.end())
                {
                    CBlock block;
                    block.ReadFromDisk((*mi).second);
                    if (inv.type == MSG_BLOCK)
                        pfrom->PushMessage("block", block);
                    else // MSG_FILTERED_BLOCK)
                    {
                        LOCK(pfrom->cs_filter);
                        if (pfrom->pfilter)
                        {
                            CMerkleBlock merkleBlock(block, *pfrom->pfilter);
                            pfrom->PushMessage("merkleblock", merkleBlock);
                            // CMerkleBlock just contains hashes, so also push any transactions in the block the client did not see
                            // This avoids hurting performance by pointlessly requiring a round-trip
                            // Note that there is currently no way for a node to request any single transactions we didnt send here -
                            // they must either disconnect and retry or request the full block.
                            // Thus, the protocol spec specified allows for us to provide duplicate txn here,
                            // however we MUST always provide at least what the remote peer needs
                            typedef std::pair<unsigned int, uint256> PairType;
                            BOOST_FOREACH(PairType& pair, merkleBlock.vMatchedTxn)
                                if (!pfrom->setInventoryKnown.count(CInv(MSG_TX, pair.second)))
                                    pfrom->PushMessage("tx", block.vtx[pair.first]);
                        }
                        // else
                            // no response
                    }

                    // Trigger them to send a getblocks request for the next batch of inventory
                    if (inv.hash == pfrom->hashContinue)
                    {
                        // Bypass PushInventory, this must send even if redundant,
                        // and we want it right after the last block so they don't
                        // wait for other stuff first.
                        vector<CInv> vInv;
                        vInv.push_back(CInv(MSG_BLOCK, hashBestChain));
                        pfrom->PushMessage("inv", vInv);
                        pfrom->hashContinue = 0;
                    }
                }
            }
            else if (inv.IsKnownType())
            {
                // Send stream from relay memory
                bool pushed = false;
                {
                    LOCK(cs_mapRelay);
                    map<CInv, CDataStream>::iterator mi = mapRelay.find(inv);
                    if (mi != mapRelay.end()) {
                        pfrom->PushMessage(inv.GetCommand(), (*mi).second);
                        pushed = true;
                    }
                }
                if (!pushed && inv.type == MSG_TX) {
                    LOCK(mempool.cs);
                    if (mempool.exists(inv.hash)) {
                        CTransaction tx = mempool.lookup(inv.hash);
                        CDataStream ss(SER_NETWORK, PROTOCOL_VERSION);
                        ss.reserve(1000);
                        ss << tx;
                        pfrom->PushMessage("tx", ss);
                        pushed = true;
                    }
                }
                if (!pushed) {
                    vNotFound.push_back(inv);
                }
            }

            // Track requests for our stuff.
            Inventory(inv.hash);
        }
    }

    pfrom->vRecvGetData.erase(pfrom->vRecvGetData.begin(), it);

    if (!vNotFound.empty()) {
        // Let the peer know that we didn't find what it asked for, so it doesn't
        // have to wait around forever. Currently only SPV clients actually care
        // about this message: it's needed when they are recursively walking the
        // dependencies of relevant unconfirmed transactions. SPV clients want to
        // do that because they want to know about (and store and rebroadcast and
        // risk analyze) the dependencies of transactions relevant to them, without
        // having to download the entire memory pool.
        pfrom->PushMessage("notfound", vNotFound);
    }
}

bool static ProcessMessage(CNode* pfrom, string strCommand, CDataStream& vRecv)
{
    RandAddSeedPerfmon();
    if (fDebug)
        printf("received: %s (%"PRIszu" bytes)\n", strCommand.c_str(), vRecv.size());
    if (mapArgs.count("-dropmessagestest") && GetRand(atoi(mapArgs["-dropmessagestest"])) == 0)
    {
        printf("dropmessagestest DROPPING RECV MESSAGE\n");
        return true;
    }





    if (strCommand == "version")
    {
        // Each connection can only send one version message
        if (pfrom->nVersion != 0)
        {
            pfrom->Misbehaving(1);
            return false;
        }

        int64 nTime;
        CAddress addrMe;
        CAddress addrFrom;
        uint64 nNonce = 1;
        vRecv >> pfrom->nVersion >> pfrom->nServices >> nTime >> addrMe;
        if (pfrom->nVersion < MIN_PROTO_VERSION)
        {
            // Since February 20, 2012, the protocol is initiated at version 209,
            // and earlier versions are no longer supported
            printf("partner %s using obsolete version %i; disconnecting\n", pfrom->addr.ToString().c_str(), pfrom->nVersion);
            pfrom->fDisconnect = true;
            return false;
        }

        if (pfrom->nVersion == 10300)
            pfrom->nVersion = 300;
        if (!vRecv.empty())
            vRecv >> addrFrom >> nNonce;
        if (!vRecv.empty())
            vRecv >> pfrom->strSubVer;
        if (!vRecv.empty())
            vRecv >> pfrom->nStartingHeight;
        if (!vRecv.empty())
            vRecv >> pfrom->fRelayTxes; // set to true after we get the first filter* message
        else
            pfrom->fRelayTxes = true;

        if (pfrom->fInbound && addrMe.IsRoutable())
        {
            pfrom->addrLocal = addrMe;
            SeenLocal(addrMe);
        }

        // Disconnect if we connected to ourself
        if (nNonce == nLocalHostNonce && nNonce > 1)
        {
            printf("connected to self at %s, disconnecting\n", pfrom->addr.ToString().c_str());
            pfrom->fDisconnect = true;
            return true;
        }

        // Be shy and don't send version until we hear
        if (pfrom->fInbound)
            pfrom->PushVersion();

        pfrom->fClient = !(pfrom->nServices & NODE_NETWORK);

        AddTimeData(pfrom->addr, nTime);

        // Change version
        pfrom->PushMessage("verack");
        pfrom->ssSend.SetVersion(min(pfrom->nVersion, PROTOCOL_VERSION));

        if (!pfrom->fInbound)
        {
            // Advertise our address
            if (!fNoListen && !IsInitialBlockDownload())
            {
                CAddress addr = GetLocalAddress(&pfrom->addr);
                if (addr.IsRoutable())
                    pfrom->PushAddress(addr);
            }

            // Get recent addresses
            if (pfrom->fOneShot || pfrom->nVersion >= CADDR_TIME_VERSION || addrman.size() < 1000)
            {
                pfrom->PushMessage("getaddr");
                pfrom->fGetAddr = true;
            }
            addrman.Good(pfrom->addr);
        } else {
            if (((CNetAddr)pfrom->addr) == (CNetAddr)addrFrom)
            {
                addrman.Add(addrFrom, addrFrom);
                addrman.Good(addrFrom);
            }
        }

        // Relay alerts
        {
            LOCK(cs_mapAlerts);
            BOOST_FOREACH(PAIRTYPE(const uint256, CAlert)& item, mapAlerts)
                item.second.RelayTo(pfrom);
        }

        pfrom->fSuccessfullyConnected = true;

        printf("receive version message: version %d, blocks=%d, us=%s, them=%s, peer=%s\n", pfrom->nVersion, pfrom->nStartingHeight, addrMe.ToString().c_str(), addrFrom.ToString().c_str(), pfrom->addr.ToString().c_str());

        cPeerBlockCounts.input(pfrom->nStartingHeight);
    }


    else if (pfrom->nVersion == 0)
    {
        // Must have a version message before anything else
        pfrom->Misbehaving(1);
        return false;
    }


    else if (strCommand == "verack")
    {
        pfrom->SetRecvVersion(min(pfrom->nVersion, PROTOCOL_VERSION));
    }


    else if (strCommand == "addr")
    {
        vector<CAddress> vAddr;
        vRecv >> vAddr;

        // Don't want addr from older versions unless seeding
        if (pfrom->nVersion < CADDR_TIME_VERSION && addrman.size() > 1000)
            return true;
        if (vAddr.size() > 1000)
        {
            pfrom->Misbehaving(20);
            return error("message addr size() = %"PRIszu"", vAddr.size());
        }

        // Store the new addresses
        vector<CAddress> vAddrOk;
        int64 nNow = GetAdjustedTime();
        int64 nSince = nNow - 10 * 60;
        BOOST_FOREACH(CAddress& addr, vAddr)
        {
            boost::this_thread::interruption_point();

            if (addr.nTime <= 100000000 || addr.nTime > nNow + 10 * 60)
                addr.nTime = nNow - 5 * 24 * 60 * 60;
            pfrom->AddAddressKnown(addr);
            bool fReachable = IsReachable(addr);
            if (addr.nTime > nSince && !pfrom->fGetAddr && vAddr.size() <= 10 && addr.IsRoutable())
            {
                // Relay to a limited number of other nodes
                {
                    LOCK(cs_vNodes);
                    // Use deterministic randomness to send to the same nodes for 24 hours
                    // at a time so the setAddrKnowns of the chosen nodes prevent repeats
                    static uint256 hashSalt;
                    if (hashSalt == 0)
                        hashSalt = GetRandHash();
                    uint64 hashAddr = addr.GetHash();
                    uint256 hashRand = hashSalt ^ (hashAddr<<32) ^ ((GetTime()+hashAddr)/(24*60*60));
                    hashRand = Hash(BEGIN(hashRand), END(hashRand));
                    multimap<uint256, CNode*> mapMix;
                    BOOST_FOREACH(CNode* pnode, vNodes)
                    {
                        if (pnode->nVersion < CADDR_TIME_VERSION)
                            continue;
                        unsigned int nPointer;
                        memcpy(&nPointer, &pnode, sizeof(nPointer));
                        uint256 hashKey = hashRand ^ nPointer;
                        hashKey = Hash(BEGIN(hashKey), END(hashKey));
                        mapMix.insert(make_pair(hashKey, pnode));
                    }
                    int nRelayNodes = fReachable ? 2 : 1; // limited relaying of addresses outside our network(s)
                    for (multimap<uint256, CNode*>::iterator mi = mapMix.begin(); mi != mapMix.end() && nRelayNodes-- > 0; ++mi)
                        ((*mi).second)->PushAddress(addr);
                }
            }
            // Do not store addresses outside our network
            if (fReachable)
                vAddrOk.push_back(addr);
        }
        addrman.Add(vAddrOk, pfrom->addr, 2 * 60 * 60);
        if (vAddr.size() < 1000)
            pfrom->fGetAddr = false;
        if (pfrom->fOneShot)
            pfrom->fDisconnect = true;
    }


    else if (strCommand == "inv")
    {
        vector<CInv> vInv;
        vRecv >> vInv;
        if (vInv.size() > MAX_INV_SZ)
        {
            pfrom->Misbehaving(20);
            return error("message inv size() = %"PRIszu"", vInv.size());
        }

        // find last block in inv vector
        unsigned int nLastBlock = (unsigned int)(-1);
        for (unsigned int nInv = 0; nInv < vInv.size(); nInv++) {
            if (vInv[vInv.size() - 1 - nInv].type == MSG_BLOCK) {
                nLastBlock = vInv.size() - 1 - nInv;
                break;
            }
        }
        for (unsigned int nInv = 0; nInv < vInv.size(); nInv++)
        {
            const CInv &inv = vInv[nInv];

            boost::this_thread::interruption_point();
            pfrom->AddInventoryKnown(inv);

            bool fAlreadyHave = AlreadyHave(inv);
            if (fDebug)
                printf("  got inventory: %s  %s\n", inv.ToString().c_str(), fAlreadyHave ? "have" : "new");

            if (!fAlreadyHave) {
                if (!fImporting && !fReindex)
                    pfrom->AskFor(inv);
            } else if (inv.type == MSG_BLOCK && mapOrphanBlocks.count(inv.hash)) {
                pfrom->PushGetBlocks(pindexBest, GetOrphanRoot(mapOrphanBlocks[inv.hash]));
            } else if (nInv == nLastBlock) {
                // In case we are on a very long side-chain, it is possible that we already have
                // the last block in an inv bundle sent in response to getblocks. Try to detect
                // this situation and push another getblocks to continue.
                pfrom->PushGetBlocks(mapBlockIndex[inv.hash], uint256(0));
                if (fDebug)
                    printf("force request: %s\n", inv.ToString().c_str());
            }

            // Track requests for our stuff
            Inventory(inv.hash);
        }
    }


    else if (strCommand == "getdata")
    {
        vector<CInv> vInv;
        vRecv >> vInv;
        if (vInv.size() > MAX_INV_SZ)
        {
            pfrom->Misbehaving(20);
            return error("message getdata size() = %"PRIszu"", vInv.size());
        }

        if (fDebugNet || (vInv.size() != 1))
            printf("received getdata (%"PRIszu" invsz)\n", vInv.size());

        if ((fDebugNet && vInv.size() > 0) || (vInv.size() == 1))
            printf("received getdata for: %s\n", vInv[0].ToString().c_str());

        pfrom->vRecvGetData.insert(pfrom->vRecvGetData.end(), vInv.begin(), vInv.end());
        ProcessGetData(pfrom);
    }


    else if (strCommand == "getblocks")
    {
        CBlockLocator locator;
        uint256 hashStop;
        vRecv >> locator >> hashStop;

        // Find the last block the caller has in the main chain
        CBlockIndex* pindex = locator.GetBlockIndex();

        // Send the rest of the chain
        if (pindex)
            pindex = pindex->pnext;
        int nLimit = 500;
        printf("getblocks %d to %s limit %d\n", (pindex ? pindex->nHeight : -1), hashStop.ToString().c_str(), nLimit);
        for (; pindex; pindex = pindex->pnext)
        {
            if (pindex->GetBlockHash() == hashStop)
            {
                printf("  getblocks stopping at %d %s\n", pindex->nHeight, pindex->GetBlockHash().ToString().c_str());
                break;
            }
            pfrom->PushInventory(CInv(MSG_BLOCK, pindex->GetBlockHash()));
            if (--nLimit <= 0)
            {
                // When this block is requested, we'll send an inv that'll make them
                // getblocks the next batch of inventory.
                printf("  getblocks stopping at limit %d %s\n", pindex->nHeight, pindex->GetBlockHash().ToString().c_str());
                pfrom->hashContinue = pindex->GetBlockHash();
                break;
            }
        }
    }


    else if (strCommand == "getheaders")
    {
        CBlockLocator locator;
        uint256 hashStop;
        vRecv >> locator >> hashStop;

        CBlockIndex* pindex = NULL;
        if (locator.IsNull())
        {
            // If locator is null, return the hashStop block
            map<uint256, CBlockIndex*>::iterator mi = mapBlockIndex.find(hashStop);
            if (mi == mapBlockIndex.end())
                return true;
            pindex = (*mi).second;
        }
        else
        {
            // Find the last block the caller has in the main chain
            pindex = locator.GetBlockIndex();
            if (pindex)
                pindex = pindex->pnext;
        }

        // we must use CBlocks, as CBlockHeaders won't include the 0x00 nTx count at the end
        vector<CBlock> vHeaders;
        int nLimit = 2000;
        printf("getheaders %d to %s\n", (pindex ? pindex->nHeight : -1), hashStop.ToString().c_str());
        for (; pindex; pindex = pindex->pnext)
        {
            vHeaders.push_back(pindex->GetBlockHeader());
            if (--nLimit <= 0 || pindex->GetBlockHash() == hashStop)
                break;
        }
        pfrom->PushMessage("headers", vHeaders);
    }


    else if (strCommand == "tx")
    {
        vector<uint256> vWorkQueue;
        vector<uint256> vEraseQueue;
        CDataStream vMsg(vRecv);
        CTransaction tx;
        vRecv >> tx;

        CInv inv(MSG_TX, tx.GetHash());
        pfrom->AddInventoryKnown(inv);

        bool fMissingInputs = false;
        CValidationState state;
        if (tx.AcceptToMemoryPool(state, true, true, &fMissingInputs))
        {
            RelayTransaction(tx, inv.hash);
            mapAlreadyAskedFor.erase(inv);
            vWorkQueue.push_back(inv.hash);
            vEraseQueue.push_back(inv.hash);

            // Recursively process any orphan transactions that depended on this one
            for (unsigned int i = 0; i < vWorkQueue.size(); i++)
            {
                uint256 hashPrev = vWorkQueue[i];
                for (set<uint256>::iterator mi = mapOrphanTransactionsByPrev[hashPrev].begin();
                     mi != mapOrphanTransactionsByPrev[hashPrev].end();
                     ++mi)
                {
                    const uint256& orphanHash = *mi;
                    const CTransaction& orphanTx = mapOrphanTransactions[orphanHash];
                    bool fMissingInputs2 = false;
                    // Use a dummy CValidationState so someone can't setup nodes to counter-DoS based on orphan
                    // resolution (that is, feeding people an invalid transaction based on LegitTxX in order to get
                    // anyone relaying LegitTxX banned)
                    CValidationState stateDummy;

                    if (tx.AcceptToMemoryPool(stateDummy, true, true, &fMissingInputs2))
                    {
                        printf("   accepted orphan tx %s\n", orphanHash.ToString().c_str());
                        RelayTransaction(orphanTx, orphanHash);
                        mapAlreadyAskedFor.erase(CInv(MSG_TX, orphanHash));
                        vWorkQueue.push_back(orphanHash);
                        vEraseQueue.push_back(orphanHash);
                    }
                    else if (!fMissingInputs2)
                    {
                        // invalid or too-little-fee orphan
                        vEraseQueue.push_back(orphanHash);
                        printf("   removed orphan tx %s\n", orphanHash.ToString().c_str());
                    }
                }
            }

            BOOST_FOREACH(uint256 hash, vEraseQueue)
                EraseOrphanTx(hash);
        }
        else if (fMissingInputs)
        {
            AddOrphanTx(tx);

            // DoS prevention: do not allow mapOrphanTransactions to grow unbounded
            unsigned int nEvicted = LimitOrphanTxSize(MAX_ORPHAN_TRANSACTIONS);
            if (nEvicted > 0)
                printf("mapOrphan overflow, removed %u tx\n", nEvicted);
        }
        int nDoS;
        if (state.IsInvalid(nDoS))
            pfrom->Misbehaving(nDoS);
    }


    else if (strCommand == "block" && !fImporting && !fReindex) // Ignore blocks received while importing
    {
        CBlock block;
        vRecv >> block;

        printf("received block %s\n", block.GetHash().ToString().c_str());
        // block.print();

        CInv inv(MSG_BLOCK, block.GetHash());
        pfrom->AddInventoryKnown(inv);

        CValidationState state;
        if (ProcessBlock(state, pfrom, &block))
            mapAlreadyAskedFor.erase(inv);
        int nDoS;
        if (state.IsInvalid(nDoS))
            pfrom->Misbehaving(nDoS);
    }


    else if (strCommand == "getaddr")
    {
        pfrom->vAddrToSend.clear();
        vector<CAddress> vAddr = addrman.GetAddr();
        BOOST_FOREACH(const CAddress &addr, vAddr)
            pfrom->PushAddress(addr);
    }


    else if (strCommand == "mempool")
    {
        std::vector<uint256> vtxid;
        LOCK2(mempool.cs, pfrom->cs_filter);
        mempool.queryHashes(vtxid);
        vector<CInv> vInv;
        BOOST_FOREACH(uint256& hash, vtxid) {
            CInv inv(MSG_TX, hash);
            if ((pfrom->pfilter && pfrom->pfilter->IsRelevantAndUpdate(mempool.lookup(hash), hash)) ||
               (!pfrom->pfilter))
                vInv.push_back(inv);
            if (vInv.size() == MAX_INV_SZ)
                break;
        }
        if (vInv.size() > 0)
            pfrom->PushMessage("inv", vInv);
    }


    else if (strCommand == "ping")
    {
        if (pfrom->nVersion > BIP0031_VERSION)
        {
            uint64 nonce = 0;
            vRecv >> nonce;
            // Echo the message back with the nonce. This allows for two useful features:
            //
            // 1) A remote node can quickly check if the connection is operational
            // 2) Remote nodes can measure the latency of the network thread. If this node
            //    is overloaded it won't respond to pings quickly and the remote node can
            //    avoid sending us more work, like chain download requests.
            //
            // The nonce stops the remote getting confused between different pings: without
            // it, if the remote node sends a ping once per second and this node takes 5
            // seconds to respond to each, the 5th ping the remote sends would appear to
            // return very quickly.
            pfrom->PushMessage("pong", nonce);
        }
    }


    else if (strCommand == "alert")
    {
        CAlert alert;
        vRecv >> alert;

        uint256 alertHash = alert.GetHash();
        if (pfrom->setKnown.count(alertHash) == 0)
        {
            if (alert.ProcessAlert())
            {
                // Relay
                pfrom->setKnown.insert(alertHash);
                {
                    LOCK(cs_vNodes);
                    BOOST_FOREACH(CNode* pnode, vNodes)
                        alert.RelayTo(pnode);
                }
            }
            else {
                // Small DoS penalty so peers that send us lots of
                // duplicate/expired/invalid-signature/whatever alerts
                // eventually get banned.
                // This isn't a Misbehaving(100) (immediate ban) because the
                // peer might be an older or different implementation with
                // a different signature key, etc.
                pfrom->Misbehaving(10);
            }
        }
    }


    else if (!fBloomFilters &&
             (strCommand == "filterload" ||
              strCommand == "filteradd" ||
              strCommand == "filterclear"))
    {
        pfrom->Misbehaving(100);
        return error("peer %s attempted to set a bloom filter even though we do not advertise that service",
                     pfrom->addr.ToString().c_str());
    }

    else if (strCommand == "filterload")
    {
        CBloomFilter filter;
        vRecv >> filter;

        if (!filter.IsWithinSizeConstraints())
            // There is no excuse for sending a too-large filter
            pfrom->Misbehaving(100);
        else
        {
            LOCK(pfrom->cs_filter);
            delete pfrom->pfilter;
            pfrom->pfilter = new CBloomFilter(filter);
            pfrom->pfilter->UpdateEmptyFull();
        }
        pfrom->fRelayTxes = true;
    }


    else if (strCommand == "filteradd")
    {
        vector<unsigned char> vData;
        vRecv >> vData;

        // Nodes must NEVER send a data item > 520 bytes (the max size for a script data object,
        // and thus, the maximum size any matched object can have) in a filteradd message
        if (vData.size() > MAX_SCRIPT_ELEMENT_SIZE)
        {
            pfrom->Misbehaving(100);
        } else {
            LOCK(pfrom->cs_filter);
            if (pfrom->pfilter)
                pfrom->pfilter->insert(vData);
            else
                pfrom->Misbehaving(100);
        }
    }


    else if (strCommand == "filterclear")
    {
        LOCK(pfrom->cs_filter);
        delete pfrom->pfilter;
        pfrom->pfilter = new CBloomFilter();
        pfrom->fRelayTxes = true;
    }


    else
    {
        // Ignore unknown commands for extensibility
    }


    // Update the last seen time for this node's address
    if (pfrom->fNetworkNode)
        if (strCommand == "version" || strCommand == "addr" || strCommand == "inv" || strCommand == "getdata" || strCommand == "ping")
            AddressCurrentlyConnected(pfrom->addr);


    return true;
}

// requires LOCK(cs_vRecvMsg)
bool ProcessMessages(CNode* pfrom)
{
    //if (fDebug)
    //    printf("ProcessMessages(%zu messages)\n", pfrom->vRecvMsg.size());

    //
    // Message format
    //  (4) message start
    //  (12) command
    //  (4) size
    //  (4) checksum
    //  (x) data
    //
    bool fOk = true;

    if (!pfrom->vRecvGetData.empty())
        ProcessGetData(pfrom);

    std::deque<CNetMessage>::iterator it = pfrom->vRecvMsg.begin();
    while (!pfrom->fDisconnect && it != pfrom->vRecvMsg.end()) {
        // Don't bother if send buffer is too full to respond anyway
        if (pfrom->nSendSize >= SendBufferSize())
            break;

        // get next message
        CNetMessage& msg = *it;

        //if (fDebug)
        //    printf("ProcessMessages(message %u msgsz, %zu bytes, complete:%s)\n",
        //            msg.hdr.nMessageSize, msg.vRecv.size(),
        //            msg.complete() ? "Y" : "N");

        // end, if an incomplete message is found
        if (!msg.complete())
            break;

        // at this point, any failure means we can delete the current message
        it++;

        // Scan for message start
        if (memcmp(msg.hdr.pchMessageStart, pchMessageStart, sizeof(pchMessageStart)) != 0) {
            printf("\n\nPROCESSMESSAGE: INVALID MESSAGESTART\n\n");
            fOk = false;
            break;
        }

        // Read header
        CMessageHeader& hdr = msg.hdr;
        if (!hdr.IsValid())
        {
            printf("\n\nPROCESSMESSAGE: ERRORS IN HEADER %s\n\n\n", hdr.GetCommand().c_str());
            continue;
        }
        string strCommand = hdr.GetCommand();

        // Message size
        unsigned int nMessageSize = hdr.nMessageSize;

        // Checksum
        CDataStream& vRecv = msg.vRecv;
        uint256 hash = Hash(vRecv.begin(), vRecv.begin() + nMessageSize);
        unsigned int nChecksum = 0;
        memcpy(&nChecksum, &hash, sizeof(nChecksum));
        if (nChecksum != hdr.nChecksum)
        {
            printf("ProcessMessages(%s, %u bytes) : CHECKSUM ERROR nChecksum=%08x hdr.nChecksum=%08x\n",
               strCommand.c_str(), nMessageSize, nChecksum, hdr.nChecksum);
            continue;
        }

        // Process message
        bool fRet = false;
        try
        {
            {
                LOCK(cs_main);
                fRet = ProcessMessage(pfrom, strCommand, vRecv);
            }
            boost::this_thread::interruption_point();
        }
        catch (std::ios_base::failure& e)
        {
            if (strstr(e.what(), "end of data"))
            {
                // Allow exceptions from under-length message on vRecv
                printf("ProcessMessages(%s, %u bytes) : Exception '%s' caught, normally caused by a message being shorter than its stated length\n", strCommand.c_str(), nMessageSize, e.what());
            }
            else if (strstr(e.what(), "size too large"))
            {
                // Allow exceptions from over-long size
                printf("ProcessMessages(%s, %u bytes) : Exception '%s' caught\n", strCommand.c_str(), nMessageSize, e.what());
            }
            else
            {
                PrintExceptionContinue(&e, "ProcessMessages()");
            }
        }
        catch (boost::thread_interrupted) {
            throw;
        }
        catch (std::exception& e) {
            PrintExceptionContinue(&e, "ProcessMessages()");
        } catch (...) {
            PrintExceptionContinue(NULL, "ProcessMessages()");
        }

        if (!fRet)
            printf("ProcessMessage(%s, %u bytes) FAILED\n", strCommand.c_str(), nMessageSize);
    }

    // In case the connection got shut down, its receive buffer was wiped
    if (!pfrom->fDisconnect)
        pfrom->vRecvMsg.erase(pfrom->vRecvMsg.begin(), it);

    return fOk;
}


bool SendMessages(CNode* pto, bool fSendTrickle)
{
    TRY_LOCK(cs_main, lockMain);
    if (lockMain) {
        // Don't send anything until we get their version message
        if (pto->nVersion == 0)
            return true;

        // Keep-alive ping. We send a nonce of zero because we don't use it anywhere
        // right now.
// FBX ping every 2 minutes, timeout 5 minutes
//        if (pto->nLastSend && GetTime() - pto->nLastSend > 30 * 60 && pto->vSendMsg.empty()) {
        if (pto->nLastSend && GetTime() - pto->nLastSend > 2 * 60 && pto->vSendMsg.empty()) {
            uint64 nonce = 0;
            if (pto->nVersion > BIP0031_VERSION)
                pto->PushMessage("ping", nonce);
            else
                pto->PushMessage("ping");
        }

        // Start block sync
        if (pto->fStartSync && !fImporting && !fReindex) {
            pto->fStartSync = false;
            pto->PushGetBlocks(pindexBest, uint256(0));
        }

        // Resend wallet transactions that haven't gotten in a block yet
        // Except during reindex, importing and IBD, when old wallet
        // transactions become unconfirmed and spams other nodes.
        if (!fReindex && !fImporting && !IsInitialBlockDownload())
        {
            ResendWalletTransactions();
        }

        // Address refresh broadcast
        static int64 nLastRebroadcast;
        if (!IsInitialBlockDownload() && (GetTime() - nLastRebroadcast > 24 * 60 * 60))
        {
            {
                LOCK(cs_vNodes);
                BOOST_FOREACH(CNode* pnode, vNodes)
                {
                    // Periodically clear setAddrKnown to allow refresh broadcasts
                    if (nLastRebroadcast)
                        pnode->setAddrKnown.clear();

                    // Rebroadcast our address
                    if (!fNoListen)
                    {
                        CAddress addr = GetLocalAddress(&pnode->addr);
                        if (addr.IsRoutable())
                            pnode->PushAddress(addr);
                    }
                }
            }
            nLastRebroadcast = GetTime();
        }

        //
        // Message: addr
        //
        if (fSendTrickle)
        {
            vector<CAddress> vAddr;
            vAddr.reserve(pto->vAddrToSend.size());
            BOOST_FOREACH(const CAddress& addr, pto->vAddrToSend)
            {
                // returns true if wasn't already contained in the set
                if (pto->setAddrKnown.insert(addr).second)
                {
                    vAddr.push_back(addr);
                    // receiver rejects addr messages larger than 1000
                    if (vAddr.size() >= 1000)
                    {
                        pto->PushMessage("addr", vAddr);
                        vAddr.clear();
                    }
                }
            }
            pto->vAddrToSend.clear();
            if (!vAddr.empty())
                pto->PushMessage("addr", vAddr);
        }


        //
        // Message: inventory
        //
        vector<CInv> vInv;
        vector<CInv> vInvWait;
        {
            LOCK(pto->cs_inventory);
            vInv.reserve(pto->vInventoryToSend.size());
            vInvWait.reserve(pto->vInventoryToSend.size());
            BOOST_FOREACH(const CInv& inv, pto->vInventoryToSend)
            {
                if (pto->setInventoryKnown.count(inv))
                    continue;

                // trickle out tx inv to protect privacy
                if (inv.type == MSG_TX && !fSendTrickle)
                {
                    // 1/4 of tx invs blast to all immediately
                    static uint256 hashSalt;
                    if (hashSalt == 0)
                        hashSalt = GetRandHash();
                    uint256 hashRand = inv.hash ^ hashSalt;
                    hashRand = Hash(BEGIN(hashRand), END(hashRand));
                    bool fTrickleWait = ((hashRand & 3) != 0);

                    // always trickle our own transactions
                    if (!fTrickleWait)
                    {
                        CWalletTx wtx;
                        if (GetTransaction(inv.hash, wtx))
                            if (wtx.fFromMe)
                                fTrickleWait = true;
                    }

                    if (fTrickleWait)
                    {
                        vInvWait.push_back(inv);
                        continue;
                    }
                }

                // returns true if wasn't already contained in the set
                if (pto->setInventoryKnown.insert(inv).second)
                {
                    vInv.push_back(inv);
                    if (vInv.size() >= 1000)
                    {
                        pto->PushMessage("inv", vInv);
                        vInv.clear();
                    }
                }
            }
            pto->vInventoryToSend = vInvWait;
        }
        if (!vInv.empty())
            pto->PushMessage("inv", vInv);


        //
        // Message: getdata
        //
        vector<CInv> vGetData;
        int64 nNow = GetTime() * 1000000;
        while (!pto->mapAskFor.empty() && (*pto->mapAskFor.begin()).first <= nNow)
        {
            const CInv& inv = (*pto->mapAskFor.begin()).second;
            if (!AlreadyHave(inv))
            {
                if (fDebugNet)
                    printf("sending getdata: %s\n", inv.ToString().c_str());
                vGetData.push_back(inv);
                if (vGetData.size() >= 1000)
                {
                    pto->PushMessage("getdata", vGetData);
                    vGetData.clear();
                }
            }
            pto->mapAskFor.erase(pto->mapAskFor.begin());
        }
        if (!vGetData.empty())
            pto->PushMessage("getdata", vGetData);

    }
    return true;
}

//////////////////////////////////////////////////////////////////////////////
//
// LitecoinMiner
//

int static FormatHashBlocks(void* pbuffer, unsigned int len)
{
    unsigned char* pdata = (unsigned char*)pbuffer;
    unsigned int blocks = 1 + ((len + 8) / 64);
    unsigned char* pend = pdata + 64 * blocks;
    memset(pdata + len, 0, 64 * blocks - len);
    pdata[len] = 0x80;
    unsigned int bits = len * 8;
    pend[-1] = (bits >> 0) & 0xff;
    pend[-2] = (bits >> 8) & 0xff;
    pend[-3] = (bits >> 16) & 0xff;
    pend[-4] = (bits >> 24) & 0xff;
    return blocks;
}

static const unsigned int pSHA256InitState[8] =
{0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};

void SHA256Transform(void* pstate, void* pinput, const void* pinit)
{
    SHA256_CTX ctx;
    unsigned char data[64];

    SHA256_Init(&ctx);

    for (int i = 0; i < 16; i++)
        ((uint32_t*)data)[i] = ByteReverse(((uint32_t*)pinput)[i]);

    for (int i = 0; i < 8; i++)
        ctx.h[i] = ((uint32_t*)pinit)[i];

    SHA256_Update(&ctx, data, sizeof(data));
    for (int i = 0; i < 8; i++)
        ((uint32_t*)pstate)[i] = ctx.h[i];
}

// Some explaining would be appreciated
class COrphan
{
public:
    CTransaction* ptx;
    set<uint256> setDependsOn;
    double dPriority;
    double dFeePerKb;

    COrphan(CTransaction* ptxIn)
    {
        ptx = ptxIn;
        dPriority = dFeePerKb = 0;
    }

    void print() const
    {
        printf("COrphan(hash=%s, dPriority=%.1f, dFeePerKb=%.1f)\n",
               ptx->GetHash().ToString().c_str(), dPriority, dFeePerKb);
        BOOST_FOREACH(uint256 hash, setDependsOn)
            printf("   setDependsOn %s\n", hash.ToString().c_str());
    }
};


uint64 nLastBlockTx = 0;
uint64 nLastBlockSize = 0;

// We want to sort transactions by priority and fee, so:
typedef boost::tuple<double, double, CTransaction*> TxPriority;
class TxPriorityCompare
{
    bool byFee;
public:
    TxPriorityCompare(bool _byFee) : byFee(_byFee) { }
    bool operator()(const TxPriority& a, const TxPriority& b)
    {
        if (byFee)
        {
            if (a.get<1>() == b.get<1>())
                return a.get<0>() < b.get<0>();
            return a.get<1>() < b.get<1>();
        }
        else
        {
            if (a.get<0>() == b.get<0>())
                return a.get<1>() < b.get<1>();
            return a.get<0>() < b.get<0>();
        }
    }
};

CBlockTemplate* CreateNewBlock(CReserveKey& reservekey)
{
    // Create new block
    auto_ptr<CBlockTemplate> pblocktemplate(new CBlockTemplate());
    if(!pblocktemplate.get())
        return NULL;
    CBlock *pblock = &pblocktemplate->block; // pointer for convenience

    // Create coinbase tx
    CTransaction txNew;
    txNew.vin.resize(1);
    txNew.vin[0].prevout.SetNull();
    txNew.vout.resize(1);
    CPubKey pubkey;
    if (!reservekey.GetReservedKey(pubkey))
        return NULL;
    txNew.vout[0].scriptPubKey << pubkey << OP_CHECKSIG;

    // Add our coinbase tx as first transaction
    pblock->vtx.push_back(txNew);
    pblocktemplate->vTxFees.push_back(-1); // updated at end
    pblocktemplate->vTxSigOps.push_back(-1); // updated at end

    // Largest block you're willing to create:
    unsigned int nBlockMaxSize = GetArg("-blockmaxsize", MAX_BLOCK_SIZE_GEN/4);
    // Limit to betweeen 1K and MAX_BLOCK_SIZE-1K for sanity:
    nBlockMaxSize = std::max((unsigned int)1000, std::min((unsigned int)(MAX_BLOCK_SIZE-1000), nBlockMaxSize));

    // How much of the block should be dedicated to high-priority transactions,
    // included regardless of the fees they pay
    unsigned int nBlockPrioritySize = GetArg("-blockprioritysize", 12000);
    nBlockPrioritySize = std::min(nBlockMaxSize, nBlockPrioritySize);

    // Minimum block size you want to create; block will be filled with free transactions
    // until there are no more or the block reaches this size:
    unsigned int nBlockMinSize = GetArg("-blockminsize", 0);
    nBlockMinSize = std::min(nBlockMaxSize, nBlockMinSize);

    // Collect memory pool transactions into the block
    int64 nFees = 0;
    {
        LOCK2(cs_main, mempool.cs);
        CBlockIndex* pindexPrev = pindexBest;
        CCoinsViewCache view(*pcoinsTip, true);

        // Priority order to process transactions
        list<COrphan> vOrphan; // list memory doesn't move
        map<uint256, vector<COrphan*> > mapDependers;
        bool fPrintPriority = GetBoolArg("-printpriority");

        // This vector will be sorted into a priority queue:
        vector<TxPriority> vecPriority;
        vecPriority.reserve(mempool.mapTx.size());
        for (map<uint256, CTransaction>::iterator mi = mempool.mapTx.begin(); mi != mempool.mapTx.end(); ++mi)
        {
            CTransaction& tx = (*mi).second;
            if (tx.IsCoinBase() || !tx.IsFinal())
                continue;

            COrphan* porphan = NULL;
            double dPriority = 0;
            int64 nTotalIn = 0;
            bool fMissingInputs = false;
            BOOST_FOREACH(const CTxIn& txin, tx.vin)
            {
                // Read prev transaction
                if (!view.HaveCoins(txin.prevout.hash))
                {
                    // This should never happen; all transactions in the memory
                    // pool should connect to either transactions in the chain
                    // or other transactions in the memory pool.
                    if (!mempool.mapTx.count(txin.prevout.hash))
                    {
                        printf("ERROR: mempool transaction missing input\n");
                        if (fDebug) assert("mempool transaction missing input" == 0);
                        fMissingInputs = true;
                        if (porphan)
                            vOrphan.pop_back();
                        break;
                    }

                    // Has to wait for dependencies
                    if (!porphan)
                    {
                        // Use list for automatic deletion
                        vOrphan.push_back(COrphan(&tx));
                        porphan = &vOrphan.back();
                    }
                    mapDependers[txin.prevout.hash].push_back(porphan);
                    porphan->setDependsOn.insert(txin.prevout.hash);
                    nTotalIn += mempool.mapTx[txin.prevout.hash].vout[txin.prevout.n].nValue;
                    continue;
                }
                const CCoins &coins = view.GetCoins(txin.prevout.hash);

                int64 nValueIn = coins.vout[txin.prevout.n].nValue;
                nTotalIn += nValueIn;

                int nConf = pindexPrev->nHeight - coins.nHeight + 1;

                dPriority += (double)nValueIn * nConf;
            }
            if (fMissingInputs) continue;

            // Priority is sum(valuein * age) / txsize
            unsigned int nTxSize = ::GetSerializeSize(tx, SER_NETWORK, PROTOCOL_VERSION);
            dPriority /= nTxSize;

            // This is a more accurate fee-per-kilobyte than is used by the client code, because the
            // client code rounds up the size to the nearest 1K. That's good, because it gives an
            // incentive to create smaller transactions.
            double dFeePerKb =  double(nTotalIn-tx.GetValueOut()) / (double(nTxSize)/1000.0);

            if (porphan)
            {
                porphan->dPriority = dPriority;
                porphan->dFeePerKb = dFeePerKb;
            }
            else
                vecPriority.push_back(TxPriority(dPriority, dFeePerKb, &(*mi).second));
        }

        // Collect transactions into block
        uint64 nBlockSize = 1000;
        uint64 nBlockTx = 0;
        int nBlockSigOps = 100;
        bool fSortedByFee = (nBlockPrioritySize <= 0);

        TxPriorityCompare comparer(fSortedByFee);
        std::make_heap(vecPriority.begin(), vecPriority.end(), comparer);

        while (!vecPriority.empty())
        {
            // Take highest priority transaction off the priority queue:
            double dPriority = vecPriority.front().get<0>();
            double dFeePerKb = vecPriority.front().get<1>();
            CTransaction& tx = *(vecPriority.front().get<2>());

            std::pop_heap(vecPriority.begin(), vecPriority.end(), comparer);
            vecPriority.pop_back();

            // Size limits
            unsigned int nTxSize = ::GetSerializeSize(tx, SER_NETWORK, PROTOCOL_VERSION);
            if (nBlockSize + nTxSize >= nBlockMaxSize)
                continue;

            // Legacy limits on sigOps:
            unsigned int nTxSigOps = tx.GetLegacySigOpCount();
            if (nBlockSigOps + nTxSigOps >= MAX_BLOCK_SIGOPS)
                continue;

            // Skip free transactions if we're past the minimum block size:
            if (fSortedByFee && (dFeePerKb < CTransaction::nMinTxFee) && (nBlockSize + nTxSize >= nBlockMinSize))
                continue;

            // Prioritize by fee once past the priority size or we run out of high-priority
            // transactions:
            if (!fSortedByFee &&
// FBX 288 blocks found a day, but priority cutoff is 144 (like bitcoin)
//                ((nBlockSize + nTxSize >= nBlockPrioritySize) || (dPriority < COIN * 576 / 250)))
                ((nBlockSize + nTxSize >= nBlockPrioritySize) || (dPriority < COIN * 144 / 250)))
            {
                fSortedByFee = true;
                comparer = TxPriorityCompare(fSortedByFee);
                std::make_heap(vecPriority.begin(), vecPriority.end(), comparer);
            }

            if (!tx.HaveInputs(view))
                continue;

            int64 nTxFees = tx.GetValueIn(view)-tx.GetValueOut();

            nTxSigOps += tx.GetP2SHSigOpCount(view);
            if (nBlockSigOps + nTxSigOps >= MAX_BLOCK_SIGOPS)
                continue;

            CValidationState state;
            if (!tx.CheckInputs(state, view, true, SCRIPT_VERIFY_P2SH))
                continue;

            CTxUndo txundo;
            uint256 hash = tx.GetHash();
            tx.UpdateCoins(state, view, txundo, pindexPrev->nHeight+1, hash);

            // Added
            pblock->vtx.push_back(tx);
            pblocktemplate->vTxFees.push_back(nTxFees);
            pblocktemplate->vTxSigOps.push_back(nTxSigOps);
            nBlockSize += nTxSize;
            ++nBlockTx;
            nBlockSigOps += nTxSigOps;
            nFees += nTxFees;

            if (fPrintPriority)
            {
                printf("priority %.1f feeperkb %.1f txid %s\n",
                       dPriority, dFeePerKb, tx.GetHash().ToString().c_str());
            }

            // Add transactions that depend on this one to the priority queue
            if (mapDependers.count(hash))
            {
                BOOST_FOREACH(COrphan* porphan, mapDependers[hash])
                {
                    if (!porphan->setDependsOn.empty())
                    {
                        porphan->setDependsOn.erase(hash);
                        if (porphan->setDependsOn.empty())
                        {
                            vecPriority.push_back(TxPriority(porphan->dPriority, porphan->dFeePerKb, porphan->ptx));
                            std::push_heap(vecPriority.begin(), vecPriority.end(), comparer);
                        }
                    }
                }
            }
        }

        nLastBlockTx = nBlockTx;
        nLastBlockSize = nBlockSize;
        printf("CreateNewBlock(): total size %"PRI64u"\n", nBlockSize);

        pblock->vtx[0].vout[0].nValue = GetBlockValue(pindexPrev->nHeight+1, nFees);
        pblocktemplate->vTxFees[0] = -nFees;

        // Fill in header
        pblock->hashPrevBlock  = pindexPrev->GetBlockHash();
        pblock->UpdateTime(pindexPrev);
        pblock->nBits          = GetNextWorkRequired(pindexPrev, pblock);
        pblock->nNonce         = 0;
        pblock->vtx[0].vin[0].scriptSig = CScript() << OP_0 << OP_0;
        pblocktemplate->vTxSigOps[0] = pblock->vtx[0].GetLegacySigOpCount();

        CBlockIndex indexDummy(*pblock);
        indexDummy.pprev = pindexPrev;
        indexDummy.nHeight = pindexPrev->nHeight + 1;
        CCoinsViewCache viewNew(*pcoinsTip, true);
        CValidationState state;
        if (!pblock->ConnectBlock(state, &indexDummy, viewNew, true))
            throw std::runtime_error("CreateNewBlock() : ConnectBlock failed");
    }

    return pblocktemplate.release();
}


void IncrementExtraNonce(CBlock* pblock, CBlockIndex* pindexPrev, unsigned int& nExtraNonce)
{
    // Update nExtraNonce
    static uint256 hashPrevBlock;
    if (hashPrevBlock != pblock->hashPrevBlock)
    {
        nExtraNonce = 0;
        hashPrevBlock = pblock->hashPrevBlock;
    }
    ++nExtraNonce;
    unsigned int nHeight = pindexPrev->nHeight+1; // Height first in coinbase required for block.version=2
    pblock->vtx[0].vin[0].scriptSig = (CScript() << nHeight << CBigNum(nExtraNonce)) + COINBASE_FLAGS;
    assert(pblock->vtx[0].vin[0].scriptSig.size() <= 100);

    pblock->hashMerkleRoot = pblock->BuildMerkleTree();
}


void FormatHashBuffers(CBlock* pblock, char* pmidstate, char* pdata, char* phash1)
{
    //
    // Pre-build hash buffers
    //
    struct
    {
        struct unnamed2
        {
            int nVersion;
            uint256 hashPrevBlock;
            uint256 hashMerkleRoot;
            unsigned int nTime;
            unsigned int nBits;
            unsigned int nNonce;
        }
        block;
        unsigned char pchPadding0[64];
        uint256 hash1;
        unsigned char pchPadding1[64];
    }
    tmp;
    memset(&tmp, 0, sizeof(tmp));

    tmp.block.nVersion       = pblock->nVersion;
    tmp.block.hashPrevBlock  = pblock->hashPrevBlock;
    tmp.block.hashMerkleRoot = pblock->hashMerkleRoot;
    tmp.block.nTime          = pblock->nTime;
    tmp.block.nBits          = pblock->nBits;
    tmp.block.nNonce         = pblock->nNonce;

    FormatHashBlocks(&tmp.block, sizeof(tmp.block));
    FormatHashBlocks(&tmp.hash1, sizeof(tmp.hash1));

    // Byte swap all the input buffer
    for (unsigned int i = 0; i < sizeof(tmp)/4; i++)
        ((unsigned int*)&tmp)[i] = ByteReverse(((unsigned int*)&tmp)[i]);

    // Precalc the first half of the first hash, which stays constant
    SHA256Transform(pmidstate, &tmp.block, pSHA256InitState);

    memcpy(pdata, &tmp.block, 128);
    memcpy(phash1, &tmp.hash1, 64);
}


bool CheckWork(CBlock* pblock, CWallet& wallet, CReserveKey& reservekey)
{
    uint256 hash = pblock->GetPoWHash();
    uint256 hashTarget = CBigNum().SetCompact(pblock->nBits).getuint256();

    if (hash > hashTarget)
        return false;

    //// debug print
    printf("Litecoin RPCMiner:\n");
    printf("proof-of-work found  \n  hash: %s  \ntarget: %s\n", hash.GetHex().c_str(), hashTarget.GetHex().c_str());
    pblock->print();
    printf("generated %s\n", FormatMoney(pblock->vtx[0].vout[0].nValue).c_str());

    // Found a solution
    {
        LOCK(cs_main);
        if (pblock->hashPrevBlock != hashBestChain)
            return error("LitecoinMiner : generated block is stale");

        // Remove key from key pool
        reservekey.KeepKey();

        // Track how many getdata requests this block gets
        {
            LOCK(wallet.cs_wallet);
            wallet.mapRequestCount[pblock->GetHash()] = 0;
        }

        // Process this block the same as if we had received it from another node
        CValidationState state;
        if (!ProcessBlock(state, NULL, pblock))
            return error("LitecoinMiner : ProcessBlock, block not accepted");
    }

    return true;
}

// Amount compression:
// * If the amount is 0, output 0
// * first, divide the amount (in base units) by the largest power of 10 possible; call the exponent e (e is max 9)
// * if e<9, the last digit of the resulting number cannot be 0; store it as d, and drop it (divide by 10)
//   * call the result n
//   * output 1 + 10*(9*n + d - 1) + e
// * if e==9, we only know the resulting number is not zero, so output 1 + 10*(n - 1) + 9
// (this is decodable, as d is in [1-9] and e is in [0-9])

uint64 CTxOutCompressor::CompressAmount(uint64 n)
{
    if (n == 0)
        return 0;
    int e = 0;
    while (((n % 10) == 0) && e < 9) {
        n /= 10;
        e++;
    }
    if (e < 9) {
        int d = (n % 10);
        assert(d >= 1 && d <= 9);
        n /= 10;
        return 1 + (n*9 + d - 1)*10 + e;
    } else {
        return 1 + (n - 1)*10 + 9;
    }
}

uint64 CTxOutCompressor::DecompressAmount(uint64 x)
{
    // x = 0  OR  x = 1+10*(9*n + d - 1) + e  OR  x = 1+10*(n - 1) + 9
    if (x == 0)
        return 0;
    x--;
    // x = 10*(9*n + d - 1) + e
    int e = x % 10;
    x /= 10;
    uint64 n = 0;
    if (e < 9) {
        // x = 9*n + d - 1
        int d = (x % 9) + 1;
        x /= 9;
        // x = n
        n = x*10 + d;
    } else {
        n = x+1;
    }
    while (e) {
        n *= 10;
        e--;
    }
    return n;
}


class CMainCleanup
{
public:
    CMainCleanup() {}
    ~CMainCleanup() {
        // block headers
        std::map<uint256, CBlockIndex*>::iterator it1 = mapBlockIndex.begin();
        for (; it1 != mapBlockIndex.end(); it1++)
            delete (*it1).second;
        mapBlockIndex.clear();

        // orphan blocks
        std::map<uint256, CBlock*>::iterator it2 = mapOrphanBlocks.begin();
        for (; it2 != mapOrphanBlocks.end(); it2++)
            delete (*it2).second;
        mapOrphanBlocks.clear();

        // orphan transactions
        mapOrphanTransactions.clear();
    }
} instance_of_cmaincleanup;
