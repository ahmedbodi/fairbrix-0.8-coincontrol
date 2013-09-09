// Copyright (c) 2010 Satoshi Nakamoto
// Copyright (c) 2009-2012 The Bitcoin developers
// Distributed under the MIT/X11 software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#include <boost/assign/list_of.hpp>

#include "base58.h"
#include "bitcoinrpc.h"
#include "db.h"
#include "init.h"
#include "main.h"
#include "net.h"
#include "wallet.h"

using namespace std;
using namespace boost;
using namespace boost::assign;
using namespace json_spirit;

//
// Utilities: convert hex-encoded Values
// (throws error if not hex).
//
uint256 ParseHashV(const Value& v, string strName)
{
    string strHex;
    if (v.type() == str_type)
        strHex = v.get_str();
    if (!IsHex(strHex)) // Note: IsHex("") is false
        throw JSONRPCError(RPC_INVALID_PARAMETER, strName+" must be hexadecimal string (not '"+strHex+"')");
    uint256 result;
    result.SetHex(strHex);
    return result;
}
uint256 ParseHashO(const Object& o, string strKey)
{
    return ParseHashV(find_value(o, strKey), strKey);
}
vector<unsigned char> ParseHexV(const Value& v, string strName)
{
    string strHex;
    if (v.type() == str_type)
        strHex = v.get_str();
    if (!IsHex(strHex))
        throw JSONRPCError(RPC_INVALID_PARAMETER, strName+" must be hexadecimal string (not '"+strHex+"')");
    return ParseHex(strHex);
}
vector<unsigned char> ParseHexO(const Object& o, string strKey)
{
    return ParseHexV(find_value(o, strKey), strKey);
}

void ScriptPubKeyToJSON(const CScript& scriptPubKey, Object& out)
{
    txnouttype type;
    vector<CTxDestination> addresses;
    int nRequired;

    out.push_back(Pair("asm", scriptPubKey.ToString()));
    out.push_back(Pair("hex", HexStr(scriptPubKey.begin(), scriptPubKey.end())));

    if (!ExtractDestinations(scriptPubKey, type, addresses, nRequired))
    {
        out.push_back(Pair("type", GetTxnOutputType(TX_NONSTANDARD)));
        return;
    }

    out.push_back(Pair("reqSigs", nRequired));
    out.push_back(Pair("type", GetTxnOutputType(type)));

    Array a;
    BOOST_FOREACH(const CTxDestination& addr, addresses)
        a.push_back(CBitcoinAddress(addr).ToString());
    out.push_back(Pair("addresses", a));
}

void TxToJSON(const CTransaction& tx, const uint256 hashBlock, Object& entry)
{
    entry.push_back(Pair("txid", tx.GetHash().GetHex()));
    entry.push_back(Pair("version", tx.nVersion));
    entry.push_back(Pair("locktime", (boost::int64_t)tx.nLockTime));
    Array vin;
    BOOST_FOREACH(const CTxIn& txin, tx.vin)
    {
        Object in;
        if (tx.IsCoinBase())
            in.push_back(Pair("coinbase", HexStr(txin.scriptSig.begin(), txin.scriptSig.end())));
        else
        {
            in.push_back(Pair("txid", txin.prevout.hash.GetHex()));
            in.push_back(Pair("vout", (boost::int64_t)txin.prevout.n));
            Object o;
            o.push_back(Pair("asm", txin.scriptSig.ToString()));
            o.push_back(Pair("hex", HexStr(txin.scriptSig.begin(), txin.scriptSig.end())));
            in.push_back(Pair("scriptSig", o));
        }
        in.push_back(Pair("sequence", (boost::int64_t)txin.nSequence));
        vin.push_back(in);
    }
    entry.push_back(Pair("vin", vin));
    Array vout;
    for (unsigned int i = 0; i < tx.vout.size(); i++)
    {
        const CTxOut& txout = tx.vout[i];
        Object out;
        out.push_back(Pair("value", ValueFromAmount(txout.nValue)));
        out.push_back(Pair("n", (boost::int64_t)i));
        Object o;
        ScriptPubKeyToJSON(txout.scriptPubKey, o);
        out.push_back(Pair("scriptPubKey", o));
        vout.push_back(out);
    }
    entry.push_back(Pair("vout", vout));

    if (hashBlock != 0)
    {
        entry.push_back(Pair("blockhash", hashBlock.GetHex()));
        map<uint256, CBlockIndex*>::iterator mi = mapBlockIndex.find(hashBlock);
        if (mi != mapBlockIndex.end() && (*mi).second)
        {
            CBlockIndex* pindex = (*mi).second;
            if (pindex->IsInMainChain())
            {
                entry.push_back(Pair("confirmations", 1 + nBestHeight - pindex->nHeight));
                entry.push_back(Pair("time", (boost::int64_t)pindex->nTime));
                entry.push_back(Pair("blocktime", (boost::int64_t)pindex->nTime));
            }
            else
                entry.push_back(Pair("confirmations", 0));
        }
    }
}

Value getrawtransaction(const Array& params, bool fHelp)
{
    if (fHelp || params.size() < 1 || params.size() > 2)
        throw runtime_error(
            "getrawtransaction <txid> [verbose=0]\n"
            "If verbose=0, returns a string that is\n"
            "serialized, hex-encoded data for <txid>.\n"
            "If verbose is non-zero, returns an Object\n"
            "with information about <txid>.");

    uint256 hash = ParseHashV(params[0], "parameter 1");

    bool fVerbose = false;
    if (params.size() > 1)
        fVerbose = (params[1].get_int() != 0);

    CTransaction tx;
    uint256 hashBlock = 0;
    if (!GetTransaction(hash, tx, hashBlock, true))
        throw JSONRPCError(RPC_INVALID_ADDRESS_OR_KEY, "No information available about transaction");

    CDataStream ssTx(SER_NETWORK, PROTOCOL_VERSION);
    ssTx << tx;
    string strHex = HexStr(ssTx.begin(), ssTx.end());

    if (!fVerbose)
        return strHex;

    Object result;
    result.push_back(Pair("hex", strHex));
    TxToJSON(tx, hashBlock, result);
    return result;
}

// FBX proof of stake voting test
int64 posvValueDiff;    // aggregate money inflow and outflow for all matching addresses
int64 posvTxValue;      // output of one tx
int posvSig;            // inflow or outflow
int posvVerbose;
int posvFastmode;
int posvWarningCount;   // more than 1 address per output
int posvErrCount;       // more than 1 address per output, vanitygen addresses used for oracle are affected

// oracle related vars start here (not used in older test functions svdebugblock and svscanblocks)
bool posv_fOracle;      // poll the oracle

int posv_nMaturity;
#define POSV_MATURITY_MAX 2
int64 posv_nMaturityBlockStart[POSV_MATURITY_MAX] = {172000, 180000}; // TODO: use block time
int64 posv_nMaturityBlockEnd[POSV_MATURITY_MAX] = {180000, 188000};
std::string posv_strMaturityVanity[POSV_MATURITY_MAX] = {"TV138", "TV139"};
std::string posv_strMaturityDesc[POSV_MATURITY_MAX] = {"20130831 22:00 EST", "20130930 22:00 EST"};

int posv_nPair;
#define POSV_PAIR_MAX 4
#define POSV_CHOICE_MAX 10                                          // IMPORTANT: this can't be >10
int64 posv_nOracleValuePerChoice[POSV_PAIR_MAX][POSV_CHOICE_MAX];   // voting results
int64 posv_OraclePairDivisor[POSV_PAIR_MAX] = {1, 10, 100, 1000};   // to get the relevant digit

// TODO: these should be p2p-listed via a 2nd set of votings, not hardcoded
string posv_strOraclePairs[POSV_PAIR_MAX] = {"BTC/USD", "BTC/XRP", "LTC/BTC", "FBX/BTC"};
std::string posv_strOracleResult[POSV_PAIR_MAX][POSV_CHOICE_MAX] =
         {{"invalid", "don't know", ">250", ">200", ">150", ">125", "<=125", "<=100", "<=75", "<=50"},
          {"invalid", "don't know", ">70000", ">40000", ">20000", ">10000", "<=10000", "<=5000", "<=3000", "<=2000"},
          {"invalid", "don't know", ">0.04", ">0.03", ">0.025", ">0.02", "<=0.02", "<=0.015", "<=0.01", "<=0.0075"},
          {"invalid", "don't know", ">0.00003", ">0.00002", ">0.00001", ">0.000005", "<=0.000005", "<=0.000003", "<=0.000002", "<=0.000001"}};
std::string posv_strOracleImplied[POSV_PAIR_MAX][POSV_CHOICE_MAX] =
         {{"invalid", "or don't care", "n/a", "<=250", "<=200", "<=150", ">100", ">75", ">50", "n/a"},
          {"invalid", "or don't care", "n/a", "<=70000", "<=40000", "<=20000", ">5000", ">3000", ">2000", "n/a"},
          {"invalid", "or don't care", "n/a", "<=0.04", "<=0.03", "<=0.025", ">0.015", ">0.01", ">0.0075", "n/a"},
          {"invalid", "or don't care", "n/a", "<=0.00003", "<=0.00002", "<=0.00001", ">0.000003", ">0.000002", ">0.000001", "n/a"}};

std::string strPosvTest;
void posvScriptPubKeyToJSON(const CScript& scriptPubKey, Object& out)
{
    txnouttype type;
    vector<CTxDestination> addresses;
    int nRequired;

    int verbose = (posvFastmode ? 0 : posvVerbose);

    if (verbose)
    {
        out.push_back(Pair("asm", scriptPubKey.ToString()));
        out.push_back(Pair("hex", HexStr(scriptPubKey.begin(), scriptPubKey.end())));
    }

// TODO: extra warning (could be a exploit)
    if (!ExtractDestinations(scriptPubKey, type, addresses, nRequired))
    {
        out.push_back(Pair("type", GetTxnOutputType(TX_NONSTANDARD)));
        return;
    }

    if (verbose)
    {
        out.push_back(Pair("reqSigs", nRequired));
        out.push_back(Pair("type", GetTxnOutputType(type)));
    }

    Array a;
    bool fmatch = false;
    int n = 0;
    BOOST_FOREACH(const CTxDestination& addr, addresses)
    {
        a.push_back(CBitcoinAddress(addr).ToString());
        string s = CBitcoinAddress(addr).ToString();

        fmatch = true;
        n++;
        unsigned int l = strPosvTest.length();
        for (unsigned int i = 0; i < l; i++)
        {
            if (s[i+1] != strPosvTest[i])
            {
                fmatch = false;
                break;
            }
        }
        if (fmatch)
        {
            posvValueDiff += (posvTxValue * posvSig);
            if (posv_fOracle)
            {
                int64 v = (posvTxValue / posv_OraclePairDivisor[posv_nPair]) % 10;
                posv_nOracleValuePerChoice[posv_nPair][v] += (posvTxValue * posvSig);
// todo: extra warning if POSV_CHOICE_MAX!=10 ?
            }
        }
    }

    if(!posvFastmode)
    {
        if (fmatch)
            out.push_back(Pair("addresses (match detected)", a));
        else
            out.push_back(Pair("addresses", a));
    }

    // notice outputs with many addresses
    if (n > 1)
    {
        posvWarningCount++;
        if (fmatch) posvErrCount++;
    }
}
Object posvblockToJSON(const CBlock& block, const CBlockIndex* blockindex)
{
    int verbose = (posvFastmode ? 0 : posvVerbose);

    Object result;
    if (verbose) result.push_back(Pair("hash", block.GetHash().GetHex()));
    CMerkleTx txGen(block.vtx[0]);
    txGen.SetMerkleBranch(&block);
    if (verbose)
    {
        result.push_back(Pair("confirmations", (int)txGen.GetDepthInMainChain()));
        result.push_back(Pair("size", (int)::GetSerializeSize(block, SER_NETWORK, PROTOCOL_VERSION)));
        result.push_back(Pair("height", blockindex->nHeight));
        result.push_back(Pair("version", block.nVersion));
        result.push_back(Pair("merkleroot", block.hashMerkleRoot.GetHex()));
    }

    int i1 = 0;
    BOOST_FOREACH(const CTransaction&tx, block.vtx)
    {
        i1++;
        if (!posvFastmode)
        {
        result.push_back(Pair("transaction#", i1));
        result.push_back(Pair("tx hash", tx.GetHash().GetHex()));
        }
        else
        {
            // we probably don't need to do this
            tx.GetHash().GetHex();
        }

        // code from TxToJSON
        BOOST_FOREACH(const CTxIn& txin, tx.vin)
        {
            if (tx.IsCoinBase())
                // always do this (error 'No information available about transaction' if not)
                if (!posvFastmode)
                    result.push_back(Pair("coinbase", HexStr(txin.scriptSig.begin(), txin.scriptSig.end())));
                else
                    HexStr(txin.scriptSig.begin(), txin.scriptSig.end());
            else
            {
                if (!posvFastmode)
                    result.push_back(Pair("input from txid", txin.prevout.hash.GetHex()));
                // we probably don't need to do this
                else
                    txin.prevout.hash.GetHex();

                int64 vout2 = (boost::int64_t)txin.prevout.n;
                if (!posvFastmode)
                    result.push_back(Pair("begin processing previous tx, vout", vout2));

                // backtrace transactions if not mined (need to know where the coins coming from)
                CTransaction tx2;
                uint256 hashBlock2 = 0;
                if (!GetTransaction(txin.prevout.hash, tx2, hashBlock2, true))
                    throw JSONRPCError(RPC_INVALID_ADDRESS_OR_KEY, "No information available about transaction");

                if (!posvFastmode)
                {
                BOOST_FOREACH(const CTxIn& txin2, tx2.vin)
                {
                        if (tx2.IsCoinBase())
                            result.push_back(Pair("previous tx: coinbase", HexStr(txin2.scriptSig.begin(), txin2.scriptSig.end())));
                        else
                            result.push_back(Pair("previous tx: input from txid", txin2.prevout.hash.GetHex()));
                }
                }

                // check outputs of each 'previous' transaction (to count coins sent *from* special addresses)
                posvTxValue = 0;
                posvSig = -1;
                for (unsigned int i = 0; i < tx2.vout.size(); i++)
                {
                    const CTxOut& txout4 = tx2.vout[i];
                    if (!posvFastmode)
                        result.push_back(Pair("previous tx output#", (int) i));
                    posvTxValue = txout4.nValue;
                    if (!posvFastmode)
                        result.push_back(Pair("value", ValueFromAmount(txout4.nValue)));

                    if (i != vout2) continue;

                    Object o4;
                    posvScriptPubKeyToJSON(txout4.scriptPubKey, o4);
                    if (!posvFastmode)
                        result.push_back(Pair("scriptPubKey previous tx", o4));
                }
                // end processing 'previous' tx

                if (verbose)
                {
                    Object o;
                    o.push_back(Pair("asm", txin.scriptSig.ToString()));
                    o.push_back(Pair("hex", HexStr(txin.scriptSig.begin(), txin.scriptSig.end())));
                    result.push_back(Pair("scriptSig", o));
                }
            }
            if (verbose) result.push_back(Pair("sequence", (boost::int64_t)txin.nSequence));
        }

        // check outputs of each transaction (to count coins sent *to* special addresses)
        posvTxValue = 0;
        posvSig = 1;
        for (unsigned int i = 0; i < tx.vout.size(); i++)
        {
            const CTxOut& txout = tx.vout[i];
            if (!posvFastmode)
                result.push_back(Pair("output#", (int) i));
            posvTxValue = txout.nValue;
            if (!posvFastmode)
                result.push_back(Pair("value", ValueFromAmount(txout.nValue)));
            Object o;
            posvScriptPubKeyToJSON(txout.scriptPubKey, o);
            if (!posvFastmode)
                result.push_back(Pair("scriptPubKey", o));
        }
    }

    if (verbose)
    {
        result.push_back(Pair("time", (boost::int64_t)block.GetBlockTime()));
        result.push_back(Pair("nonce", (boost::uint64_t)block.nNonce));
        result.push_back(Pair("bits", HexBits(block.nBits)));
        result.push_back(Pair("difficulty", GetDifficulty(blockindex)));

        if (blockindex->pprev)
            result.push_back(Pair("previousblockhash", blockindex->pprev->GetBlockHash().GetHex()));
        if (blockindex->pnext)
            result.push_back(Pair("nextblockhash", blockindex->pnext->GetBlockHash().GetHex()));
    }

    return result;
}
Value svdebugblock(const Array& params, bool fHelp)
{
    if (fHelp || params.size() != 3)
        throw runtime_error(
            "svdebugblock <index> <verbose 0|1> <vanity string>\n"
            "Scans block <index> and returns info about balance change\n"
            "of addresses containing <vanity string>.");

    posvVerbose = params[1].get_int();
    posvFastmode = 0;
    posvValueDiff = 0;
    posvWarningCount = 0;
    posvErrCount = 0;
    posv_fOracle = false;
    posv_nPair = posv_nMaturity = 0;

    int nHeight = params[0].get_int();
    if (nHeight < 0 || nHeight > nBestHeight)
        throw runtime_error("Block number out of range.");

    strPosvTest = params[2].get_str();
    if (strPosvTest.length() > 10)
        throw runtime_error("vanity string too long");

    CBlockIndex* pblockindex = FindBlockByHeight(nHeight);
    CBlock block;
    block.ReadFromDisk(pblockindex);

    Object result = posvblockToJSON(block, pblockindex);

    result.push_back(Pair("total balance change", (double)posvValueDiff/(double)COIN));
    if (posvWarningCount) result.push_back(Pair("outputs with more than 1 addr", posvWarningCount));
    if (posvErrCount) result.push_back(Pair("matching outputs with more than 1 addr", posvErrCount));

    return result;
}
Value svscanblocks(const Array& params, bool fHelp)
{
    if (fHelp || params.size() != 3)
        throw runtime_error(
            "svscanblocks <index start> <index end> <vanity string>\n"
            "Scans all blocks between <index start> and <index end>.\n"
            "Returns info about aggregate balance change of addresses\n"
            "containing <vanity string>.");

    posvVerbose = 0;
    posvFastmode = 1;
    posvValueDiff = 0;
    posvWarningCount = 0;
    posvErrCount = 0;
    posv_fOracle = false;
    posv_nPair = posv_nMaturity = 0;

    int nHeight0 = params[0].get_int();
    int nHeight1 = params[1].get_int();
    if (nHeight1 < nHeight0 || nHeight0 < 0 || nHeight1 > nBestHeight)
        throw runtime_error("Block number(s) out of range.");

    strPosvTest = params[2].get_str();
    if (strPosvTest.length() > 10)
        throw runtime_error("vanity string too long");

    for (int nHeight = nHeight0; nHeight <= nHeight1; nHeight++)
    {
        CBlockIndex* pblockindex = FindBlockByHeight(nHeight);
        CBlock block;
        block.ReadFromDisk(pblockindex);

        posvblockToJSON(block, pblockindex);
    }

    Object result2;
    result2.push_back(Pair("total balance change", (double)posvValueDiff/(double)COIN));
    if (posvWarningCount) result2.push_back(Pair("outputs with more than 1 addr", posvWarningCount));
    if (posvErrCount) result2.push_back(Pair("matching outputs with more than 1 addr", posvErrCount));

    return result2;
}
Value svgetoracle(const Array& params, bool fHelp)
{
    if (fHelp || params.size() < 1 || params.size() > 2)
        throw runtime_error(
            "svgetoracle <pair index> [voting index=0]\n"
            "(use 'svlistpairs' to list available currency pairs,\n"
            " and 'svlistvotings' to list available votings/maturities)");

    posvVerbose = 0;
    posvFastmode = 1;
    posvValueDiff = 0;
    posvWarningCount = 0;
    posvErrCount = 0;
    posv_fOracle = true;

    // clear all votes
    for (unsigned int i = 0; i < POSV_PAIR_MAX; i++)
        for (unsigned int j = 0; j < 10; j++)
            posv_nOracleValuePerChoice[i][j] = 0;

    posv_nPair = params[0].get_int();
    if ((posv_nPair < 0) || (posv_nPair >= POSV_PAIR_MAX))
        throw runtime_error("Pair index out of range.");

    if (params.size() >= 2)
        posv_nMaturity = params[1].get_int();
    else
        posv_nMaturity = 0;
    if ((posv_nMaturity < 0) || (posv_nMaturity >= POSV_MATURITY_MAX))
        throw runtime_error("Maturity number out of range.");

    strPosvTest = posv_strMaturityVanity[posv_nMaturity];
    int nHeight0 = posv_nMaturityBlockStart[posv_nMaturity];
    int nHeight1 = posv_nMaturityBlockEnd[posv_nMaturity];

    // voting in progress
    if (nBestHeight < nHeight1)
        nHeight1 = nBestHeight;

    if (nHeight1 < nHeight0 || nHeight0 < 0 || nHeight1 > nBestHeight)
        throw runtime_error("Block number(s) out of range.");
    if (strPosvTest.length() > 10)
        throw runtime_error("vanity string too long");

    for (int nHeight = nHeight0; nHeight <= nHeight1; nHeight++)
    {
        CBlockIndex* pblockindex = FindBlockByHeight(nHeight);
        CBlock block;
        block.ReadFromDisk(pblockindex);

        posvblockToJSON(block, pblockindex);
    }

    Object result2;
    result2.push_back(Pair("total voting coins", (double)posvValueDiff/(double)COIN));
    if (posvWarningCount) result2.push_back(Pair("outputs with more than 1 addr", posvWarningCount));
    if (posvErrCount) result2.push_back(Pair("matching outputs with more than 1 addr", posvErrCount));

    result2.push_back(Pair(posv_strOraclePairs[posv_nPair], posv_strMaturityDesc[posv_nMaturity]));

    // skip 'invalid' choice (i==0)
    for (unsigned int i = 1; i < POSV_CHOICE_MAX; i++)
    {
        int64 v = posv_nOracleValuePerChoice[posv_nPair][i];
        if (v > 0)
        {
            result2.push_back(Pair(" ", " "));
            result2.push_back(Pair("result", posv_strOracleResult[posv_nPair][i]));
            result2.push_back(Pair("implied", posv_strOracleImplied[posv_nPair][i]));
            result2.push_back(Pair("voting coins", (double)v/(double)COIN));
        }
    }

    return result2;
}
Value svlistvotings(const Array& params, bool fHelp)
{
    if (fHelp || params.size() > 1)
        throw runtime_error(
            "svlistvotings [verbose=0]\n"
            "Lists all votings/maturity dates and 'special addresses'\n"
            "available for the proof of stake voting system.\n"
            "(verbose=0 to list only votings currently in progress)");

    bool fVerbose = false;
    if (params.size() > 0)
        fVerbose = (params[0].get_int() != 0);

    Object result2;
    result2.push_back(Pair("list of available maturities, total count", int(POSV_MATURITY_MAX)));

    for (unsigned int i = 0; i < POSV_MATURITY_MAX; i++)
    {
        result2.push_back(Pair(" ", " "));
        result2.push_back(Pair("index", int(i)));
        result2.push_back(Pair("date and time", posv_strMaturityDesc[i]));

        int state = 0;
        if (nBestHeight < posv_nMaturityBlockStart[i])
        {
            result2.push_back(Pair("status", "voting not yet started"));
        }
        else if (nBestHeight < posv_nMaturityBlockEnd[i])
        {
            result2.push_back(Pair("status", "voting in progress"));
            state = 1;
        }
        else
        {
            result2.push_back(Pair("status", "voting closed"));
            state = 2;
        }
        if (!fVerbose && state != 1) continue;

        if (state == 1)
            result2.push_back(Pair("basics", "vote by sending coins to an address generated with vanitygen"));
        string s = "vanitygen -X95 f" + posv_strMaturityVanity[i];
        result2.push_back(Pair("vanitygen command", s));
        result2.push_back(Pair("voting start (block height)", posv_nMaturityBlockStart[i]));
        result2.push_back(Pair("voting start (block time)", "n/a"));
        result2.push_back(Pair("voting end (block height)", posv_nMaturityBlockEnd[i]));
        result2.push_back(Pair("voting end (block time)", "n/a"));
    }

    return result2;
}
Value svlistpairs(const Array& params, bool fHelp)
{
    if (fHelp || params.size() > 1)
        throw runtime_error(
            "svlistpairs [verbose=1]\n"
            "Lists all currency pairs available\n"
            "for the proof of stake voting system,\n"
            "and voting instructions for each pair.");

    bool fVerbose = true;
    if (params.size() > 0)
        fVerbose = (params[0].get_int() != 0);

    Object result2;
    result2.push_back(Pair("list of available pairs, total count", int(POSV_PAIR_MAX)));
    if (fVerbose)
        result2.push_back(Pair("basics", "least significant digits of amount sent are used for voting"));

    for (unsigned int i = 0; i < POSV_PAIR_MAX; i++)
    {
        result2.push_back(Pair(" ", " "));
        result2.push_back(Pair("index", int(i)));
        result2.push_back(Pair("name", posv_strOraclePairs[i]));

        if (fVerbose)
        {
        if (i == 0)
            result2.push_back(Pair("digit for voting", "8th after the decimal point (single satoshis)"));
        else if (i == 1)
            result2.push_back(Pair("digit for voting", "7th after the decimal point (satoshis * 10)"));
        else if (i == 2)
            result2.push_back(Pair("digit for voting", "6th after the decimal point (satoshis * 100)"));
        else if (i == 3)
            result2.push_back(Pair("digit for voting", "5th after the decimal point (satoshis * 1000)"));
        else if (i == 4)
            result2.push_back(Pair("digit for voting", "4th after the decimal point (satoshis * 10000)"));
        else if (i == 5)
            result2.push_back(Pair("digit for voting", "3rd after the decimal point (1/1000 coins)"));
        else if (i == 6)
            result2.push_back(Pair("digit for voting", "2nd after the decimal point (1/100 coins)"));
        else if (i == 7)
            result2.push_back(Pair("digit for voting", "1st after the decimal point (1/10 coins)"));
        }

        result2.push_back(Pair("available choices", "1..9"));
        for (unsigned int i2 = 1; i2 < POSV_CHOICE_MAX; i2++)
        {
            string s = posv_strOracleResult[i][i2];
            if (i2 == 2)
                s = s + ", no upper limit";
            else if (i2 == 9)
                s = s + ", no lower limit";
            else
                s = s + ", " + posv_strOracleImplied[i][i2];
            result2.push_back(Pair(s, int(i2)));
        }
    }
    return result2;
}

Value listunspent(const Array& params, bool fHelp)
{
    if (fHelp || params.size() > 3)
        throw runtime_error(
            "listunspent [minconf=1] [maxconf=9999999]  [\"address\",...]\n"
            "Returns array of unspent transaction outputs\n"
            "with between minconf and maxconf (inclusive) confirmations.\n"
            "Optionally filtered to only include txouts paid to specified addresses.\n"
            "Results are an array of Objects, each of which has:\n"
            "{txid, vout, scriptPubKey, amount, confirmations}");

    RPCTypeCheck(params, list_of(int_type)(int_type)(array_type));

    int nMinDepth = 1;
    if (params.size() > 0)
        nMinDepth = params[0].get_int();

    int nMaxDepth = 9999999;
    if (params.size() > 1)
        nMaxDepth = params[1].get_int();

    set<CBitcoinAddress> setAddress;
    if (params.size() > 2)
    {
        Array inputs = params[2].get_array();
        BOOST_FOREACH(Value& input, inputs)
        {
            CBitcoinAddress address(input.get_str());
            if (!address.IsValid())
                throw JSONRPCError(RPC_INVALID_ADDRESS_OR_KEY, string("Invalid Fairbrix address: ")+input.get_str());
            if (setAddress.count(address))
                throw JSONRPCError(RPC_INVALID_PARAMETER, string("Invalid parameter, duplicated address: ")+input.get_str());
           setAddress.insert(address);
        }
    }

    Array results;
    vector<COutput> vecOutputs;
    pwalletMain->AvailableCoins(vecOutputs, false);
    BOOST_FOREACH(const COutput& out, vecOutputs)
    {
        if (out.nDepth < nMinDepth || out.nDepth > nMaxDepth)
            continue;

        if (setAddress.size())
        {
            CTxDestination address;
            if (!ExtractDestination(out.tx->vout[out.i].scriptPubKey, address))
                continue;

            if (!setAddress.count(address))
                continue;
        }

        int64 nValue = out.tx->vout[out.i].nValue;
        const CScript& pk = out.tx->vout[out.i].scriptPubKey;
        Object entry;
        entry.push_back(Pair("txid", out.tx->GetHash().GetHex()));
        entry.push_back(Pair("vout", out.i));
        CTxDestination address;
        if (ExtractDestination(out.tx->vout[out.i].scriptPubKey, address))
        {
            entry.push_back(Pair("address", CBitcoinAddress(address).ToString()));
            if (pwalletMain->mapAddressBook.count(address))
                entry.push_back(Pair("account", pwalletMain->mapAddressBook[address]));
        }
        entry.push_back(Pair("scriptPubKey", HexStr(pk.begin(), pk.end())));
        if (pk.IsPayToScriptHash())
        {
            CTxDestination address;
            if (ExtractDestination(pk, address))
            {
                const CScriptID& hash = boost::get<const CScriptID&>(address);
                CScript redeemScript;
                if (pwalletMain->GetCScript(hash, redeemScript))
                    entry.push_back(Pair("redeemScript", HexStr(redeemScript.begin(), redeemScript.end())));
            }
        }
        entry.push_back(Pair("amount",ValueFromAmount(nValue)));
        entry.push_back(Pair("confirmations",out.nDepth));
        results.push_back(entry);
    }

    return results;
}

Value createrawtransaction(const Array& params, bool fHelp)
{
    if (fHelp || params.size() != 2)
        throw runtime_error(
            "createrawtransaction [{\"txid\":txid,\"vout\":n},...] {address:amount,...}\n"
            "Create a transaction spending given inputs\n"
            "(array of objects containing transaction id and output number),\n"
            "sending to given address(es).\n"
            "Returns hex-encoded raw transaction.\n"
            "Note that the transaction's inputs are not signed, and\n"
            "it is not stored in the wallet or transmitted to the network.");

    RPCTypeCheck(params, list_of(array_type)(obj_type));

    Array inputs = params[0].get_array();
    Object sendTo = params[1].get_obj();

    CTransaction rawTx;

    BOOST_FOREACH(const Value& input, inputs)
    {
        const Object& o = input.get_obj();

        uint256 txid = ParseHashO(o, "txid");

        const Value& vout_v = find_value(o, "vout");
        if (vout_v.type() != int_type)
            throw JSONRPCError(RPC_INVALID_PARAMETER, "Invalid parameter, missing vout key");
        int nOutput = vout_v.get_int();
        if (nOutput < 0)
            throw JSONRPCError(RPC_INVALID_PARAMETER, "Invalid parameter, vout must be positive");

        CTxIn in(COutPoint(txid, nOutput));
        rawTx.vin.push_back(in);
    }

    set<CBitcoinAddress> setAddress;
    BOOST_FOREACH(const Pair& s, sendTo)
    {
        CBitcoinAddress address(s.name_);
        if (!address.IsValid())
            throw JSONRPCError(RPC_INVALID_ADDRESS_OR_KEY, string("Invalid Fairbrix address: ")+s.name_);

        if (setAddress.count(address))
            throw JSONRPCError(RPC_INVALID_PARAMETER, string("Invalid parameter, duplicated address: ")+s.name_);
        setAddress.insert(address);

        CScript scriptPubKey;
        scriptPubKey.SetDestination(address.Get());
        int64 nAmount = AmountFromValue(s.value_);

        CTxOut out(nAmount, scriptPubKey);
        rawTx.vout.push_back(out);
    }

    CDataStream ss(SER_NETWORK, PROTOCOL_VERSION);
    ss << rawTx;
    return HexStr(ss.begin(), ss.end());
}

Value decoderawtransaction(const Array& params, bool fHelp)
{
    if (fHelp || params.size() != 1)
        throw runtime_error(
            "decoderawtransaction <hex string>\n"
            "Return a JSON object representing the serialized, hex-encoded transaction.");

    vector<unsigned char> txData(ParseHexV(params[0], "argument"));
    CDataStream ssData(txData, SER_NETWORK, PROTOCOL_VERSION);
    CTransaction tx;
    try {
        ssData >> tx;
    }
    catch (std::exception &e) {
        throw JSONRPCError(RPC_DESERIALIZATION_ERROR, "TX decode failed");
    }

    Object result;
    TxToJSON(tx, 0, result);

    return result;
}

Value signrawtransaction(const Array& params, bool fHelp)
{
    if (fHelp || params.size() < 1 || params.size() > 4)
        throw runtime_error(
            "signrawtransaction <hex string> [{\"txid\":txid,\"vout\":n,\"scriptPubKey\":hex,\"redeemScript\":hex},...] [<privatekey1>,...] [sighashtype=\"ALL\"]\n"
            "Sign inputs for raw transaction (serialized, hex-encoded).\n"
            "Second optional argument (may be null) is an array of previous transaction outputs that\n"
            "this transaction depends on but may not yet be in the block chain.\n"
            "Third optional argument (may be null) is an array of base58-encoded private\n"
            "keys that, if given, will be the only keys used to sign the transaction.\n"
            "Fourth optional argument is a string that is one of six values; ALL, NONE, SINGLE or\n"
            "ALL|ANYONECANPAY, NONE|ANYONECANPAY, SINGLE|ANYONECANPAY.\n"
            "Returns json object with keys:\n"
            "  hex : raw transaction with signature(s) (hex-encoded string)\n"
            "  complete : 1 if transaction has a complete set of signature (0 if not)"
            + HelpRequiringPassphrase());

    RPCTypeCheck(params, list_of(str_type)(array_type)(array_type)(str_type), true);

    vector<unsigned char> txData(ParseHexV(params[0], "argument 1"));
    CDataStream ssData(txData, SER_NETWORK, PROTOCOL_VERSION);
    vector<CTransaction> txVariants;
    while (!ssData.empty())
    {
        try {
            CTransaction tx;
            ssData >> tx;
            txVariants.push_back(tx);
        }
        catch (std::exception &e) {
            throw JSONRPCError(RPC_DESERIALIZATION_ERROR, "TX decode failed");
        }
    }

    if (txVariants.empty())
        throw JSONRPCError(RPC_DESERIALIZATION_ERROR, "Missing transaction");

    // mergedTx will end up with all the signatures; it
    // starts as a clone of the rawtx:
    CTransaction mergedTx(txVariants[0]);
    bool fComplete = true;

    // Fetch previous transactions (inputs):
    CCoinsView viewDummy;
    CCoinsViewCache view(viewDummy);
    {
        LOCK(mempool.cs);
        CCoinsViewCache &viewChain = *pcoinsTip;
        CCoinsViewMemPool viewMempool(viewChain, mempool);
        view.SetBackend(viewMempool); // temporarily switch cache backend to db+mempool view

        BOOST_FOREACH(const CTxIn& txin, mergedTx.vin) {
            const uint256& prevHash = txin.prevout.hash;
            CCoins coins;
            view.GetCoins(prevHash, coins); // this is certainly allowed to fail
        }

        view.SetBackend(viewDummy); // switch back to avoid locking mempool for too long
    }

    bool fGivenKeys = false;
    CBasicKeyStore tempKeystore;
    if (params.size() > 2 && params[2].type() != null_type)
    {
        fGivenKeys = true;
        Array keys = params[2].get_array();
        BOOST_FOREACH(Value k, keys)
        {
            CBitcoinSecret vchSecret;
            bool fGood = vchSecret.SetString(k.get_str());
            if (!fGood)
                throw JSONRPCError(RPC_INVALID_ADDRESS_OR_KEY, "Invalid private key");
            CKey key = vchSecret.GetKey();
            tempKeystore.AddKey(key);
        }
    }
    else
        EnsureWalletIsUnlocked();

    // Add previous txouts given in the RPC call:
    if (params.size() > 1 && params[1].type() != null_type)
    {
        Array prevTxs = params[1].get_array();
        BOOST_FOREACH(Value& p, prevTxs)
        {
            if (p.type() != obj_type)
                throw JSONRPCError(RPC_DESERIALIZATION_ERROR, "expected object with {\"txid'\",\"vout\",\"scriptPubKey\"}");

            Object prevOut = p.get_obj();

            RPCTypeCheck(prevOut, map_list_of("txid", str_type)("vout", int_type)("scriptPubKey", str_type));

            uint256 txid = ParseHashO(prevOut, "txid");

            int nOut = find_value(prevOut, "vout").get_int();
            if (nOut < 0)
                throw JSONRPCError(RPC_DESERIALIZATION_ERROR, "vout must be positive");

            vector<unsigned char> pkData(ParseHexO(prevOut, "scriptPubKey"));
            CScript scriptPubKey(pkData.begin(), pkData.end());

            CCoins coins;
            if (view.GetCoins(txid, coins)) {
                if (coins.IsAvailable(nOut) && coins.vout[nOut].scriptPubKey != scriptPubKey) {
                    string err("Previous output scriptPubKey mismatch:\n");
                    err = err + coins.vout[nOut].scriptPubKey.ToString() + "\nvs:\n"+
                        scriptPubKey.ToString();
                    throw JSONRPCError(RPC_DESERIALIZATION_ERROR, err);
                }
                // what todo if txid is known, but the actual output isn't?
            }
            if ((unsigned int)nOut >= coins.vout.size())
                coins.vout.resize(nOut+1);
            coins.vout[nOut].scriptPubKey = scriptPubKey;
            coins.vout[nOut].nValue = 0; // we don't know the actual output value
            view.SetCoins(txid, coins);

            // if redeemScript given and not using the local wallet (private keys
            // given), add redeemScript to the tempKeystore so it can be signed:
            if (fGivenKeys && scriptPubKey.IsPayToScriptHash())
            {
                RPCTypeCheck(prevOut, map_list_of("txid", str_type)("vout", int_type)("scriptPubKey", str_type)("redeemScript",str_type));
                Value v = find_value(prevOut, "redeemScript");
                if (!(v == Value::null))
                {
                    vector<unsigned char> rsData(ParseHexV(v, "redeemScript"));
                    CScript redeemScript(rsData.begin(), rsData.end());
                    tempKeystore.AddCScript(redeemScript);
                }
            }
        }
    }

    const CKeyStore& keystore = (fGivenKeys ? tempKeystore : *pwalletMain);

    int nHashType = SIGHASH_ALL;
    if (params.size() > 3 && params[3].type() != null_type)
    {
        static map<string, int> mapSigHashValues =
            boost::assign::map_list_of
            (string("ALL"), int(SIGHASH_ALL))
            (string("ALL|ANYONECANPAY"), int(SIGHASH_ALL|SIGHASH_ANYONECANPAY))
            (string("NONE"), int(SIGHASH_NONE))
            (string("NONE|ANYONECANPAY"), int(SIGHASH_NONE|SIGHASH_ANYONECANPAY))
            (string("SINGLE"), int(SIGHASH_SINGLE))
            (string("SINGLE|ANYONECANPAY"), int(SIGHASH_SINGLE|SIGHASH_ANYONECANPAY))
            ;
        string strHashType = params[3].get_str();
        if (mapSigHashValues.count(strHashType))
            nHashType = mapSigHashValues[strHashType];
        else
            throw JSONRPCError(RPC_INVALID_PARAMETER, "Invalid sighash param");
    }

    bool fHashSingle = ((nHashType & ~SIGHASH_ANYONECANPAY) == SIGHASH_SINGLE);

    // Sign what we can:
    for (unsigned int i = 0; i < mergedTx.vin.size(); i++)
    {
        CTxIn& txin = mergedTx.vin[i];
        CCoins coins;
        if (!view.GetCoins(txin.prevout.hash, coins) || !coins.IsAvailable(txin.prevout.n))
        {
            fComplete = false;
            continue;
        }
        const CScript& prevPubKey = coins.vout[txin.prevout.n].scriptPubKey;

        txin.scriptSig.clear();
        // Only sign SIGHASH_SINGLE if there's a corresponding output:
        if (!fHashSingle || (i < mergedTx.vout.size()))
            SignSignature(keystore, prevPubKey, mergedTx, i, nHashType);

        // ... and merge in other signatures:
        BOOST_FOREACH(const CTransaction& txv, txVariants)
        {
            txin.scriptSig = CombineSignatures(prevPubKey, mergedTx, i, txin.scriptSig, txv.vin[i].scriptSig);
        }
        if (!VerifyScript(txin.scriptSig, prevPubKey, mergedTx, i, SCRIPT_VERIFY_P2SH | SCRIPT_VERIFY_STRICTENC, 0))
            fComplete = false;
    }

    Object result;
    CDataStream ssTx(SER_NETWORK, PROTOCOL_VERSION);
    ssTx << mergedTx;
    result.push_back(Pair("hex", HexStr(ssTx.begin(), ssTx.end())));
    result.push_back(Pair("complete", fComplete));

    return result;
}

Value sendrawtransaction(const Array& params, bool fHelp)
{
    if (fHelp || params.size() < 1 || params.size() > 1)
        throw runtime_error(
            "sendrawtransaction <hex string>\n"
            "Submits raw transaction (serialized, hex-encoded) to local node and network.");

    // parse hex string from parameter
    vector<unsigned char> txData(ParseHexV(params[0], "parameter"));
    CDataStream ssData(txData, SER_NETWORK, PROTOCOL_VERSION);
    CTransaction tx;

    // deserialize binary data stream
    try {
        ssData >> tx;
    }
    catch (std::exception &e) {
        throw JSONRPCError(RPC_DESERIALIZATION_ERROR, "TX decode failed");
    }
    uint256 hashTx = tx.GetHash();

    bool fHave = false;
    CCoinsViewCache &view = *pcoinsTip;
    CCoins existingCoins;
    {
        fHave = view.GetCoins(hashTx, existingCoins);
        if (!fHave) {
            // push to local node
            CValidationState state;
            if (!tx.AcceptToMemoryPool(state, true, false))
                throw JSONRPCError(RPC_DESERIALIZATION_ERROR, "TX rejected"); // TODO: report validation state
        }
    }
    if (fHave) {
        if (existingCoins.nHeight < 1000000000)
            throw JSONRPCError(RPC_INVALID_ADDRESS_OR_KEY, "transaction already in block chain");
        // Not in block, but already in the memory pool; will drop
        // through to re-relay it.
    } else {
        SyncWithWallets(hashTx, tx, NULL, true);
    }
    RelayTransaction(tx, hashTx);

    return hashTx.GetHex();
}
