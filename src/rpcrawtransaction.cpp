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
int64 posvValueDiff;    // aggregate money inflow and outflow for all addresses matching the vanity test string
int64 posvTxValue;      // output of one tx
int posvSig;            // inflow or outflow
int posvVerbose;
int posvFastmode;
int posvWarningCount;   // more than 1 address per output
int posvErrCount;       // more than 1 address per output, vanitygen addresses used for oracle are affected
std::string strPosvTest;// vanity test string, now generated from hash of starting block
                        // (so it's not known in advance)

// oracle related vars start here (not used in older test functions svdebugblock and svscanblocks)
bool posv_fOracle;      // poll the oracle

int posv_nMaturity;
#define POSV_MATURITY_MAX 4
// don't give exact time because daily OHLC data should be "good enough"
std::string posv_strMaturityDesc[POSV_MATURITY_MAX] = {"close of Aug 31, 2013", "close of Sep 30, 2013", "close of Oct 31, 2013", "close of Dec 31, 2013"};
int64 posv_nMaturityTime[POSV_MATURITY_MAX] = {1377979200, 1380571200, 1383253200, 1388523600};
// start block should be several weeks before maturity date even if it costs performance
int64 posv_nMaturityBlockStart[POSV_MATURITY_MAX] = {168000, 173110, 178000, 181000};
// this is a raw estimate, block time stamp is used for actual voting end
int64 posv_nMaturityBlockEnd[POSV_MATURITY_MAX] = {180000, 188000, 196000, 200000};

int posv_nPair;
#define POSV_PAIR_MAX 4
#define POSV_CHOICE_MAX 10                                          // this shouldn't be !=10
int64 posv_nOracleValuePerChoice[POSV_PAIR_MAX][POSV_CHOICE_MAX];   // voting results

// base58.h included here, but not in rpcblockchain.cpp
//static const char* pszBase58 = "123456789ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz";
int posvSanCount = 0; // failed sanity test
int posvCritCount = 0; // other critical error

// "checkpoint" for listing codes -- amounts of satoshis to be parked at a special vanitygen address
// to list anything (not yet implemented)
// ('svcalc LTCBTC' etc to get the numbers)
int64 posv_nPairName[POSV_MATURITY_MAX][POSV_PAIR_MAX] = {{20613063, 20615077, 0, 0},
                                                          {20613063, 20615077, 139414004, 0},
                                                          {20613063, 20615077, 139414004, 0},
                                                          {20613063, 20615077, 139414004, 0}};

// the following 5 arrays are only used in posv_strDescOracleResult
int posv_nPairScale[POSV_MATURITY_MAX][POSV_PAIR_MAX] = {{5, 2, 0, 0},
                                                         {5, 2, 4, 0},
                                                         {5, 2, 4, 0},
                                                         {5, 3, 4, 0}};
int posv_nPairMidpoint[POSV_MATURITY_MAX][POSV_PAIR_MAX] = {{11, 15,  0, 0},
                                                            {11, 16, 46, 0},
                                                            {11, 16, 46, 0},
                                                            {11,  8, 46, 0}};
bool posv_fPairIsAlt[POSV_MATURITY_MAX][POSV_PAIR_MAX] = {{0, 0, 0, 0},
                                                          {0, 0, 0, 0},
                                                          {0, 0, 0, 0},
                                                          {0, 0, 0, 0}};
#define POSV_SCALE_MAX 7 // up to 10 scales
// alternative scales
const double posv_AltScales[POSV_SCALE_MAX][10] =
{{ 0, 0, 0, 1, 2, 3, 4, 5, 6, 7 },
 { 0, 0, 0, 1, 2, 4, 8, 16, 32, 64},
 { 0, 0, 0, 1, 3, 10, 30, 100, 300, 1000},
 { 0, 0, 0, 1, 10, 100, 1000, 10000, 100000, 1000000},
 { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // reserved
 { 0, 0, 0, 4, 5, 6, 7, 8, 9, 10 },
 { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } // reserved
};
// logarithmic scales (for currency pairs)
const double posv_Strikes[POSV_SCALE_MAX][64] =
{{0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,   100000000, 10000000, 1000000, 100000, 10000, 1000, 100, 10,
  1, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001,  0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0 },

 {0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 3000000, 1000000, 300000, 100000, 30000, 10000,   3000, 1000, 300, 100, 30, 10, 3,
  1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001, 0.0003,   0.0001, 0.00003, 0.00001, 0.000003, 0.000001, 0.0000003, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0 },

// 4 strikes per order of magnitude, max strike diff 1/2  (e.g. 100/200, 500/1000)
 {0, 0, 0, 0, 0, 0, 0, 2000000,   1000000, 500000, 300000, 200000, 100000, 50000, 30000, 20000,   10000, 5000, 3000, 2000, 1000, 500, 300, 200,   100, 50, 30, 20, 10, 5, 3, 2,
  1, 0.5, 0.3, 0.2, 0.1, 0.05, 0.03, 0.02,    0.01, 0.005, 0.003, 0.002, 0.001, 0.0005, 0.0003, 0.0002,   0.0001, 0.00005, 0.00003, 0.00002, 0.00001, 0.000005, 0.000003, 0.000002,    0.000001, 0, 0, 0, 0, 0, 0, 0 },

 {0, 0, 100000, 75000, 50000, 30000, 20000, 15000,    10000, 7500, 5000, 3000, 2000, 1500, 1000, 750,    500, 300, 200, 150, 100, 75, 50, 30,    20, 15, 10, 7.5, 5, 3, 2, 1.5,
  1, 0.75, 0.5, 0.3, 0.2, 0.15, 0.1, 0.05,    0.03, 0.02, 0.015, 0.01, 0.0075, 0.005, 0.003, 0.002,     0.0015, 0.001, 0.00075, 0.0005, 0.0003, 0.0002, 0.00015, 0.0001,   0.000075, 0.00005, 0.00003, 0.00002, 0, 0, 0, 0 },

// 8 strikes per order of magnitude, max strike diff 2/3 (e.g. 50/75, 100/150)
 {10000, 7500, 5000, 4000, 3000, 2500, 2000, 1500,   1000, 750, 500, 400, 300, 250, 200, 150,    100, 75, 50, 40, 30, 25, 20, 15,    10, 7.5, 5, 4, 3, 2.5, 2, 1.5,
  1, 0.75, 0.5, 0.4, 0.3, 0.25, 0.2, 0.15,   0.1, 0.075, 0.05, 0.04, 0.03, 0.025, 0.02, 0.015,   0.01, 0.0075, 0.005, 0.004, 0.003, 0.0025, 0.002, 0.0015,   0.001, 0.00075, 0.0005, 0.0004, 0.0003, 0.00025, 0.0002, 0.00015 },

// 10 strikes per order of magnitude, max strike diff 3/4 (e.g. 750/1000, 300/400, 150/200)
 {1500, 1250, 1000, 750, 600, 500, 400, 300,    250, 200, 150, 125, 100, 75, 60, 50,    40, 30, 25, 20, 15, 12.5, 10, 7.5,    6, 5, 4, 3, 2.5, 2, 1.5, 1.25,
  1, 0.75, 0.6, 0.5, 0.4, 0.3, 0.25, 0.2,    0.15, 0.125, 0.1, 0.075, 0.06, 0.05, 0.04, 0.03,    0.025, 0.02, 0.015, 0.0125, 0.01, 0.0075, 0.006, 0.005,    0.004, 0.003, 0.0025, 0.002, 0.0015, 0.00125, 0.001, 0.00075 },

// 13 strikes per order of magnitude, max strike diff 4/5 (e.g. 40/50, 80/100, 100/125, 200/250)
 {300, 250, 200, 175, 150, 125, 100, 80,    70, 60, 50, 40, 35, 30, 25, 20,    17.5, 15, 12.5, 10, 8, 7, 6, 5,    4, 3.5, 3, 2.5, 2, 1.75, 1.5, 1.25,
  1, 0.8, 0.7, 0.6, 0.5, 0.4, 0.35, 0.3,    0.25, 0.2, 0.175, 0.15, 0.125, 0.1, 0.08, 0.07,    0.06, 0.05, 0.04, 0.035, 0.03, 0.025, 0.02, 0.0175,    0.015, 0.0125, 0.01, 0.008, 0.007, 0.006, 0.005, 0.004 }
};

static bool posv_IsListed(int ipair, int ivoting)
{
    if (ipair < 0 || ipair >= POSV_PAIR_MAX || ivoting < 0 || ivoting >= POSV_MATURITY_MAX)
        return false;

    if (!posv_nPairName[ivoting][ipair])
        return false;

    return true;
}
static int svhextoint(char c)
{
    if ((c >= 'a') && (c <= 'f')) return (int)(c - 'a' + 10);
    if ((c >= 'A') && (c <= 'F')) return (int)(c - 'A' + 10);
    if ((c >= '0') && (c <= '9')) return (int)(c - '0');
    posvSanCount++;
    return 0;
}
static bool svblockhashtovanity(int nHeight)
{
    if (nHeight < 0 || nHeight > nBestHeight) return false;

    CBlock block;
    CBlockIndex* pblockindex = FindBlockByHeight(nHeight);
    block.ReadFromDisk(pblockindex);

    string s = block.GetHash().GetHex();
    int ls = s.length();
    strPosvTest = "ERROR";
    for (unsigned int k = 0; k < 5; k++)
    {
        int kk = 2 + (k * 3);
        if (nHeight > 200000) kk = (ls / 2) + (k * 3); // fix: don't use 1st half of chars to avoid leading zeros
        int j = ( (svhextoint(s[kk]) * 256) + (svhextoint(s[kk+1]) * 16) + svhextoint(s[kk+2]) ) % 58;

        // vanitygen v0.22 can only generate addresses from fE7... to fdS...
        // pszBase58[13] == 'E'
        // pszBase58[36] == 'd'
        if (k == 0)
        {
            if (j >= 36)
                j = (j - 36) + 14;
            else if (j <= 13)
                j += 14;
        }

        if (j < 0 || j > 57)
            return false;

        strPosvTest[k] = pszBase58[j];
    }
    return true;
}
// result should be 6 uppercase chars (like "BTCUSD" etc)
static string posv_strDynPairName(int i, int ivoting)
{
    string sd = "??????";
    if (i < 0 || i >= POSV_PAIR_MAX || ivoting < 0 || ivoting > POSV_MATURITY_MAX)
        return sd;

    int64 id = posv_nPairName[ivoting][i];
    for (int kd = 5; kd >= 0; kd--)
    {
        sd[kd] = id % 26 + 'A';
        id /= 26;
    }
    return sd;
}
static string posv_strDescOracleResult(int ivoting, int i, int i2, bool show_result, bool show_implied)
{
    string sd = "error";
    if (ivoting < 0 || ivoting >= POSV_MATURITY_MAX || i < 0 || i >= POSV_PAIR_MAX || i2 < 0 || i2 >= POSV_CHOICE_MAX)
        return sd;

    if (i2 < 2)
    {
        sd = " ";
        if (show_result)
        {
          if (i2 == 0)
              sd = "invalid";
          else //if (i2 == 1)
              sd = "don't know or don't care";
        }
        return sd;
    }

    bool show_both = (show_result&&show_implied ? true : false);
    int y = posv_nPairScale [ivoting] [i];
    if (y < 0 || y >= POSV_SCALE_MAX)
        return sd;

    // alternative scales
    if (posv_fPairIsAlt [ivoting] [i])
    {
        double d = posv_AltScales [y] [i2];
        double dnext = (i2 < 9) ? posv_AltScales [y] [i2 + 1] : 0.0;

        if (show_result)
        {
            if ((i2 == 9) || ((i2 < 9) && (dnext - d > 1.5)))
                sd = ">=" + to_string(d);
            else
                sd = "==" + to_string(d);
        }
        if (show_implied)
        {
            if (show_both)
                sd = sd + ", ";
            else
                sd = "";

            if ((i2 < 9) && (dnext - d > 1.5))
                sd = sd + "<" + to_string(dnext);
        }

        return sd;
    }

    int x = posv_nPairMidpoint [ivoting] [i];
    if (x < 4 || x >= 63-4)                    // fixme: is 3 ok?, use constant
        return sd;

    if (show_result)
    {
        if (i2 <= 5)
            sd = ">" + to_string(posv_Strikes [y] [x - 5 + i2]);
        else // if (i2 >= 6)
            sd = "<=" + to_string(posv_Strikes [y] [x - 6 + i2]);
    }

    if (show_implied)
    {
        if (show_both)
            sd = sd + ", ";
        else
            sd = "";

        if (i2 == 2)
            sd = sd + "no upper limit";
        else if (i2 == 9)
            sd = sd + "no lower limit";
        else if (i2 <= 5)
            sd = sd + "<=" + to_string(posv_Strikes [y] [x - 5 + i2 - 1]);
        else // if (i2 >= 6)
            sd = sd + ">" + to_string(posv_Strikes [y] [x - 6 + i2 + 1]);
    }

    return sd;
}
static int posv_nMaturityBlockEndEx(int i)
{
    if (i < 0 || i >= POSV_MATURITY_MAX)
        return 0;

    // voting closed 30 days after maturity date
    int nHeight = posv_nMaturityBlockEnd[i];
    int64 t = posv_nMaturityTime[i] + (30 * 24 * 3600);

    if (nBestHeight > nHeight)
    {
        CBlockIndex* pblockindex = FindBlockByHeight(nHeight);

        // if more blocks than estimated have been produced
        while (pblockindex->GetBlockTime() < t && nHeight < nBestHeight)
        {
            nHeight++;
            pblockindex = pblockindex->pnext;
        }
    }

    if (GetTime() > t)
    {
        CBlockIndex* pblockindex = nBestHeight<nHeight ? FindBlockByHeight(nBestHeight) : FindBlockByHeight(nHeight);
        if (pblockindex->GetBlockTime() > t)
        {
            if (nBestHeight<nHeight) nHeight = nBestHeight;

            // if less blocks than estimated have been produced
            while (pblockindex->GetBlockTime() > t && nHeight > 0)
            {
                nHeight--;
                pblockindex = pblockindex->pprev;
            }
        }
    }

    return nHeight;
}
#define POSV_STATE_NOTYETOPEN 0
#define POSV_STATE_WARMUP 1
#define POSV_STATE_INPROGRESS 2
#define POSV_STATE_CLOSED 3
static int posv_nGetVotingState(int i)
{
    if (i < 0)
        return POSV_STATE_CLOSED;
    if (i >= POSV_MATURITY_MAX)
        return POSV_STATE_NOTYETOPEN;

    if (nBestHeight < posv_nMaturityBlockStart[i])
        return POSV_STATE_NOTYETOPEN;

    if (nBestHeight <= posv_nMaturityBlockEndEx(i))
    {
        CBlockIndex* pblockindex = FindBlockByHeight(nBestHeight);
        if (pblockindex->GetBlockTime() < posv_nMaturityTime[i])
            return POSV_STATE_WARMUP;
        else
            return POSV_STATE_INPROGRESS;
    }

    return POSV_STATE_CLOSED;
}
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

    if (!ExtractDestinations(scriptPubKey, type, addresses, nRequired))
    {
        posvCritCount++; // possible exploit

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

                int divisor = 1;
                for (unsigned int k = 0; k < POSV_PAIR_MAX; k++)
                {
                    // get oracle answer for all pairs
                    int64 v = (posvTxValue / divisor) % 10;
                    posv_nOracleValuePerChoice[k][v] += (posvTxValue * posvSig);
                    divisor *= 10;
                }
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
    posvWarningCount = posvErrCount = posvCritCount = 0;
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
    if (posvCritCount) result.push_back(Pair("illegible tx's", posvCritCount));

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
    posvWarningCount = posvErrCount = posvCritCount = 0;
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
    if (posvCritCount) result2.push_back(Pair("illegible tx's", posvCritCount));

    return result2;
}
Value svgetoracle(const Array& params, bool fHelp)
{
    if (fHelp || params.size() > 2)
        throw runtime_error(
            "svgetoracle [voting index] [pair index]\n"
            "Returns the oracle's answer for a given maturity,\n"
            "as determined by the proof of stake voting system.\n"
            "If no voting index is given, returns a concise answer\n"
            "for last closed voting.\n"
            "(Use 'svlistpairs' to list available currency pairs.\n"
            " Use 'svlistvotings' to list available votings/maturities)");

    posvVerbose = 0;
    posvFastmode = 1;
    posvValueDiff = 0;
    posvWarningCount = posvErrCount = posvCritCount = posvSanCount = 0;
    posv_fOracle = true;

    // clear all votes
    for (unsigned int i = 0; i < POSV_PAIR_MAX; i++)
        for (unsigned int j = 0; j < 10; j++)
            posv_nOracleValuePerChoice[i][j] = 0;

    // get voting index
    posv_nMaturity = 0;
    if (params.size() >= 1)
    {
        posv_nMaturity = params[0].get_int();
        if ((posv_nMaturity < 0) || (posv_nMaturity >= POSV_MATURITY_MAX))
            throw runtime_error("Maturity number out of range.");
    }
    // no user input -- use last "closed"
    else
    {
        for (unsigned int j = 0; j < POSV_MATURITY_MAX; j++)
            if (posv_nGetVotingState(j) >= POSV_STATE_CLOSED)
                posv_nMaturity = j;
    }

    posv_nPair = 0;
    unsigned int nPairMin = 0;
    unsigned int nPairMax = POSV_PAIR_MAX-1;
    if (params.size() > 1)
    {
        posv_nPair = params[1].get_int();
        nPairMin = nPairMax = posv_nPair;
        if ((posv_nPair < 0) || (posv_nPair >= POSV_PAIR_MAX))
            throw runtime_error("Pair index out of range.");

        if (!posv_IsListed(posv_nPair, posv_nMaturity))
            throw runtime_error("Nothing listed.");
    }

    bool fVerbose = params.size() < 1 ? false : true;

    if (!svblockhashtovanity(posv_nMaturityBlockStart[posv_nMaturity])) // calculates strPosvTest
        throw runtime_error("Couldn't get vanitygen string.");

    int nHeight0 = posv_nMaturityBlockStart[posv_nMaturity];
    int nHeight1 = posv_nMaturityBlockEndEx(posv_nMaturity);

    // voting in progress
    if (nBestHeight < nHeight1)
        nHeight1 = nBestHeight;

    if (nHeight1 < nHeight0 || nHeight0 < 0 || nHeight1 > nBestHeight)
        throw runtime_error("Block number(s) out of range.");
    if (strPosvTest.length() > 10)
        throw runtime_error("vanity string too long");

    // collect the oracle votes (for 1 maturity date, but all pairs)
    for (int nHeight = nHeight0; nHeight <= nHeight1; nHeight++)
    {
        CBlockIndex* pblockindex = FindBlockByHeight(nHeight);
        CBlock block;
        block.ReadFromDisk(pblockindex);

        posvblockToJSON(block, pblockindex);
    }

    Object result2;

    result2.push_back(Pair("maturity", posv_strMaturityDesc[posv_nMaturity]));
    result2.push_back(Pair("total voting coins", (double)posvValueDiff/(double)COIN));
    if (posv_nGetVotingState(posv_nMaturity) >= POSV_STATE_CLOSED)
        result2.push_back(Pair("voting closed", "the results are final"));

    if (posvWarningCount) result2.push_back(Pair("outputs with more than 1 addr", posvWarningCount));
    if (posvErrCount) result2.push_back(Pair("matching outputs with more than 1 addr", posvErrCount));
    if (posvCritCount) result2.push_back(Pair("illegible tx's", posvCritCount));
    if (posvSanCount) result2.push_back(Pair("sanity violations", posvSanCount));

    for (unsigned int j = nPairMin; j <= nPairMax; j++)
    {
        if (!posv_IsListed(j, posv_nMaturity)) break;
        posv_nPair = j; // is this needed?

        result2.push_back(Pair(" ", " "));
        result2.push_back(Pair("name", posv_strDynPairName(posv_nPair, posv_nMaturity)));

        // determine winner for concise list
        // (skip 'invalid' choice i==0)
        unsigned int imax = 1;
        int64 vmax = 0;
        if (!fVerbose)
        {
            for (unsigned int i = 1; i < POSV_CHOICE_MAX; i++)
            {
                int64 v = posv_nOracleValuePerChoice[posv_nPair][i];
                if (v > vmax)
                {
                    vmax = v;
                    imax = i;
                }
            }
        }

        int nResults = 0;
        for (unsigned int i = 1; i < POSV_CHOICE_MAX; i++)
        {
            int64 v = posv_nOracleValuePerChoice[posv_nPair][i];
            if ((v > 0) && (fVerbose || i == imax))
            {
                result2.push_back(Pair("result", posv_strDescOracleResult(posv_nMaturity, posv_nPair, i, true, true)));
                result2.push_back(Pair("voting coins", (double)v/(double)COIN));
                nResults++;
            }
        }
        if (!nResults)
            result2.push_back(Pair("no result", "no votes found"));
    }

    return result2;
}
Value svlistvotings(const Array& params, bool fHelp)
{
    if (fHelp || params.size() > 1)
        throw runtime_error(
            "svlistvotings [verbose=0]\n"
            "Lists all votings/maturity dates and info about\n"
            "'special addresses' for the proof of stake voting system.\n"
            "If verbose=1, returns info even if a voting is closed.");

    bool fVerbose = false;
    if (params.size() > 0)
        fVerbose = (params[0].get_int() != 0);

    Object result2;
    result2.push_back(Pair("available maturities", int(POSV_MATURITY_MAX)));

    for (unsigned int i = 0; i < POSV_MATURITY_MAX; i++)
    {
        result2.push_back(Pair(" ", " "));
        result2.push_back(Pair("index", int(i)));
        result2.push_back(Pair("maturity date", posv_strMaturityDesc[i]));

        int state = posv_nGetVotingState(i);

        if (state == POSV_STATE_NOTYETOPEN)
            result2.push_back(Pair("voting status", "not yet started"));
        else if (state == POSV_STATE_WARMUP)
            result2.push_back(Pair("voting status", "open for testing"));
        else if (state == POSV_STATE_INPROGRESS)
            result2.push_back(Pair("voting status", "in progress"));
        else if (state == POSV_STATE_CLOSED)
            result2.push_back(Pair("voting status", "closed"));

        if (!fVerbose && state != POSV_STATE_INPROGRESS)
            continue;

        if (state == POSV_STATE_INPROGRESS || (state == POSV_STATE_WARMUP || fVerbose))
        {
            result2.push_back(Pair("voting method", "coins in vanitygen address, when voting ends"));

            string s = "?????";
            if (svblockhashtovanity(posv_nMaturityBlockStart[i]))
                s = "vanitygen -X95 f" + strPosvTest;
            result2.push_back(Pair("vanitygen command", s));

            result2.push_back(Pair("command for help", "svlistpairs " + to_string(i)));
        }

        result2.push_back(Pair("voting start (fixed block height)", posv_nMaturityBlockStart[i]));
        if (state < POSV_STATE_CLOSED)
        {
            result2.push_back(Pair("voting end (est. block height)", posv_nMaturityBlockEnd[i]));
            result2.push_back(Pair("voting end (date)", "maturity date + 30 days"));
        }
        else
            result2.push_back(Pair("voting end (final block height)", posv_nMaturityBlockEndEx(i)));
    }

    return result2;
}
Value svlistpairs(const Array& params, bool fHelp)
{
    if (fHelp || params.size() < 1 || params.size() > 2)
        throw runtime_error(
            "svlistpairs <voting index> [verbose=1]\n"
            "Lists all pairs available for the given\n"
            "proof of stake voting. (pairs are usually\n"
            "currency pairs like 'LTCBTC')\n"
            "If verbose=1, returns voting instructions per pair.\n"
            "(Use 'svlistvotings' for a list of available\n"
            "votings/maturity dates.)");

    bool fVerbose = true;
    if (params.size() > 1)
        fVerbose = (params[1].get_int() != 0);

    int nVoting = params[0].get_int();
    if (nVoting < 0 || nVoting >= POSV_MATURITY_MAX)
        throw runtime_error("Maturity date index out of range.");

    Object result2;
    result2.push_back(Pair("list of available pairs for maturity date", posv_strMaturityDesc[nVoting]));

    // how many pairs are listed, not just array size
    unsigned int nPairs = 0;
    for (unsigned int i = 0; i < POSV_PAIR_MAX; i++)
    {
        if (!posv_nPairName[nVoting][i])
            break;

        nPairs++;
    }
    if (fVerbose)
        result2.push_back(Pair("voting method", to_string(nPairs) + " least significant digits of amount sent are used"));

    for (unsigned int i = 0; i < nPairs; i++)
    {
        result2.push_back(Pair(" ", " "));

        result2.push_back(Pair("index", int(i)));

        result2.push_back(Pair("name", posv_strDynPairName(i, nVoting)));

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

        result2.push_back(Pair("available choices", "1..9"));
        for (unsigned int i2 = 1; i2 < POSV_CHOICE_MAX; i2++)
        {
            string s = posv_strDescOracleResult(nVoting, i, i2, true, true);
            result2.push_back(Pair(s, int(i2)));
        }
        }
    }
    return result2;
}
Value svcalc(const Array& params, bool fHelp)
{
    if (fHelp || params.size() < 1 || params.size() > 2)
        throw runtime_error(
            "svcalc <string of 6 uppercase letters>\n"
            "returns the listing code for a given pair name.\n"
            "svcalc 1 <voting start block height>\n"
            "returns the vanitygen string for the voting.");

    Object result2;
    string s = params[0].get_str();

    if (params.size() == 1)
    {
        if (s.length() != 6)
            throw runtime_error("Need string length of 6.");

        int64 m = 1;
        int64 i = 0;
        for (int k = 5; k >=0; k--)
        {
            i += (int)((s[k]-'A') * m);
            m *= 26;
        }
        result2.push_back(Pair(s, i));
    }
    else if (s.length() == 1 && s[0] == '1')
    {
        int nHeight = params[1].get_int();
        if (nHeight < 0 || nHeight > nBestHeight)
            throw runtime_error("Block number out of range.");

        strPosvTest = "?????";
        if (!svblockhashtovanity(nHeight))
            throw runtime_error("Could not make vanitygen string from block hash hex.");
        if (posvSanCount)
            result2.push_back(Pair("sanity check failed", posvSanCount));
        result2.push_back(Pair("vanitygen string", strPosvTest));
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
