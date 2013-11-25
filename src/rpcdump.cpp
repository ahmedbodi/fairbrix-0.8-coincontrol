// Copyright (c) 2009-2012 Bitcoin Developers
// Distributed under the MIT/X11 software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#include "init.h" // for pwalletMain
#include "bitcoinrpc.h"
#include "ui_interface.h"
#include "base58.h"

#include <boost/lexical_cast.hpp>

#define printf OutputDebugStringF

using namespace json_spirit;
using namespace std;

class CTxDump
{
public:
    CBlockIndex *pindex;
    int64 nValue;
    bool fSpent;
    CWalletTx* ptx;
    int nOut;
    CTxDump(CWalletTx* ptx = NULL, int nOut = -1)
    {
        pindex = NULL;
        nValue = 0;
        fSpent = false;
        this->ptx = ptx;
        this->nOut = nOut;
    }
};

// FBX proof of stake voting test
Value svimportprivkeys(const Array& params, bool fHelp)
{
    if (fHelp || params.size() != 0)
        throw runtime_error(
            "svimportprivkeys\n"
            "Add private keys from file 'svimportprivkeys.txt' to your wallet.\n"
            "Format:\n"
            "<privkey> <label>\n"
            "<privkey> <label>\n"
            "etc.");

    if (!GetBoolArg("-xmode", false))
        throw runtime_error(
            "Disabled unless xmode=true.");

    FILE *fp;
    fp = fopen("svimportprivkeys.txt", "r");
    if (fp == NULL)
        throw runtime_error(
            "Can't read svimportprivkeys.txt.");

    char sk[100], sl[100];

    for (unsigned int i = 0; i < 3192; i++)
    {
        if (fscanf(fp, "%52s ", sk) < 1)
        {
            fclose(fp);
            return "End of file (please close and restart fairbrix)";
        }

        if (fscanf(fp, "%99s ", sl) < 1)
        {
            fclose(fp);
            return "Error: couldn't read label";
        }

        string strSecret;
        string strLabel;
        strSecret.assign(sk);
        strLabel.assign(sl);

        CBitcoinSecret vchSecret;
        bool fGood = vchSecret.SetString(strSecret);

        if (!fGood) throw JSONRPCError(RPC_INVALID_ADDRESS_OR_KEY, "Invalid private key");

        CKey key = vchSecret.GetKey();
        CPubKey pubkey = key.GetPubKey();
        CKeyID vchAddress = pubkey.GetID();
        {
            LOCK2(cs_main, pwalletMain->cs_wallet);

            pwalletMain->MarkDirty();
            pwalletMain->SetAddressBookName(vchAddress, strLabel);

            if (!pwalletMain->AddKeyPubKey(key, pubkey))
                throw JSONRPCError(RPC_WALLET_ERROR, "Error adding key to wallet");

        }

        MilliSleep(10);

    }

    fclose(fp);
    return "Imported 3192 keys (please close and restart fairbrix)";
}

Value importprivkey(const Array& params, bool fHelp)
{
    if (fHelp || params.size() < 1 || params.size() > 3)
        throw runtime_error(
            "importprivkey <fairbrixprivkey> [label] [rescan=true]\n"
            "Adds a private key (as returned by dumpprivkey) to your wallet.");

    string strSecret = params[0].get_str();
    string strLabel = "";
    if (params.size() > 1)
        strLabel = params[1].get_str();

    // Whether to perform rescan after import
    bool fRescan = true;
    if (params.size() > 2)
        fRescan = params[2].get_bool();

    CBitcoinSecret vchSecret;
    bool fGood = vchSecret.SetString(strSecret);

    if (!fGood) throw JSONRPCError(RPC_INVALID_ADDRESS_OR_KEY, "Invalid private key");

    CKey key = vchSecret.GetKey();
    CPubKey pubkey = key.GetPubKey();
    CKeyID vchAddress = pubkey.GetID();
    {
        LOCK2(cs_main, pwalletMain->cs_wallet);

        pwalletMain->MarkDirty();
        pwalletMain->SetAddressBookName(vchAddress, strLabel);

        if (!pwalletMain->AddKeyPubKey(key, pubkey))
            throw JSONRPCError(RPC_WALLET_ERROR, "Error adding key to wallet");

        if (fRescan) {
            pwalletMain->ScanForWalletTransactions(pindexGenesisBlock, true);
            pwalletMain->ReacceptWalletTransactions();
        }
    }

    return Value::null;
}

Value dumpprivkey(const Array& params, bool fHelp)
{
    if (fHelp || params.size() != 1)
        throw runtime_error(
            "dumpprivkey <fairbrixaddress>\n"
            "Reveals the private key corresponding to <fairbrixaddress>.");

    string strAddress = params[0].get_str();
    CBitcoinAddress address;
    if (!address.SetString(strAddress))
        throw JSONRPCError(RPC_INVALID_ADDRESS_OR_KEY, "Invalid Fairbrix address");
    CKeyID keyID;
    if (!address.GetKeyID(keyID))
        throw JSONRPCError(RPC_TYPE_ERROR, "Address does not refer to a key");
    CKey vchSecret;
    if (!pwalletMain->GetKey(keyID, vchSecret))
        throw JSONRPCError(RPC_WALLET_ERROR, "Private key for address " + strAddress + " is not known");
    return CBitcoinSecret(vchSecret).ToString();
}
