// Copyright (c) 2009-2010 Satoshi Nakamoto
// Copyright (c) 2009-2012 The Bitcoin developers
// Distributed under the MIT/X11 software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.
#ifndef BITCOIN_INIT_H
#define BITCOIN_INIT_H

#include "wallet.h"

extern CWallet* pwalletMain;

void StartShutdown();
bool ShutdownRequested();
void Shutdown();
bool AppInit2(boost::thread_group& threadGroup);
std::string HelpMessage();

// FBX proof of stake voting test
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
extern std::string strPosxAddresses[POSX_ROUND_MAX][POSX_PAIR_MAX][POSX_ANSWER_MAX][POSX_ASK_MAX][POSX_PROB_MAX];

#endif
