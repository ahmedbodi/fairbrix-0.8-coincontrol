Fairbrix is the second oldest scrypt based cryptocoin.
 - scrypt as a proof of work scheme
 - 5 minute block targets
 - 25 coins per block (constant forever)
 - 2016 blocks (1 week) to retarget difficulty

This version (Fairbrix 0.8.3.7beta) is based on Litecoin version 'exp-0.8.3.7-cc'.

Litecoin is a lite version of Bitcoin using scrypt as a proof of work scheme.
 - 2.5 minute block targets
 - subsidy halves in 840k blocks (~4 years)
 - ~84 million total coins
The rest is the same as bitcoin.
 - 50 coins per block
 - 2016 blocks to retarget difficulty

Bitcoin is a free open source peer-to-peer electronic cash system that is
completely decentralized, without the need for a central server or trusted
parties.  Users hold the crypto keys to their own money and transact directly
with each other, with the help of a P2P network to check for double-spending.


Build instructions 
===================

Debian, Ubuntu, Mint
--------------------

Install Qt Creator.

Install libminiupnpc-dev.

Make sure that the required packages for Qt4 development of your
distribution are installed, for Debian and Ubuntu these are:

    apt-get install qt4-qmake libqt4-dev build-essential libboost-dev libboost-system-dev \
        libboost-filesystem-dev libboost-program-options-dev libboost-thread-dev \
        libssl-dev libdb++-dev

then execute the following:

    qmake
    make

Alternatively, use Qt Creator and open the `bitcoin-qt.pro` file.

An executable named `fairbrix-qt` will be built.

To compile the daemon with STATIC option:

    cd src
    make -f makefile.unix STATIC=1


Windows
--------

see https://bitcointalk.org/index.php?topic=149479.0
("Building headless Bitcoin and Bitcoin-qt on Windows")

 - "Qt 4.8.5 command prompt" means: run C:\Qt\4.8.5\bin\qtvars.bat from Windows command prompt (as administrator)

 - Steps 3.2 and 4.2 are already done.

 - An executable (in the \release folder) will be built.

Notes:

 - You will also need to distribute mingwm10.dll along with the executable(s).
 - Only use libgcc_s_dw2-1.dll, libstdc++-6.dll and mingwm10.dll from C:\MinGW\bin folder,
   and QtCore4.dll, QtGui4.dll and QtNetwork4.dll from C:\Qt\4.8.5\bin folder.
   (DLLs with same name from other folders crash the executable)
 - Keep Qt 4.8.5 installed. (optional but recommended)


License
========

Fairbrix is released under the terms of the MIT license. See `COPYING` for more
information or see http://opensource.org/licenses/MIT.

