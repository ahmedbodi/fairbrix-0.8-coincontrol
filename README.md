Fairbrix is the second oldest scrypt based cryptocoin (after Tenebrix).
 - scrypt as a proof of work scheme
 - 5 minute block targets
 - 25 coins per block (constant forever)
 - 2016 blocks (1 week) to retarget difficulty

This version (Fairbrix v0.8.5.1beta) is based on Litecoin version 'exp-0.8.3.7-cc',
Bitcoin v0.8.1 (for IRC node discovery) and Litecoin v0.8.5.1 (for security fixes).

Q&D list of planned features: https://github.com/wiggi/fairbrix-0.8-coincontrol/blob/master/README.md


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

To compile using Qt 5.1.1 see https://bitcointalk.org/index.php?topic=149479.0
("Building headless Bitcoin and Bitcoin-qt on Windows")

To compile using Qt 4.8.5

1 Prepare your build system. I suggest setting up a clean virtual machine (Windows7 SP1 64bit)
via Virtualbox or similar.

1.1 Install MinGW: https://sourceforge.net/downloads/mingw
Make sure to select prepackaged repository catalogues (use gcc 4.6.2 as 4.7.2 won't work).
A minimal setup will require "C Compiler", "C++ Compiler" and "MSYS Basic System" to be installed.

1.2.a Install ActivePerl Community Edition: http://www.activestate.com/activeperl/downloads
Tested with ActivePerl-5.16.3.1603-MSWin32-x64, but newest x86 should work just fine.

1.2.b From a MinGw shell (MSYS) run the following:

    mingw-get install msys-perl

1.3 Add MinGW bin folder to your PATH environment variable (C:\MinGW\bin if you used installer defaults).

2 Download, unpack and build required dependencies. (save them in c:\deps folder)

2.1 OpenSSL: http://www.openssl.org/source/openssl-1.0.1e.tar.gz
From a MinGw shell (MSYS), unpack the source archive with tar (this will avoid symlink issues) then configure and make:

    cd /c/deps/
    tar xvfz openssl-1.0.1e.tar.gz
    cd openssl-1.0.1e
    ./config
    make

2.2 Berkeley DB: http://download.oracle.com/berkeley-db/db-4.8.30.NC.tar.gz
From a MinGW shell:

    cd /c/deps/
    tar xvfz db-4.8.30.NC.tar.gz
    cd db-4.8.30.NC/build_unix
    ../dist/configure --disable-replication --enable-mingw --enable-cxx

Edit your build_unix/db.h by replacing line 113 `typedef pthread_t db_threadid_t;`
with `typedef u_int32_t db_threadid_t;` and:

    make

2.3 Boost: http://sourceforge.net/projects/boost/files/boost/1.54.0/
Unzip boost inside your C:\deps folder, then bootstrap and compile from a Windows command prompt (as administrator):

    cd C:\deps\boost_1_54_0\
    bootstrap.bat mingw
    b2 --build-type=complete --with-chrono --with-filesystem --with-program_options --with-system --with-thread toolset=gcc stage

3 Compile leveldb, then compile fairbrix.

3.1 Extract fairbrix (for example to C:\fairbrix) then start MinGW shell and change into leveldb folder:

    cd /C/fairbrix/src/leveldb
    TARGET_OS=NATIVE_WINDOWS make libleveldb.a libmemenv.a

3.2 From a Windows command prompt (as administrator) run:

    cd C:\fairbrix\src
    mingw32-make -f makefile.mingw
    strip fairbrixd.exe

4 Setup Qt 4.8.5 and compile Fairbrix-qt

4.1 Install Qt 4.8.5 http://download.qt-project.org/official_releases/qt/4.8/4.8.5/qt-win-opensource-4.8.5-mingw.exe
Setup will probably complain about MinGw being an unsupported version, accept anyway.

4.2 From "Qt 4.8.5 command prompt" configure then make:

    cd C:\fairbrix
    qmake "USE_UPNP=-" bitcoin-qt.pro
    mingw32-make -f Makefile.Release

Notes:

 - "Qt 4.8.5 command prompt" means: run C:\Qt\4.8.5\bin\qtvars.bat from Windows command prompt (as administrator)

 - Distribute libgcc_s_dw2-1.dll, libstdc++-6.dll and mingwm10.dll from C:\MinGW\bin folder,
   and QtCore4.dll, QtGui4.dll and QtNetwork4.dll from C:\Qt\4.8.5\bin folder along with the executable(s).
   (DLLs with same name from other folders crash the executable)

 - GUI animation (update spinner) doesn't play without having Qt 4.8.5 installed.


Precompiled packages
=====================

Windows Qt-client and daemon (fairbrix-0.8-coincontrol-20131217.zip):

    https://mega.co.nz/#!XUlzjSjR!LgCbOE0xfaEcVeWAuxgEpwNzx4tNl56eb2tK5O4RUW4


Linux daemon, statically linked.
To run in a new empty server (e.g. Ubuntu 32bit Amazon EC2 micro instance):

    wget http://www.xmail.net/wiggi/fairbrixd-0.8-20130904.zip
    unzip fairbrixd-0.8-20130904.zip
    mkdir .fairbrix
    cp -p fairbrixd-0.8-20130904/fbx.conf .fairbrix/
    cd fairbrixd-0.8-20130904
    ./fairbrixd -daemon
    ./fairbrixd getinfo


License
========

Fairbrix is released under the terms of the MIT license. See `COPYING` for more
information or see http://opensource.org/licenses/MIT.

