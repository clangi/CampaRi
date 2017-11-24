#!/bin/bash
# $Id: run_test.sh,v 1.1 2008/05/06 18:09:24 matt Exp $
if [ `uname -i` == "i386" ]; then
  ABSINTH_BIN=_ABSINTH_HOME_/bin/i386/absinth
else
  ABSINTH_BIN=_ABSINTH_HOME_/bin/x86_64/absinth
fi


${ABSINTH_BIN} -k t3p.key > test.log;

