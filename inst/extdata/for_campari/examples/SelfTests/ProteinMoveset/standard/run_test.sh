#!/bin/bash
if [ `uname -i` == "i386" ]; then
  ABSINTH_BIN=/packages/absinth/bin/i386/absinth
else
  ABSINTH_BIN=/packages/absinth/bin/x86_64/absinth
fi


${ABSINTH_BIN} -k gprotu.key > test.log;
