#!/bin/bash

if [ $# -lt 1 ] ; then
  echo "usage: $0 <mscript.m>"
  echo "       Calls matlab or octave to execute the mscript."
  echo "       Tries matlab first and then octave."
fi
mscript=$1

if [ -x "$(command -v matlab)" ]   ; then
  matlab -batch $mscript
elif [ -x "$(command -v octave)" ] ; then
  echo "$mscript" | octave --no-gui
else
  echo "error: install matlab or octave to run the mscript."
fi

