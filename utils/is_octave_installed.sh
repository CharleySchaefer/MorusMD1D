#!/bin/bash
exe="octave"
if [ -x "$(command -v $exe)" ]
then
  echo "1"
else
  echo "0"
fi

