#!/bin/bash
executable=$1
if [ -x "$(command -v $executable)" ]
then
  echo "1"
else
  echo "0"
fi

