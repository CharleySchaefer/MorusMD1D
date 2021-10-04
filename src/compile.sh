#!/bin/bash

executable="MorusMD"
debug=0
dir=".."
src=$dir/src
utils=$dir/utils
doc=$dir/doc

# CHECK FILE EXISTENCE
echo "COMPILE SIMULATION CODE" 2>&1 | tee compilation.err
if [ ! -d "$src" ] ; then
  echo "  Error:     directory $src is missing" 2>&1 | tee -a compilation.err
else
  echo "  Success:   directory $src found" 2>&1 | tee -a compilation.err
fi
if [ ! -d "$utils" ] ; then
  echo "  Error:     directory $utils is missing" 2>&1 | tee -a compilation.err
else
  echo "  Success:   directory $utils found" 2>&1 | tee -a compilation.err
fi

if [ ! `$utils/is_installed.sh gcc` -eq 1 ] ; then
  echo "  Error:     gcc is not installed" 2>&1 | tee -a compilation.err
else
  echo "  Success:   `gcc --version | head -n 1` is installed" 2>&1 | tee -a compilation.err
fi

echo "  Compiling: $executable" 2>&1 | tee -a compilation.err
if [ $debug -eq 1 ] ; then
  gcc -o $executable $src/main.c $src/printHelp.c  -lm -g 2>&1 | tee -a compilation.err
else
  gcc -o $executable $src/main.c $src/printHelp.c -lm 2>&1 | tee -a compilation.err
fi


if [ -f "$executable" ] ; then
  echo "  Success:   `./$executable --version | head -n 1` is installed" 2>&1 | tee -a compilation.err
else
  echo "  Error:     $executable is not installed" 2>&1 | tee -a compilation.err
fi




# CHECK FILE EXISTENCE
echo "CHECKING POSTPROCESSING TOOLS" 2>&1 | tee -a compilation.err
if [ ! `$utils/is_installed.sh gnuplot` -eq 1 ] ; then
  echo "  Warning:     gnuplot is not installed." 2>&1 | tee -a compilation.err
  echo "             install: sudo apt-get install gnuplot" 2>&1 | tee -a compilation.err
else
  echo "  Success:   `gnuplot --version | head -n 1` is installed" 2>&1 | tee -a compilation.err
fi

if [ ! `$utils/is_installed.sh matlab` -eq 1 ] ; then
  echo "  Warning:   Matlab is not installed." 2>&1 | tee -a compilation.err
  echo "             install Matlab OR GNU Octave" 2>&1 | tee -a compilation.err
else
  echo "  Success:   `matlab --version | head -n 1` is installed" 2>&1 | tee -a compilation.err
fi

if [ ! `$utils/is_installed.sh octave` -eq 1 ] ; then
  echo "  Warning:     gnuplot is not installed." 2>&1 | tee -a compilation.err
  echo "               install: sudo apt-get install octave" 2>&1 | tee -a compilation.err
else
  echo "  Success:   `octave --version | head -n 1` is installed" 2>&1 | tee -a compilation.err
fi


# MANUAL
echo "USAGE:" 2>&1 | tee -a compilation.err
echo "  See $doc/Manual.md (Markdown format)" 2>&1 | tee -a compilation.err


echo "DONE." 2>&1 | tee -a compilation.err
