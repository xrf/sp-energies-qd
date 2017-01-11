#!/bin/sh
set -eu
dir=`dirname "$0"`
cd "$dir"
i=${PLOTFLAGS-} # use -i for interactive

./plot_fit.py $i --num-filled=1 --freq=1.0 --fit-start=6 ground
./plot_fit.py $i --num-filled=1 --freq=0.28 --fit-start=6 ground
./plot_fit.py $i --num-filled=2 --freq=1.0 --fit-start=6 ground
./plot_fit.py $i --num-filled=2 --freq=0.28 --fit-start=6 ground
./plot_fit.py $i --num-filled=3 --freq=1.0 --fit-start=9 ground
./plot_fit.py $i --num-filled=3 --freq=0.28 --fit-start=9 ground
./plot_fit.py $i --num-filled=4 --freq=1.0 --fit-start=10 ground
#needs more data:
#./plot_fit.py $i --num-filled=4 --freq=0.28 --fit-start=¿ ground
./plot_fit.py $i --num-filled=5 --freq=1.0 --fit-start=15 ground
./plot_fit.py $i --num-filled=5 --freq=0.28 --fit-start=17 ground
./plot_fit.py $i --num-filled=6 --freq=1.0 --fit-start=16 ground

./plot_fit.py $i --num-filled=1 --freq=1.0 --fit-start=7 add
./plot_fit.py $i --num-filled=1 --freq=0.28 --fit-start=7 add
./plot_fit.py $i --num-filled=2 --freq=1.0 --fit-start=10 add
./plot_fit.py $i --num-filled=2 --freq=0.28 --fit-start=10 add
./plot_fit.py $i --num-filled=3 --freq=1.0 --fit-start=13 add
./plot_fit.py $i --num-filled=3 --freq=0.28 --fit-start=13 add
./plot_fit.py $i --num-filled=4 --freq=1.0 --fit-start=13 add

./plot_fit.py $i --num-filled=1 --freq=1.0 --fit-start=7 rm
./plot_fit.py $i --num-filled=1 --freq=0.28 --fit-start=7 rm
./plot_fit.py $i --num-filled=2 --freq=1.0 --fit-start=10 rm
./plot_fit.py $i --num-filled=2 --freq=0.28 --fit-start=10 rm
./plot_fit.py $i --num-filled=3 --freq=1.0 --fit-start=13 rm
./plot_fit.py $i --num-filled=3 --freq=0.28 --fit-start=13 rm
./plot_fit.py $i --num-filled=4 --freq=1.0 --fit-start=13 rm
