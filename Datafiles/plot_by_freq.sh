#!/bin/sh
set -eu
dir=`dirname "$0"`
cd "$dir"

./plot_by_freq.py ground
./plot_by_freq.py add
./plot_by_freq.py rm
