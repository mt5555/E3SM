#!/bin/bash
# Build an animated GIF of a plotxz.py field across all output timesteps.
#
# Usage:
#   ./mtest-animation.sh <netcdf_file> <varname> [contour_min,max]
#
# Example:
#   ./mtest-animation.sh /pscratch/sd/t/taylorm/homme/m2/movies-L60-nu_top0/mtest21.nc w -1.0,1.0

set -e

NCFILE=$1
VARNAME=$2
CRANGE=${3:--1.0,1.0}

if [[ -z $NCFILE || -z $VARNAME ]]; then
    echo "usage: $0 <netcdf_file> <varname> [contour_min,max]"
    exit 1
fi

PLOTXZ=$HOME/codes/nclscript/contour/plotxz.py
OUTDIR=$(dirname "$NCFILE")
BASE=$(basename "$NCFILE" .nc)
FRAMES=$OUTDIR/${BASE}_${VARNAME}_frames

# plotxz.py writes frames next to its input file, so symlink the netcdf
# file into a scratch working directory to keep the frame pngs isolated
mkdir -p "$FRAMES"
ln -sf "$(realpath "$NCFILE")" "$FRAMES/$(basename "$NCFILE")"

# generate one plot png per timestep
python3 "$PLOTXZ" -i "$FRAMES/$(basename "$NCFILE")" -c "$CRANGE" -t -1 "$VARNAME"

# combine frames into a looping animated gif (ffmpeg handles large/many
# frames far more efficiently than ImageMagick's convert, which tends to
# exhaust its pixel-cache memory on a login node), then remove the
# working directory of per-frame pngs
GIF=$OUTDIR/${BASE}_${VARNAME}_animation.gif
ffmpeg -nostdin -y -framerate 5 -i "$FRAMES/${BASE}.${VARNAME}_t%03d.png" \
    -vf "split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" "$GIF"
rm -rf "$FRAMES"

echo "wrote $GIF"
