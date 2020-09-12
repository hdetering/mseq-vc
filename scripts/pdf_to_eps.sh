#!/bin/bash
# vim: syntax=sh tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert PDF file to EPS with specified geometry.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-09-12
#------------------------------------------------------------------------------

PDF=$1
X_MAX_PX="2250"
Y_MAX_PX="2625"
FORMAT="eps" # output format

# calculate width/height ratio of max size
max_dim_xratio=$(echo "${X_MAX_PX}/${Y_MAX_PX}" | bc -l)
# get PDF dimensions
pdf_dim_xratio=$(identify ${PDF} | awk '{split($3,a,"x");print a[1]/a[2]}')

if (( $(echo "${pdf_dim_xratio} > ${max_dim_xratio}" | bc -l) )); then
  # determine target geometry based on width
  x_px=${X_MAX_PX}
  y_px=$(awk -v x=${x_px} -v r=${pdf_dim_xratio} 'BEGIN{printf "%.0f",x/r}')  
else
  # determine target geometry based on height
  y_px=${Y_MAX_PX}
  x_px=$(awk -v y=${y_px} -v r=${pdf_dim_xratio} 'BEGIN{printf "%.0f",y*r}')
fi

# convert from PDF to TIFF
convert -resize ${x_px}x${y_px} ${PDF} ${PDF%.pdf}.${FORMAT}
