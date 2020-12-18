#!/bin/sh

echo
echo "SHELL SCRIPT FOR MOVIE USING FFMPEG"
echo



# WRITE YOUR FILE NAME IN BETWEEN ""
imageDIR="aniImagesPara"
movieDIR="."
fps=10

################## SEARCH AND STORE ########################

ffmpeg -r 5 -i "$imageDIR/ani.%04d.png" -c:v libx264 -vf fps=$fps -pix_fmt yuv420p "$movieDIR/animate.mp4"

echo "DONE!"
