#!/bin/bash
mogrify -format png dumps/ag/*.bmp

cd dumps/ag/
mencoder mf://*.png -mf w=180:h=180:fps=10:type=png -ovc lavc -lavcopts vcodec=mpeg4:vbitrate=12000:vhq:vpass=1 -o test.avi
mencoder mf://*.png -mf w=180:h=180:fps=10:type=png -ovc lavc -lavcopts vcodec=mpeg4:vbitrate=12000:vhq:vpass=2 -o test.avi
cd ../../
rm dumps/ag/*.bmp
mv dumps/ag/test.avi agab.avi
