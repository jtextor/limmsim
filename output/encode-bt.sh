#!/bin/bash
mogrify -format png dumps/b/*.bmp

for file in dumps/b/*.png 
do
	convert $file -channel blue -fx '0' $file
done

mogrify -format tga dumps/b/*.png

cd dumps/b/
mencoder mf://*.tga -mf w=180:h=180:fps=10:type=tga -ovc lavc -lavcopts vcodec=mjpeg:vbitrate=12000:vhq:vpass=1 -o test.avi
mencoder mf://*.tga -mf w=180:h=180:fps=10:type=tga -ovc lavc -lavcopts vcodec=mjpeg:vbitrate=12000:vhq:vpass=2 -o test.avi
cd ../../
rm -f dumps/b/*.bmp
rm -f dumps/b/*.tga
mv dumps/b/test.avi b-t.avi
