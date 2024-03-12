#!/usr/bin/gnuplot

set data style lines

plot 'data/b.dat' using 1:3 title 'B internalisiert',\
'data/b.dat' using 1:4 title 'B praesentiert AG',\
'data/b.dat' using 1:5 title 'B teilt sich',\
'data/b.dat' using 1:6 title 'B Plasma'

pause -1 

plot 'data/t.dat' using 1:2 title 'T aktiv',\
'data/t.dat' using 1:3 title 'T teilt sich',\
'data/t.dat' using 1:4 title 'T bewaffnet'

pause -1 

plot 'data/b.dat' using 1:($2+$3+$4+$5+$6) title 'B gesamt',\
'data/t.dat' using 1:($2+$3+$4) title 'T gesamt'
#,\

pause -1

plot 'data/ag.dat' using 1:2 title 'Antigen',\
'data/ab.dat' using 1:($2/10) title 'Antikoerper / 10'

pause -1

