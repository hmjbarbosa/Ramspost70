'open catt-brams_g2.ctl'
'set mpdset ./mres'
'enable print topo'
'set gxout shaded'
'd topo'
'draw title Topography'
'run ./cbarn'
'print'
'disable print'
'!./gxgif -r -i topo'
'clear'
'set mpdset ./mres'
'enable print temp'
'set gxout contour'
'd tempc'
'draw title Temperature 00H(GMT) 25/02/2007'
'print'
'disable print'
'!./gxgif -r -i temp'
'clear'
'set mpdset ./mres'
'enable print co'
'set gxout contour'
'set t last'
'd CO'
'draw title CO 00H(GMT) 25/02/2007'
'print'
'disable print'
'!./gxgif -r -i co'
'clear'
'set mpdset ./mres'
'enable print pm25'
'set gxout contour'
'd pm25'
'set t last'
'draw title Particulate Material 00H(GMT) 25/02/2007'
'print'
'disable print'
'!./gxgif -r -i pm25'
'clear'
'set mpdset ./mres'
'enable print src1'
'set gxout contour'
'd src1'
'set t last'
'draw title Source-1 00H(GMT) 25/02/2007'
'print'
'disable print'
'!./gxgif -r -i src1'
'quit'
