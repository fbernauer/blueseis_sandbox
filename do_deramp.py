#! /usr/bin/env python

import os

doy = ['246', '247', '248', '249', '250', '251', '252', '253', '254']

sta = 'BS3'

for d in doy:
    print('processind day '+d)
    os.system("./process_blueseis.py -F /bay_mobil/stromboli2018/BlueSeis/archive/2018/XX/"+sta+"/HJ*.D/XX."+sta+"..HJ*.D.2018."+d+" -R /bay_mobil/stromboli2018/BlueSeis/archive/2018/XX/"+sta+"/YR*.D/XX."+sta+"..YR*.D.2018."+d+" -M /home/fbernauer/Experiments/Huddle_BS123/rotation_matrix_eudemo.txt -O /bay_mobil/Stromboli2018/processed/ -l 3600000 -m 10 -a mean -o 1 -p 0")
