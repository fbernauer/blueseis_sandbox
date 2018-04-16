# blueseis_sandbox

Version from Apr. 16. 2018

################################################################

The script "process_ramp.py" can be used to deramp BlueSeis data recorded before Feb. 20. 2018.
It processes a single channel. Day-files (or shorter) in miniseed-format can be processed.

command to execute process_ramp.py:
./process_ramp.py test_data/FOG.mseed test_data/RAMP.mseed test_data/deramped/ 7500 0 10
 
#################################################################
 
The script "process_blueseis.py" can be used to process BlueSeis data recorded after Feb. 20. 2018.
It processes three channel.
First, the traces are rotated using the rotation matrix file.
Second, the deramping is done for every single trace.
Third, the traces are rotated back.
 
command to execute process_blueseis.py:
./process_blueseis.py -F test_data/test_HJ1.mseed test_data/test_HJ2.mseed test_data/test_HJ3.mseed -R test_data/test_YR1.mseed test_dat    a/test_YR2.mseed test_data/test_YR3.mseed -O test_data/deramped/ -l 7500 -o 0 -m 10 -M test_data/rotation_matrix.txt

##################################################################
 
The directory "notebooks" contains a Jupyter notebook that demonstrates how the deramping works.
 
##################################################################
 
The directory "test_data" contains data examples (Fog-data and ramp-data) that is needed to run the scripts and the notebook, as well as     a rotation matrix file.
 
##################################################################
 
The directory "plots" contains plots demonstrating the output of "process_blueseis.py". The plots can be reproduced running the script "    plot_data.py".
