# lm2blocksino5d
C code is saved in  ./src
Matlab code is saved in ./matlab_version
Please find the other information in '20220207-SIMSET-listmode2-5Dsino.pptx'
During compile, if you have the following issue:
'cc1plus: error: unrecognized command line option...'
Solution: change the cmake file,you'll have to use -std=c++0x, for more recent versions you can use -std=c++11.
demo script for runing the convertion
./bin/lm2blocksino_5d your_listmode_directory your_sinogram_saving_name ./include/index_blockpairs_trans_int16_2x77x60.raw
