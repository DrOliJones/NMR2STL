# NMR2STL
rough matlab and Python code to turn 2D NMR files to STL format
Instructions for NMR2STL code

This is an initial version of the program for testing. At this time, I have left the files as raw matlab code as compiling is currently breaking some features.


The application still has some limitations. Stability is not promised.


With these limitations, you must therefore do the following:

   -Your TOPSPIN raw data should be exported as a text file and named: NMR_TOPSPIN.txt (other names won’t work)

   -Your TOPSPIN raw data should be the full spectra, do not zoom in and/or cut section before you export it. Do not delete any of the data in the file before use.

   -You must have python installed. (We use up to v3.7.3, and haven't tested this with other versions. We can't see why it wouldn't work – but never say never)

   -Extract the zip to a folder, and put your raw data in this same folder. (I'll make it so it works in any directory at a later date if I get time)

   -Double click start_nmr2stl.m to open Matlab. Click Run. (Big green play button)


NB: At present this code only works with Windows and with Bruker (TOPSPIN) data.
