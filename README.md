#Brittany
========

plsr approach to extracting water quality data from high frequency spectrophotometric observations

This project uses R, and is focused on learning how to make rational decisions for taking samples to calibrate PLSR models. The metrics we are using to evaluate the quality of the PLSR model are based on Nash-Sutcliffe goodness of fit, and a tool developed by Ritter and Muñoz for bootstrapping an assessment of the model's goodness of prediction. 

Some functions have been written to compartmentalize the reading of data, and several stages of analysis. Nevertheless, the code is messy and unweildly. It would be a good idea to think a little about how to smooth out the data management process. Maybe adding a paths file...

One major anoying thing is the initiation of the Chemographic table. I am opening the table at the beginning of the main script, but the correct number of columns is not know until subsetting is done. This should be solved sooner than later.

Also, the main file keeps changing a lot. This is fine, but I feel like the functions and accessory files are not general enough to make these changes robust among different versions. There should be a way of centralizing some definitions and function needs.

Also, I am doing a poor job of resource management. Any one prediction requires the use of a suite of files, but the needed files may be a subset of about 16-18 files. It seemed to make sense to pre-load all of these, as re-loading them repeatedly is not very time responsible, but at the same time, the initial load takes quite a long time, and comits a large amound of data to memory. This is not ideal, but hasn't caused a system crash yet. 


