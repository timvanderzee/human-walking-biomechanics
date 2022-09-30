# Welcome to the 'human-walking-biomechanics' repository!

This GitHub repository includes MATLAB scripts and functions that can be used to **select, visualize, and analyze** human walking biomechanics data from an open access data repository (LINK HERE). The data is described in more detail in a corresponding data descriptor paper (LINK HERE). 

Data exists on three different levels:
1. Raw motions and forces (.c3d)
2. Processed motion and forces in Visual3D software (.cmo)
3. Exported inverse dynamics data, including biomechanical variables (.mat)

This GitHub repository exclusively deals with level 3 of the data. Specifically, it can be used to:
1. **Select** 5 good quality strides of data from one of 33 different 90 s walking trials
2. **Visualize** the 5 good quality strides
3. **Analyze** the 5 good quality strides, e.g., computes mechanical work and average force

## main.m ##

Running the main.m script results in recreating both all strides files and 5 strides files from exported data. If the provided exported data is the same as in the data repository (LINK HERE), then the recreated all strides files and 5 strides files should also be the same as the ones in the data repository (LINK HERE). 

The main.m script therefore shows how all strides files and 5 strides files are created from exported files.
