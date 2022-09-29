# Welcome to the 'human-walking-biomechanics' repository!

This GitHub repository includes MATLAB scripts and functions that can be used to **select, visualize, and analyze** human walking biomechanics data from an open access data repository (LINK HERE). The data is described in more detail in a corresponding data descriptor paper (LINK HERE). 

Data exists on three different levels:
1. Raw motions and forces (.c3d)
2. Processed motion and forces in Visual3D software (.cmo)
3. Exported inverse dynamics data, including biomechanical variables (.mat)

This GitHub repository exclusively deals with level 3 of the data. Specifically, it:
1. **Selects** 5 good quality strides of data from 60 s walking trials in exported inverse dynamics data
2. **Visualizes** the 5 good quality strides
3. **Analyzes** the 5 good quality strides, e.g., computes mechanical work done per stride
