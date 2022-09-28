# Selection #

Data includes 3D marker positions, joint angles, moments, forces, powers from healthy individuals (N = 10) walking at various combinations of speed, step length, step frequency and step length.
For each of 26 combinations of speed, step length etc., 60 s of recorded data is available. In addition, we have selected 5 strides of good quality data for each trial.

This folders includes code that allows the users to select a different 5 strides of good quality data, using 3 scripts:

1. select_5steps.m
2. process_5steps.m

## 1. select_5steps.m ##
This script allows the user to select a different 5 strides for further analysis, and saves the indices defining the new interval in 5steps_heelstrikes.mat. 

## 2. process_5steps.m ##
This script creates new 5 strides data files using the provided 5steps_heelstrikes.mat file and the full (60 s) trial data. 


