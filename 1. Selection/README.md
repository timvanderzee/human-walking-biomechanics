# Selection #

Data includes 3D marker positions, joint angles, moments, forces, powers from healthy individuals (N = 10) walking at various combinations of speed, step length, step frequency and step length. For each of 33 combinations of speed, step length etc., 60 s of recorded data is available. In addition, we have selected 5 strides of good quality data for each trial.

This folders includes code that allows the users to select a different 5 strides of good quality data, using 4 scripts:

1. combine_trials.m
2. select_5strides.m
3. process_5strides.m
4. compare_versions.m

![picture](dataflow.png)

## combine_trials.m ##
Data from each of 33 trials for each of 10 subjects has been exported from Visual3D (V3D) into a separate MAT file (pXexport_TX.mat). combine_trials.m combines all 33 trials of one subject into a single MAT file for convenience, called pX_AllStridesData.mat. 

## select_5strides.m ##
This script loads the pX_AllStridesData.mat file and allows the user to select a different 5 strides for further analysis. The indices defining the new interval are saved in 5strides_heelstrikes.mat. 

## process_5strides.m ##
This script loads the pX_AllStridesData.mat file and the 5strides_heelstrikes.mat file to create files containing only the 5 strides of good quality data, called pX_5StridesData.mat

## compare_versions.m ##
The pX_5StridesData.mat file depends on the chosen interval. The user may yield different versions of 5 strides files. The compare_versions.m script facilitates comparing such versions.
