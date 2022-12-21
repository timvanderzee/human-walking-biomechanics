%% Init
clear 
close all
clc

%% Ask for User Input - Data Storage Directory
% ------------------------------------- Folder where Level 3 data is stored
msg1 = ['Select folder where Level 3 data is stored'...
                            ' (should end with \Level 3 - MATLAB files\)']; 
disp(msg1);
datafolder = uigetdir(cd, msg1);

% ----------------------------------- Folder where GitHub repo is cloned to
msg2 = ['Select folder where GitHub repo is cloned to'...
                         ' (should end with \human-walking-biomechanics)'];
disp(msg2);
codefolder = uigetdir(cd, msg2);

addpath(genpath(codefolder))

%% Combine Exported Files
disp('Combining Trials')
disp('')
combine_trials

%% Create 5 stride files
disp('Creating 5 Stride Files')
disp('')
process_5strides

%% Analyze 5 strides files
disp('Analyzing 5 strides files')
disp('')
analyse_biomechanics_script

%% Analyze summary data
disp('Analyzing summary data')
disp('')
analyse_soft_tissue_work