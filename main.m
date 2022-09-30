clear all; close all; clc

% folder where Level 3 data is stored (should end with \Level 3 - MATLAB files\)
datafolder = uigetdir;

% folder where GitHub repo is cloned to (should end with \human-walking-biomechanics\)
codefolder = uigetdir;

%% Combine
disp('Combining')
disp('')
import_folder = [datafolder,'V3D exported data'];
export_folder = [datafolder,'All Strides Data files'];
combine_trials

%% Create 5 stride files
disp('Creating 5 stride files')
disp('')
import_folder = export_folder;
export_folder = [datafolder,'5 Strides Data files'];
process_5steps

%% Analyze 5 strides files
disp('Analyzing 5 strides files')
disp('')
import_folder = export_folder;
export_folder = [codefolder,'3. Analysis'];
analyse_biomechanics_script

%% Analyze summary data
disp('Analyzing summary data')
disp('')
folder = export_folder;
analyse_soft_tissue_work