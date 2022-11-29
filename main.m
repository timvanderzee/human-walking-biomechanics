clear all; close all; clc
% loads exported files, combines and saves, selects 5 strides and saves, analyzing 5 strides and saves summary data 

% folder where Level 3 data is stored (should end with \Level 3 - MATLAB files\)
fprintf('Select folder where Level 3 data is stored (should end with \Level 3 - MATLAB files)');
datafolder = uigetdir(cd, 'Select folder where Level 3 data is stored (should end with \Level 3 - MATLAB files)');

% folder where GitHub repo is cloned to (should end with \human-walking-biomechanics\)
fprintf('Select folder where GitHub repo is cloned to (should end with \human-walking-biomechanics)')
codefolder = uigetdir(cd, 'Select folder where GitHub repo is cloned to (should end with \human-walking-biomechanics)');

%% Combine exported files
disp('Combining')
disp('')
import_folder = fullfile(datafolder,'V3D exported data');
export_folder = fullfile(datafolder,'All Strides Data files');
combine_trials

%% Create 5 stride files
disp('Creating 5 stride files')
disp('')
import_folder = fullfile(datafolder,'All Strides Data files');
export_folder = fullfile(datafolder,'5 Strides Data files');
process_5strides

%% Analyze 5 strides files
disp('Analyzing 5 strides files')
disp('')
import_folder = export_folder;
export_folder = import_folder;
analyse_biomechanics_script

%% Analyze summary data
disp('Analyzing summary data')
disp('')
folder = export_folder;
analyse_soft_tissue_work