%% ------------------------------------------------------------------- Init
clear
close all 
clc

%% Ask for User Input to Data Directory
msg = ['Select folder where the files are that have been exported '...
                                'from Visual3D (Level 3 - MATLAB files).'];
disp(msg);                                                    
folder = uigetdir(cd, msg);
addpath(genpath(folder));

%% Choose Variables

variable_type = 'power'; % 'angle', 'moment', or 'power'
leg           = 'left';  % 'right' or 'left'
joint         = 'ankle'; % 'ankle', 'knee', or 'hip'

experiment    = 'preferred walking'; % 'constant step length',
                                     % 'constant step frequency', 
                                     % 'constant speed', or 
                                     % 'preferred walking'

%% Look-up variable names and trials
[variablename, typename, trials] = ...
               lookup_variable_name(variable_type, leg, joint, experiment);

%% Subjects' Means for Several Trials
col   = variableplot(typename, variablename, trials);

% ------------------ trial that will be plotted in the following 2 sections
trial = input(['What trial are you interested in? '...
                            '(1-33 excluding step width trials (26-30))']);

%% Subjects' Data for a Single Trial
visualizejointkin(trial, joint, col)

% subject that will be plotted in the following section/follwoing 3 figures
subj = input('Which subject (from 1-9) are you interested in?');

%% Variable of interest (5 strides for one subject, one trial)
plot5strides(typename, variablename, trial, subj, col)

%% Ground reaction force (5 strides for one subject, one trial)
GRFplot(subj, trial)

%% 3x3 plot (5 strides for one subject, one trial)
plot3x3(subj, trial)

