clear all; close all; clc

% add code folder to datapath
scriptname = 'examplescript.m';
scriptfile = which(scriptname);
scriptfolder = scriptfile(1:(end-length(scriptname)));
addpath(genpath(scriptfolder));
cd(scriptfolder)

% ALSO MAKE SURE YOU ADD THE DATAFOLDER TO PATH!

%% choose variable
% in this section of code define the following 4 variables

variable_type = 'power';
% variable type = 'angle', 'moment', or 'power'

leg = 'left';
%leg = 'right' or 'left'

joint = 'ankle';
% joint = 'ankle', 'knee', or 'hip'

experiment = 'preferred walking';
% experiment= 'constant step length', 'constant step frequency', 'constant speed', or 'preferred walking'

%% look-up variable and trials. 
% Variable name is a string continaing the name of the variable which is
% found in the data files
% typename is a sting containing the name of the field where the variable
% can be found
% trials is a variable that contains the trials numbers for trials
% matching the experiment defined above. These trials are returned from
% lowest speed to highest speed

[variablename, typename, trials] = lookup_variable_name(variable_type, leg, joint, experiment);


%% calls on functions to look at subject means for a number of trials
% col contains the column of data that will be plotted in the following 2 
% sections
col=variableplot(typename, variablename, trials);
trial=input('What trial are you interested in?  [1-33]');
% trial is a single number for the trial that will be plotted in the
% following 2 sections


%% calls on function to look at each subjects' data for a single trial

visualizejointkin(trial, joint, col)
subj=input('Which subject (from 1-9) are you interested in?');

% subj is a single number for the subject that will be plotted in the
% following section


%% look at specific subject
% one subject, one trial 
plot5steps(typename, variablename, trial, subj, col)
GRFplot(datafolder, codefolder, subj, trial)

