function plot5steps(field, variable, trial, subject, direction)

% Plots 5 steps and the mean of 5 steps for a given motion capture variable
% in a specified plane of motion

% INPUTS: field- a string that contains the name of the field that contains
% the variable
% variable- a string that contains the name of the variable of interest
% subjects can be defined as a single number from 1-9 
% trial can be defined as a single trial from 1-33 excluding trials 26-30
% direction can be defined as 1 for data in the sagittal plane, 2 for data in
% the frontal plane, or 3 for data in the transverse plane

%OUTPUT: a single plot. For variables of joint angles, moments, or power for
%the ankle, knee, or hip joints there will be yaxis labels for other
%variables yaxis labels will be blank

load(['p',num2str(subject),'_5stepsdata.mat'],'data')

%% extract heel strikes 
% Extract grf data from data file
grfl = data(trial).grf.force1;
grfr = data(trial).grf.force2;       

% Get heelstrikes
[hsl, tol, hsr, tor] = invDynGrid_getHS_TO(grfl,grfr, 20);
hsr(6) = length(grfl); % often not found with the function

hsl = unique(hsl); hsl(find(diff(hsl)<5)) = [];
hsr = unique(hsr); hsr(find(diff(hsr)<5)) = [];

hs = [hsl hsr];

% divide by ten since mocap was collected at a frequency 10x less
% than grfs were collected
hs_mo = ceil(hs/10); 

%% Define string variables depending on the input variable

if isempty(strfind(variable, 'r_'))==0
    side='Right';
    index=2;
elseif isempty(strfind(variable, 'l_'))==0
    index=1;
    side='Left';
end

if isempty(strfind(variable, 'ank'))==0
    joint='Ankle';
elseif isempty(strfind(variable, 'kne'))==0
    joint='Knee';
elseif isempty(strfind(variable, 'hip'))==0
    joint='Hip';
end

if isempty(strfind(variable, 'angle'))==0
    type='Angle';
    varunit='(Deg)';
elseif isempty(strfind(variable, 'moment'))==0
    type='Moment';
    varunit='(Nm/kg)';
elseif isempty(strfind(variable, 'power'))==0
    type='Power';
    varunit='(W/kg)';
else
    type= ' ';
    varunits= ' ';
end

%% Interpolate the data for the variable of interest to 101 data points
var=interpolate_to_percgaitcycle(data(trial).(field).(variable)(:,direction),hs_mo(:,index),101);

%% Plot each 5 steps and mean of 5 steps for the variable 
figure('name', (strcat(['5 steps and Mean for Subject ', num2str(subject), ' Trial ', num2str(trial),' --- ', side, ' ', joint, ' ', type])));
plot(mean(var,2), 'linewidth', 5, 'color', 'k'); hold on
plot(var, 'linewidth', 2)
xlabel ('% gait cycle')
ylabel (strcat([side, ' ', joint, ' ', type, ' ', varunit]))
title (strcat(['Subject ', num2str(subject), ' Trial ', num2str(trial),': ', side, ' ', joint, ' ', type]));
legend('Mean', 'Step 1', 'Step 2', 'Step 3', 'Step 4', 'Step 5', 'location','best')
grid on
axis tight

end

