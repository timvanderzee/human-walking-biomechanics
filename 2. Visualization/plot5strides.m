function [] = plot5strides(typename, variablename, trial, subject, direction)
% -------------------------------------------------------------------------
% [] = plot5strides(typename, variablename, trial, subject, direction)
% 
% plots 5 strides and the mean of 5 strides for a given motion capture 
% variable in a specified plane of motion
% 
% INPUTS: 
%   * typename - string with name of the field that contains the variable
%   * variable - string that contains the name of the variable of interest
%   * trial    _ tiral number from 1-33 excluding trials 26-30
%   * subject  - subjects number from 1-9 
%   * direction - plane to be plotted: 1 for data in the sagittal plane
%                                      2 for data in the frontal plane
%                                      3 for data in the transverse plane
%
% OUTPUT: single plot for variables of joint angles, moments, or power for
%         the ankle, knee, or hip joints there will be yaxis labels for 
%         other variables yaxis labels will be blank
%
%--------------------------------------------------------------------------

load(['p',num2str(subject),'_5StridesData.mat'],'data')

% extract heel strikes 
% Extract grf data from data file
grfl = data(trial).Force.force1;
grfr = data(trial).Force.force2;       

% Get heelstrikes
[hsl, ~, hsr, ~] = invDynGrid_getHS_TO(grfl,grfr, 20);
hsr(6)           = length(grfl); % often not found with the function

hsl              = unique(hsl); 
hsl(diff(hsl)<5) = [];

hsr              = unique(hsr); 
hsr(diff(hsr)<5) = [];

hs = [hsl hsr];

% divide by ten since mocap was collected at a frequency 10x less
% than grfs were collected
hs_mo = ceil(hs/10); 

%% Define string variables depending on the input variable
if contains(variablename, 'r_') 
    side  = 'Right';
    index = 2;
    
elseif contains(variablename, 'l_') 
    side  = 'Left';
    index = 1;
end

if contains(variablename, 'ank') 
    joint = 'Ankle';
    
elseif contains(variablename, 'kne') 
    joint = 'Knee';
    
elseif contains(variablename, 'hip') 
    joint = 'Hip';
end

if contains(variablename, 'angle') 
    type    = 'Angle';
    varunit = '(Deg)';
    
elseif contains(variablename, 'moment')
    type    = 'Moment';
    varunit = '(Nm/kg)';
    
elseif contains(variablename, 'power')
    type    = 'Power';
    varunit = '(W/kg)';
    
else
    type    = ' ';
    varunit = ' ';
end

% Interpolate the data for the variable of interest to 101 data points
var = interpolate_to_percgaitcycle(...
        data(trial).(typename).(variablename)(:,direction),hs_mo(:,index),...
                                                                      101);

%% Plot figure 
%Plot each 5 strides and mean of 5 strides for the variable 
figure('name', (strcat(['5 strides and Mean for Subject ', num2str(subject),...
        ' Trial ', num2str(trial),' --- ', side, ' ', joint, ' ', type])));
    
hold on   
plot(mean(var,2), 'linewidth', 5, 'color', 'k'); 
plot(var, 'linewidth', 2)

xlabel ('% gait cycle')
ylabel (strcat([side, ' ', joint, ' ', type, ' ', varunit]))

title (strcat(['Subject ', num2str(subject), ' Trial ',...
                       num2str(trial),': ', side, ' ', joint, ' ', type]));

legend('Mean', 'stride 1', 'stride 2', 'stride 3', 'stride 4', 'stride 5',...
                                                         'location','best')
grid on
axis tight

end
