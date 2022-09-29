clear all; close all; clc

% folder where the files are that have been exported from Visual3D
folder = uigetdir;
% import_folder = '';

addpath(genpath(folder))

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
% 
% typename is a string containing the name of the field where the variable
% can be found
% 
% trials is a variable that contains the trial numbers for trials
% matching the experiment defined above. These trials are returned from
% lowest speed to highest speed

[variablename, typename, trials] = lookup_variable_name(variable_type, leg, joint, experiment);

%% calls on functions to look at subjects' means for several trials
% col contains the column of data that will be plotted in the following 2 
% sections
col=variableplot(typename, variablename, trials);
trial=input('What trial are you interested in?  (1-33 excluding step width trials(26-30))');

% trial is a single number for the trial that will be plotted in the
% following 2 sections

%% calls on function to look at each subjects' data for a single trial

visualizejointkin(trial, joint, col)
subj=input('Which subject (from 1-9) are you interested in?');

% subj is a single number for the subject that will be plotted in the
% following section/follwoing 3 figures

%% Variable of interest (5 strides for one subject, one trial)
plot5steps(typename, variablename, trial, subj, col)

%% Ground reaction force (5 strides for one subject, one trial)
GRFplot(subj, trial)

%% 3x3 plot (5 strides for one subject, one trial)
plot3x3(subj, trial)

%%% FUNCTIONS %%%
%% lookup_variable_name Function
function[variablename,typename, trials] = lookup_variable_name(variable_type, leg, joint, experiment)
%This function returns the variable name for Link_Model_Based Angle,
%Moment, and Power variables
%It also returns a list of the trials for a specific experiment type

typename = 'Link_Model_Based';

if strcmp(joint, 'ankle')
    jointn='ank';
elseif strcmp (joint, 'knee')
    jointn='kne';
elseif strcmp(joint, 'hip')
    jointn='hip';
else
    disp 'unable to identify joint name'
end

variablename= strcat(leg(1), '_', jointn, '_', variable_type);

if strcmp(experiment, 'constant step length')
    trials= [1, 4, 7, 33, 12, 15];
elseif strcmp(experiment, 'constant step frequency')
    trials= [3, 6, 9, 31, 10, 13];
elseif strcmp(experiment, 'constant speed')
    trials= [17:19, 23:25]; %20-22 are preferred but at this speed, could include one of them in this array
elseif strcmp(experiment, 'preferred walking')
    trials= [2, 5, 8, 20, 32, 11, 14, 16];
end

end

            
%% variableplot Function
function col=variableplot(field,variable, trials)

%This function works for mocap variables that are angles, moments, or power
%for the ankle, knee, or hip joints. It calculates the mean over subjects
%for each trial specified as an input
%INPUTS: 
%field should equal a string that is the name of the field that contains
%the variable. Ex. 'Link_Model_Based'
%varibale: A string that is equal to the variable name  ex. 'r_kne_power'. 
%It is important that the variable is written the same way that it is in 
%the data file so that this function can call on it.
% trials is equal to any combination of trials from 1-33 that you are
% interested in plotting. Below the trials contained within 1 of the 4
% experiments are listed

% OUTPUT: A plot of the mean over subjects for each trial.

% trials for exp. 1-constant step length [1, 4, 7, 33, 12, 15]   
% exp. 2 constant step freq [3, 6, 9, 10, 13, 31]
% exp. 3 constant speed [17:25]
% exp. 4 perferred speed [2, 5, 8, 11, 14, 16, 20:22, 32]
%datafolder='/Users/emilymundinger/Desktop/Inverse Dynamics Grid/'
%codefolder='/Users/emilymundinger/Desktop/Lab Repos/humanwalkingbiomechanicsdata/'

%put trials in order from low to high speed
trialsvel=[1:9,17:33, 10:16]; %trials ordered from lowest speed to highest
for t=1:length(trials)
    index(1,t)=find(trialsvel==trials(t));
    index(2,t)=trials(t);
end
trialsort=sortrows(index');
trials=trialsort(:,2);

%Load data and define column of interest within the variable of interest

for subject= 1:9 %Loops over all subjects that have 5 strides data files
    data= [];
    load ((strcat('p', num2str(subject), '_5StridesData.mat')))
    
    if subject==1
        varsize=size(data(trials(1)).(field).(variable));
        if varsize(2)>1
            % Asks user which column of data they are interested in
            col=input('Which plane? [1,2,3] (Type 1 for Sagittal plane, 2 for Frontal plane, and 3 for Transverse plane)');
            if col==1
                varplane= 'Sagittal plane';
            elseif col==2
                varplane= 'Frontal plane';
            elseif col==3
                varplane= 'Transverse plane';
            end
            disp(['You chose: ', varplane])
        else
            disp('Invalid choice, defaulted to Sagittal plane')
            col=1;
        end
    end
    
    disp(['Retrieving data for subject: ', num2str(subject)])
   
%Extract grfs for each trial

    for trial=1:length(trials) %Loops over trials
        
        if isempty(data(trials(trial)).Force)==1 
            continue
        else 
            GRFL=[data(trials(trial)).Force.force1(:,1), data(trials(trial)).Force.force1(:,2), data(trials(trial)).Force.force1(:,3)];
            GRFR=[data(trials(trial)).Force.force2(:,1), data(trials(trial)).Force.force2(:,2), data(trials(trial)).Force.force2(:,3)];
        end   
        % Get heelstrikes for the left and right (hsl, hsr)
        [hsl, tol, hsr, tor] = invDynGrid_getHS_TO(GRFL, GRFR,40);
        hsr(6) = length(GRFL); % often not found with the function

        % divide by ten since mocap was collected at a frequency 10x less
        % than grfs were collected
        hsl=ceil(hsl/10);
        hsr=ceil(hsr/10);
        
        hsl = unique(hsl); hsl(find(diff(hsl)<5)) = [];
        hsr = unique(hsr); hsr(find(diff(hsr)<5)) = [];
        
        % determines if we are looking at a varibale on the left or right
        if isempty(strfind(variable, 'l_'))==0
            HS=hsl;
            sidename= 'Left';
        elseif isempty(strfind(variable, 'r_'))==0
            HS=hsr;
            sidename= 'Right';
        end
        %var1 contains the mean of 5 steps for each subject and variable
        var1(:, trial,subject)=mean(interpolate_to_percgaitcycle(data(trials(trial)).(field).(variable)(:,col), HS,101),2,'omitnan');
    end
end

%Take the average over subjects 

for trial=1:length(trials)
    %subjmean takes the mean of var1 over subjects and therefore only
    %contains mean data over subjects for each trial
    subjmean(:,trial)= mean(var1(:,trial,:), 3);
end

%Define variables for producing the Figure

% Define names and axis labels depending on the type of variable
if isempty(strfind(variable, 'angle'))==0
    yname= 'Angle (Degrees)' ;
    tname='Angle';
elseif isempty(strfind(variable, 'moment'))==0
    yname= 'Moment (Nm/kg)' ;
    tname= 'Moment';
elseif isempty(strfind(variable, 'power'))==0
    yname= 'Power (W/kg)' ;
    tname= 'Power';
end

if isempty(strfind(variable, 'ank'))==0
    jname= 'Ankle';
elseif isempty(strfind(variable, 'kne'))==0
    jname= 'Knee';
elseif isempty(strfind(variable, 'hip'))==0
    jname= 'Hip';
end

% get trialsname based on trials
trialsname=lookup_trial_names(trials); 

%Define colours for plotting
colour=jet(length(trials));
c=[1:length(trials)];

%Plotting 
figure('name', strcat ([sidename, ' ', jname, ' ', tname, ' in the ' varplane, ': subject means for ', num2str(length(trials)),' ', 'trials']))
for t=1:length(trials)
    plot (subjmean(:,t), 'LineWidth', 2, 'Color', colour(c(t),:)); hold on;
    xlabel ('% gait cycle');
    ylabel (yname)
    axis tight;
    title ([sidename,' ', jname, ' ', tname, ': subject means']);
    triallegend{t}= strcat(['trial', num2str(trials(t)), ': ', trialsname{t}]);
end
legend (triallegend, 'location', 'best')
xlim([1,100])
grid on

end

%% visualizejointin Function

function[] = visualizejointkin(trial, joint, direction)
% This function produces a plot of each subjects data for a single variable
% for a single trial
%INPUTS: trial- single number for the trial of interest
% joint- can be defined as 'ankle', 'knee', or 'hip'
% direction- can be defined as 1 for data in the sagittal plane, 2 for data in
% the frontal plane, or 3 for data in the transverse plane

% OUTPUT: a figure with 6 subplots showing the joint angles, moment, and
% power for the left leg and the right leg. Each subplot have the mean over
% subjects and each subjects' individual data

% Set current directory (cd) to data folder
subjects = [1:9];
trials = 1:33;
npoints=101; %change this value to change points for interpolation

%Define variables
Ajoint_pc = nan(npoints, subjects(end), 2);
Mjoint_pc = nan(npoints, subjects(end), 2);
Pjoint_pc = nan(npoints, subjects(end), 2);

%Loop over subjects and load subject data 
for subj = subjects
    data=[];
    disp(strcat(['Loading data for Subject:', ' ', num2str(subj)]))
    load(['p',num2str(subj),'_5Stridesdata.mat'],'data')

        
        if isempty(data(trials(trial)).Force)==1
            continue
        else        
%Define joint Moment, Angle, and Power for the given joint
            if strcmp(joint,'ankle')
                Mjoint = [data(trial).Link_Model_Based.l_ank_moment(:,direction) data(trial).Link_Model_Based.r_ank_moment(:,direction)];
                Ajoint = [data(trial).Link_Model_Based.l_ank_angle(:,direction) data(trial).Link_Model_Based.r_ank_angle(:,direction)];
                Pjoint = [data(trial).Link_Model_Based.l_ank_power(:,direction) data(trial).Link_Model_Based.r_ank_power(:,direction)];
                
            elseif strcmp(joint,'knee')
                Mjoint = [data(trial).Link_Model_Based.l_kne_moment(:,direction) data(trial).Link_Model_Based.r_kne_moment(:,direction)];
                Ajoint = [data(trial).Link_Model_Based.l_kne_angle(:,direction) data(trial).Link_Model_Based.r_kne_angle(:,direction)];
                Pjoint = [data(trial).Link_Model_Based.l_kne_power(:,direction) data(trial).Link_Model_Based.r_kne_power(:,direction)];
                                
            elseif strcmp(joint,'hip')
                Mjoint = [data(trial).Link_Model_Based.l_hip_moment(:,direction) data(trial).Link_Model_Based.r_hip_moment(:,direction)];
                Ajoint = [data(trial).Link_Model_Based.l_hip_angle(:,direction) data(trial).Link_Model_Based.r_hip_angle(:,direction)];
                Pjoint = [data(trial).Link_Model_Based.l_hip_power(:,direction) data(trial).Link_Model_Based.r_hip_power(:,direction)];
            else
                disp('no joint defined')
            end
            
            %% extract heel strikes
            % Extract GRF variable
            GRFL=[data(trials(trial)).Force.force1(:,1), data(trials(trial)).Force.force1(:,2), data(trials(trial)).Force.force1(:,3)];
            GRFR=[data(trials(trial)).Force.force2(:,1), data(trials(trial)).Force.force2(:,2), data(trials(trial)).Force.force2(:,3)];

            % Get heelstrikes
            [hsl, tol, hsr, tor] = invDynGrid_getHS_TO(GRFL, GRFR,40);
            hsr(6) = length(GRFL); % often not found with the function
            
            % divide by ten since mocap was collected at a frequency 10x less
            % than grfs were collected
            hsl=ceil(hsl/10);
            hsr=ceil(hsr/10);
            
            hsl = unique(hsl); hsl(find(diff(hsl)<5)) = [];
            hsr = unique(hsr); hsr(find(diff(hsr)<5)) = [];
            
            hs = [hsl hsr];
            %% interpolate to percentage gait cycle (101 points)
            if sum(isnan(Mjoint(:)))<1 && sum(isnan(Pjoint(:)))<1
                for L = 1:2
                    Ajoint_pc(:,subj,L) = mean(interpolate_to_percgaitcycle(Ajoint(:,L),hs(:,L),npoints),2,'omitnan');
                    Mjoint_pc(:,subj,L) = mean(interpolate_to_percgaitcycle(Mjoint(:,L),hs(:,L),npoints),2,'omitnan');
                    Pjoint_pc(:,subj,L) = mean(interpolate_to_percgaitcycle(Pjoint(:,L),hs(:,L),npoints),2,'omitnan');
                end
            end
        end
% Take the mean over the subjects
Ajoint_pcmean=mean(Ajoint_pc, 2, 'omitnan');
Mjoint_pcmean=mean(Mjoint_pc, 2, 'omitnan');
Pjoint_pcmean=mean(Pjoint_pc, 2, 'omitnan');
        
end

%Variables for Plotting
% colours for plotting different subjects' data
colour=jet(11);
c=[2,4:12];

% npoints can be changed to change the nunber of points for interpolation.
% xa is set below so that the x axis is still from 1-100% gait cycle
xa=linspace(0,npoints,101);

% define string for the plane of data being plotted
if direction==1
    varplane= 'Sagittal plane';
elseif direction==2
    varplane= 'Frontal plane';
elseif direction==3
    varplane= 'Transverse plane';
end

% Making the Plot
figure('name',['Rotational Joint kin in the ' varplane, ' for trial: ', num2str(trial)],'units','normalized','outerposition',[0 0 .5 1])

% Plot the means of the subjects
subplot(321)
plot(xa, Ajoint_pcmean(:,1), 'linewidth', 2, 'color', 'k'); hold on
subplot(322)
plot(xa, Ajoint_pcmean(:,2), 'linewidth', 2, 'color', 'k'); hold on
subplot(323)
plot(xa, Mjoint_pcmean(:,1), 'linewidth', 2, 'color', 'k'); hold on
subplot(324)
plot(xa, Mjoint_pcmean(:,2), 'linewidth', 2, 'color', 'k'); hold on
subplot(325)
plot(xa, Pjoint_pcmean(:,1), 'linewidth', 2, 'color', 'k'); hold on
subplot(326)
plot(xa, Pjoint_pcmean(:,2), 'linewidth', 2, 'color', 'k'); hold on


%Plot data for each individual subject
for s=1:length(subjects);
    if isempty(Mjoint_pc(:,s,1))==1
        continue
    else

        % left joint angle
        subplot(321)
        plot(xa, Ajoint_pc(:,s,1), 'linewidth', 1.5, 'color', colour(c(s),:)); hold on
        title(strcat(['Left ', joint, ' Angle']))
        axis tight; yl1(1,:) = get(gca,'ylim');
        ylabel('Angle (deg)');

        % right joint angle
        subplot(322)
        plot(xa, Ajoint_pc(:,s,2), 'linewidth', 1.5,'color', colour(c(s),:)); hold on
        title(strcat(['Right ', joint, ' Angle']))
        axis tight; yl1(2,:) = get(gca,'ylim');
        ylabel(' Angle (deg)');
        
        % left joint moment
        subplot(323)
        plot(xa, Mjoint_pc(:,s, 1),'linewidth', 1.5,'color', colour(c(s),:)); hold on
        title(strcat(['Left ', joint, ' Moment']))
        axis tight; yl2(1,:) = get(gca,'ylim');
        ylabel('Moment (N-m/kg)');

        % right joint moment
        subplot(324)
        plot(xa, Mjoint_pc(:,s,2), 'linewidth', 1.5, 'color', colour(c(s),:)); hold on
        title(strcat(['Right ', joint, ' Moment']))
        axis tight; yl2(2,:) = get(gca,'ylim');
        ylabel('Moment (N-m/kg)');

        % left joint power
        subplot(325)
        plot(xa, Pjoint_pc(:,s,1), 'linewidth', 1.5, 'color', colour(c(s),:)); hold on
        title(strcat(['Left ', joint, ' Power']))
        axis tight; yl3(1,:) = get(gca,'ylim');
        ylabel('Power (W/kg)')

        % right joint power
        subplot(326)
        plot(xa, Pjoint_pc(:,s,2), 'linewidth', 1.5, 'color', colour(c(s),:)); hold on
        title(strcat(['Right ', joint, ' Power']))
        axis tight; yl3(2,:) = get(gca,'ylim');
        ylabel('Power (W/kg)');
        
    end
end
        

% make nice
        for i = 1:6
            subplot(3,2,i)
            xlim([0 npoints-1])
            xlabel('Percentage gait cycle');
            grid on
        end
        
        for i=1:2
            subplot(3,2,i)
            ylim([min(yl1(:,1)), max(yl1(:,2))]);
        end
        
        for i=3:4
            subplot(3,2,i)
            ylim([min(yl2(:,1)), max(yl2(:,2))]);
        end
        
        for i=5:6
            subplot(3,2,i)
            ylim([min(yl3(:,1)), max(yl3(:,2))]);
        end
        
        for i=1:4
            subplot(3,2,i)
            yyaxis right
            set(gca, 'YTick', [], 'YColor', 'k')
            ylabel('<---- Ext. ------------ Flex. ---->');
        end

legend('Mean','1','2','3','4','5','6','7','8','9','location','best')


end


%% plot5steps Function

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

data=[];
load(['p',num2str(subject),'_5StridesData.mat'],'data')

% extract heel strikes 
% Extract grf data from data file
grfl = data(trial).Force.force1;
grfr = data(trial).Force.force2;       

% Get heelstrikes
[hsl, tol, hsr, tor] = invDynGrid_getHS_TO(grfl,grfr, 20);
hsr(6) = length(grfl); % often not found with the function

hsl = unique(hsl); hsl(find(diff(hsl)<5)) = [];
hsr = unique(hsr); hsr(find(diff(hsr)<5)) = [];

hs = [hsl hsr];

% divide by ten since mocap was collected at a frequency 10x less
% than grfs were collected
hs_mo = ceil(hs/10); 

% Define string variables depending on the input variable

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

% Interpolate the data for the variable of interest to 101 data points
var=interpolate_to_percgaitcycle(data(trial).(field).(variable)(:,direction),hs_mo(:,index),101);

%Plot each 5 steps and mean of 5 steps for the variable 
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

%% GRFplot Function

function GRFplot(subjects, trials)

% Figure has 3 subplots for GRF. It plots the left and right
% leg one on top of the other and the 3 subplots are the GRF in the
% medial/lateral, anterior/posterior, and vertical directions.
 
% Subjects can be a single number from 1 to 9 or an array including any
% combinations of subjects from 1 to 9
% Trials can be a single trial number from 1:33 or an array of trials
% Ouputs: For each subject/trial combination, a single plot will appear,


% 5 steps analysis

%loops over subject and loads subject's data and then loops over trials
for subj = subjects
    data=[];
    load(['p',num2str(subj),'_5StridesData'])
    for trial = trials
% Figure: 5 step figure: ground reaction forces

        figure('name', (['GRF for Subject ', num2str(subj),', Trial: ', num2str(trial)]))
        
        title (['GRF for Subject ', num2str(subj),', trial ', num2str(trial)])
        if isnan(data(trial).Force.force1)==1
            subplot(3,1,1)
            plot (0,0)
            subplot(3,1,2)
            subplot (3,1,3)
        elseif isnan(data(trial).Force.force1)==0
            subplot(3,1,1) %plot for medial/lateral GRFs
            plot(data(trial).Force.force1(:,1)); hold on;
            plot(data(trial).Force.force2(:,1)); hold on;
            ylabel ('Force (N)')
            title('Medial/Lateral')
            
            subplot(3,1,2) %plot for anterior/posterior GRFs
            plot(data(trial).Force.force1(:,2)); hold on;
            plot(data(trial).Force.force2(:,2)); hold on;
            ylabel ('Force (N)')
            title('Anterior/Posterior')
            
            subplot(3,1,3) %plot for vertical GRFs
            plot(data(trial).Force.force1(:,3)); hold on;
            plot(data(trial).Force.force2(:,3)); hold on;
            ylabel ('Force (N)')
            title ('Vertical')
        end
        for k=1:3
            subplot(3,1,k)
            xlabel('Frame Number')
        end
        legend ('Left Leg', 'Right Leg', 'location', 'best')
    end
end
        
end

%% lookup_trial_names Function

function trialnames=lookup_trial_names(trials)
%returns a cell array with the trial descriptions for each trial in the
%input variable: 'trials'

names={'0.7 m/s constant step length'
'0.7 m/s preferred'
'0.7 m/s constant step frequency'
'0.9 m/s constant step length'
'0.9 m/s preferred'
'0.9 m/s constant step frequency'
'1.1 m/s constant step length'
'1.1 m/s preferred'
'1.1 m/s constant step frequency'
'1.6 m/s constant step frequency'
'1.6 m/s preferred'
'1.6 m/s constant step length'
'1.8 m/s constant step frequency'
'1.8 m/s preferred'
'1.8 m/s constant step length'
'2.0 m/s preferred'
'1.25 m/s lowest step frequency'
'1.25 m/s lower step frequency'
'1.25 m/s low step frequency'
'1.25 m/s preferred'
'1.25 m/s preferred'
'1.25 m/s preferred'
'1.25 m/s high step frequency'
'1.25 m/s higher step frequency'
'1.25 m/s highest step frequency'
'1.25 m/s zero step width'
'1.25 m/s 10 cm step width'
'1.25 m/s 20 cm step width'
'1.25 m/s 30 cm step width'
'1.25 m/s 40 cm step width'
'1.40 m/s constant step length'
'1.40 m/s preferred'
'1.40 m/s constant step frequency'};

for i=1:length(trials)
    trialnames{i} = names{trials(i)};
end
end

%% plot3x3 Function

function plot3x3(s,t)

%INPUTS: s- subject, t- trial
% input single number not an array of numbers
%OUTPUT: a 3x3 figure of the ankle, knee, and hip angle, moment, and power
%for the left (blue) and right(red) legs avereaged over 5 steps

trial=t;
subj=s;

data=[];
load(strcat('p', num2str(subj), '_5StridesData.mat'))
    
    %% Extract grfs for each trial
  
if isempty(data(trial).Force)==1 
   disp('data.Force is empty')
   keyboard
else
    GRFL=[data(trial).Force.force1(:,1), data(trial).Force.force1(:,2), data(trial).Force.force1(:,3)];
    GRFR=[data(trial).Force.force2(:,1), data(trial).Force.force2(:,2), data(trial).Force.force2(:,3)];
end   
% Get heelstrikes for the left and right (hsl, hsr)
[hsl, tol, hsr, tor] = invDynGrid_getHS_TO(GRFL, GRFR,40);
hsr(6) = length(GRFL); % often not found with the function

% divide by ten since mocap was collected at a frequency 10x less
% than grfs were collected
hsl=ceil(hsl/10);
hsr=ceil(hsr/10);

hsl = unique(hsl); hsl(find(diff(hsl)<5)) = [];
hsr = unique(hsr); hsr(find(diff(hsr)<5)) = [];

%Angle and Moment data (x) power (sum of x, y, z) 
%Ajoint(i,j) i=1 ankle, i=2 knee, i=3 hip, j=1 left, j=2 right
Ajoint(1,1,:)=-(mean(interpolate_to_percgaitcycle(data(trial).Link_Model_Based.l_ank_angle(:,1), hsl, 201),2,'omitnan'));
Ajoint(1,2,:)=-(mean(interpolate_to_percgaitcycle(data(trial).Link_Model_Based.r_ank_angle(:,1), hsr, 201),2,'omitnan'));
Mjoint(1,1,:)=-(mean(interpolate_to_percgaitcycle(data(trial).Link_Model_Based.l_ank_moment(:,1), hsl, 201),2,'omitnan'));
Mjoint(1,2,:)=-(mean(interpolate_to_percgaitcycle(data(trial).Link_Model_Based.r_ank_moment(:,1), hsr, 201),2,'omitnan'));
Pjoint(1,1,:)=mean(interpolate_to_percgaitcycle(sum(data(trial).Link_Model_Based.l_ank_power,2), hsl, 201),2,'omitnan');
Pjoint(1,2,:)=mean(interpolate_to_percgaitcycle(sum(data(trial).Link_Model_Based.r_ank_power,2), hsr, 201),2,'omitnan');

Ajoint(2,1,:)=mean(interpolate_to_percgaitcycle(data(trial).Link_Model_Based.l_kne_angle(:,1), hsl, 201),2,'omitnan');
Ajoint(2,2,:)=mean(interpolate_to_percgaitcycle(data(trial).Link_Model_Based.r_kne_angle(:,1), hsr, 201),2,'omitnan');
Mjoint(2,1,:)=mean(interpolate_to_percgaitcycle(data(trial).Link_Model_Based.l_kne_moment(:,1), hsl, 201),2,'omitnan');
Mjoint(2,2,:)=mean(interpolate_to_percgaitcycle(data(trial).Link_Model_Based.r_kne_moment(:,1), hsr, 201),2,'omitnan');
Pjoint(2,1,:)=mean(interpolate_to_percgaitcycle(sum(data(trial).Link_Model_Based.l_kne_power,2), hsl, 201),2,'omitnan');
Pjoint(2,2,:)=mean(interpolate_to_percgaitcycle(sum(data(trial).Link_Model_Based.r_kne_power,2), hsr, 201),2,'omitnan');

Ajoint(3,1,:)=-(mean(interpolate_to_percgaitcycle(data(trial).Link_Model_Based.l_hip_angle(:,1), hsl, 201),2,'omitnan'));
Ajoint(3,2,:)=-(mean(interpolate_to_percgaitcycle(data(trial).Link_Model_Based.r_hip_angle(:,1), hsr, 201),2,'omitnan'));
Mjoint(3,1,:)=-(mean(interpolate_to_percgaitcycle(data(trial).Link_Model_Based.l_hip_moment(:,1), hsl, 201),2,'omitnan'));
Mjoint(3,2,:)=-(mean(interpolate_to_percgaitcycle(data(trial).Link_Model_Based.r_hip_moment(:,1), hsr, 201),2,'omitnan'));
Pjoint(3,1,:)=mean(interpolate_to_percgaitcycle(sum(data(trial).Link_Model_Based.l_hip_power,2), hsl, 201),2,'omitnan');
Pjoint(3,2,:)=mean(interpolate_to_percgaitcycle(sum(data(trial).Link_Model_Based.r_hip_power,2), hsr, 201),2,'omitnan');

figure('name', (['3x3 plot for Subject ', num2str(s),', Trial: ', num2str(t)]))

for i = 1:9
    subplot(3,3,i);
    plot([0 100], [0 0], 'k-'); hold on
    xlabel('% Gait cycle'); box off
end

%Change y limits

for i=1:3
subplot(3,3,i);
plot(linspace(0,100,201), squeeze(Ajoint(i,:,:)),'linewidth',2);

ylim([-70 25]); ylabel('Angle (deg)')

subplot(3,3,3+i);
plot(linspace(0,100,201), squeeze(Mjoint(i,:,:)),'linewidth',2);

ylim([-1 2]); ylabel('Moment (N-m/kg)')

subplot(3,3,6+i);
plot(linspace(0,100,201), squeeze(Pjoint(i,:,:)),'linewidth',2);

ylim([-2 3]); ylabel('Power (W/kg)')
end

subplot(331); title('Ankle');
subplot(332); title('Knee');
subplot(333); title('Hip');

end
