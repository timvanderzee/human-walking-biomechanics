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

%% put trials in order from low to high speed
trialsvel=[1:9,17:33, 10:16]; %trials ordered from lowest speed to highest
for t=1:length(trials)
    index(1,t)=find(trialsvel==trials(t));
    index(2,t)=trials(t);
end
trialsort=sortrows(index');
trials=trialsort(:,2);

%% Load data and define column of interest within the variable of interest

for subject= 1:9 %Loops over all subjects that have 5 step data files
    load ((strcat('p', num2str(subject), '_5StepsData.mat')))
    
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
   
%% Extract grfs for each trial

    for trial=1:length(trials) %Loops over trials
        
        if isempty(data(trials(trial)).grf.force1)==1 
            continue
        elseif isempty(data(trials(trial)).grf.force2)==1 
            continue
        else 
            GRFL=[data(trials(trial)).grf.force1(:,1), data(trials(trial)).grf.force1(:,2), data(trials(trial)).grf.force1(:,3)];
            GRFR=[data(trials(trial)).grf.force2(:,1), data(trials(trial)).grf.force2(:,2), data(trials(trial)).grf.force2(:,3)];
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

%% Take the average iver subjects 

for trial=1:length(trials)
    %subjmean takes the mean of var1 over subjects and therefore only
    %contains mean data over subjects for each trial
    subjmean(:,trial)= mean(var1(:,trial,:), 3);
end

%% Define variables for producing the Figure

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
    janme= 'Hip';
end

% get trialsname based on trials
trialsname=lookup_trial_names(trials); 

%Define colours for plotting
colour=jet(length(trials));
c=[1:length(trials)];

%% Plotting 
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

