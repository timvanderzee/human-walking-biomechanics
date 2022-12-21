function col = variableplot(field, variable, trials)
% -------------------------------------------------------------------------
% col = variableplot(field, variable, trials)
%
% works for mocap variables that are angles, moments, or power for ankle, 
% knee, or hip joints. It calculates the mean over subjects for each trial 
% specified as an input
% 
% INPUTS: 
%   - field: string of the name of the field that contains the variable.
%                                               Example: 'Link_Model_Based'
%
%   - varibale: A string that is equal to the variable name. 
%                                                    Example: 'r_kne_power' 
%               Note: the variable has to be written the same way that it
%                 is in the data file so that this function can call on it.
%
%   - trials: any combination of trials from 1-33 you are interested in 
%             plotting. Below the trials contained within 1 of the 4
%             experiments are listed.
%
% OUTPUT: 
%   - col: plane that is plotted [1,2,3] with 1 = sagittal plane
%                                             2 = frontal plane
%                                             3 = transverse plane
% 
%   - Plot of the mean over subjects for each trial.
% 
% TRIAL LIST:
%            * trials for exp. 1-constant step length [1, 4, 7, 33, 12, 15]   
%            * exp. 2 constant step freq [3, 6, 9, 10, 13, 31]
%            * exp. 3 constant speed [17:25]
%            * exp. 4 perferred speed [2, 5, 8, 11, 14, 16, 20:22, 32]
%
% -------------------------------------------------------------------------

%% Init vars
no_subjects      = 9;   % all subjects that have 5 strides data files
no_interp_points = 101;

index       = zeros(2, length(trials)); 
var1        = zeros(no_interp_points, length(trials), no_subjects);
subjmean    = zeros(length(var1),length(trials));
triallegend = cell(1, length(trials));

%% Extract Data
% ------------------------------ Put trials in order from low to high speed
trialsvel = [1:9,17:33, 10:16]; 

for t = 1:length(trials)
    
    index(1,t) = find(trialsvel==trials(t));
    index(2,t) = trials(t);
    
end
trialsort = sortrows(index');
trials    = trialsort(:,2);


% ------------------------- For all subjects that have 5 strides data files
for subject = 1:no_subjects 
    
    % --------- Load data and define column of interest within the variable
    load ((strcat('p', num2str(subject), '_5StridesData.mat')),'data')
    
    % For first subj: Asks user which column of data they are interested in
    if subject == 1 
        varsize = size(data(trials(1)).(field).(variable));
        
        if varsize(2) > 1
            col = input(['Which plane? do you want to visualize? [1,2,3]'...
                         '(1 = Sagittal plane,'...
                         ' 2 = Frontal plane,'...
                         ' 3 = Transverse plane)']);       
         switch col
             case 1
                 varplane = 'Sagittal plane';

             case 2
                 varplane = 'Frontal plane';

             case 3
                 varplane = 'Transverse plane';

             otherwise
                warning(['variableplot: Invalid plane choice.' ...
                                    ' Set to default 1 = sagittal plane.'])
                col = 1;    
         end
         disp(['Your choice: ', varplane])
        end
    end
    
    disp(['Retrieving data for subject: ', num2str(subject)])
   
% --------------------------------------------- Extract grfs for all trials
    for trial = 1:length(trials)
        
        if isempty(data(trials(trial)).Force) == 1 
            continue
            
        else 
            GRFL = [data(trials(trial)).Force.force1(:,1), ...
                        data(trials(trial)).Force.force1(:,2), ...
                            data(trials(trial)).Force.force1(:,3)];
                
            GRFR = [data(trials(trial)).Force.force2(:,1),...
                        data(trials(trial)).Force.force2(:,2),...
                            data(trials(trial)).Force.force2(:,3)];
        end   
        
        % Get heelstrikes for the left and right (hsl, hsr)
        [hsl, ~, hsr, ~] = invDynGrid_getHS_TO(GRFL, GRFR,40);
        hsr(6) = length(GRFL); % often not found with the function

        % divide by ten since mocap was collected at a frequency 10x less
        % than grfs were collected
        hsl = ceil(hsl/10);
        hsr = ceil(hsr/10);
        
        hsl              = unique(hsl);
        hsl(diff(hsl)<5) = [];
        
        hsr              = unique(hsr); 
        hsr(diff(hsr)<5) = [];
        
        % determines if we are looking at left or right foot
        if contains(variable, 'l_') 
            
            HS       = hsl;
            sidename = 'Left';
            
        elseif contains(variable, 'r_') 
            HS       = hsr;
            sidename = 'Right';
        end
        
        %var1 contains the mean of 5 strides for each subject and variable
        interp = interpolate_to_percgaitcycle(...
                        data(trials(trial)).(field).(variable)(:,col), ...
                                                      HS ,no_interp_points);
        var1(:, trial, subject) = mean(interp,2,'omitnan');
    end
    
end

% Take the average over subjects 
for trial = 1:length(trials)
    % subjmean takes the mean of var1 over subjects and therefore only
    % contains mean data over subjects for each trial
    subjmean(:,trial) = mean(var1(:,trial,:), 3);
end

%% Define figure variables

% Define names and axis labels depending on the type of variable
if contains(variable, 'angle') 
    yname = 'Angle (Degrees)' ;
    tname = 'Angle';
    
    elseif contains(variable, 'moment') 
        yname = 'Moment (Nm/kg)' ;
        tname = 'Moment';
    
        elseif contains(variable, 'power') 
            yname = 'Power (W/kg)' ;
            tname = 'Power';
end

if contains(variable, 'ank') 
    jname= 'Ankle';
    
    elseif contains(variable, 'kne') 
        jname= 'Knee';
        
        elseif contains(variable, 'hip') 
            jname= 'Hip';
end

% get trialsname
trialsname = lookup_trial_names(trials); 

% define colours for plotting
colour  = jet(length(trials));
c       = 1:length(trials);

%% Plotting 
figure('name', strcat ([sidename, ' ', jname, ' ', tname, ' in the ' ...
                         varplane, ': subject means for ',...
                                  num2str(length(trials)),' ', 'trials']));

for t = 1:length(trials)
    
    plot (subjmean(:,t), 'LineWidth', 2, 'Color', colour(c(t),:))
    hold on
    title ([sidename,' ', jname, ' ', tname, ': subject means'])
    
    xlabel ('% gait cycle')
    ylabel (yname)
    axis tight
    
 triallegend{t} = strcat(['trial',num2str(trials(t)), ': ',trialsname{t}]);
 
end

legend (triallegend, 'location', 'best')
xlim([1,100])
grid on
end
