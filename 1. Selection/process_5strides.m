% -------------------------------------------------------------------------
% process_5strides.m
%                       Creates .mat files with 5 strides of good quality 
%                       data from exported data
% -------------------------------------------------------------------------

%% Set Import and Export Directories 
import_folder = fullfile(datafolder,'All Strides Data files');
export_folder = fullfile(datafolder,'5 Strides Data files');

%% Settings
subjects  = 1:9;
no_trials = 33;

if ~exist('import_folder','var')
    import_folder = uigetdir;
end

if ~exist('export_folder','var')
    export_folder = uigetdir;
end

%% Change Directory to Export Path and Load Data for Processing
cd(export_folder)

% ----------------------------------------------- loading heel strikes file
if exist('5strides_heelstrikes.mat', 'file')
    load('5strides_heelstrikes.mat','hsl_grf','hsr_grf')
else
    disp('Did not find 5strides_heelstrikes.mat file')
    return
end

% --------------------- convert to mocap (10 times less sampling frequency)
hsl_grf_mocap = round(hsl_grf/10) + 1;
hsr_grf_mocap = round(hsr_grf/10) + 1;

for subj = subjects
    disp(['Subject: ', num2str(subj)])
    
    % ---- skip if 5 stride data already exists, else load input data array 
    cd(export_folder)
    if exist(['p',num2str(subj),'_5StridesData.mat'],'file')
        load(['p',num2str(subj),'_5StridesData'])
        disp(['5 stride data for subject ' num2str(subj) ...
                                  ' already exists. Skipping processing.'])
        continue 
    end
    
    cd(import_folder);
    disp('Loading data ...')
    load(['p',num2str(subj),'_AllStridesData'], 'data')
    disp('Loading complete')
    
    % -------------- reverse for loop for preallocation of data vector size
    for trial = no_trials:-1:1
        
        % --------------- check if selected trial is excluded from analysis 
        if (subj == 6 && trial == 21) || (subj == 6 && trial == 31) || ...
                                                 (subj == 7 && trial == 24)
            disp(['Trial number: ', num2str(trial), ...
                     ' - note: trial is missing (see Supplementary File)'])
            continue
        elseif (subj == 3 && trial == 4) || (subj == 9 && trial == 14) ||...
                                                  (subj == 4 && trial == 1)
            disp(['Trial number: ', num2str(trial), ...
                    ' - note: trial is excluded (see Supplementary File)'])
            continue
        elseif (trial > 25 && trial < 31) 
            disp(['Trial number: ', num2str(trial), ...
                    ' - note: trial is excluded (see Supplementary File)'])
            continue
        else
            disp(['Trial number: ', num2str(trial)])
        end
       
        %% Select 5 strides
        fn = fieldnames(data(trial));
        selected_data(trial) = data(trial);
        % ----------------------------------------- loop through all fields

        for i = 1:numel(fn)
            fn1 = fieldnames(data(trial).(fn{i}));

            if strcmp(fn{i}, 'Analog') || strcmp(fn{i}, 'Force')
                hsl = hsl_grf;
                hsr = hsr_grf;
            else
                hsl = hsl_grf_mocap;
                hsr = hsr_grf_mocap;
            end

            for j = 1:numel(fn1)
                if isstruct(data(trial).(fn{i}).(fn1{j}))
                    fn2=fieldnames(data(trial).(fn{i}).(fn1{j}));
                    for k=1: numel(fn2)
                        % ------------------------------------- access data
                        selected_data(trial).(fn{i}).(fn1{j}).(fn2{k}) = ...
                        data(trial).(fn{i}).(fn1{j}).(fn2{k})(hsl(subj,trial):hsr(subj,trial),:);     
                    end 
                else
                    if ~isempty(data(trial).(fn{i}).(fn1{j})) && ...
                                (length(data(trial).(fn{i}).(fn1{j})) > (hsr(subj,trial)-hsl(subj,trial)))
                        selected_data(trial).(fn{i}).(fn1{j}) = ...
                            data(trial).(fn{i}).(fn1{j})(hsl(subj,trial):hsr(subj,trial),:);  
                    else
                        selected_data(trial).(fn{i}).(fn1{j}) = data(trial).(fn{i}).(fn1{j});
                    end
                end
            end
        end 
          
    end
    
    %% Save data
    data = selected_data;
    cd(export_folder)
    disp('Saving data ...')
    save(['p',num2str(subj),'_5StridesData'], 'data')
end



