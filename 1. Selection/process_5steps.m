% Creates .mat files with 5 steps of good quality data from exported data

%% Settings
subjects = 1:9;
trials = 1:33;

if ~exist('import_folder','var')
    import_folder = uigetdir;
end

if ~exist('export_folder','var')
    export_folder = uigetdir;
end

%% Code
cd(export_folder)

% Loading heel strikes file
if exist('5steps_heelstrikes.mat', 'file')
    load('5steps_heelstrikes.mat','hsl_grf','hsr_grf')
else
    disp('Did not find 5steps_heelstrikes.mat file')
    return
end

% convert to mocap
hsl_grf_mocap = round(hsl_grf/10) + 1;
hsr_grf_mocap = round(hsr_grf/10) + 1;

for subj = subjects
    disp(['Subject: ', num2str(subj)])
    
    % start with what you have (if something exists), else preallocate (required for newer MATLAB versions)    
    cd(export_folder)
    if exist(['p',num2str(subj),'_5StridesData'], 'file')
        load(['p',num2str(subj),'_5StridesData'])
    else, data=[];
    end
    
    % cd to import folder
    cd(import_folder);
    
    disp('Loading data ...')
    load(['p',num2str(subj),'_AllStridesData'], 'data')
    disp('Loading complete')
    
    for trial = trials
        
        if (subj == 6 && trial == 21) || (subj == 6 && trial == 31) || (subj == 7 && trial == 24)
            disp(['Trial number: ', num2str(trial), ' - note: trial is missing (see Supplementary File)'])
            continue
        elseif (subj == 3 && trial == 4) || (subj == 9 && trial == 14) || (subj == 4 && trial == 1)
            disp(['Trial number: ', num2str(trial), ' - note: trial is excluded (see Supplementary File)'])
            continue
        elseif (trial > 25 && trial < 31) 
            disp(['Trial number: ', num2str(trial), ' - note: trial is excluded (see Supplementary File)'])
            continue
        else
            disp(['Trial number: ', num2str(trial)])
        end
       
        %% Select 5 strides
        fn=fieldnames(data(trial));
        selected_data(trial) = data(trial);
        %loop through the fields

        for i=1: numel(fn)
            fn1=fieldnames(data(trial).(fn{i}));

            if strcmp(fn{i}, 'Analog') || strcmp(fn{i}, 'Force')
                hsl = hsl_grf;
                hsr = hsr_grf;
            else
                hsl = hsl_grf_mocap;
                hsr = hsr_grf_mocap;
            end

            for j=1: numel(fn1)
                if isstruct(data(trial).(fn{i}).(fn1{j}))
                    fn2=fieldnames(data(trial).(fn{i}).(fn1{j}));
                    for k=1: numel(fn2)
                        %access the data
                        selected_data(trial).(fn{i}).(fn1{j}).(fn2{k}) = data(trial).(fn{i}).(fn1{j}).(fn2{k})(hsl(subj,trial):hsr(subj,trial),:);     
                    end 
                else
                    if ~isempty(data(trial).(fn{i}).(fn1{j})) && (length(data(trial).(fn{i}).(fn1{j})) > (hsr(subj,trial)-hsl(subj,trial)))
                        selected_data(trial).(fn{i}).(fn1{j}) = data(trial).(fn{i}).(fn1{j})(hsl(subj,trial):hsr(subj,trial),:);  
                    else
                        selected_data(trial).(fn{i}).(fn1{j}) = data(trial).(fn{i}).(fn1{j});
                    end
                end
            end
        end 
          
    end
    
    %% Save
    data = selected_data;
    cd(export_folder)
    disp('Saving data ...')
    save(['p',num2str(subj),'_5StridesData'], 'data')
end



