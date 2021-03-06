function [] = process_5steps(subjects, trials)
%Add datafolder with the PnExported file folders to path before running funtion
%INPUTS: 
    %subjects: specify which subjects to make all strides file, 
    %ex. subjects=[1:3]
    % trials: specify which trials to include in all strides files, exclude
    % missing trials and step width trials
    % ex. trials= [1:25, 31:33] will exculde step width trials
    
%When running function for subject 9, in the target section of the function
%remove lines that involve upper body markers (LAC, RAC, LEP, REP, LWR,
%RWR)

% Loads:
% 1. 5 steps heelstikes
% 2. Raw data that is exported from visual3D
% Saves: 5 steps data stored for each participant in a struct
% Make sure the data and the heelstrike files are in the path and that the
% cd contains the subjects' raw data file folders

load('5steps_heelstrikes.mat','hsl_grf','hsr_grf')
subjnames = {'1' ,'2', '3', '4', '5', '6', '7', '8', '9', '10'};

% convert to mocap
hsl_grf_mocap = round(hsl_grf/10) + 1;
hsr_grf_mocap = round(hsr_grf/10) + 1;

for subj = subjects
    disp(subj)
    data=[];
    % start with what you have
    load(['p',num2str(subj),'_5steps_data'], 'data')
    
    folder=which(['p', num2str(subjnames{subj}), 'export_T03.mat']);
    folder=[folder(1:(end-17))];
    cd(folder)
    ind=strfind(folder, 'R');
    datafolder=folder(1:ind)
    
    for trial = trials
        disp(trial)
        
        if trial<10
        files=dir(['*','_T0', num2str(trial),'.mat']);
        elseif trial>9
        files=dir(['*','export_T', num2str(trial),'.mat']);
        end

        if ~isempty(files)
            load(files.name)
        else, continue
        end

        if ~isempty(files)
            load(files.name)
        else, continue
        end

        
        if ~isnan(hsl_grf(subj,trial)) && ~isempty(force1{1,1})

            % if hip is empty
            if isempty(l_hip_power{1})
                l_hip_power{1} = nan(size(l_kne_power{1}));
                r_hip_power{1} = nan(size(r_kne_power{1}));
            end
           

            
            %% Store Force Platform Data
            data(trial).Platform.ForcePlatformOrigin = [ForcePlatformOrigin{1}];
            data(trial).Platform.ForcePlatformCorners = [ForcePlatformCorners{1}];
            data(trial).Platform.ForcePlatformCalibration = [ForcePlatformCalibration{1}];
            data(trial).Platform.ForcePlatformZero = [ForcePlatformZero{1}];
            data(trial).Platform.ForcePlatformZeros = [ForcePlatformZeros{1}];
            data(trial).Platform.MatrixStore = [MatrixStore{1}];
            
             %% Store Analog Data 
            data(trial).Analog.f1original= [f1xoriginal{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:), f1yoriginal{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:), f1zoriginal{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:)];
            data(trial).Analog.m1original= [f2xoriginal{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:), f2yoriginal{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:), f2zoriginal{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:)];
            data(trial).Analog.f2original= [m1xoriginal{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:), m1yoriginal{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:), m1zoriginal{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:)];
            data(trial).Analog.m2original= [m2xoriginal{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:), m2yoriginal{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:), m2zoriginal{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:)];
            data(trial).Analog.f1processed= [f1xprocessed{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:), f1yprocessed{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:), f1zprocessed{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:)];
            data(trial).Analog.m1processed= [f2xprocessed{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:), f2yprocessed{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:), f2zprocessed{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:)];
            data(trial).Analog.f2processed= [m1xprocessed{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:), m1yprocessed{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:), m1zprocessed{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:)];
            data(trial).Analog.m2processed= [m2xprocessed{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:), m2yprocessed{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:), m2zprocessed{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:)];
            
            
            %% Store Target Data
        
            data(trial).TargetData.L5TH_pos= [L5TH_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
           data(trial).TargetData.LAC_pos= [LAC_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LASI_pos= [LASI_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LCAL_pos= [LCAL_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
           data(trial).TargetData.LEP_pos= [LEP_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LGTR_pos= [LGTR_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LLEP_pos= [LLEP_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LLML_pos= [LLML_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LMEP_pos= [LMEP_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LMML_pos= [LMML_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LSH1_pos= [LSH1_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LSH2_pos= [LSH2_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LSH3_pos= [LSH3_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LTH1_pos= [LTH1_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LTH2_pos= [LTH2_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LTH3_pos= [LTH3_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
           data(trial).TargetData.LWR_pos= [LWR_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            
            data(trial).TargetData.R5TH_pos= [R5TH_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
           data(trial).TargetData.RAC_pos= [RAC_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RASI_pos= [RASI_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RCAL_pos= [RCAL_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
           data(trial).TargetData.REP_pos= [REP_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RGTR_pos= [RGTR_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RLEP_pos= [RLEP_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RLML_pos= [RLML_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RMEP_pos= [RMEP_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RMML_pos= [RMML_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RSH1_pos= [RSH1_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RSH2_pos= [RSH2_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RSH3_pos= [RSH3_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RTH1_pos= [RTH1_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RTH2_pos= [RTH2_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RTH3_pos= [RTH3_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
           data(trial).TargetData.RWR_pos= [RWR_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.SACR_pos= [SACR_pos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            
            data(trial).TargetData.L5TH_pos_proc= [L5TH_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
           data(trial).TargetData.LAC_pos_proc= [LAC_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LASI_pos_proc= [LASI_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LCAL_pos_proc= [LCAL_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
           data(trial).TargetData.LEP_pos_proc= [LEP_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LGTR_pos_proc= [LGTR_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LLEP_pos_proc= [LLEP_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LLML_pos_proc= [LLML_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LMEP_pos_proc= [LMEP_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LMML_pos_proc= [LMML_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LSH1_pos_proc= [LSH1_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LSH2_pos_proc= [LSH2_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LSH3_pos_proc= [LSH3_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LTH1_pos_proc= [LTH1_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LTH2_pos_proc= [LTH2_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.LTH3_pos_proc= [LTH3_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
           data(trial).TargetData.LWR_pos_proc= [LWR_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            
            data(trial).TargetData.R5TH_pos_proc= [R5TH_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
           data(trial).TargetData.RAC_pos_proc= [RAC_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RASI_pos_proc= [RASI_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RCAL_pos_proc= [RCAL_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
           data(trial).TargetData.REP_pos_proc= [REP_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RGTR_pos_proc= [RGTR_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RLEP_pos_proc= [RLEP_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RLML_pos_proc= [RLML_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RMEP_pos_proc= [RMEP_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RMML_pos_proc= [RMML_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RSH1_pos_proc= [RSH1_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RSH2_pos_proc= [RSH2_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RSH3_pos_proc= [RSH3_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RTH1_pos_proc= [RTH1_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RTH2_pos_proc= [RTH2_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.RTH3_pos_proc= [RTH3_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
           data(trial).TargetData.RWR_pos_proc= [RWR_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).TargetData.SACR_pos_proc= [SACR_pos_proc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            
            %% Store Landmark Data
            data(trial).Landmark.HHL= [HHL{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Landmark.HHR= [HHR{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            
            %% Store Time Data 
            data(trial).Time.AnalogTime= [AnalogTime{1}];
            data(trial).Time.FRAMES= [FRAMES{1}];
            data(trial).Time.TIME= [TIME{1}];
            
            %% Store Ground Reaction Forces
            data(trial).Force.force1 = [force1{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:)]; 
            data(trial).Force.force2 = [force2{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:)];
            
            data(trial).Force.cop1 = [cop1{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:)];
            data(trial).Force.cop2 = [cop2{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:)];
            
            data(trial).Force.freemoment1 = [freemoment1{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:)];
            data(trial).Force.freemoment2 = [freemoment2{1}(hsl_grf(subj,trial):hsr_grf(subj,trial),:)];
            
            
            %% Store Kinetic_Kinematic
            data(trial).Kinetic_Kinematic.lFtAngAcc = [lFtAngAcc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.lShAngAcc = [lShAngAcc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.lThAngAcc = [lThAngAcc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rFtAngAcc = [rFtAngAcc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rSkAngAcc = [rSkAngAcc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rThAngAcc = [rThAngAcc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rPvAngAcc = [rPvAngAcc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            
            data(trial).Kinetic_Kinematic.lFtAngVel = [lFtAngVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.lSkAngVel = [lShAngVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.lThAngVel = [lThAngVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rFtAngVel = [rFtAngVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rSkAngVel = [rSkAngVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rThAngVel = [rThAngVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rPvAngVel = [rPvAngVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            
            data(trial).Kinetic_Kinematic.lFtCGAcc = [lFtCGAcc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.lSkCGAcc = [lSkCGAcc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.lThCGAcc = [lThCGAcc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rFtCGAcc = [rFtCGAcc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rSkCGAcc = [rSkCGAcc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rThCGAcc = [rThCGAcc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rPvCGAcc = [rPvCGAcc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            
            data(trial).Kinetic_Kinematic.lFtCGPos = [lFtCGPos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.lSkCGPos = [lSkCGPos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.lThCGPos = [lThCGPos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rFtCGPos = [rFtCGPos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rSkCGPos = [rSkCGPos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rThCGPos = [rThCGPos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rPvCGPos = [rPvCGPos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            
            data(trial).Kinetic_Kinematic.lFtCGVel = [lFtCGVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.lSkCGVel = [lSkCGVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.lThCGVel = [lThCGVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rFtCGVel = [rFtCGVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rSkCGVel = [rSkCGVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rThCGVel = [rThCGVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rPvCGVel = [rPvCGVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            
            data(trial).Kinetic_Kinematic.lFtDistEndPos = [lFtDistEndPos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.lSkDistEndPos = [lShDistEndPos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.lThDistEndPos = [lThDistEndPos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rFtDistEndPos = [rFtDistEndPos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rShDistEndPos = [rShDistEndPos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rThDistEndPos = [rThDistEndPos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rPvDistEndPos = [rPvDistEndPos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            
            data(trial).Kinetic_Kinematic.lFtDistEndVel = [lFtDistEndVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.lSkDistEndVel = [lShDistEndVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.lThDistEndVel = [lThDistEndVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rFtDistEndVel = [rFtDistEndVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rShDistEndVel = [rShDistEndVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rThDistEndVel = [rThDistEndVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rPvDistEndVel = [rPvDistEndVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            
            data(trial).Kinetic_Kinematic.lFtProxEndForce = [lFtProxEndForce{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.lSkProxEndForce = [lShProxEndForce{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.lThProxEndForce = [lThProxEndForce{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rFtProxEndForce = [rFtProxEndForce{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rShProxEndForce = [rShProxEndForce{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rThProxEndForce = [rThProxEndForce{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rPvProxEndForce = [rPvProxEndForce{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            
            data(trial).Kinetic_Kinematic.lFtProxEndPos = [lFtProxEndPos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.lSkProxEndPos = [lShProxEndPos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.lThProxEndPos = [lThProxEndPos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rFtProxEndPos = [rFtProxEndPos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rShProxEndPos = [rShProxEndPos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rThProxEndPos = [rThProxEndPos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rPvProxEndPos = [rPvProxEndPos{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            
            data(trial).Kinetic_Kinematic.lFtProxEndTorque = [lFtProxEndTorque{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.lSkProxEndTorque = [lShProxEndTorque{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.lThProxEndTorque = [lThProxEndTorque{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rFtProxEndTorque = [rFtProxEndTorque{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rShProxEndTorque = [rShProxEndTorque{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rThProxEndTorque = [rThProxEndTorque{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rPvProxEndTorque = [rPvProxEndTorque{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            
            data(trial).Kinetic_Kinematic.lFtProxEndVel = [lFtProxEndVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.lSkProxEndVel = [lShProxEndVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.lThProxEndVel = [lThProxEndVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rFtProxEndVel = [rFtProxEndVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rShProxEndVel = [rShProxEndVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rThProxEndVel = [rThProxEndVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rPvProxEndVel = [rPvProxEndVel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            
            data(trial).Kinetic_Kinematic.lFtSegResidual = [lFtSegResidual{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.lSkSegResidual = [lShSegResidual{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.lThSegResidual = [lThSegResidual{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rFtSegResidual = [rFtSegResidual{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rShSegResidual = [rShSegResidual{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rThSegResidual = [rThSegResidual{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Kinetic_Kinematic.rPvSegResidual = [rPvSegResidual{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            
            %% Store Link_Model_Based
 
            data(trial).Link_Model_Based.l_ank_angle = [l_ank_angle{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.l_hip_angle = [l_hip_angle{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.l_kne_angle = [l_kne_angle{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.r_ank_angle = [r_ank_angle{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.r_hip_angle = [r_hip_angle{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.r_kne_angle = [r_kne_angle{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];

            data(trial).Link_Model_Based.l_ank_moment = [l_ank_moment{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.l_hip_moment = [l_hip_moment{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.l_kne_moment = [l_kne_moment{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.r_ank_moment = [r_ank_moment{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.r_hip_moment = [r_hip_moment{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.r_kne_moment = [r_kne_moment{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];

            data(trial).Link_Model_Based.l_ank_power = [l_ank_power{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.l_hip_power = [l_hip_power{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.l_kne_power = [l_kne_power{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.r_ank_power = [r_ank_power{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.r_hip_power = [r_hip_power{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.r_kne_power = [r_kne_power{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];

            data(trial).Link_Model_Based.l_ank_force_loc = [l_ank_force_loc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.l_hip_force_loc = [l_hip_force_loc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.l_kne_force_loc = [l_kne_force_loc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.r_ank_force_loc = [r_ank_force_loc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.r_hip_force_loc = [r_hip_force_loc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.r_kne_force_loc = [r_kne_force_loc{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];

            data(trial).Link_Model_Based.lft_rotenergy = [lft_rotenergy{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.rft_rotenergy = [rft_rotenergy{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.lsk_rotenergy = [lsk_rotenergy{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.rsk_rotenergy = [rsk_rotenergy{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.lth_rotenergy = [lth_rotenergy{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.rth_rotenergy = [rth_rotenergy{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];

            data(trial).Link_Model_Based.l_ank_vel = [l_ank_vel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.l_kne_vel = [l_kne_vel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.l_hip_vel = [l_hip_vel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.r_ank_vel = [r_ank_vel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.r_kne_vel = [r_kne_vel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.r_hip_vel = [r_hip_vel{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            
            data(trial).Link_Model_Based.lsk_wrt_lank = [lsk_wrt_lank{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.rsk_wrt_rank = [rsk_wrt_rank{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.lth_wrt_lkne = [lth_wrt_lkne{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.rth_wrt_rkne = [rth_wrt_rkne{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.rpv_wrt_lhip = [rpv_wrt_lhip{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];
            data(trial).Link_Model_Based.rpv_wrt_rhip = [rpv_wrt_rhip{1}(hsl_grf_mocap(subj,trial):hsr_grf_mocap(subj,trial),:)];

        else
            data(trial).name = nan;
            data(trial).LandmarkData = nan;
            data(trial).TargetData = nan;
            data(trial).AnalogData = nan;
            data(trial).Force = nan;
            data(trial).Kinetic_Kinematic = nan;
            data(trial).Link_Model_Based = nan;
            

        end
    end

    cd(datafolder)
    save(['p',num2str(subj),'_5StridesData'], 'data')
end
end

