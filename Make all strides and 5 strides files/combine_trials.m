function combine_trials(datafolder, subjects, trials)

addpath(genpath(datafolder))
subjnames = {'1' ,'2', '3', '4', '5', '6', '7', '8', '9', '10'};


for subj = subjects
    disp(subj)
    
    % go to subject folder
    subjname = char(subjnames(subj));
    cd(datafolder)
    
    
    folder=(strcat(datafolder,'\P ', subjnames(subj),'exportedfiles'))
    cd(folder{1,1})
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
        
        if  ~isempty(force1)
        
            %% Store Force Platform Data
            data(trial).Platform.ForcePlatformOrigin = [ForcePlatformOrigin{1}];
            data(trial).Platform.ForcePlatformCorners = [ForcePlatformCorners{1}];
            data(trial).Platform.ForcePlatformCalibration = [ForcePlatformCalibration{1}];
            data(trial).Platform.ForcePlatformZero = [ForcePlatformZero{1}];
            data(trial).Platform.ForcePlatformZeros = [ForcePlatformZeros{1}];
            data(trial).Platform.MatrixStore = [MatrixStore{1}];
          
            
            %% Store Analog Data 
            data(trial).Analog.f1original= [f1xoriginal{1}, f1yoriginal{1}, f1zoriginal{1}];
            data(trial).Analog.m1original= [f2xoriginal{1}, f2yoriginal{1}, f2zoriginal{1}];
            data(trial).Analog.f2original= [m1xoriginal{1}, m1yoriginal{1}, m1zoriginal{1}];
            data(trial).Analog.m2original= [m2xoriginal{1}, m2yoriginal{1}, m2zoriginal{1}];
            
            data(trial).Analog.f1processed= [f1xprocessed{1}, f1yprocessed{1}, f1zprocessed{1}];
            data(trial).Analog.m1processed= [f2xprocessed{1}, f2yprocessed{1}, f2zprocessed{1}];
            data(trial).Analog.f2processed= [m1xprocessed{1}, m1yprocessed{1}, m1zprocessed{1}];
            data(trial).Analog.m2processed= [m2xprocessed{1}, m2yprocessed{1}, m2zprocessed{1}];
            
            
            %% Store Target Data
            data(trial).TargetData.L5TH_pos= [L5TH_pos{1}];
           data(trial).TargetData.LAC_pos= [LAC_pos{1}];
            data(trial).TargetData.LASI_pos= [LASI_pos{1}];
            data(trial).TargetData.LCAL_pos= [LCAL_pos{1}];
           data(trial).TargetData.LEP_pos= [LEP_pos{1}];
            data(trial).TargetData.LGTR_pos= [LGTR_pos{1}];
            data(trial).TargetData.LLEP_pos= [LLEP_pos{1}];
            data(trial).TargetData.LLML_pos= [LLML_pos{1}];
            data(trial).TargetData.LMEP_pos= [LMEP_pos{1}];
            data(trial).TargetData.LMML_pos= [LMML_pos{1}];
            data(trial).TargetData.LSH1_pos= [LSH1_pos{1}];
            data(trial).TargetData.LSH2_pos= [LSH2_pos{1}];
            data(trial).TargetData.LSH3_pos= [LSH3_pos{1}];
            data(trial).TargetData.LTH1_pos= [LTH1_pos{1}];
            data(trial).TargetData.LTH2_pos= [LTH2_pos{1}];
            data(trial).TargetData.LTH3_pos= [LTH3_pos{1}];
           data(trial).TargetData.LWR_pos= [LWR_pos{1}];
            
            data(trial).TargetData.R5TH_pos= [R5TH_pos{1}];
           data(trial).TargetData.RAC_pos= [RAC_pos{1}];
            data(trial).TargetData.RASI_pos= [RASI_pos{1}];
            data(trial).TargetData.RCAL_pos= [RCAL_pos{1}];
           data(trial).TargetData.REP_pos= [REP_pos{1}];
            data(trial).TargetData.RGTR_pos= [RGTR_pos{1}];
            data(trial).TargetData.RLEP_pos= [RLEP_pos{1}];
            data(trial).TargetData.RLML_pos= [RLML_pos{1}];
            data(trial).TargetData.RMEP_pos= [RMEP_pos{1}];
            data(trial).TargetData.RMML_pos= [RMML_pos{1}];
            data(trial).TargetData.RSH1_pos= [RSH1_pos{1}];
            data(trial).TargetData.RSH2_pos= [RSH2_pos{1}];
            data(trial).TargetData.RSH3_pos= [RSH3_pos{1}];
            data(trial).TargetData.RTH1_pos= [RTH1_pos{1}];
            data(trial).TargetData.RTH2_pos= [RTH2_pos{1}];
            data(trial).TargetData.RTH3_pos= [RTH3_pos{1}];
           data(trial).TargetData.RWR_pos= [RWR_pos{1}];
            data(trial).TargetData.SACR_pos= [SACR_pos{1}];
                        
            data(trial).TargetData.L5TH_pos_proc= [L5TH_pos_proc{1}];
           data(trial).TargetData.LAC_pos_proc= [LAC_pos_proc{1}];
            data(trial).TargetData.LASI_pos_proc= [LASI_pos_proc{1}];
            data(trial).TargetData.LCAL_pos_proc= [LCAL_pos_proc{1}];
           data(trial).TargetData.LEP_pos_proc= [LEP_pos_proc{1}];
            data(trial).TargetData.LGTR_pos_proc= [LGTR_pos_proc{1}];
            data(trial).TargetData.LLEP_pos_proc= [LLEP_pos_proc{1}];
            data(trial).TargetData.LLML_pos_proc= [LLML_pos_proc{1}];
            data(trial).TargetData.LMEP_pos_proc= [LMEP_pos_proc{1}];
            data(trial).TargetData.LMML_pos_proc= [LMML_pos_proc{1}];
            data(trial).TargetData.LSH1_pos_proc= [LSH1_pos_proc{1}];
            data(trial).TargetData.LSH2_pos_proc= [LSH2_pos_proc{1}];
            data(trial).TargetData.LSH3_pos_proc= [LSH3_pos_proc{1}];
            data(trial).TargetData.LTH1_pos_proc= [LTH1_pos_proc{1}];
            data(trial).TargetData.LTH2_pos_proc= [LTH2_pos_proc{1}];
            data(trial).TargetData.LTH3_pos_proc= [LTH3_pos_proc{1}];
           data(trial).TargetData.LWR_pos_proc= [LWR_pos_proc{1}];
            
            data(trial).TargetData.R5TH_pos_proc= [R5TH_pos_proc{1}];
           data(trial).TargetData.RAC_pos_proc= [RAC_pos_proc{1}];
            data(trial).TargetData.RASI_pos_proc= [RASI_pos_proc{1}];
            data(trial).TargetData.RCAL_pos_proc= [RCAL_pos_proc{1}];
           data(trial).TargetData.REP_pos_proc= [REP_pos_proc{1}];
            data(trial).TargetData.RGTR_pos_proc= [RGTR_pos_proc{1}];
            data(trial).TargetData.RLEP_pos_proc= [RLEP_pos_proc{1}];
            data(trial).TargetData.RLML_pos_proc= [RLML_pos_proc{1}];
            data(trial).TargetData.RMEP_pos_proc= [RMEP_pos_proc{1}];
            data(trial).TargetData.RMML_pos_proc= [RMML_pos_proc{1}];
            data(trial).TargetData.RSH1_pos_proc= [RSH1_pos_proc{1}];
            data(trial).TargetData.RSH2_pos_proc= [RSH2_pos_proc{1}];
            data(trial).TargetData.RSH3_pos_proc= [RSH3_pos_proc{1}];
            data(trial).TargetData.RTH1_pos_proc= [RTH1_pos_proc{1}];
            data(trial).TargetData.RTH2_pos_proc= [RTH2_pos_proc{1}];
            data(trial).TargetData.RTH3_pos_proc= [RTH3_pos_proc{1}];
           data(trial).TargetData.RWR_pos_proc= [RWR_pos_proc{1}];
            data(trial).TargetData.SACR_pos_proc= [SACR_pos_proc{1}];
            
            %% Store Landmark Data
            data(trial).Landmark.HHL= [HHL{1}];
            data(trial).Landmark.HHR= [HHR{1}];
            
            %% Store Time Data 
            data(trial).Time.AnalogTime= [AnalogTime{1}];
            data(trial).Time.FRAMES= [FRAMES{1}];
            data(trial).Time.TIME= [TIME{1}];
           
            
            %% Store Ground Reaction Forces
            data(trial).Force.force1 = [force1{1}]; 
            data(trial).Force.force2 = [force2{1}];
            
            data(trial).Force.cop1 = [cop1{1}];
            data(trial).Force.cop2 = [cop2{1}];
            
            data(trial).Force.freemoment1 = [freemoment1{1}];
            data(trial).Force.freemoment2 = [freemoment2{1}];
            
            %% Store Kinetic_Kinematic
            data(trial).Kinetic_Kinematic.lFtAngAcc = [lFtAngAcc{1}];
            data(trial).Kinetic_Kinematic.lSkAngAcc = [lShAngAcc{1}];
            data(trial).Kinetic_Kinematic.lThAngAcc = [lThAngAcc{1}];
            data(trial).Kinetic_Kinematic.rFtAngAcc = [rFtAngAcc{1}];
            data(trial).Kinetic_Kinematic.rSkAngAcc = [rSkAngAcc{1}];
            data(trial).Kinetic_Kinematic.rThAngAcc = [rThAngAcc{1}];
            data(trial).Kinetic_Kinematic.rPvAngAcc = [rPvAngAcc{1}];
            
            data(trial).Kinetic_Kinematic.lFtAngVel = [lFtAngVel{1}];
            data(trial).Kinetic_Kinematic.lSkAngVel = [lShAngVel{1}];
            data(trial).Kinetic_Kinematic.lThAngVel = [lThAngVel{1}];
            data(trial).Kinetic_Kinematic.rFtAngVel = [rFtAngVel{1}];
            data(trial).Kinetic_Kinematic.rSkAngVel = [rSkAngVel{1}];
            data(trial).Kinetic_Kinematic.rThAngVel = [rThAngVel{1}];
            data(trial).Kinetic_Kinematic.rPvAngVel = [rPvAngVel{1}];
            
            data(trial).Kinetic_Kinematic.lFtCGAcc = [lFtCGAcc{1}];
            data(trial).Kinetic_Kinematic.lSkCGAcc = [lSkCGAcc{1}];
            data(trial).Kinetic_Kinematic.lThCGAcc = [lThCGAcc{1}];
            data(trial).Kinetic_Kinematic.rFtCGAcc = [rFtCGAcc{1}];
            data(trial).Kinetic_Kinematic.rSkCGAcc = [rSkCGAcc{1}];
            data(trial).Kinetic_Kinematic.rThCGAcc = [rThCGAcc{1}];
            data(trial).Kinetic_Kinematic.rPvCGAcc = [rPvCGAcc{1}];
            
            data(trial).Kinetic_Kinematic.lFtCGPos = [lFtCGPos{1}];
            data(trial).Kinetic_Kinematic.lSkCGPos = [lSkCGPos{1}];
            data(trial).Kinetic_Kinematic.lThCGPos = [lThCGPos{1}];
            data(trial).Kinetic_Kinematic.rFtCGPos = [rFtCGPos{1}];
            data(trial).Kinetic_Kinematic.rSkCGPos = [rSkCGPos{1}];
            data(trial).Kinetic_Kinematic.rThCGPos = [rThCGPos{1}];
            data(trial).Kinetic_Kinematic.rPvCGPos = [rPvCGPos{1}];
            
            data(trial).Kinetic_Kinematic.lFtCGVel = [lFtCGVel{1}];
            data(trial).Kinetic_Kinematic.lSkCGVel = [lSkCGVel{1}];
            data(trial).Kinetic_Kinematic.lThCGVel = [lThCGVel{1}];
            data(trial).Kinetic_Kinematic.rFtCGVel = [rFtCGVel{1}];
            data(trial).Kinetic_Kinematic.rSkCGVel = [rSkCGVel{1}];
            data(trial).Kinetic_Kinematic.rThCGVel = [rThCGVel{1}];
            data(trial).Kinetic_Kinematic.rPvCGVel = [rPvCGVel{1}];
            
            data(trial).Kinetic_Kinematic.lFtDistEndPos = [lFtDistEndPos{1}];
            data(trial).Kinetic_Kinematic.lSkDistEndPos = [lShDistEndPos{1}];
            data(trial).Kinetic_Kinematic.lThDistEndPos = [lThDistEndPos{1}];
            data(trial).Kinetic_Kinematic.rFtDistEndPos = [rFtDistEndPos{1}];
            data(trial).Kinetic_Kinematic.rSkDistEndPos = [rShDistEndPos{1}];
            data(trial).Kinetic_Kinematic.rThDistEndPos = [rThDistEndPos{1}];
            data(trial).Kinetic_Kinematic.rPvDistEndPos = [rPvDistEndPos{1}];
            
            data(trial).Kinetic_Kinematic.lFtDistEndVel = [lFtDistEndVel{1}];
            data(trial).Kinetic_Kinematic.lSkDistEndVel = [lShDistEndVel{1}];
            data(trial).Kinetic_Kinematic.lThDistEndVel = [lThDistEndVel{1}];
            data(trial).Kinetic_Kinematic.rFtDistEndVel = [rFtDistEndVel{1}];
            data(trial).Kinetic_Kinematic.rSkDistEndVel = [rShDistEndVel{1}];
            data(trial).Kinetic_Kinematic.rThDistEndVel = [rThDistEndVel{1}];
            data(trial).Kinetic_Kinematic.rPvDistEndVel = [rPvDistEndVel{1}];
            
            data(trial).Kinetic_Kinematic.lFtProxEndForce = [lFtProxEndForce{1}];
            data(trial).Kinetic_Kinematic.lSkProxEndForce = [lShProxEndForce{1}];
            data(trial).Kinetic_Kinematic.lThProxEndForce = [lThProxEndForce{1}];
            data(trial).Kinetic_Kinematic.rFtProxEndForce = [rFtProxEndForce{1}];
            data(trial).Kinetic_Kinematic.rSkProxEndForce = [rShProxEndForce{1}];
            data(trial).Kinetic_Kinematic.rThProxEndForce = [rThProxEndForce{1}];
            data(trial).Kinetic_Kinematic.rPvProxEndForce = [rPvProxEndForce{1}];
           
            data(trial).Kinetic_Kinematic.lFtProxEndPos = [lFtProxEndPos{1}];
            data(trial).Kinetic_Kinematic.lSkProxEndPos = [lShProxEndPos{1}];
            data(trial).Kinetic_Kinematic.lThProxEndPos = [lThProxEndPos{1}];
            data(trial).Kinetic_Kinematic.rFtProxEndPos = [rFtProxEndPos{1}];
            data(trial).Kinetic_Kinematic.rSkProxEndPos = [rShProxEndPos{1}];
            data(trial).Kinetic_Kinematic.rThProxEndPos = [rThProxEndPos{1}];
            data(trial).Kinetic_Kinematic.rPvProxEndPos = [rPvProxEndPos{1}];
            
            data(trial).Kinetic_Kinematic.lFtProxEndTorque = [lFtProxEndTorque{1}];
            data(trial).Kinetic_Kinematic.lSkProxEndTorque = [lShProxEndTorque{1}];
            data(trial).Kinetic_Kinematic.lThProxEndTorque = [lThProxEndTorque{1}];
            data(trial).Kinetic_Kinematic.rFtProxEndTorque = [rFtProxEndTorque{1}];
            data(trial).Kinetic_Kinematic.rSkProxEndTorque = [rShProxEndTorque{1}];
            data(trial).Kinetic_Kinematic.rThProxEndTorque = [rThProxEndTorque{1}];
            data(trial).Kinetic_Kinematic.rPvProxEndTorque = [rPvProxEndTorque{1}];
            
            data(trial).Kinetic_Kinematic.lFtProxEndVel = [lFtProxEndVel{1}];
            data(trial).Kinetic_Kinematic.lSkProxEndVel = [lShProxEndVel{1}];
            data(trial).Kinetic_Kinematic.lThProxEndVel = [lThProxEndVel{1}];
            data(trial).Kinetic_Kinematic.rFtProxEndVel = [rFtProxEndVel{1}];
            data(trial).Kinetic_Kinematic.rSkProxEndVel = [rShProxEndVel{1}];
            data(trial).Kinetic_Kinematic.rThProxEndVel = [rThProxEndVel{1}];
            data(trial).Kinetic_Kinematic.rPvProxEndVel = [rPvProxEndVel{1}];
            
            data(trial).Kinetic_Kinematic.lFtSegResidual = [lFtSegResidual{1}];
            data(trial).Kinetic_Kinematic.lSkSegResidual = [lShSegResidual{1}];
            data(trial).Kinetic_Kinematic.lThSegResidual = [lThSegResidual{1}];
            data(trial).Kinetic_Kinematic.rFtSegResidual = [rFtSegResidual{1}];
            data(trial).Kinetic_Kinematic.rSkSegResidual = [rShSegResidual{1}];
            data(trial).Kinetic_Kinematic.rThSegResidual = [rThSegResidual{1}];
            data(trial).Kinetic_Kinematic.rPvSegResidual = [rPvSegResidual{1}];
            
            %% Store Link_Model_Based
 
            data(trial).Link_Model_Based.l_ank_angle = [l_ank_angle{1}];
            data(trial).Link_Model_Based.l_hip_angle = [l_hip_angle{1}];
            data(trial).Link_Model_Based.l_kne_angle = [l_kne_angle{1}];
            data(trial).Link_Model_Based.r_ank_angle = [r_ank_angle{1}];
            data(trial).Link_Model_Based.r_hip_angle = [r_hip_angle{1}];
            data(trial).Link_Model_Based.r_kne_angle = [r_kne_angle{1}];

            data(trial).Link_Model_Based.l_ank_moment = [l_ank_moment{1}];
            data(trial).Link_Model_Based.l_hip_moment = [l_hip_moment{1}];
            data(trial).Link_Model_Based.l_kne_moment = [l_kne_moment{1}];
            data(trial).Link_Model_Based.r_ank_moment = [r_ank_moment{1}];
            data(trial).Link_Model_Based.r_hip_moment = [r_hip_moment{1}];
            data(trial).Link_Model_Based.r_kne_moment = [r_kne_moment{1}];

            data(trial).Link_Model_Based.l_ank_power = [l_ank_power{1}];
            data(trial).Link_Model_Based.l_hip_power = [l_hip_power{1}];
            data(trial).Link_Model_Based.l_kne_power = [l_kne_power{1}];
            data(trial).Link_Model_Based.r_ank_power = [r_ank_power{1}];
            data(trial).Link_Model_Based.r_hip_power = [r_hip_power{1}];
            data(trial).Link_Model_Based.r_kne_power = [r_kne_power{1}];

            data(trial).Link_Model_Based.l_ank_force_loc = [l_ank_force_loc{1}];
            data(trial).Link_Model_Based.l_hip_force_loc = [l_hip_force_loc{1}];
            data(trial).Link_Model_Based.l_kne_force_loc = [l_kne_force_loc{1}];
            data(trial).Link_Model_Based.r_ank_force_loc = [r_ank_force_loc{1}];
            data(trial).Link_Model_Based.r_hip_force_loc = [r_hip_force_loc{1}];
            data(trial).Link_Model_Based.r_kne_force_loc = [r_kne_force_loc{1}];

            data(trial).Link_Model_Based.lft_rotenergy = [lft_rotenergy{1}];
            data(trial).Link_Model_Based.rft_rotenergy = [rft_rotenergy{1}];
            data(trial).Link_Model_Based.lsk_rotenergy = [lsk_rotenergy{1}];
            data(trial).Link_Model_Based.rsk_rotenergy = [rsk_rotenergy{1}];
            data(trial).Link_Model_Based.lth_rotenergy = [lth_rotenergy{1}];
            data(trial).Link_Model_Based.rth_rotenergy = [rth_rotenergy{1}];

            data(trial).Link_Model_Based.l_ank_vel = [l_ank_vel{1}];
            data(trial).Link_Model_Based.l_kne_vel = [l_kne_vel{1}];
            data(trial).Link_Model_Based.l_hip_vel = [l_hip_vel{1}];
            data(trial).Link_Model_Based.r_ank_vel = [r_ank_vel{1}];
            data(trial).Link_Model_Based.r_kne_vel = [r_kne_vel{1}];
            data(trial).Link_Model_Based.r_hip_vel = [r_hip_vel{1}];
            
            data(trial).Link_Model_Based.lsk_wrt_lank = [lsk_wrt_lank{1}];
            data(trial).Link_Model_Based.rsk_wrt_rank = [rsk_wrt_rank{1}];
            data(trial).Link_Model_Based.lth_wrt_lkne = [lth_wrt_lkne{1}];
            data(trial).Link_Model_Based.rth_wrt_rkne = [rth_wrt_rkne{1}];
            data(trial).Link_Model_Based.rpv_wrt_lhip = [rpv_wrt_lhip{1}];
            data(trial).Link_Model_Based.rpv_wrt_rhip = [rpv_wrt_rhip{1}];

            %% Store Subject Variables
            data(trial).Subject.Height = [subjectHeight{1}];
            data(trial).Subject.Mass = [subjectMass{1}];
            data(trial).Subject.centerofmass = [centerofmass{1}];

        else
           keyboard
            
        end
    end

    cd(datafolder)
    save(['p',num2str(subj),'_AllStridesData.mat'], 'data', '-v7.3')
    disp(data)
end
end



