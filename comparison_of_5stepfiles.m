%% Load new and old data

cd('/Users/emilymundinger/Desktop/human-walking-biomechanics')
load('p4_5StridesData.mat');

cd('/Users/emilymundinger/Desktop/human-walking-biomechanics')
dataold= load('p4_5StridesData0403.mat', 'data');

trials=[2:25,32:33]

for t=1:length(trials)
trial=trials(t)

%% make variables
force1_new=data(trial).Force.force1;
force2_new=data(trial).Force.force2;

LSH2_proc_new= data(trial).TargetData.LSH2_pos_proc(:,1:3);
RTH1_proc_new= data(trial).TargetData.RTH1_pos_proc(:,1:3);

lAnkAng_new= data(trial).Link_Model_Based.l_ank_angle;
rKneMom_new= data(trial).Link_Model_Based.r_kne_moment;
lHipPow_new= data(trial).Link_Model_Based.l_hip_power;
rHipVel_new= data(trial).Link_Model_Based.r_hip_vel;


force1_old= dataold.data(trial).Force.force1;
force2_old= dataold.data(trial).Force.force2;

LSH2_proc_old= dataold.data(trial).TargetData.LSH2_pos_proc(:,1:3);
RTH1_proc_old= dataold.data(trial).TargetData.RTH1_pos_proc(:,1:3);

lAnkAng_old= dataold.data(trial).Link_Model_Based.l_ank_angle;
rKneMom_old= dataold.data(trial).Link_Model_Based.r_kne_moment;
lHipPow_old= dataold.data(trial).Link_Model_Based.l_hip_power;
rHipVel_old= dataold.data(trial).Link_Model_Based.r_hip_vel;

%% Differences

diffforce1= force1_new-force1_old;
diffforce2= force2_new-force2_old;

diffLSH2= LSH2_proc_new-LSH2_proc_old;
diffRTH1= RTH1_proc_new-RTH1_proc_old;

difflAnkAng= lAnkAng_new-lAnkAng_old;
diffrKneMom= rKneMom_new-rKneMom_old;
difflHipPow= lHipPow_new-lHipPow_old;
diffrHipVel= rHipVel_new-rHipVel_old;

sum(diffforce1, 'omitnan')
sum(diffforce2, 'omitnan')
sum(diffLSH2, 'omitnan')
sum(diffRTH1, 'omitnan')
sum(difflAnkAng, 'omitnan')
sum(diffrKneMom, 'omitnan')
sum(difflHipPow, 'omitnan')
sum(diffrHipVel, 'omitnan')


if sum(diffforce1>0, 'omitnan')
    keyboard
elseif sum(diffforce2>0, 'omitnan')
    keyboard
elseif sum(diffLSH2, 'omitnan')
    keyboard
elseif sum(diffRTH1, 'omitnan')
    keyboard
elseif sum(difflAnkAng, 'omitnan')
    keyboard
elseif sum(diffrKneMom, 'omitnan')
    keyboard
elseif sum(difflHipPow, 'omitnan')
    keyboard
elseif sum(diffrHipVel, 'omitnan')
    keyboard
end
    



end   
    
    
%% Plot to compare new and old variables

% figure('name', 'Force data Difference')
% 
% plot(force1_new(:,1), force1_old(:,1)); hold on
% plot(force2_new(:,1), force2_old(:,1)); hold on
% plot(force1_new(:,2), force1_old(:,2)); hold on
% plot(force2_new(:,2), force2_old(:,2)); hold on
% plot(force1_new(:,3), force1_old(:,3)); hold on
% plot(force2_new(:,3), force2_old(:,3)); hold on
% 
% figure('name', 'Processed target data Difference')
% 
% plot(LSH2_proc_new(:,1)); hold on; plot(LSH2_proc_old(:,1)); hold on
% plot(RTH1_proc_new(:,1)); hold on; plot(RTH1_proc_old(:,1)); hold on
% plot(LSH2_proc_new(:,2)); hold on; plot(LSH2_proc_old(:,2)); hold on
% plot(RTH1_proc_new(:,2)); hold on; plot(RTH1_proc_old(:,2)); hold on
% plot(LSH2_proc_new(:,3)); hold on; plot(LSH2_proc_old(:,3)); hold on
% plot(RTH1_proc_new(:,3)); hold on; plot(RTH1_proc_old(:,3)); hold on
% 
% figure('name', 'Link_Model_Based Items Difference')
% 
% plot(lAnkAng_new, lAnkAng_old); hold on
% plot(RTH1_proc_new, RTH1_proc_old); 
