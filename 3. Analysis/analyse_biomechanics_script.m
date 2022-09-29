clear all; close all; clc
% Analyzes biomechanical data in 5 strides files and create a summary file (Wsoft.mat)

%% Settings
subjects = 1:9;
trials = 1:33;

% folder where the files are that have been exported from Visual3D
import_folder = uigetdir;
% import_folder = '';

% folder where to save the summary file
export_folder = uigetdir;
% export_folder = '';

%% Participant and data collection parameters
vwalks = [.7 .7 .7 .9 .9 .9 1.1 1.1 1.1 1.6 1.6 1.6 1.8 1.8 1.8 2.0...
    1.25 1.25 1.25 1.25 1.25 1.25 1.25 1.25 1.25 1.25 1.25 1.25 1.25 1.25 1.4 1.4 1.4]; % m/s

% participant body mass
mass = [81.8 57.3 97.5 57 56.7 72.6 86.2 88.6 77]; % [kg]

% assumed relative segment mass
relsegmass = [0.0145 0.0465 .1 0.0145 0.0465 .1 .142]; % according to visual 3D
% https://www.c-motion.com/v3dwiki/index.php/Segment_Mass

% data collection sample frequencies
fsgrf = 1200; % [Hz]
fsmoc = 120;  % [Hz]

% corresponding delta(t) values
dtgrf = 1/fsgrf;
dtmoc = 1/fsmoc;

nt = 33; ns = 9; nl = 2;

%% preallocate variables of interest
Wsoft = nan(nt,ns,nl);
Wankle = nan(nt,ns,nl);
Wknee = nan(nt,ns,nl);
Whip = nan(nt,ns,nl);
Wperi = nan(nt,ns,nl);
Wcom = nan(nt,ns,nl);
Wbody = nan(nt,ns,nl);
Wbody_neg = nan(nt,ns,nl);
Wbody_pos = nan(nt,ns,nl);
Wjoint = nan(nt,ns,nl);
Wjoint_pos = nan(nt,ns,nl);
Wjoint_neg = nan(nt,ns,nl);
Wbodycoll = nan(nt,ns,nl);
WCOMcoll = nan(nt,ns,nl);
Wsoftcoll = nan(nt,ns,nl);
Wjointcoll = nan(nt,ns,nl);
Fpeak = nan(nt,ns,nl);
Tstride = nan(nt,ns,nl);
vcom_hs = nan(nt,ns,nl);
vcom_hs_alt = nan(nt,ns,nl);
deltav = nan(nt,ns,nl);
nsteps = nan(nt,ns,nl);
Wjointtranscoll = nan(nt,ns,nl);
Wjoint_trans = nan(nt,ns,nl);
Wankle_trans = nan(nt,ns,nl);
Wjoint_rot = nan(nt,ns,nl);
Wperi_neg = nan(nt,ns,nl);
FTI = nan(nt,ns,nl);

%% Main code
for subj = subjects
    
    % preallocate (required for newer MATLAB versions)  
    data = [];
    
    % lood 5 strides data
    cd(import_folder)
    
    if exist(['p',num2str(subj),'_5StridesData.mat'], 'file')
        load(['p',num2str(subj),'_5StridesData.mat'])    
 
    else
        disp('No .mat file found')
    end
    
    disp(['subj:', num2str(subj)])
    disp('trial:')
    
%for trial = [1]    
for trial = 1:33
        
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
        
    if (~isempty(data(trial).Force))
        if  isfield(data(trial).Force,'force1')

          
    %% Force platform: GRF, vCOM and PCOM
    grfl = data(trial).Force.force1(:,1:3);
    grfr = data(trial).Force.force2(:,1:3);
    grf = grfl + grfr;
    
    % peak grf
    Fpeak(trial,subj) = max(data(trial).Force.force1,[],'all');

    % make time series
    ngrf = length(data(trial).Force.force1);
    nmoc = length(data(trial).Link_Model_Based.rpv_wrt_rhip);
    
    tgrf = 0:dtgrf:(dtgrf*(ngrf-1));
    tmoc = 0:dtmoc:(dtmoc*(nmoc-1));
    
    %% Gait events
    if ~(subj == 7 && trial == 27), threshold = 30;
    else, threshold = 50;
    end
        
    [hsl, tol, hsr, tor] = invDynGrid_getHS_TO(grfl, grfr, threshold);
    
    hsr(6) = length(grfl); % often not found with the function
    
    hsl = unique(hsl); hsl(find(diff(hsl)<5)) = [];
    hsr = unique(hsr); hsr(find(diff(hsr)<5)) = [];
    tol = unique(tol); tol(find(diff(tol)<5)) = [];
    tor = unique(tor); tor(find(diff(tor)<5)) = [];

    nsteps(trial,subj,1) = length(hsl);
    nsteps(trial,subj,2) = length(hsr);
    
    actual_mass = mean(grf(hsl(1):hsl(end),3))/9.81;

    FTI(trial,subj,1) = mean(grfl(hsl(1):hsl(end),3));
    FTI(trial,subj,2) = mean(grfr(hsr(1):hsr(end),3));
    
%% Center of mass velocity
    vcom_raw = []; vcom = []; Pcom = [];
    vcom_raw(:,1) = cumtrapz(tgrf,(grf(:,1))/actual_mass); % Donelan et al. 2002, med-lat
    vcom_raw(:,2) = cumtrapz(tgrf,(grf(:,2))/actual_mass) - vwalks(trial); % Donelan et al. 2002, for-aft
    vcom_raw(:,3) = cumtrapz(tgrf,(grf(:,3)-actual_mass*9.81)/actual_mass); % Donelan et al. 2002, vertical
    
    p = nan(3,2);
    for i = 1:3
        p(i,:) = polyfit(hsl(1):hsl(end), vcom_raw(hsl(1):hsl(end),i)', 1);
        vcom(:,i) = vcom_raw(:,i) - polyval(p(i,:), 1:length(tgrf))';
    end
    
    vcom(:,2) = vcom(:,2) - vwalks(trial);
    
    % center of mass power
    Pcom(:,1) = sum(grfl .* vcom,2);
    Pcom(:,2) = sum(grfr .* vcom,2);
    
    % flip for-aft to make it positive
    Vcom = [vcom(:,1) -vcom(:,2) vcom(:,3)];
    
    %% Velocity events   
    % COM veloicty events
    hstohstohs = [hsl(1:5) tor(1:5) hsr(1:5) tol(1:5)]';
    out = comvelocityevents(Vcom,[hstohstohs(1:numel(hstohstohs)) hsl(end)]);

    % v-up events for left and right leg
    vupl = out.vcomevents(4:12:end);
    vupr = out.vcomevents(10:12:end); 
    
    % v-down events for left and right leg
    vdownl = out.vcomevents(6:12:end);
    vdownr = out.vcomevents(12:12:end);
    
    % convert to motion capture frames
    hs = [hsl hsr];
    hs_mo = ceil(hs/10);
    if hs_mo(end,2) > length(tmoc), hs_mo(end,2) = length(tmoc);
    end
    vup_mo = ceil([vupl vupr]/10);
    
    % velocity at heelstrike
    vcom_hs(trial,subj,1) = norm(mean(Vcom(vdownl,:)));
    vcom_hs(trial,subj,2) = norm(mean(Vcom(vdownr,:)));
    
    vcom_hs_alt(trial,subj,1) = norm(mean(Vcom(hsl-1,:)));
    vcom_hs_alt(trial,subj,2) = norm(mean(Vcom(hsr-1,:)));
    
%     if trial == 20
%         figure(1);
%         plot(Vcom(:,2),Vcom(:,3),'k'); hold on
%         
%         plot(Vcom(vdownl,2), Vcom(vdownl,3),'ro');
%         plot(Vcom(vdownr,2), Vcom(vdownr,3),'bo');
%         
%         plot(Vcom(vupl,2), Vcom(vupl,3),'rx');
%         plot(Vcom(vupr,2), Vcom(vupr,3),'bx');
%         
%         
%         figure
%         plot(Vcom(:,3));
%         hold on
%         plot(vupl,Vcom(vupl,3),'rx');
%        plot(vdownr,Vcom(vdownr,3),'bo');
%        
%        plot(vupr,Vcom(vupr,3),'bx');
%        plot(vdownl,Vcom(vdownl,3),'ro');
%         
% %         keyboard
%     end
%         if subj == 5
%             saveV = Vcom;
%         end
%     end

    
    % delta velocity 
    deltav(trial,subj,1) = mean(atan2d(Vcom(vupr,3), Vcom(vupr,2)) - atan2d(Vcom(vdownl,3), Vcom(vdownl,2)));
    deltav(trial,subj,2) = mean(atan2d(Vcom(vupl(2:end),3), Vcom(vupl(2:end),2)) - atan2d(Vcom(vdownr(1:4),3), Vcom(vdownr(1:4),2)));
    
    %% Stride time
    % stride time (take length(grf) because things end with right
    % heelstrike, but this last heelstrike is not included in hsr)
    Tstride(trial,subj,1) = ((hsl(end)-hsl(1))*dtgrf)/5;
    Tstride(trial,subj,2) = ((hsr(end)-hsr(1))*dtgrf)/5;
%     
    %% COM power and velocity in motion capture time
    Pcom_mo = [];
    Pcom_mo(:,1) = interp1(tgrf,Pcom(:,1),tmoc);
    Pcom_mo(:,2) = interp1(tgrf,Pcom(:,2),tmoc);

    vcom_mo = interp1(tgrf, vcom, tmoc);
    vcom_mo(:,2) = vcom_mo(:,2) + vwalks(trial); % relative to the lab frame

    % alternative COM power: from pelvis velocity
    Pcom_mo_alt = []; vcom_mo_alt = []; grf_mo = [];
    grf_mo(:,1:3) = interp1(tgrf,grfl,tmoc); grf_mo(:,4:6) = interp1(tgrf,grfr,tmoc);
    
    %% Peripheral work 
    segmentmass = mass(subj) * relsegmass;
    
    % segment_velocities: CG velocity of segment, must be Nx21 for (1) foot, (2) shank, (3) thigh and
    % last 3 columns are pelvis, first 9 columns for left leg, last 9 for
    % right leg. global coordinate system
    segment_velocities = [];
    segment_velocities(:,1:3) = data(trial).Kinetic_Kinematic.lFtCGVel;
    segment_velocities(:,4:6) = data(trial).Kinetic_Kinematic.lSkCGVel;
    segment_velocities(:,7:9) = data(trial).Kinetic_Kinematic.lThCGVel;
    segment_velocities(:,10:12) = data(trial).Kinetic_Kinematic.rFtCGVel;
    segment_velocities(:,13:15) = data(trial).Kinetic_Kinematic.rSkCGVel;
    segment_velocities(:,16:18) = data(trial).Kinetic_Kinematic.rThCGVel;
    segment_velocities(:,19:21) = data(trial).Kinetic_Kinematic.rPvCGVel;
    
    % rotational_energies: rotational energy of each segment, must be Nx6
    % or (1) foot, (2) shank, (3) thigh, first 3 columns for left leg, last
    % 3 for right leg
    rotational_energies = [data(trial).Link_Model_Based.lft_rotenergy data(trial).Link_Model_Based.lsk_rotenergy data(trial).Link_Model_Based.lth_rotenergy...
                           data(trial).Link_Model_Based.rft_rotenergy data(trial).Link_Model_Based.rsk_rotenergy data(trial).Link_Model_Based.rth_rotenergy];
            
    [Pper] = CalcPeripheralPower(segment_velocities, rotational_energies, segmentmass, vcom_mo, fsmoc);
    
    %% Translational joint power
    % local_CG_position is the position of segments CG wrt to the distal joint (proximal point on
    % distal segment) in proximal segment coordinate system(?)
    local_CG_position = [data(trial).Link_Model_Based.lsk_wrt_lank data(trial).Link_Model_Based.lth_wrt_lkne data(trial).Link_Model_Based.rpv_wrt_lhip...
                         data(trial).Link_Model_Based.rsk_wrt_rank data(trial).Link_Model_Based.rth_wrt_rkne data(trial).Link_Model_Based.rpv_wrt_rhip];
    
    % local_joint_force is the force of the proximal segment on the distal
    % segment, expressed in proximal segment coordinate
    local_joint_force = [data(trial).Link_Model_Based.l_ank_force_loc data(trial).Link_Model_Based.l_kne_force_loc data(trial).Link_Model_Based.l_hip_force_loc...
                         data(trial).Link_Model_Based.r_ank_force_loc data(trial).Link_Model_Based.r_kne_force_loc data(trial).Link_Model_Based.r_hip_force_loc];
    
    transjointpower = CalcTransJointPower_Locally(local_CG_position, local_joint_force, fsmoc);
%     
    %% Joint powers
    % ankle power
    Pank = nan(nmoc, 2,2);
    Pank(:,:,1) = [sum(data(trial).Link_Model_Based.l_ank_power,2) sum(transjointpower(:,1:3),2)] * mass(subj);
    Pank(:,:,2) = [sum(data(trial).Link_Model_Based.r_ank_power,2) sum(transjointpower(:,10:12),2)] * mass(subj);
%     
    % knee power
    Pkne = nan(nmoc, 2,2);
    Pkne(:,:,1) = [sum(data(trial).Link_Model_Based.l_kne_power,2) sum(transjointpower(:,4:6),2)] * mass(subj);
    Pkne(:,:,2) = [sum(data(trial).Link_Model_Based.r_kne_power,2) sum(transjointpower(:,13:15),2)] * mass(subj);
%     
    % hip power
    Phip = nan(nmoc, 2,2);
    Phip(:,:,1) = [sum(data(trial).Link_Model_Based.l_hip_power,2) sum(transjointpower(:,7:9),2)] * mass(subj);
    Phip(:,:,2) = [sum(data(trial).Link_Model_Based.r_hip_power,2) sum(transjointpower(:,16:18),2)] * mass(subj);
%     
    %% Summed joint and soft tissue
    Pjoint_rot = Pank(:,1,:)  + Pkne(:,1,:) + Phip(:,1,:); % summed rotational power
    Pjoint_trans = Pank(:,2,:)  + Pkne(:,2,:) + Phip(:,2,:); % summed translational power
    Pjoint = squeeze(Pjoint_rot + Pjoint_trans);
    Pbody = Pcom_mo + Pper; % whole-body power
    Psoft = Pbody - Pjoint; % soft tissue power
%     
    %% Work terms: full stride
    for i = 1:2
        % individual joint work
        Wankle(trial,subj,i) = fintrapz(tmoc, sum(Pank(:,:,i),2),[],[hs_mo(1,i) hs_mo(end,i)]) / 5;
        Wknee(trial,subj,i) = fintrapz(tmoc, sum(Pkne(:,:,i),2),[],[hs_mo(1,i) hs_mo(end,i)]) / 5;
        Whip(trial,subj,i) = fintrapz(tmoc, sum(Phip(:,:,i),2),[],[hs_mo(1,i) hs_mo(end,i)]) / 5;
        
        % summed joint
        Wjoint(trial,subj,i) = fintrapz(tmoc, Pjoint(:,i),[],[hs_mo(1,i) hs_mo(end,i)]) / 5;
        Wjoint_pos(trial,subj,i) = fintrapz(tmoc, Pjoint(:,i),'+',[hs_mo(1,i) hs_mo(end,i)]) / 5;
        Wjoint_neg(trial,subj,i) = fintrapz(tmoc, Pjoint(:,i),'-',[hs_mo(1,i) hs_mo(end,i)]) / 5;
        Wjoint_trans(trial,subj,i) = fintrapz(tmoc, Pjoint_trans(:,i),[],[hs_mo(1,i) hs_mo(end,i)]) / 5;
        Wankle_trans(trial,subj,i) = fintrapz(tmoc, Pank(:,2,i),[],[hs_mo(1,i) hs_mo(end,i)]) / 5;
        Wjoint_rot(trial,subj,i) = fintrapz(tmoc, Pjoint_rot(:,i),[],[hs_mo(1,i) hs_mo(end,i)]) / 5;

%         whole-body work
        Wperi(trial,subj,i) = fintrapz(tmoc, Pper(:,i),[],[hs_mo(1,i) hs_mo(end,i)]) / 5;
        Wcom(trial,subj,i) = fintrapz(tgrf, Pcom(:,i),[], [hs(1,i) hs(end,i)]) / 5;
        Wbody(trial,subj,i) = fintrapz(tmoc, Pbody(:,i),[],[hs_mo(1,i) hs_mo(end,i)]) / 5;
        Wbody_neg(trial,subj,i) = fintrapz(tmoc, Pbody(:,i),'-',[hs_mo(1,i) hs_mo(end,i)]) / 5;
        Wbody_pos(trial,subj,i) = fintrapz(tmoc, Pbody(:,i),'+',[hs_mo(1,i) hs_mo(end,i)]) / 5;
        Wperi_neg(trial,subj,i) = fintrapz(tmoc, Pper(:,i),'-',[hs_mo(1,i) hs_mo(end,i)], Pbody(:,i)) / 5;
        
        % soft tissue work
        Wsoft(trial,subj,i) = fintrapz(tmoc, Psoft(:,i),[],[hs_mo(1,i) hs_mo(end,i)]) / 5;
    end
    
    %% Work terms: collision phase
    % collision work
    Wbodycol_step = nan(5,2); WCOMcol_step = nan(5,2);  Wsoftcol_step = nan(5,2);   Wjointcol_step = nan(5,2); Wjointtrans_col_step = nan(5,2);
    
    for i = 1:2
        for istep = 1:5       
                Wbodycol_step(istep,i) = fintrapz(tmoc, Pbody(:,i), '-', [hs_mo(istep,i) vup_mo(istep,i)]);
                WCOMcol_step(istep,i) = fintrapz(tmoc, Pcom_mo(:,i), '-', [hs_mo(istep,i) vup_mo(istep,i)]);
                Wsoftcol_step(istep,i) = fintrapz(tmoc, Psoft(:,i), '-', [hs_mo(istep,i) vup_mo(istep,i)], Pbody(:,i));  
                Wjointcol_step(istep,i) = fintrapz(tmoc, Pjoint(:,i), '-', [hs_mo(istep,i) vup_mo(istep,i)], Pbody(:,i));  
                Wjointtrans_col_step(istep,i) = fintrapz(tmoc, Pjoint_trans(:,i), '-', [hs_mo(istep,i) vup_mo(istep,i)], Pbody(:,i));  
        end
        
        % average over steps
        WCOMcoll(trial,subj,i) = mean(WCOMcol_step(:,i));
        Wbodycoll(trial,subj,i) = mean(Wbodycol_step(:,i));
        Wsoftcoll(trial,subj,i) = mean(Wsoftcol_step(:,i));
        Wjointcoll(trial,subj,i) = mean(Wjointcol_step(:,i)); 
        Wjointtranscoll(trial,subj,i) = mean(Wjointtrans_col_step(:,i)); 
    end   
    
        else
            disp(['p',num2str(subj),'_5StepsData : trial', num2str(trial), ' is NaN'])
        end
    else
        disp(['p',num2str(subj),'_5StepsData : trial', num2str(trial), ' is empty'])
    end
end

end

% missing trials
[r,c] = find(isnan(nsteps(:,:,1)));
missing = [c r];

%% Saving
cd(export_folder)
save('Wsoft_old')

%% Functions
function [Pper] = CalcPeripheralPower(segmentvelocity, rotenergy, segmenmass, vcom_mo, fsmo)

Etrans = nan(length(vcom_mo),6);
k = 0;
for s = 1:3:size(segmentvelocity,2)
    k = k+1;
    vrel = segmentvelocity(:,s:s+2) - vcom_mo;
    Etrans(:,k) = .5 * segmenmass(k) * dot(vrel', vrel')';
end

% sum over ankle, knee and hip per leg
Eperl = sum(Etrans(:,1:3),2) + sum(rotenergy(:,1:3),2);
Eperr = sum(Etrans(:,4:6),2) + sum(rotenergy(:,4:6),2);

% take derivative
Pper = nan(length(Eperl),2);
Pper(isfinite(Eperl),1) = grad5(Eperl(isfinite(Eperl)), 1/fsmo);
Pper(isfinite(Eperr),2) = grad5(Eperr(isfinite(Eperr)), 1/fsmo);

%% optional check: compare pelvis velocity with COM velocity, should be fairly similar
% figure;
% color = get(gca,'colororder');
% for i = 1:3
%     plot(segmentvelocity(:,end-3+i),'color',color(i,:)); hold on
%     plot(vcom_mo(:,i), '--','color',color(i,:))
% end
%     
end

function[transjointpower] = CalcTransJointPower_Locally(proxsegmentcom_loc_wrt_jointpos, fjoint_loc, fs)
% [transjointpower] = CalcTransJointPower_Locally(proxsegmentcom_loc_wrt_jointpos, fjoint_loc, fs)
% is the translational joint power expressed in proximal segment local
% coordinate system (Nx3) 
%
% jointpos_wrt_proxsegmentcom_loc is the segment center of gravity (CG) position, 
% computed with the V3D fucntion SEG_CGPOSITION
% (https://c-motion.com/v3dwiki/index.php?title=Model_Based_Item:_SEG_CGPOSITION)
% using the proximal segment as segment, the distal segment as reference
% and the proximal segment as resolution 
%
% fjoint_loc is the joint force, computed with the V3D function JOINT_FORCE
% (https://c-motion.com/v3dwiki/index.php?title=Model_Based_Item:_JOINT_FORCE)
% using the joint as joint/segment, and the proximal segment as resolution
% fs is the sample frequency (Hz) 
%
% Tim van der Zee, University of Calgary, 2020

% We want the opposite of proxsegmentcom_loc_wrt_jointpos, namely
% jointpos_wrt_proxsegmentcom_loc:
jointpos_wrt_proxsegmentcom_loc = -proxsegmentcom_loc_wrt_jointpos;

% Now take the time derivative of that:
jointvel_wrt_proxsegmentcom_loc = grad5(jointpos_wrt_proxsegmentcom_loc, 1/fs);

% And calculate power as the element-wise product with force
transjointpower = fjoint_loc .* jointvel_wrt_proxsegmentcom_loc;
end

%% comvelocityevents Function

function out = comvelocityevents(Vcom,hstohstohs)
% out = comvelocityevents(Vcom,hstohstohs);
% A function to return the indices, Velocities at, Symbols for different
% events of interest in a hodograph. Output is a structure.  
% 
% Vcom = COM velocity in 3D, where X is Lateral, Y is Forward and Z is Up. 
% 
% hstohstohs = a vector of 2*nsteps+1 indices indicating Heel Strikes and
% Toe-Offs within the Vcom signal. "nsteps" is the number of STEPS in the
% Vcom signal. The first hs is the initial heel strike,  and the last hs is
% the "closing" heel strike, usually coincident with the last sample of Vcom

n = size(Vcom,1);

% Determine the number of STEPS in the Vcom signal: 
% length(hstohstohs) = 2*nsteps+1. 
% CURRENTLY LIMITED to 2 steps! (never tried with more) 2009-11-23 PGA
% if length(hstohstohs) == 3
%     nsteps = 1;
% elseif length(hstohstohs) == 5; 
%     nsteps = 2;
% else
%     error('Wrong number of steps for "comvelocitymetrics"')
% end

% Try to make it work with more steps 2020-06-25 TJZ
nsteps = (length(hstohstohs)-1)/2;

% create "vdup", a duplicate of Vcom that will have a bunch of
% regions-of-no-interest blanked out with NaN at different times.
vdup = Vcom;

% parse heel strikes and toe-offs 
hs = hstohstohs(1:2:end);
to = hstohstohs(2:2:end); 
    to(end+1) = to(1)+hs(end)-1; % add an after-step toe off, just to give a programmatic index for NaNing data for vdn

vcomevents = hs(1);  % indices of events of interest: first is the initial heel strike
vcomeventsyms = {'s'}; % symbols used to mark these events on a plot. HS gets a square ('s')

for ii = 1:nsteps  % loop over steps
    
    %% Mid-Double-Support
    vdup(1:(hs(ii)-1),:) = NaN;  % set anything before the current step to NaN
    if vdup(hs(ii),3)>0   % check if V(hs,3) is greater than zero. If so, hs occurred with rising com velocity. 
        imidds(ii) = hs(ii);
    else
        imidds(ii) = find(vdup(:,3)>0,1,'first'); % otherwise look for mid-stance as an upward zero-crossing in vcom.
    end
    
    %% V-Up
    vdup(1:(imidds(ii)),:) = NaN; % cancel out times in the past
    vdup(round(mean([to(ii) hs(ii+1)])):end,:) = NaN; % cancel out times after roughly mid-step (halfway between toe off and heel strike)
    [bla, blah] = max(atan2(vdup(:,3),abs(vdup(:,2)))); % find "vup", the steepest upward velocity of the COM
    ivup(ii) = blah;
        vdup = Vcom; % reset vdup
    
    %% Mid-Single Support
    vdup(1:(ivup(ii)-1),:) = NaN; % cancel out times in the past
    vdup(hs(ii+1):end,:) = NaN; % cancel out times in the next step
%     if ii == 7, keyboard
%     end
    imidss(ii) = find(vdup(:,3)<0,1,'first'); % Mid-single--support is a downward zero-crossing in Vcom
    
    %% V-Down
%     vdup(1:(imidss(ii)-1),:) = NaN;
%     vdup(hs(ii+1):end,:) = NaN; % don't look too far this time...
    vdup(vdup(:,3)>0,:) = NaN;  % if the COM is rising, this can't be the steepest fall
    vdup(to(ii+1):end,:) = NaN; % don't look too far (end at next toe off)
    [bla blah] = min( atan2(vdup(:,3),abs(vdup(:,2))) );  % find "vdn", the steepest downward COM velocity
    ivdn(ii) = blah ;
%     if ~isempty(blah)
%         ivdn(ii) = blah ;
%     else % special case for a "vdn" that has not happened before the stride ends. Wrap around to the beginning. 
%         vdup = Vcom; % reset vdup so we have the beginning again
%         vdup(vdup(:,3)>0) = NaN;  % if the COM is rising, this can't be the steepest fall
%         vdup (to(1)+1:end,:) = NaN; % since we are only looking at the beginning, get rid of the later part. 
%         blah = find( (diff(atan2(vdup(:,3),abs(vdup(:,2)))) < 0), 1, 'last'); 
%         if ~isempty(blah)
%             ivdn(ii) = blah; % use the wrapped index for vdn, if it was found
%         else
%             ivdn(ii) = hs(end) ; % use the last heel strike as the index of vdn. 
%         end
%     end
        vdup = Vcom; % Reset vdup for next time around;
    
    % Concatenate all these events as: [previousevents, mid-double-support, toe-off, vup, mid-single-support, vdn, next-heel-strike]
    % First time, previous is [hs], so this yields [hs midds to vup midss vdn hs]
    vcomevents = [vcomevents', imidds(ii), to(ii), ivup(ii), imidss(ii), ivdn(ii), hs(ii+1)]';
    vcomeventsyms = {vcomeventsyms{:}, '.', 'd', '^', 'x', 'v', 's' }';
end

vcomeventVs = Vcom(vcomevents,:); % get the actual COM velocity at these events

% Make the output Structure: events, symbols, velocities
out.vcomevents = vcomevents;
out.vcomeventsyms = cell2mat(vcomeventsyms);
out.vcomeventVs = vcomeventVs;

end


%% fintrapz Function

function [I] = fintrapz(x,y,method,range,y2)

if nargin > 2
    if nargin < 5
        y2 = y;
    end
    
    if strcmp(method,'+') == 1, y(y2<0) = 0;
    elseif strcmp(method,'-') == 1, y(y2>0) = 0;
    end
end

if nargin > 3
    r = range;
else
    r = [1 length(y)];
end

Y = y(r(1):r(2));
X = x(r(1):r(2));

i = isfinite(Y);
I = trapz(X(i),Y(i));

end

%% grad5 Function

function yp=grad5(y,dx)
% function yp=grad5(y,dx)
% 031091 KvS
% purpose: calculation of 5 point derivative
% inputs : y : vector containing signal as a function of x
%          dx: change in x between steps
% output : yp: derivative of y with respect to x
 
if nargin~=2
  disp('OOPS: grad5 expects 2 input arguments')
  disp('      hit any key to continue');pause
  return
end

%y=y(:);
[nrow,ncol]=size(y);


yp(1,1:ncol)=zeros(1,ncol);
yp(2,1:ncol )=zeros(1,ncol);
yp(nrow-1,1:ncol)=zeros(1,ncol);
yp(nrow,1:ncol)=zeros(1,ncol);

yp(1,1:ncol)=(y(2,:)-y(1,:))/dx;
yp(2,1:ncol)=(y(3,:)-y(1,:))/(2*dx);
yp(nrow-1,1:ncol)=(y(nrow,:)-y(nrow-2,:))/(2*dx);
yp(nrow,1:ncol)=(y(nrow,:)-y(nrow-1,:))/dx;

coef=[1 -8 0 8 -1];
for i=3:nrow-2;
  yp(i,:)=(coef*y(i-2:i+2,:))/(12*dx);
end;

%if flip; y=y';yp=yp';end;
end
