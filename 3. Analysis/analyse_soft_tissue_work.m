clear all; close all;

% folder where the files are that have been exported from Visual3D
% folder = uigetdir;
folder = 'C:\Users\timvd\Documents\Inverse dynamics\Level 3 - MATLAB files\Level 3 - MATLAB files - reproduced\5 Strides Data files from process_5steps_new';
cd(folder);

%% Code
% load('Wsoft_29-Sep-2022.mat')
load('Wsoft.mat')
terms = {'Collision', 'Soft tissue', 'Negative'};
old_order = 0;

% in van der Zee & Kuo (2021)
coefficients_in_paper = [4.512, 2.847, 5.083];
SE_in_paper = [0.079, 0.059, 0.108];

for j = 1:3
    term = terms{j};

if strcmp(term, 'Soft tissue')
    Wterm = Wsoftcoll;
    modelString = 'Work ~ -1 + GaitPar';
    i = 1;
elseif strcmp(term, 'Collision')
    Wterm = Wbodycoll;
    modelString = 'Work ~ -1 + GaitPar';
    i = 1;
else
    Wterm = Wbody_neg;
    modelString = 'Work ~ GaitPar';
    i = 2;
end

if old_order
% parameters
parms.leglength = [0.8941 0.9398 1.04  0.8636 0.8636 0.9398 0.9906 0.99 0.9398 0.876]; % m
parms.mass = [81.8000 57.3000 97.5000 60.3000 56.7000 72.6000 86.2000 88.6000 77 57]; % kg

else
% 1-9
parms.leglength = [0.8941 0.9398 1.04 0.876 0.8636 0.9398 0.9906 0.99 0.9398]; % m
parms.mass = [81.8000 57.3000 97.5000 57 56.7000 72.6000 86.2000 88.6000 77]; % kg
end

parms.g = 9.81; % m/s2
parms.dimensionless = 1;

% reorganize
[~, SPs, Wmat, SPs_av, dVs] = get_Wmat_v2(Wterm, parms,'velocity-based&heelstrike',old_order);

 % dVs is total redirection, and we assume collision takes up half
Delta = tan(dVs * pi/180 / 2);

% Eq. 4-6 in van der Zee & Kuo (2021)
if old_order
    subjs = [1:3, 5:10];
else
subjs = 1:9;

end
X = .5 .* SPs(:,subjs).^2 .* Delta(:,subjs).^2;
Y = Wmat(:,subjs);
S = repmat((1:size(X,2))',1, size(X,1))';

% lump together in one big table
T = array2table([(Y(1:numel(Y)))' (X(1:numel(X)))' (S(1:numel(S)))']);
T.Properties.VariableNames = {'Work' 'GaitPar' 'Subject'};

% do linear regression!
LM = fitlm(T,modelString);

% display coefficients
disp(['Coefficient for ',term,': ', num2str(round(LM.Coefficients.Estimate(i),3)), ' +- ',  num2str(round(LM.Coefficients.SE(i),3))])
disp(['Coef. in paper for ',term,': ', num2str(coefficients_in_paper(j)), ' +- ',  num2str(SE_in_paper(j))])
end

function[SLs, SPs, Wmat, SPs_av, dVs] = get_Wmat_v2(Wterm, parms,type, old_order)

% cd(datafolder)
load('Wsoft.mat','Tstride','vcom_hs','deltav','vcom_hs_alt')

if nargin == 2
    type = [];
end

if strcmp(type,'velocity-based&heelstrike')
    vcom_hs = vcom_hs_alt;
end

if old_order
% % exclude subject 4 and trial 4 of subject 3
Wterm(:,4,:) = nan;
% Wterm(4,3,:) = nan;
Tstride(:,4,:) = nan;
end
% Tstride(4,3,:) = nan;

%% Reorganize based on freq and SL and speed
if old_order
    N = 10; % number of subjects
else
    N = 9;
end
M = 10; % numbers of trials

% define variables
vCSLs = nan(M,N);
vCSPs = nan(M,N);
vCFRs = nan(M,N);
vPRFs = nan(M,N);
vCSPs_hs = nan(M,N);

dVCSPs = nan(M,N);
dVCSLs = nan(M,N);
dVCFRs = nan(M,N);
dVPRFs = nan(M,N);

TstepPRF = nan(M,N);
TstepcSP = nan(M,N);
TstepcSL = nan(M,N);
TstepcFR = nan(M,N);

WcSL = nan(M,N,2);
WcFR = nan(M,N,2);
WPRF = nan(M,N,2);
WcSP = nan(M,N,2);

% trials
icSL = [1 4 7 12 15 20 33];
icFR = [3 6 9 10 13 20 31];
iPREF = [2 5 8 11 14 16 20 32];
icSP = [17:20, 23:25];

% corresponding speeds / step lengths
vCSLs(1:7,:) = repmat([.7 .9 1.1 1.6 1.8 1.25 1.4]',1,N);
vCFRs(1:7,:) = repmat([.7 .9 1.1 1.6 1.8 1.25 1.4]',1,N);
vPRFs(1:8,:) = repmat([.7 .9 1.1 1.6 1.8 2.0 1.25 1.4]',1,N);

Tstep = Tstride / 2; % stride - > step

vCSPs = 1.25 * ones(size(vCSPs));

% heelstrike speeds
vCSPs_hs(1:length(icSP),:) = mean(vcom_hs(icSP,:,:),3);
vCSLs_hs(1:length(icSL),:) = mean(vcom_hs(icSL,:,:),3);
vCFRs_hs(1:length(icFR),:) = mean(vcom_hs(icFR,:,:),3);
vPRFs_hs(1:length(iPREF),:) = mean(vcom_hs(iPREF,:,:),3);

% delta velocities
dVCSPs(1:length(icSP),:) = mean(deltav(icSP,:,:),3);
dVCSLs(1:length(icSL),:) = mean(deltav(icSL,:,:),3);
dVCFRs(1:length(icFR),:) = mean(deltav(icFR,:,:),3);
dVPRFs(1:length(iPREF),:) = mean(deltav(iPREF,:,:),3);

% step times
TstepcSL(1:7,:) = mean(Tstep(icSL,:,:),3);
TstepcFR(1:7,:) = mean(Tstep(icFR,:,:),3);
TstepcSP(1:7,:) = mean(Tstep(icSP,:,:),3);
TstepPRF(1:8,:) = mean(Tstep(iPREF,:,:),3);

% step lenghts
sCSLs = vCSLs .* TstepcSL;
sCFRs = vCFRs .* TstepcFR;
sCSPs = vCSPs .* TstepcSP;
SLPRF = vPRFs .* TstepPRF;

% assign work term to conditions
WcSL(1:7,:,:) = Wterm(icSL,:,:);
WcFR(1:7,:,:) = Wterm(icFR,:,:);
WPRF(1:8,:,:) = Wterm(iPREF,:,:);
WcSP(1:7,:,:) = Wterm(icSP,:,:);

% put everything in one big matrix
SLs = [SLPRF(1:8,:); sCFRs(1:7,:); sCSPs(1:7,:); sCSLs(1:7,:)];
SPs_av = [vPRFs(1:8,:); vCFRs(1:7,:); vCSPs(1:7,:); vCSLs(1:7,:)];
SPs = [vPRFs_hs(1:8,:); vCFRs_hs(1:7,:); vCSPs_hs(1:7,:); vCSLs_hs(1:7,:)];
Wmat = -[mean(WPRF(1:8,:,:),3); mean(WcFR(1:7,:,:),3); mean(WcSP(1:7,:,:),3); mean(WcSL(1:7,:,:),3)];
dVs =  [dVPRFs(1:8,:); dVCFRs(1:7,:); dVCSPs(1:7,:); dVCSLs(1:7,:)];

%% Normalize
Mg = parms.mass .* parms.g;
MgL = Mg .* parms.leglength;

%% make gait parameters dimensionless (optional)

if parms.dimensionless == 1
    SLs = SLs ./ repmat(parms.leglength, length(SLs),1);
    SPs = SPs ./ repmat(sqrt(parms.g*parms.leglength), length(SPs),1);  
    SPs_av = SPs_av ./ repmat(sqrt(parms.g*parms.leglength), length(SPs),1);
    Wmat = Wmat ./ repmat(MgL,length(Wmat),1);
end

end

function[SLs, SPs, Wmat, SPs_av, dVs] = get_Wmat(datafolder,Wterm, parms,type)

cd(datafolder)
load('Wsoft.mat','Tstride','vcom_hs','deltav','vcom_hs_alt')

if nargin == 2
    type = [];
end

if strcmp(type,'velocity-based&heelstrike')
    vcom_hs = vcom_hs_alt;
end

% % exclude subject 4 and trial 4 of subject 3
% Wterm(:,4,:) = nan;
% Wterm(4,3,:) = nan;
% Tstride(:,4,:) = nan;
% Tstride(4,3,:) = nan;

%% Reorganize based on freq and SL and speed
N = 9; % number of subjects
M = 10; % numbers of trials

% define variables
vCSLs = nan(M,N);
vCSPs = nan(M,N);
vCFRs = nan(M,N);
vPRFs = nan(M,N);
vCSPs_hs = nan(M,N);

dVCSPs = nan(M,N);
dVCSLs = nan(M,N);
dVCFRs = nan(M,N);
dVPRFs = nan(M,N);

TstepPRF = nan(M,N);
TstepcSP = nan(M,N);
TstepcSL = nan(M,N);
TstepcFR = nan(M,N);

WcSL = nan(M,N,2);
WcFR = nan(M,N,2);
WPRF = nan(M,N,2);
WcSP = nan(M,N,2);

% trials
icSL = [1 4 7 12 15 20 33];
icFR = [3 6 9 10 13 20 31];
iPREF = [2 5 8 11 14 16 20 32];
icSP = [17:20, 23:25];

% corresponding speeds / step lengths
vCSLs(1:7,:) = repmat([.7 .9 1.1 1.6 1.8 1.25 1.4]',1,N);
vCFRs(1:7,:) = repmat([.7 .9 1.1 1.6 1.8 1.25 1.4]',1,N);
vPRFs(1:8,:) = repmat([.7 .9 1.1 1.6 1.8 2.0 1.25 1.4]',1,N);

Tstep = Tstride / 2; % stride - > step

vCSPs = 1.25 * ones(size(vCSPs));

% heelstrike speeds
vCSPs_hs(1:length(icSP),:) = mean(vcom_hs(icSP,:,:),3);
vCSLs_hs(1:length(icSL),:) = mean(vcom_hs(icSL,:,:),3);
vCFRs_hs(1:length(icFR),:) = mean(vcom_hs(icFR,:,:),3);
vPRFs_hs(1:length(iPREF),:) = mean(vcom_hs(iPREF,:,:),3);

% delta velocities
dVCSPs(1:length(icSP),:) = mean(deltav(icSP,:,:),3);
dVCSLs(1:length(icSL),:) = mean(deltav(icSL,:,:),3);
dVCFRs(1:length(icFR),:) = mean(deltav(icFR,:,:),3);
dVPRFs(1:length(iPREF),:) = mean(deltav(iPREF,:,:),3);

% step times
TstepcSL(1:7,:) = mean(Tstep(icSL,:,:),3);
TstepcFR(1:7,:) = mean(Tstep(icFR,:,:),3);
TstepcSP(1:7,:) = mean(Tstep(icSP,:,:),3);
TstepPRF(1:8,:) = mean(Tstep(iPREF,:,:),3);

% step lenghts
sCSLs = vCSLs .* TstepcSL;
sCFRs = vCFRs .* TstepcFR;
sCSPs = vCSPs .* TstepcSP;
SLPRF = vPRFs .* TstepPRF;

% assign work term to conditions
WcSL(1:7,:,:) = Wterm(icSL,:,:);
WcFR(1:7,:,:) = Wterm(icFR,:,:);
WPRF(1:8,:,:) = Wterm(iPREF,:,:);
WcSP(1:7,:,:) = Wterm(icSP,:,:);

% put everything in one big matrix
SLs = [SLPRF(1:8,:); sCFRs(1:7,:); sCSPs(1:7,:); sCSLs(1:7,:)];
SPs_av = [vPRFs(1:8,:); vCFRs(1:7,:); vCSPs(1:7,:); vCSLs(1:7,:)];
SPs = [vPRFs_hs(1:8,:); vCFRs_hs(1:7,:); vCSPs_hs(1:7,:); vCSLs_hs(1:7,:)];
Wmat = -[mean(WPRF(1:8,:,:),3); mean(WcFR(1:7,:,:),3); mean(WcSP(1:7,:,:),3); mean(WcSL(1:7,:,:),3)];
dVs =  [dVPRFs(1:8,:); dVCFRs(1:7,:); dVCSPs(1:7,:); dVCSLs(1:7,:)];

%% Normalize
Mg = parms.mass .* parms.g;
MgL = Mg .* parms.leglength;

%% make gait parameters dimensionless (optional)

if parms.dimensionless == 1
    SLs = SLs ./ repmat(parms.leglength, length(SLs),1);
    SPs = SPs ./ repmat(sqrt(parms.g*parms.leglength), length(SPs),1);  
    SPs_av = SPs_av ./ repmat(sqrt(parms.g*parms.leglength), length(SPs),1);
    Wmat = Wmat ./ repmat(MgL,length(Wmat),1);
end

end

