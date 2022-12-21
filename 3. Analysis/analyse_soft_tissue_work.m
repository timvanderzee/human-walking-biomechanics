% -------------------------------------------------------------------------
% analyse_soft_tissue_work.m 
% -------------------------------------------------------------------------

folder = fullfile(datafolder,'All Strides Data files');

% folder where the files are that have been exported from Visual3D
% % folder = uigetdir;
% if ~exist('folder','var')
%     folder = uigetdir;
% end
% cd(folder);

if ishandle(1), close(1); end; figure(1)

%% Parameters
cd(folder)
load(['Wsoft_',date,'.mat'])
terms = {'Collision', 'Soft tissue', 'Negative'};
old_order = 0;

% in van der Zee & Kuo (2021)
coefficients_in_paper = [4.512, 2.847, 5.083];
SE_in_paper = [0.079, 0.059, 0.108];

% parameters
parms.leglength = [0.8941 0.9398 1.04 0.876 0.8636 0.9398 0.9906 0.99 0.9398]; % m
parms.mass = [81.8000 57.3000 97.5000 57 56.7000 72.6000 86.2000 88.6000 77]; % kg

parms.idx = [1 9 16 23 30];
parms.g = 9.81; % m/s2
parms.dimensionless = 1;
subjs = 1:9;

%% Loop over work terms
for j = 3:-1:1
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

% reorganize
[~, SPs, Wmat, SPs_av, dVs] = get_Wmat(Wterm, parms,'velocity-based&heelstrike');

 % dVs is total redirection, and we assume collision takes up half
Delta = tan(dVs * pi/180 / 2);

% Eq. 4-6 in van der Zee & Kuo (2021)
X = .5 .* SPs(:,subjs).^2 .* Delta(:,subjs).^2;
Y = Wmat(:,subjs);
S = repmat((1:size(X,2))',1, size(X,1))';

% lump together in one big table
T = array2table([(Y(1:numel(Y)))' (X(1:numel(X)))' (S(1:numel(S)))']);
T.Properties.VariableNames = {'Work' 'GaitPar' 'Subject'};

% do linear regression!
LM = fitlm(T,modelString);

% display coefficients
disp(['Coefficient for ',term,': ', ...
                num2str(round(LM.Coefficients.Estimate(i),3)), ' +- ',...
                                  num2str(round(LM.Coefficients.SE(i),3))])
                              
disp(['Coef. in paper for ',term,': ', num2str(coefficients_in_paper(j)),...
                                         ' +- ',  num2str(SE_in_paper(j))])

%% 3D plots
figure(1)
color = get(gca,'colororder');
subplot(1,3,j);
[~, RMSE_R2(j,:), ~] = make_3D_plot(Wterm, parms,'Dofit', color,...
                                          'velocity-based&heelstrike', LM);

titles = {'Whole-body collision', 'Soft tissue collision', ...
                                                      'Whole-body stride'};

for i = 1:3
    subplot(1,3,j);
    set(gca,'Xtick', 0:.3:1.2,'Ytick',0:.2:.8);
    xlabel('Step length');ylabel('Speed'); zlabel('Dissipation');
    axis([0 .45 0 .8 0 .3]); view(-45,20); title(titles{j});
    text(0, 0.7, 0.25, ['R^2 = ', num2str(round(RMSE_R2(j,2),2))])
end

set(gcf,'units','normalized','position', [.1 .3 .9 .4])

end
    

%% Functions
function[Wmat, RMSE_R2, LM] = make_3D_plot(Wterm, parms, optional_fit, color, type, LM)
 
if nargin == 2
    optional_fit = nan;
end

if nargin < 5
    type = 'impact';
end

if isfield(parms,'colors')
    colors = parms.colors;
else
    figure;
    colors = get(gca,'colororder');
    close
end

%% make 3D plot
idx = parms.idx;
subjs = 1:9;
xgrid = (0:.05:1.2)'; 
% xfine = (0:.01:1.2)';

[SP_grid,SL_grid] = meshgrid(xgrid); % Generate x and y data
SW = @(C,SP,SL) C*SL.^2.*SP.^2;
[SLs, SPs, Wmat,SPs_av, dVs] = get_Wmat(Wterm, parms,type);

if strcmp(type(1:7), 'average') == 1
    SPs = SPs_av;
    c = 1/8;
end

if strcmp(type, 'average&exclude') == 1
    Wmat(idx(3):idx(4)-1,:) = nan;
end

if strcmp(type, 'velocity-based') == 1
    SLs = dVs /180*pi;
end

if strcmp(type, 'velocity-based&heelstrike') == 1
    SLs = tan(dVs /180*pi/2);
    c = .5;
end

stats.coef = LM.Coefficients.Estimate;
coef = stats.coef(end) * c;

if length(stats.coef) == 1
    offset = 0;
else
    offset = stats.coef(1);
end

%% Plotting
symb = {'o','s','d','v'};
fillcolor = colors.^2;

if strcmp(optional_fit,'Dofit')
   hold on; grid on
   
   for i = 1:4
       SLsvec = SLs(idx(i):idx(i+1)-1,:); 
       SPsvec = SPs(idx(i):idx(i+1)-1,:); 
       Wvec = Wmat(idx(i):idx(i+1)-1,:);
       
       h = plot3(SLsvec(:), SPsvec(:), Wvec(:), symb{i});
       set(h,'marker',symb{i},'color',color(i,:),'markersize',7,...
                                         'markerfacecolor',fillcolor(i,:)); 
       hold on
       
   end
    
    % zero plane
    patch([1 1 -1 -1],[1 -1 -1 1], [0 0 0 0], 'LineStyle', 'None', ...
                                'FaceColor', 'k','FaceAlpha', .05); 
    hold on
    
    % plot surface
    surf(SL_grid,SP_grid,SW(coef, SP_grid,SL_grid) + offset,color(1,:),...
          'FaceColor', color(1,:),'FaceAlpha', .2, 'EdgeColor', color(1,:))

    % make nice
    ylabel('Speed'); 
    xlabel('Step length'); 
    zlabel('Soft tissue dissipation')
end

%% RMSE & R2
prediction = SW(coef, SPs,SLs) + offset;
actual = Wmat(:,subjs);
pred = prediction(:,subjs);

RMSE_R2(1) = rms(actual(:)-pred(:),'omitnan');
RES = sum((actual(:)-pred(:)).^2,'omitnan');    
TOT = sum((actual(:)-mean(actual(:),'omitnan')).^2,'omitnan');
RMSE_R2(2) = 1 - RES/TOT; 
end

function[SLs, SPs, Wmat, SPs_av, dVs] = get_Wmat(Wterm, parms,type)

% cd(datafolder)
load(['Wsoft_',date,'.mat'],'Tstride','vcom_hs','deltav','vcom_hs_alt')

if nargin == 2
    type = [];
end

if strcmp(type,'velocity-based&heelstrike')
    vcom_hs = vcom_hs_alt;
end

%% Reorganize based on freq and SL and speed
M = 10; % numbers of trials
N = 9;  % number of subjects

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


