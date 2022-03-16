function[] = Wsoft_to_coef(datafolder)
cd(datafolder)
load('Wsoft.mat')
terms = {'Collision', 'Soft tissue', 'Negative'};

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

% parameters
parms.leglength = [0.8941 0.9398 1.04  0.8636 0.8636 0.9398 0.9906 0.99 0.9398 0.876]; % m
parms.mass = [81.8000 57.3000 97.5000 60.3000 56.7000 72.6000 86.2000 88.6000 77 57]; % kg

% 1-9
parms.leglength = [0.8941 0.9398 1.04 0.876 0.8636 0.9398 0.9906 0.99 0.9398]; % m
parms.mass = [81.8000 57.3000 97.5000 57 56.7000 72.6000 86.2000 88.6000 77]; % kg

parms.g = 9.81; % m/s2
parms.dimensionless = 1;

% reorganize
[~, SPs, Wmat, SPs_av, dVs] = get_Wmat(datafolder,Wterm, parms,'velocity-based&heelstrike');

 % dVs is total redirection, and we assume collision takes up half
Delta = tan(dVs * pi/180 / 2);

% Eq. 4-6 in van der Zee & Kuo (2021)
subjs = 1:9;
X = .5 .* SPs(:,subjs).^2 .* Delta(:,subjs).^2;
Y = Wmat(:,subjs);
S = repmat((1:size(X,2))',1, size(X,1))';

% lump together in one big table
T = array2table([(Y(1:numel(Y)))' (X(1:numel(X)))' (S(1:numel(S)))']);
T.Properties.VariableNames = {'Work' 'GaitPar' 'Subject'};

% do linear regression!
LM = fitlm(T,modelString);

% display coefficients
disp(['Coefficient for ',term,': ', num2str(LM.Coefficients.Estimate(i)), ' +- ',  num2str(LM.Coefficients.SE(i))])
end


end