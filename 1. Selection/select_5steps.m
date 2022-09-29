clear all; close all; clc
% Creates .mat files with 5 steps of good quality data from exported data

%% Settings
subjects = 1;
trials = 1:33;

% if you want to redefine start indices (not recommended, but possible)
new_indices = 0;

% folder where the files are that have been exported from Visual3D
% import_folder = uigetdir;
import_folder = 'C:\Users\timvd\Documents\Inverse dynamics\Level 3 - MATLAB files\Level 3 - MATLAB files\V3D exported data';

% folder where to save the 5 strides data
% export_folder = uigetdir;
export_folder = 'C:\Users\timvd\Documents\human-walking-biomechanics\1. Selection';

%% Code
%set new_indices to 1 if you want to redefine the range of indices for the
%entire trial, if not set to 0

% start with what we have already and add to that
cd(export_folder)
if exist('5steps_heelstrikes.mat', 'file')
    load('5steps_heelstrikes.mat','hsl_grf','hsr_grf','hsl','start')
else
    disp('No existing file found, starting from scratch')
end

fs_grf = 1200; % [Hz]

%% Loop
for subj = subjects    
    cd(import_folder)
    if exist(['p',num2str(subj),'_AllStridesData.mat'],'file')
        load(['p',num2str(subj),'_AllStridesData.mat'],'data')
    else 
        disp('No AllStridesData.mat file found')
    end
    
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
        
        % start
        netshift = 0;
        shift = 10;
    
        if isempty(data(trials(trial)).Force)
            disp('No data found')
            continue
        end
        
        % Forces
        f1 = data(trials(trial)).Force.force1;
        f2 = data(trials(trial)).Force.force2;
        
        % start of the walking trial (after jump)
        if new_indices == 1 || ~exist('start','var')
           idx = select_indices(f1,f2);
        else
            idx = start(subj, trials(trial));
        end
        
        % Mocap at half the sample frequency
        idx_mo = round(idx/10);
      
        % extract ground reaction force
        grfl = f1(idx(1):end,:);
        grfr = f2(idx(1):end,:);

        % Threshold at 20 N
        grflt = grfl;
        grfrt = grfr;
        grflt(grflt(:,3)<20,:) = 0;
        grfrt(grfrt(:,3)<20,:) = 0;  
        
        % Joint powers
        pjl = [sum(data(trials(trial)).Link_Model_Based.l_ank_power(idx_mo(1):end,:),2)...
               sum(data(trials(trial)).Link_Model_Based.l_kne_power(idx_mo(1):end,:),2)...
               sum(data(trials(trial)).Link_Model_Based.l_hip_power(idx_mo(1):end,:),2)];
        pjr = [sum(data(trials(trial)).Link_Model_Based.r_ank_power(idx_mo(1):end,:),2)...
               sum(data(trials(trial)).Link_Model_Based.r_kne_power(idx_mo(1):end,:),2)...
               sum(data(trials(trial)).Link_Model_Based.r_hip_power(idx_mo(1):end,:),2)];
           
        % extract heelstrikes and toe-offs from grf
        [LHS, ~, RHS, ~] = invDynGrid_getHS_TO(grflt, grfrt, 5);
        
        % all unique heelstrikes
        hsl_all = unique(LHS);
        hsl_all_mo = unique(round(hsl_all/10) + 1); % heelstrikes for mocop (fs_mocop = .1*fs_grf)
        
        % if heel strike exists, back-calculated the net shift from midpoint
        n = length(hsl_all);
        if isfinite(hsl(subj, trials(trial)))
            [~,i] = min(abs(hsl(subj, trials(trial))-hsl_all));
                netshift = (i)-(ceil(n/2));
        else
            disp('NaN detected')
            netshift=0;
        end

        while shift ~= 0
           
           % determine the chosen heel strikes (hsl_mid)
           if (netshift) > abs(ceil(n/2))
               hsl_mid = hsl_all(1:5);
               keyboard
            elseif (netshift) < (ceil(n/2)*(-1))
                hsl_mid=hsl_all(n-5:n);
                keyboard
           else % normal case
                hsl_mid = hsl_all((ceil(n/2)+netshift):(ceil(n/2)+5+netshift)); 
            end

            hsl_mo_mid = round(hsl_mid/10) + 1;

            %% left leg
            if ishandle(1), close(1); end; figure(1)
            
            %plot GRF of entire stance phase and GRF of 5 steps
            subplot(241);
            plot(pjl); hold on
            xline(hsl_mo_mid(1),'k--')
            xline(hsl_mo_mid(end),'k--')
            ylabel ('Joint power (W)'); title('Entire trial')
                      
            subplot(242);
            plotperstride(pjl, hsl_all_mo); hold on
            plotperstride(pjl, hsl_mo_mid)
            ylabel ('Joint power (W)'); title('Entire trial')
            legend('Ankle','Knee','Hip','location', 'best'); legend boxoff
            
            subplot(243);
            plotperstride(pjl, hsl_mo_mid)
            ylabel ('Joint power (W)'); title('Selected interval')

            subplot(244);
            plotperstride(grfl, hsl_mid)
            ylabel ('Ground reaction force(N)'); title('Selected interval')

            %% right leg
            subplot(245);
            plot(pjr); hold on
            xline(hsl_mo_mid(1),'k--')
            xline(hsl_mo_mid(end),'k--')
            ylabel ('Joint power (W)'); title('Entire trial')

            subplot(246);
            plotperstride(pjr, hsl_all_mo); hold on
            plotperstride(pjr, hsl_mo_mid)
            ylabel ('Joint power (W)'); title('Entire trial')

            subplot(247);
            plotperstride(pjr, hsl_mo_mid)
            ylabel ('Joint power (W)'); title('Selected interval')

            subplot(248);
            plotperstride(grfr, hsl_mid)
            ylabel ('Ground reaction force(N)'); title('Selected interval')
            legend('X','Y','Z','location', 'best'); legend boxoff
            
            for i = 1:8
                subplot (2,4,i)
                xlabel('Frame number')
                axis tight
            end

            set(gcf,'Units','normalized', 'position', [0 .5 1 .5])

            %Ask user if they want to shift, the number entered shifts them
            %that many heelstrikes forward or backward from the current 5 step
            %range (not the first 5 step range)
            shift = input('Press a number to shift (positive for forward, 0 to stop)');
%             close(h)

            % update net shift
            netshift = netshift + shift;

        end

        hsl(subj,trial) = hsl_mid(1);
        start(subj,trial) = idx(1);
  
        %% Find new heelstrike
        % first new left heel strike: choose the closest to the first left heel strike
        [~,closest] = min(abs(LHS-hsl(subj,trials(trial))));
        hsl_close = [LHS(closest) LHS(closest+5)];

        % only consider right heel strikes after the identified left
        % heel strike and pick the first one
        RHS_future = RHS(RHS>hsl_close(2));
        hsr_first = RHS_future(1);

        %% Find 5 steps

        figure ('name', strcat(['subject: ', num2str(subj), ', trial: ', num2str(trials(trial))]));
        plot(grflt(:,3)); hold on
        plot(hsl(subj,trial), 0, 'bo')
        plot(LHS, zeros(size(LHS)), 'bx')
        plot(hsl_close, [0 0], 'g.', 'markersize', 7)
        
        xline(hsl(subj,trial), 'linewidth', 2)
        xline(hsr_first, 'linewidth', 2)
        
        plot(grfrt(:,3)); hold on
        plot(RHS, zeros(size(RHS)), 'rx')
        plot(hsr_first, 0, 'g.', 'markersize', 7)

        % ground reaction force heelstrikes
        hsl_grf(subj,trials(trial)) = hsl_close(1) +  start(subj,trials(trial)) -1;
        hsr_grf(subj,trials(trial)) = hsr_first +  start(subj,trials(trial)) -1;
    end
end  

% save('5steps_heelstrikes.mat','hsl_grf','hsr_grf','hsl','start')

function idx = select_indices(f1,f2)
figure();
plot(f1); hold on; plot(f2)
xlabel ('Force')
ylabel ('Frames')
[idx,~] = ginput(2);
close
idx = round(idx);
end

function [] = plotperstride(x, hsl)

color = get(gca,'colororder');

for k = 1:(length(hsl)-1)
    for j = 1:size(x,2)
        plot(x(hsl(k):hsl(k+1),j), 'color', color(j,:)); hold on
    end
end

end
