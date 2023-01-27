function [] = visualizejointkin(trial, joint, direction)
% -------------------------------------------------------------------------
% [] = visualizejointkin(trial, joint, direction)
%
% creates a plot of each subjects data for a single variable and trial
% 
% INPUTS: 
%   * trial     - single number for the trial of interest
%   * joint     - can be defined as 'ankle', 'knee', or 'hip'
%   * direction - plane to look at: 1 for data in the sagittal plane
%                                   2 for data in the frontal plane
%                                   3 for data in the transverse plane
% 
% OUTPUT: a figure with 6 subplots showing the joint angles, moment, and
%         power for the left leg and the right leg. Each subplot have the 
%         mean over subjects and each subjects' individual data
%
% -------------------------------------------------------------------------

% Set current directory (cd) to data folder
subjects = 1:9;
trials   = 1:33;
npoints  = 101; % change this value to change points for interpolation

% Define variables
Ajoint_pc = nan(npoints, subjects(end), 2);
Mjoint_pc = nan(npoints, subjects(end), 2);
Pjoint_pc = nan(npoints, subjects(end), 2);

% Loop over subjects and load subject data 
for subj = subjects

    disp(strcat(['Loading data for Subject:', ' ', num2str(subj)]))
    load(['p',num2str(subj),'_5Stridesdata.mat'],'data')

        
    if isempty(data(trials(trial)).Force)==1
        continue
    else        
%Define joint Moment, Angle, and Power for the given joint
        if strcmp(joint,'ankle')
        Mjoint = [data(trial).Link_Model_Based.l_ank_moment(:,direction)...
                  data(trial).Link_Model_Based.r_ank_moment(:,direction)];

        Ajoint = [data(trial).Link_Model_Based.l_ank_angle(:,direction)...
                  data(trial).Link_Model_Based.r_ank_angle(:,direction)];

        Pjoint = [data(trial).Link_Model_Based.l_ank_power(:,direction)...
                  data(trial).Link_Model_Based.r_ank_power(:,direction)];

        elseif strcmp(joint,'knee')
        Mjoint = [data(trial).Link_Model_Based.l_kne_moment(:,direction)...
                  data(trial).Link_Model_Based.r_kne_moment(:,direction)];

        Ajoint = [data(trial).Link_Model_Based.l_kne_angle(:,direction)...
                  data(trial).Link_Model_Based.r_kne_angle(:,direction)];

        Pjoint = [data(trial).Link_Model_Based.l_kne_power(:,direction)...
                  data(trial).Link_Model_Based.r_kne_power(:,direction)];

        elseif strcmp(joint,'hip')
        Mjoint = [data(trial).Link_Model_Based.l_hip_moment(:,direction)...
               data(trial).Link_Model_Based.r_hip_moment(:,direction)];

        Ajoint = [data(trial).Link_Model_Based.l_hip_angle(:,direction)...
                data(trial).Link_Model_Based.r_hip_angle(:,direction)];

        Pjoint = [data(trial).Link_Model_Based.l_hip_power(:,direction)...
                data(trial).Link_Model_Based.r_hip_power(:,direction)];

        else
            disp('no joint defined')
        end
            
    %% extract heel strikes
    % Extract GRF variable
    GRFL=[data(trials(trial)).Force.force1(:,1),...
                data(trials(trial)).Force.force1(:,2),...
                        data(trials(trial)).Force.force1(:,3)];

    GRFR=[data(trials(trial)).Force.force2(:,1),...
                data(trials(trial)).Force.force2(:,2),...
                        data(trials(trial)).Force.force2(:,3)];

    % Get heelstrikes
    [hsl, ~, hsr, ~] = invDynGrid_getHS_TO(GRFL, GRFR,40);
    hsr(6) = length(GRFL); % often not found with the function

    % divide by ten since mocap was collected at a frequency 10x less
    % than grfs were collected
    hsl=ceil(hsl/10);
    hsr=ceil(hsr/10);

    hsl = unique(hsl); 
    hsl(diff(hsl)<5) = [];

    hsr = unique(hsr); 
    hsr(diff(hsr)<5) = [];

    hs = [hsl hsr];
            
        %% interpolate to percentage gait cycle (101 points)
        if sum(isnan(Mjoint(:)))<1 && sum(isnan(Pjoint(:)))<1
            for L = 1:2
                Ajoint_pc(:,subj,L) = mean(interpolate_to_percgaitcycle(...
                                 Ajoint(:,L),hs(:,L),npoints),2,'omitnan');
                             
                Mjoint_pc(:,subj,L) = mean(interpolate_to_percgaitcycle(...
                                 Mjoint(:,L),hs(:,L),npoints),2,'omitnan');
                             
                Pjoint_pc(:,subj,L) = mean(interpolate_to_percgaitcycle(...
                                 Pjoint(:,L),hs(:,L),npoints),2,'omitnan');
            end
        end
    end
% Take the mean over the subjects
Ajoint_pcmean = mean(Ajoint_pc, 2, 'omitnan');
Mjoint_pcmean = mean(Mjoint_pc, 2, 'omitnan');
Pjoint_pcmean = mean(Pjoint_pc, 2, 'omitnan');
        
end

% Variables for Plotting
% colours for plotting different subjects' data
colour = jet(11);
c      = [2,4:12];

% npoints can be changed to change the nunber of points for interpolation.
% xa is set below so that the x axis is still from 1-100% gait cycle
xa = linspace(0,npoints,101);

% define string for the plane of data being plotted
switch direction
    case 1
        varplane = 'Sagittal plane';
    case 2
        varplane = 'Frontal plane';
    case 3
        varplane = 'Transverse plane';
end

%% Plot
figure('name',['Rotational Joint kin in the ' varplane, ' for trial: ',...
           num2str(trial)],'units','normalized','outerposition',[0 0 .5 1])

% ------------------------------------------ Plot the means of the subjects
subplot(321)
plot(xa, Ajoint_pcmean(:,1), 'linewidth', 2, 'color', 'k'); 
hold on

subplot(322)
plot(xa, Ajoint_pcmean(:,2), 'linewidth', 2, 'color', 'k'); 
hold on

subplot(323)
plot(xa, Mjoint_pcmean(:,1), 'linewidth', 2, 'color', 'k'); 
hold on

subplot(324)
plot(xa, Mjoint_pcmean(:,2), 'linewidth', 2, 'color', 'k'); 
hold on

subplot(325)
plot(xa, Pjoint_pcmean(:,1), 'linewidth', 2, 'color', 'k'); 
hold on

subplot(326)
plot(xa, Pjoint_pcmean(:,2), 'linewidth', 2, 'color', 'k'); 
hold on


% ----------------------------------- Plot data for each individual subject
for s=1:length(subjects)
    if isempty(Mjoint_pc(:,s,1))==1
        continue
    else

        % left joint angle
        subplot(321)
        plot(xa, Ajoint_pc(:,s,1), 'linewidth', 1.5, 'color', colour(c(s),:)); 
        hold on
        title(strcat(['Left ', joint, ' Angle']))
        axis tight; yl1(1,:) = get(gca,'ylim');
        ylabel('Angle (deg)');

        % right joint angle
        subplot(322)
        plot(xa, Ajoint_pc(:,s,2), 'linewidth', 1.5,'color', colour(c(s),:)); 
        hold on
        title(strcat(['Right ', joint, ' Angle']))
        axis tight; yl1(2,:) = get(gca,'ylim');
        ylabel(' Angle (deg)');
        
        % left joint moment
        subplot(323)
        plot(xa, Mjoint_pc(:,s, 1),'linewidth', 1.5,'color', colour(c(s),:)); 
        hold on
        title(strcat(['Left ', joint, ' Moment']))
        axis tight; yl2(1,:) = get(gca,'ylim');
        ylabel('Moment (N-m/kg)');

        % right joint moment
        subplot(324)
        plot(xa, Mjoint_pc(:,s,2), 'linewidth', 1.5, 'color', colour(c(s),:)); 
        hold on
        title(strcat(['Right ', joint, ' Moment']))
        axis tight; yl2(2,:) = get(gca,'ylim');
        ylabel('Moment (N-m/kg)');

        % left joint power
        subplot(325)
        plot(xa, Pjoint_pc(:,s,1), 'linewidth', 1.5, 'color', colour(c(s),:));
        hold on
        title(strcat(['Left ', joint, ' Power']))
        axis tight; yl3(1,:) = get(gca,'ylim');
        ylabel('Power (W/kg)')

        % right joint power
        subplot(326)
        plot(xa, Pjoint_pc(:,s,2), 'linewidth', 1.5, 'color', colour(c(s),:)); 
        hold on
        title(strcat(['Right ', joint, ' Power']))
        axis tight; yl3(2,:) = get(gca,'ylim');
        ylabel('Power (W/kg)');
        
    end
end
        

% ----------------------------------------------------- Limit and name axes
        for i = 1:6
            subplot(3,2,i)
            xlim([0 npoints-1])
            xlabel('Percentage gait cycle');
            grid on
        end
        
        for i=1:2
            subplot(3,2,i)
            ylim([min(yl1(:,1)), max(yl1(:,2))]);
        end
        
        for i=3:4
            subplot(3,2,i)
            ylim([min(yl2(:,1)), max(yl2(:,2))]);
        end
        
        for i=5:6
            subplot(3,2,i)
            ylim([min(yl3(:,1)), max(yl3(:,2))]);
        end
        
        for i=1:4
            subplot(3,2,i)
            yyaxis right
            set(gca, 'YTick', [], 'YColor', 'k')
            ylabel('<---- Ext. ------------ Flex. ---->');
        end

legend('Mean','1','2','3','4','5','6','7','8','9','location','best')


end

