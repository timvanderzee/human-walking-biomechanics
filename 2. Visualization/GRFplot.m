function [] = GRFplot(subjects, trials)
% -------------------------------------------------------------------------
% [] = GRFplot(subjects, trials)
%
% creates a figure with 3 subplots for GRF. It plots the left and right
% leg one on top of the other and the 3 subplots are the GRF in the
% medial/lateral, anterior/posterior, and vertical directions.
% 
% Input:
%   * subjects - single number from 1 to 9 or an array including any
%                combinations of subjects from 1 to 9
%   * trials   - single trial number from 1:33 or an array of trials
%
% Ouputs: For each subject/trial combination, a single plot will appear
%
% -------------------------------------------------------------------------



% loops over subject and loads subject's data and then loops over trials
for subj = subjects
    load(['p',num2str(subj),'_5StridesData'],'data')
    
    for trial = trials
% Figure: 5 stride figure: ground reaction forces

        figure('name', (['GRF for Subject ', num2str(subj),', Trial: ', num2str(trial)]))
        
        title (['GRF for Subject ', num2str(subj),', trial ', num2str(trial)])
        if isnan(data(trial).Force.force1)
            subplot(3,1,1)
            plot (0,0)
            subplot(3,1,2)
            subplot (3,1,3)
            
        elseif ~isnan(data(trial).Force.force1)
            subplot(3,1,1) %plot for medial/lateral GRFs
            plot(data(trial).Force.force1(:,1)); hold on;
            plot(data(trial).Force.force2(:,1)); hold on;
            ylabel ('Force (N)')
            title('Medial/Lateral')
            
            subplot(3,1,2) %plot for anterior/posterior GRFs
            plot(data(trial).Force.force1(:,2)); hold on;
            plot(data(trial).Force.force2(:,2)); hold on;
            ylabel ('Force (N)')
            title('Anterior/Posterior')
            
            subplot(3,1,3) %plot for vertical GRFs
            plot(data(trial).Force.force1(:,3)); hold on;
            plot(data(trial).Force.force2(:,3)); hold on;
            ylabel ('Force (N)')
            title ('Vertical')
        end
        for k=1:3
            subplot(3,1,k)
            xlabel('Frame Number')
        end
        legend ('Left Leg', 'Right Leg', 'location', 'best')
    end
end
        
end

