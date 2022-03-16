function [] = find_5steps(datafolder, subjects, trials)

addpath(genpath(datafolder))

load('5steps_indices.mat','hsl','start')
load('5steps_heelstrikes.mat','hsl_grf','hsr_grf')

% hsl_grf = nan(10,33);
% hsr_grf = nan(10,33);

for subj = subjects
    disp(subj)
    
    for trial = 1:length(trials)
        disp(trials(trial))

        if exist(['p',num2str(subj),'_RawData.mat'],'file')
           load(['p',num2str(subj),'_RawData.mat'],'data')
        else
            continue
        end
        
        if ~isnan(hsl(subj,trials(trial))) && ~isempty(data(trials(trial)).grf.force1)
                         
            %% Find new heelstrike
            grfl = data(trials(trial)).grf.force1(start(subj,trials(trial)):end,1:3);
            grfr = data(trials(trial)).grf.force2(start(subj,trials(trial)):end,1:3);
            
            % Threshold at 20 N
            grfl(grfl(:,3)<20,:) = 0;
            grfr(grfr(:,3)<20,:) = 0;  
            
            % find hsl with different method            
            [LHS, ~, RHS, ~] = invDynGrid_getHS_TO(grfl, grfr, 5);
            
            % make sure they are unique
            LHS = unique(LHS);
            
            % first new left heel strike: choose the closest to the first left heel strike
            [~,closest] = min(abs(LHS-hsl(subj,trials(trial))));
            hsl_close = [LHS(closest) LHS(closest+5)];

            % only consider right heel strikes after the identified left
            % heel strike and pick the first one
            RHS_future = RHS(RHS>hsl_close(2));
            hsr_first = RHS_future(1);
            
            figure ('name', strcat(['subject: ', num2str(subj), ', trial: ', num2str(trials(trial))]));
            plot(grfl(:,3)); hold on
            plot(hsl(subj,trial), 0, 'bo')
            plot(LHS, zeros(size(LHS)), 'bx')
            plot(hsl_close, [0 0], 'b+')
            
            plot(grfr(:,3)); hold on
%             plot(hsl(subj,trial), 0, 'o')
            plot(RHS, zeros(size(RHS)), 'rx')
            plot(hsr_first, [0 0], 'r+')

            % ground reaction force heelstrikes
            hsl_grf(subj,trials(trial)) = hsl_close(1) +  start(subj,trials(trial)) -1;
            hsr_grf(subj,trials(trial)) = hsr_first +  start(subj,trials(trial)) -1;
           
        end
    end
end

cd(datafolder)
save('5steps_heelstrikes.mat','hsl_grf','hsr_grf')

end
