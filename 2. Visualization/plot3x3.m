clear all; close all; clc

% add code folder to datapath
scriptname = 'plot3x3.m';
scriptfile = which(scriptname);
scriptfolder = scriptfile(1:(end-length(scriptname)));
addpath(genpath(scriptfolder));
cd(scriptfolder)


trials=27;

for subj=1:9
    disp(strcat('Subj:', num2str(subj)))
    load(strcat('p', num2str(subj), '_5StepsData.mat'))
    
    %% Extract grfs for each trial

    for trial=1:length(trials) %Loops over trials
        disp(trials(trial))
        
        if isempty(data(trials(trial)).grf)==1 
            continue
        elseif isempty(data(trials(trial)).grf)==1 
            continue
        elseif isnan(data(trials(trial)).grf.force1)==1 
            continue
        elseif isnan(data(trials(trial)).grf.force2)==1 
            continue
        else 
            GRFL=[data(trials(trial)).grf.force1(:,1), data(trials(trial)).grf.force1(:,2), data(trials(trial)).grf.force1(:,3)];
            GRFR=[data(trials(trial)).grf.force2(:,1), data(trials(trial)).grf.force2(:,2), data(trials(trial)).grf.force2(:,3)];
        end   
        % Get heelstrikes for the left and right (hsl, hsr)
        [hsl, tol, hsr, tor] = invDynGrid_getHS_TO(GRFL, GRFR,40);
        hsr(6) = length(GRFL); % often not found with the function

        % divide by ten since mocap was collected at a frequency 10x less
        % than grfs were collected
        hsl=ceil(hsl/10);
        hsr=ceil(hsr/10);
        
        hsl = unique(hsl); hsl(find(diff(hsl)<5)) = [];
        hsr = unique(hsr); hsr(find(diff(hsr)<5)) = [];
        
        %var1(:, trial,subject)=mean(interpolate_to_percgaitcycle(data(trials(trial)).(field).(variable)(:,col), HS,101),2,'omitnan');
    
        for direction=1:3
            r_ank_angle(subj, trial, direction,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.r_ank_angle(:,direction), hsr, 101),2,'omitnan');
            r_kne_angle(subj, trial, direction,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.r_kne_angle(:,direction), hsr, 101),2,'omitnan');
            r_hip_angle(subj, trial, direction,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.r_hip_angle(:,direction), hsr, 101),2,'omitnan');
            r_ank_moment(subj, trial, direction,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.r_ank_moment(:,direction), hsr, 101),2,'omitnan');
            r_kne_moment(subj, trial, direction,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.r_kne_moment(:,direction), hsr, 101),2,'omitnan');
            r_hip_moment(subj, trial, direction,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.r_hip_moment(:,direction), hsr, 101),2,'omitnan');
            r_ank_power(subj, trial, direction,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.r_ank_power(:,direction), hsr, 101),2,'omitnan');
            r_kne_power(subj, trial, direction,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.r_kne_power(:,direction), hsr, 101),2,'omitnan');
            r_hip_power(subj, trial, direction,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.r_hip_power(:,direction), hsr, 101),2,'omitnan');

            l_ank_angle(subj, trial, direction,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.l_ank_angle(:,direction), hsl, 101),2,'omitnan');
            l_kne_angle(subj, trial, direction,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.l_kne_angle(:,direction), hsl, 101),2,'omitnan');
            l_hip_angle(subj, trial, direction,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.l_hip_angle(:,direction), hsl, 101),2,'omitnan');
            l_ank_moment(subj, trial, direction,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.l_ank_moment(:,direction), hsl, 101),2,'omitnan');
            l_kne_moment(subj, trial, direction,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.l_kne_moment(:,direction), hsl, 101),2,'omitnan');
            l_hip_moment(subj, trial, direction,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.l_hip_moment(:,direction), hsl, 101),2,'omitnan');
            l_ank_power(subj, trial, direction,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.l_ank_power(:,direction), hsl, 101),2,'omitnan');
            l_kne_power(subj, trial, direction,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.l_kne_power(:,direction), hsl, 101),2,'omitnan');
            l_hip_power(subj, trial, direction,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.l_hip_power(:,direction), hsl, 101),2,'omitnan');
        end
    end
end

%% Sum planes and find mean of left/right and then for subjects
keyboard
r_ank_angle_all(subj, trial,:)=sum(r_ank_angle, 3)
r_kne_angle_all(subj, trial,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.r_kne_angle(:,direction), hsr, 101),2,'omitnan');
r_hip_angle_all(subj, trial,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.r_hip_angle(:,direction), hsr, 101),2,'omitnan');
r_ank_moment_all(subj, trial,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.r_ank_moment(:,direction), hsr, 101),2,'omitnan');
r_kne_moment_all(subj, trial,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.r_kne_moment(:,direction), hsr, 101),2,'omitnan');
r_hip_moment_all(subj, trial,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.r_hip_moment(:,direction), hsr, 101),2,'omitnan');
r_ank_power_all(subj, trial,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.r_ank_power(:,direction), hsr, 101),2,'omitnan');
r_kne_power_all(subj, trial,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.r_kne_power(:,direction), hsr, 101),2,'omitnan');
r_hip_power_all(subj, trial,:)=mean(interpolate_to_percgaitcycle(data(trials(trial)).Link_Model_Based.r_hip_power(:,direction), hsr, 101),2,'omitnan');

