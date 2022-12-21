function [] = plot3x3(subj,trial)
% -------------------------------------------------------------------------
% [] = plot3x3(subj,trial)
% 
% INPUTS: 
%       * subj  - subject, input single number not an array of numbers
%       * trial - input single number not an array of numbers
% 
% OUTPUT: a 3x3 figure of the ankle, knee, and hip angle, moment, and power
%         for the left (blue) and right(red) legs avereaged over 5 strides
% 
% -------------------------------------------------------------------------

%% Load and process data
load(strcat('p', num2str(subj), '_5StridesData.mat'),'data');
    
% Extract grfs for each trial
if isempty(data(trial).Force)==1 
   disp('data.Force is empty')
   keyboard
else
    GRFL=[data(trial).Force.force1(:,1),...
                data(trial).Force.force1(:,2),...
                        data(trial).Force.force1(:,3)];
                    
    GRFR=[data(trial).Force.force2(:,1),...
                    data(trial).Force.force2(:,2),...
                            data(trial).Force.force2(:,3)];
end   

% Get heelstrikes for the left and right (hsl, hsr)
[hsl, ~, hsr, ~] = invDynGrid_getHS_TO(GRFL, GRFR,40);
hsr(6) = length(GRFL); % often not found with the function

% divide by ten since mocap was collected at a frequency 10x less
% than grfs were collected
hsl = ceil(hsl/10);
hsr = ceil(hsr/10);

hsl              = unique(hsl); 
hsl(diff(hsl)<5) = [];

hsr              = unique(hsr); 
hsr(diff(hsr)<5) = [];

% Angle and Moment data (x) power (sum of x, y, z) 
% Ajoint(i,j) i=1 ankle, i=2 knee, i=3 hip, j=1 left, j=2 right
Ajoint(1,1,:) = -(mean(interpolate_to_percgaitcycle(    data(trial).Link_Model_Based.l_ank_angle(:,1),  hsl, 201),2,'omitnan'));
Ajoint(1,2,:) = -(mean(interpolate_to_percgaitcycle(    data(trial).Link_Model_Based.r_ank_angle(:,1),  hsr, 201),2,'omitnan'));
Mjoint(1,1,:) = -(mean(interpolate_to_percgaitcycle(    data(trial).Link_Model_Based.l_ank_moment(:,1), hsl, 201),2,'omitnan'));
Mjoint(1,2,:) = -(mean(interpolate_to_percgaitcycle(    data(trial).Link_Model_Based.r_ank_moment(:,1), hsr, 201),2,'omitnan'));
Pjoint(1,1,:) =   mean(interpolate_to_percgaitcycle(sum(data(trial).Link_Model_Based.l_ank_power,2),    hsl, 201),2,'omitnan');
Pjoint(1,2,:) =   mean(interpolate_to_percgaitcycle(sum(data(trial).Link_Model_Based.r_ank_power,2),    hsr, 201),2,'omitnan');

Ajoint(2,1,:) =   mean(interpolate_to_percgaitcycle(    data(trial).Link_Model_Based.l_kne_angle(:,1),  hsl, 201),2,'omitnan');
Ajoint(2,2,:) =   mean(interpolate_to_percgaitcycle(    data(trial).Link_Model_Based.r_kne_angle(:,1),  hsr, 201),2,'omitnan');
Mjoint(2,1,:) =   mean(interpolate_to_percgaitcycle(    data(trial).Link_Model_Based.l_kne_moment(:,1), hsl, 201),2,'omitnan');
Mjoint(2,2,:) =   mean(interpolate_to_percgaitcycle(    data(trial).Link_Model_Based.r_kne_moment(:,1), hsr, 201),2,'omitnan');
Pjoint(2,1,:) =   mean(interpolate_to_percgaitcycle(sum(data(trial).Link_Model_Based.l_kne_power,2),    hsl, 201),2,'omitnan');
Pjoint(2,2,:) =   mean(interpolate_to_percgaitcycle(sum(data(trial).Link_Model_Based.r_kne_power,2),    hsr, 201),2,'omitnan');

Ajoint(3,1,:) = -(mean(interpolate_to_percgaitcycle(    data(trial).Link_Model_Based.l_hip_angle(:,1),  hsl, 201),2,'omitnan'));
Ajoint(3,2,:) = -(mean(interpolate_to_percgaitcycle(    data(trial).Link_Model_Based.r_hip_angle(:,1),  hsr, 201),2,'omitnan'));
Mjoint(3,1,:) = -(mean(interpolate_to_percgaitcycle(    data(trial).Link_Model_Based.l_hip_moment(:,1), hsl, 201),2,'omitnan'));
Mjoint(3,2,:) = -(mean(interpolate_to_percgaitcycle(    data(trial).Link_Model_Based.r_hip_moment(:,1), hsr, 201),2,'omitnan'));
Pjoint(3,1,:) =   mean(interpolate_to_percgaitcycle(sum(data(trial).Link_Model_Based.l_hip_power,2),    hsl, 201),2,'omitnan');
Pjoint(3,2,:) =   mean(interpolate_to_percgaitcycle(sum(data(trial).Link_Model_Based.r_hip_power,2),    hsr, 201),2,'omitnan');

%% Plot 
figure('name', (['3x3 plot for Subject ', num2str(subj),', Trial: ', num2str(trial)]))

for i = 1:9
    subplot(3,3,i);
    plot([0 100], [0 0], 'k-'); hold on
    xlabel('% Gait cycle'); box off
end

% --------------------------------------------------------- Change y limits

for i=1:3
subplot(3,3,i);
plot(linspace(0,100,201), squeeze(Ajoint(i,:,:)),'linewidth',2);

ylim([-70 25]); ylabel('Angle (deg)')

subplot(3,3,3+i);
plot(linspace(0,100,201), squeeze(Mjoint(i,:,:)),'linewidth',2);

ylim([-1 2]); ylabel('Moment (N-m/kg)')

subplot(3,3,6+i);
plot(linspace(0,100,201), squeeze(Pjoint(i,:,:)),'linewidth',2);

ylim([-2 3]); ylabel('Power (W/kg)')
end

subplot(331); title('Ankle');
subplot(332); title('Knee');
subplot(333); title('Hip');

end

