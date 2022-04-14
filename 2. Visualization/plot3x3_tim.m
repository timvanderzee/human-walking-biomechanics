close all; clear all;

cd('C:\Users\tim.vanderzee\OneDrive - University of Calgary\2. Soft tissue - JEB\processed_data')

%%
close all

figure(1);
for i = 1:9
    subplot(3,3,i);
    plot([0 100], [0 0], 'k-'); hold on
    xlabel('% Gait cycle'); box off
end

for i = 1:3
    if i == 1, load('Ankle_PG.mat')
    elseif i == 2,  load('Knee_PG.mat')
    else,  load('Hip_PG.mat')
    end
    
    if i == 1 || i == 3
        Ajoint_pc = -Ajoint_pc;
        Mjoint_pc = -Mjoint_pc;
    end
    

figure(1);

subplot(3,3,i);
plot(linspace(0,100,201), squeeze(Ajoint_pc(:,1,20,:)),'linewidth',2);

ylim([-70 25]); ylabel('Angle (deg)')

subplot(3,3,3+i);
plot(linspace(0,100,201), squeeze(Mjoint_pc(:,1,20,:)),'linewidth',2);

ylim([-1 2]); ylabel('Moment (N-m/kg)')

subplot(3,3,6+i);
plot(linspace(0,100,201), squeeze(Pjoint_pc(:,1,20,:)),'linewidth',2);

ylim([-2 3]); ylabel('Power (W/kg)')
end

subplot(331); title('Ankle');
subplot(332); title('Knee');
subplot(333); title('Hip');
