% script to compare heelstrikes.mat outputs

load('5steps_heelstrikes_all4.mat')
hsl_old=hsl;
hsl_grf_old=hsl_grf;
hsr_grf_old=hsr_grf;
start_old=start;

select_5steps(1,1:5)
pause(10)

load('5steps_heelstrikes_NEW.mat')

delta1=hsl_old(1,1:10)-hsl(1,1:10);
delta2=hsl_grf_old(1,1:10)-hsl_grf(1,1:10);
delta3=hsr_grf_old(1,1:10)-hsr_grf(1,1:10);
delta4=start_old(1,1:10)-start(1,1:10);

figure();
plot(delta1); hold on
plot(delta2, 'mo'); hold on
plot(delta3, 'g--'); hold on
plot(delta4, 'c+'); hold on

