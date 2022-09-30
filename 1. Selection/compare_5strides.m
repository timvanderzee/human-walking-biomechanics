clear all; close all; clc

version1folder = 'C:\Users\timvd\Documents\Inverse dynamics\Level 3 - MATLAB files\Level 3 - MATLAB files - reproduced\5 Strides Data files from process_5steps_new';
version2folder = 'C:\Users\timvd\Documents\Inverse dynamics\Level 3 - MATLAB files\Level 3 - MATLAB files - reproduced\5 Strides Data files from corrected process_5steps';

cd(version1folder)
load('p1_5StridesData.mat')
data_v1 = data;

cd(version2folder)
load('p1_5StridesData.mat')
data_v2 = data;


%%
trial = 1;

fn=fieldnames(data(trial));


for i=1: numel(fn)
    fn1=fieldnames(data(trial).(fn{i}));

    for j=1: numel(fn1)

        if isstruct(data(trial).(fn{i}).(fn1{j}))
            fn2=fieldnames(data(trial).(fn{i}).(fn1{j}));
            for k=1: numel(fn2)
                %access the data
                
                figure(1); 
                set(gcf,'name', [fn{i},fn1{j}, fn2{k}]);
                subplot(121); plot(data_v1(trial).(fn{i}).(fn1{j}).(fn2{k}))
                subplot(122); plot(data_v2(trial).(fn{i}).(fn1{j}).(fn2{k}))
                
                pause
            end 
            
        else
            if ~isempty(data(trial).(fn{i}).(fn1{j})) && size(data(trial).(fn{i}).(fn1{j}),3) < 2 && isnumeric(data(trial).(fn{i}).(fn1{j})(1))
                figure(1); 
                set(gcf,'name', [fn{i},fn1{j}]);
                subplot(121); plot(data_v1(trial).(fn{i}).(fn1{j}))
                subplot(122); plot(data_v2(trial).(fn{i}).(fn1{j}))
                
                pause
               
            end
        end
    end
end 