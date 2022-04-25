function trialnames=lookup_trial_names(trials)
%returns a cell array with the trial descriptions for each trial in the
%input variable: 'trials'

names={'0.7 m/s constant step length',
'0.7 m/s preferred',
'0.7 m/s constant step frequency',
'0.9 m/s constant step length',
'0.9 m/s preferred',
'0.9 m/s constant step frequency',
'1.1 m/s constant step length',
'1.1 m/s preferred',
'1.1 m/s constant step frequency',
'1.6 m/s constant step frequency',
'1.6 m/s preferred',
'1.6 m/s constant step length',
'1.8 m/s constant step frequency',
'1.8 m/s preferred',
'1.8 m/s constant step length',
'2.0 m/s preferred',
'1.25 m/s lowest step frequency',
'1.25 m/s lower step frequency',
'1.25 m/s low step frequency',
'1.25 m/s preferred',
'1.25 m/s preferred',
'1.25 m/s preferred',
'1.25 m/s high step frequency',
'1.25 m/s higher step frequency',
'1.25 m/s highest step frequency',
'1.25 m/s zero step width',
'1.25 m/s 10 cm step width',
'1.25 m/s 20 cm step width',
'1.25 m/s 30 cm step width',
'1.25 m/s 40 cm step width',
'1.40 m/s constant step length',
'1.40 m/s preferred',
'1.40 m/s constant step frequency'};

for i=1:length(trials)
    trialnames{i} = names{trials(i)};
end
end