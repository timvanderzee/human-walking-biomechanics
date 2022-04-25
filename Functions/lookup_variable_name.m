function[variablename,typename, trials] = lookup_variable_name(variable_type, leg, joint, experiment)
%This function returns the variable name for Link_Model_Based Angle,
%Moment, and Power variables
%It also returns a list of the trials for a specific experiment type

typename = 'Link_Model_Based';

if strcmp(joint, 'ankle')
    jointn='ank';
elseif strcmp (joint, 'knee')
    jointn='kne';
elseif strcmp(joint, 'hip')
    jointn='hip';
else
    disp 'unable to identify joint name'
end

variablename= strcat(leg(1), '_', jointn, '_', variable_type);

if strcmp(experiment, 'constant step length')
    trials= [1, 4, 7, 33, 12, 15];
elseif strcmp(experiment, 'constant step frequency')
    trials= [3, 6, 9, 31, 10, 13];
elseif strcmp(experiment, 'constant speed')
    trials= [17:19, 23:25]; %20-22 are preferred but at this speed, could include one of them in this array
elseif strcmp(experiment, 'preferred walking')
    trials= [2, 5, 8, 20, 32, 11, 14, 16];
end

end

            