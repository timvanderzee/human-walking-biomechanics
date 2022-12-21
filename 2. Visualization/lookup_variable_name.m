function[variablename, typename, trials] = ...
                lookup_variable_name(variable_type, leg, joint, experiment)
% -------------------------------------------------------------------------
% [variablename, typename, trials] = ...
%               lookup_variable_name(variable_type, leg, joint, experiment)
% 
% returns the variable name for Link_Model_Based Angle, Moment, and Power 
% variables and a list of the trials for a specific experiment type
% 
% Inputs:  - variable_type: 'angle', 'moment', 'power'
%          - leg:           'right' or 'left'
%          - joint:         'ankle', 'knee', or 'hip'
%
%          - experiment:    'constant step length',
%                           'constant step frequency', 
%                           'constant speed',
%                           'preferred walking'
% 
% Output:  - variablename: string continaing the name of the variable which
%                          is found in the data files
%          - typename:     string containing the name of the field where 
%                          the variable can be found
%          - trials:       variable that contains the trial numbers for 
%                          trials matching the experiment defined above. 
%                          The trials are sorted from lowest speed to 
%                          highest speed
% 
% -------------------------------------------------------------------------

typename = 'Link_Model_Based';

switch joint
    case 'ankle' 
        jointn='ank';
        
    case 'knee'
        jointn='kne';
        
    case 'hip'
        jointn='hip';
        
    otherwise
        warning('unable to identify joint name');
end


variablename = strcat(leg(1), '_', jointn, '_', variable_type);

switch experiment
    case 'constant step length'
        trials= [1, 4, 7, 33, 12, 15];
        
    case 'constant step frequency'
        trials= [3, 6, 9, 31, 10, 13];
        
    case 'constant speed'
        trials= [17:19, 23:25]; 
        % 20-22 are preferred but at this speed, could include one of them 
        % in this array
        
    case 'preferred walking'
        trials= [2, 5, 8, 20, 32, 11, 14, 16];
        
    otherwise 
        warning('unable to identify experiment number');
end
end

        
