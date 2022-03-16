function out = comvelocityevents(Vcom,hstohstohs)
% out = comvelocityevents(Vcom,hstohstohs);
% A function to return the indices, Velocities at, Symbols for different
% events of interest in a hodograph. Output is a structure.  
% 
% Vcom = COM velocity in 3D, where X is Lateral, Y is Forward and Z is Up. 
% 
% hstohstohs = a vector of 2*nsteps+1 indices indicating Heel Strikes and
% Toe-Offs within the Vcom signal. "nsteps" is the number of STEPS in the
% Vcom signal. The first hs is the initial heel strike,  and the last hs is
% the "closing" heel strike, usually coincident with the last sample of Vcom

n = size(Vcom,1);

% Determine the number of STEPS in the Vcom signal: 
% length(hstohstohs) = 2*nsteps+1. 
% CURRENTLY LIMITED to 2 steps! (never tried with more) 2009-11-23 PGA
% if length(hstohstohs) == 3
%     nsteps = 1;
% elseif length(hstohstohs) == 5; 
%     nsteps = 2;
% else
%     error('Wrong number of steps for "comvelocitymetrics"')
% end

% Try to make it work with more steps 2020-06-25 TJZ
nsteps = (length(hstohstohs)-1)/2;

% create "vdup", a duplicate of Vcom that will have a bunch of
% regions-of-no-interest blanked out with NaN at different times.
vdup = Vcom;

% parse heel strikes and toe-offs 
hs = hstohstohs(1:2:end);
to = hstohstohs(2:2:end); 
    to(end+1) = to(1)+hs(end)-1; % add an after-step toe off, just to give a programmatic index for NaNing data for vdn

vcomevents = hs(1);  % indices of events of interest: first is the initial heel strike
vcomeventsyms = {'s'}; % symbols used to mark these events on a plot. HS gets a square ('s')

for ii = 1:nsteps  % loop over steps
    
    %% Mid-Double-Support
    vdup(1:(hs(ii)-1),:) = NaN;  % set anything before the current step to NaN
    if vdup(hs(ii),3)>0   % check if V(hs,3) is greater than zero. If so, hs occurred with rising com velocity. 
        imidds(ii) = hs(ii);
    else
        imidds(ii) = find(vdup(:,3)>0,1,'first'); % otherwise look for mid-stance as an upward zero-crossing in vcom.
    end
    
    %% V-Up
    vdup(1:(imidds(ii)),:) = NaN; % cancel out times in the past
    vdup(round(mean([to(ii) hs(ii+1)])):end,:) = NaN; % cancel out times after roughly mid-step (halfway between toe off and heel strike)
    [bla, blah] = max(atan2(vdup(:,3),abs(vdup(:,2)))); % find "vup", the steepest upward velocity of the COM
    ivup(ii) = blah;
        vdup = Vcom; % reset vdup
    
    %% Mid-Single Support
    vdup(1:(ivup(ii)-1),:) = NaN; % cancel out times in the past
    vdup(hs(ii+1):end,:) = NaN; % cancel out times in the next step
%     if ii == 7, keyboard
%     end
    imidss(ii) = find(vdup(:,3)<0,1,'first'); % Mid-single--support is a downward zero-crossing in Vcom
    
    %% V-Down
%     vdup(1:(imidss(ii)-1),:) = NaN;
%     vdup(hs(ii+1):end,:) = NaN; % don't look too far this time...
    vdup(vdup(:,3)>0,:) = NaN;  % if the COM is rising, this can't be the steepest fall
    vdup(to(ii+1):end,:) = NaN; % don't look too far (end at next toe off)
    [bla blah] = min( atan2(vdup(:,3),abs(vdup(:,2))) );  % find "vdn", the steepest downward COM velocity
    ivdn(ii) = blah ;
%     if ~isempty(blah)
%         ivdn(ii) = blah ;
%     else % special case for a "vdn" that has not happened before the stride ends. Wrap around to the beginning. 
%         vdup = Vcom; % reset vdup so we have the beginning again
%         vdup(vdup(:,3)>0) = NaN;  % if the COM is rising, this can't be the steepest fall
%         vdup (to(1)+1:end,:) = NaN; % since we are only looking at the beginning, get rid of the later part. 
%         blah = find( (diff(atan2(vdup(:,3),abs(vdup(:,2)))) < 0), 1, 'last'); 
%         if ~isempty(blah)
%             ivdn(ii) = blah; % use the wrapped index for vdn, if it was found
%         else
%             ivdn(ii) = hs(end) ; % use the last heel strike as the index of vdn. 
%         end
%     end
        vdup = Vcom; % Reset vdup for next time around;
    
    % Concatenate all these events as: [previousevents, mid-double-support, toe-off, vup, mid-single-support, vdn, next-heel-strike]
    % First time, previous is [hs], so this yields [hs midds to vup midss vdn hs]
    vcomevents = [vcomevents', imidds(ii), to(ii), ivup(ii), imidss(ii), ivdn(ii), hs(ii+1)]';
    vcomeventsyms = {vcomeventsyms{:}, '.', 'd', '^', 'x', 'v', 's' }';
end

vcomeventVs = Vcom(vcomevents,:); % get the actual COM velocity at these events

% Make the output Structure: events, symbols, velocities
out.vcomevents = vcomevents;
out.vcomeventsyms = cell2mat(vcomeventsyms);
out.vcomeventVs = vcomeventVs;

end