function [hsl, tol, hsr, tor] = invDynGrid_getHS_TO(GRFl, GRFr, thres)

% Takes grfl and grfr from InverseDynamicsGrid data, 
% generates left foot heel strike, left foor toe off, right HS, right TO
% ! GRF vertical threshold is set to be 40. Need to be revised. 
% ! This code may cut out 1 HS or TO from each leg for consistency. 
% by Hansol

% Adapted from Hansol to make threshold an input

GRFl(GRFl(:,3)<thres,:) = 0;
GRFr(GRFr(:,3)<thres,:) = 0;

cutoff = 6;
fs = 1200;
[b_low,a_low] = butter(2,cutoff/(fs/2),'low');

% GRFl = grfl{1};
temp = filtfilt(b_low,a_low,GRFl(:,3));

idx1 = find(temp<=250);
idx2 = find(temp>250);
ref_idx1 = intersect(idx1+1, idx2); % first point greater than 200
ref_idx2 = intersect(idx1-1, idx2); % last point greater than 200
clearvars idx1 idx2

% detect HS and TO %
hsl = nan(length(ref_idx1),1);
for q=1:length(ref_idx1)
  hs = find(GRFl(1:ref_idx1(q),3)<thres, 1, 'last');
  if(q==1&&numel(hs)==0), hsl(q) = nan;
  else
    hsl(q) = hs;
  end
end

if(isnan(hsl(1))), hsl(1) = []; end

tol = nan(length(ref_idx2),1);
for q=1:length(ref_idx2)
  to = find(GRFl(ref_idx2(q):end,3)<thres, 1, 'first');
  if(q==length(ref_idx2)&&numel(to)==0), tol(q) = nan;
  else
    tol(q) = to + ref_idx2(q) - 1;
  end
end
if(isnan(tol(end))), tol(end) = []; end

if(tol(1)<hsl(1)), tol(1)=[]; end
% if(hsl(end)>tol(end)), hsl(end)=[]; end;



% GRFr = grfr{1};
temp = filtfilt(b_low,a_low,GRFr(:,3));

idx1 = find(temp<=250);
idx2 = find(temp>250);
ref_idx1 = intersect(idx1+1, idx2); % first point greater than 200
ref_idx2 = intersect(idx1-1, idx2); % last point greater than 200
clearvars idx1 idx2

% detect HS and TO %
hsr = nan(length(ref_idx1),1);
for q=1:length(ref_idx1)
  hs = find(GRFr(1:ref_idx1(q),3)<thres, 1, 'last');
  if(q==1&&numel(hs)==0), hsr(q) = nan;
  else
    hsr(q) = hs;
  end
end

if(isnan(hsr(1))), hsr(1) = []; end

tor = nan(length(ref_idx2),1);
for q=1:length(ref_idx2)
  to = find(GRFr(ref_idx2(q):end,3)<thres, 1, 'first');
  if(q==length(ref_idx2)&&numel(to)==0), tor(q) = nan;
  else
    tor(q) = to + ref_idx2(q) - 1;
  end
end
if(isnan(tor(end))), tor(end) = []; end

% if(tor(1)<hsr(1)), tor(1)=[]; end;
if(hsr(end)>tor(end)), hsr(end)=[]; end