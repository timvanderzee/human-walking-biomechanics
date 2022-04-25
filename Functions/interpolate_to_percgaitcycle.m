function[xi] = interpolate_to_percgaitcycle(x,segm,npoints)

xi = nan(npoints,length(segm)-1);

for j = 1:(length(segm)-1)
    
    % heelstrike to heelstrike, finite only
    xseg = x(segm(j):segm(j+1));
    xfin = xseg(isfinite(xseg)); xfin = xfin(:);
    
    % make artifical time vector
    tart = linspace(0, 100, length(xfin)); tart = tart(:);
    
    % interpolate
    xi(:,j) = interp1(tart, xfin, linspace(0, 100, npoints));
end

end