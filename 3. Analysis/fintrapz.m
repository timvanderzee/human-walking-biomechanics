function [I] = fintrapz(x,y,method,range,y2)

if nargin > 2
    if nargin < 5
        y2 = y;
    end
    
    if strcmp(method,'+') == 1, y(y2<0) = 0;
    elseif strcmp(method,'-') == 1, y(y2>0) = 0;
    end
end

if nargin > 3
    r = range;
else
    r = [1 length(y)];
end

Y = y(r(1):r(2));
X = x(r(1):r(2));

i = isfinite(Y);
I = trapz(X(i),Y(i));

end
