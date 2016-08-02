function varargout = minandmax2(varargin)
%MINANDMAX2   Find global minimum and maximum of a DISKFUN.
%   M = minandmax2(F) returns the minimum and maximum values of a DISKFUN. 
%   M is a vector of length 2 such that 
%   M(1) = min(f(lambda,theta)) and M(2) = max(f(lambda,theta)).
%
%   [M, LOC] = minandmax2(F) also returns the positions of the minimum and 
%   maximum. For example,
%
%       F(LOC(1,1),LOC(1,2)) = M(1)  and  F(LOC(2,1),LOC(2,2)) = M(2)
%
% See also DISKFUN/MAX2, DISKFUN/MIN2, DISKFUN/NORM.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% we want to search over polar coords for the max and min: 
iscart = diskfun.coordsetting(varargin{:}); 
f = varargin{1}; 
f.coords = 'polar';
[out{1:nargout}] = minandmax2@separableApprox(f);
if numel(out) == 1
    varargout = out;
else
    if iscart %return loc var as Cartesian coords
        [val, loc] = out{:};
        t = loc(1, :); 
        loc(1, :) = loc(2, :).*cos(t); 
        loc(2, :) = loc(2, :).*sin(t); 
        varargout = {val, loc};
    else
        varargout = out; 
    end
end

end
