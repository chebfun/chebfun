function varargout = max2(varargin)
%MAX2   Global maximum of a DISKFUN.
%   Y = MAX2(F) returns the global maximum of F over its domain. 
%   
%   [M, LOC] = MAX2(F) returns the global maximum in M and its location in
%   LOC (as a default, LOC is returned in Cartesian coordinates).
%   
%   If the flag 'polar' is included or the diskfun is set to evaluate in
%   polar coords, LOC is returned in polar coordinates.
%
%  This command may be faster if the OPTIMIZATION TOOLBOX is installed.
% 
% See also DISKFUN/MIN2, DISKFUN/MINANDMAX2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% will search over polar coords for the max: 
iscart = diskfun.coordsetting(varargin{:}); 
f = varargin{1}; 
f.coords = 'polar';
[out{1:nargout}] = max2@separableApprox(f);
if numel(out) == 1
    varargout = out;
else
    if iscart %return loc var as Cartesian coords
        [val, loc] = out{:};
        t = loc(1); 
        loc(1) = loc(2).*cos(t); 
        loc(2) = loc(2).*sin(t); 
        varargout = {val, loc};
    else
        varargout = out; 
    end
end

    


end
