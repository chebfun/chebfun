function varargout = max2(varargin)
%MAX2   Global maximum of a DISKFUN.
%   Y = MAX2(F) returns the global maximum of F over its domain. 
%   
%   [M, LOC] = MAX2(F) returns the global maximum in M and its location in
%   LOC.
%
%  This command may be faster if the OPTIMIZATION TOOLBOX is installed.
% 
% See also DISKFUN/MIN2, DISKFUN/MINANDMAX2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = varargin{1}; 
if ( isempty(f) )
    varargout = [];
    return
end

% Give CDR to separableApprox so it evaluates in polar.

f = cart2pol(f, 'cdr');

[out{1:nargout}] = max2@separableApprox(f);
if numel(out) == 1
    varargout = out;
else
    % Return loc in Cartesian coords.
    [val, loc] = out{:};
     t = loc(1); 
     loc(1) = loc(2).*cos(t); 
     loc(2) = loc(2).*sin(t); 
     varargout = {val, loc};
end

end