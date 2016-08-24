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

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = varargin{1}; 
if ( isempty(f) )
    varargout = [];
    return
end

%give CDR to separable approx so it evaluates in polar. 
[c, d, r] = cdr(f); 
f = c*d*r.'; 
[out{1:nargout}] = max2@separableApprox(f);
if numel(out) == 1
    varargout = out;
else
    %return loc var as Cartesian coords
    [val, loc] = out{:};
     t = loc(1); 
     loc(1) = loc(2).*cos(t); 
     loc(2) = loc(2).*sin(t); 
     varargout = {val, loc};
end

end
    



