function varargout = max2(varargin)
%MAX2   Global maximum of a SPHEREFUN.
%   Y = MAX2(F) returns the global maximum of F over its domain. 
%   
%   [M, LOC] = MAX2(F) returns the global maximum in M and its location in
%   LOC (in longitude and latitude, respectively).
%
%  This command may be faster if the OPTIMIZATION TOOLBOX is installed.
% 
% See also SPHEREFUN/MIN2, SPHEREFUN/MINANDMAX2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = max2@separableApprox(varargin{:});

end
