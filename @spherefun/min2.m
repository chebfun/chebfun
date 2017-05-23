function varargout = min2(varargin)
%MIN2   Global minimum of a SPHEREFUN.
%   M = MIN2(F) returns the global minimum of F over its domain. 
%   
%   [M, LOC] = MIN2(F) returns the global minimum in M and its location in
%   LOC (in longitude and latitude, respectively).
%
%  This command may be faster if the OPTIMIZATION TOOLBOX is installed.
%
% See also SPHEREFUN/MAX2, SPHEREFUN/MINANDMAX2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = min2@separableApprox(varargin{:});

end