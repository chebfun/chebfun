function varargout = lap(varargin)
%LAP   Laplacian of a SPHEREFUN.
%   L = LAP(F) returns a SPHEREFUN representing the Laplacian of F. 
%
%   This is shorthand for LAPLACIAN(F).
%
% See also SPHEREFUN/LAPLACIAN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = laplacian(varargin{:});

end
