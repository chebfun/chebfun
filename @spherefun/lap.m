function varargout = lap(varargin)
%LAP   Laplacian of a SPHEREFUN2.
%   L = LAP(F) returns a SPHEREFUN2 representing the Laplacian of F. 
%
%   This is shorthand for LAPLACIAN(F).
%
% See also SPHEREFUN2/LAPLACIAN.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = lap@separableApprox(varargin{:});
end
