function varargout = surfacearea(varargin)
%SURFACEAREA    Surface area of a SPHEREFUN.
%
%   SURFACEAREA(F) computes the surface area of the SPHEREFUN in the domain of F.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = surfacearea@separableApprox(varargin{:});

end
