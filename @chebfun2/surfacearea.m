function varargout = surfacearea(varargin)
%SURFACEAREA    Surface area of a CHEBFUN2.
%
%   SURFACEAREA(F) computes the surface area of the CHEBFUN2 in the domain of F.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = surfacearea@separableApprox(varargin{:});

end
