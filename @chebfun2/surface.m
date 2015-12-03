function varargout = surface(varargin)
%SURFACE  Plot surface of a CHEBFUN2.
%   SURFACE(X, Y, Z, C) adds the surface in X,Y,Z,C to the current axes.
%
%   SURFACE(X, Y, Z) uses C = Z, so color is proportional to surface height.
%
%   See SURF for a complete description of the various forms that X,Y,Z,C can
%   take.
%
% See also SURF.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = surface@separableApprox(varargin{:});

end
