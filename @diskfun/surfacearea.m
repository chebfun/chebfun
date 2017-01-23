function varargout = surfacearea(varargin)
%SURFACEAREA   Surface area of a DISKFUN.
%   SURFACEAREA(F) computes the surface area of the DISKFUN in the domain 
%   of F.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = surfacearea@separableApprox(varargin{:});

end