function varargout = laplacian(varargin)
%LAPLACIAN   Laplacian of a CHEBFUN2.
%   L = LAPLACIAN(F) returns a CHEBFUN2 representing the Laplacian of F.
%
% See also CHEBFUN2/LAP.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call SEPARABLEAPPROX/LAPLACIAN:
[varargout{1:nargout}] = laplacian@separableApprox(varargin{:});

end
