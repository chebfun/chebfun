function varargout = laplacian(varargin)
%LAPLACIAN   Laplacian of a CHEBFUN2.
%
%   L = LAPLACIAN(F) returns a CHEBFUN2 representing the Laplacian of F.
%
% See also LAP.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = laplacian@separableApprox(varargin{:});

end
