function varargout = mesh(varargin)
%MESH   3-D mesh surface of a CHEBFUN2.
%   MESH is not supported for CHEBFUN2 objects, and so returns an error.
%
% See also SURF.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = mesh@separableApprox(varargin{:});

end
