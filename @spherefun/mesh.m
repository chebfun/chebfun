function varargout = mesh(varargin)
%MESH   3-D mesh surface of a SPHEREFUN.
%   MESH is not supported for SPHEREFUN objects, and so returns an error.
%
% See also SPHEREFUN/SURF.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = mesh@separableApprox(varargin{:});

end
