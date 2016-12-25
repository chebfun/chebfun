function varargout = mesh(varargin)
%MESH   Mesh surface of a DISKFUN.
%   MESH is not supported for DISKFUN objects, and so returns an error.
%
% See also DISKFUN/SURF.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = mesh@separableApprox(varargin{:});

end