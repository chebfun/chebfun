function varargout = mesh( varargin ) %#ok<STOUT>
%MESH   3-D mesh surface of a SEPARABLEAPPROX.
%   MESH is not supported for SEPARABLEAPPROX objects, and so returns an error.
%
% See also SURF.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:SEPARABLEAPPROX:mesh:noSupport', ...
                  'MESH for SEPARABLEAPPROX objects is not supported.'); 

end