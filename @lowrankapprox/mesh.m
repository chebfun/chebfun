function varargout = mesh( varargin ) %#ok<STOUT>
%MESH   3-D mesh surface of a LOWRANKAPPROX.
%   MESH is not supported for LOWRANKAPPROX objects, and so returns an error.
%
% See also SURF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:LOWRANKAPPROX:mesh:noSupport', ...
                  'MESH for LOWRANKAPPROX objects is not supported.'); 

end