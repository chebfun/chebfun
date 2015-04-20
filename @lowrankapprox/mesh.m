function varargout = mesh( varargin ) %#ok<STOUT>
%MESH   3-D mesh surface of a CHEBFUN2.
%   MESH is not supported for CHEBFUN2 objects, and so returns an error.
%
% See also SURF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:CHEBFUN2:mesh:noSupport', 'MESH is not supported by Chebfun2.'); 

end