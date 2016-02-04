function varargout = biharm(varargin)
%BIHARM   Biharmonic operator of a CHEBFUN2.
%   B = BIHARM(F) returns a CHEBFUN2 representing the biharmonic operator 
%   applied to F.
%
%   This is shorthand for BIHARMONIC(F).
%
% See also CHEBFUN2/BIHARMONIC.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call SEPARABLEAPPROX/BIHARM:
[varargout{1:nargout}] = biharm@separableApprox(varargin{:});

end
