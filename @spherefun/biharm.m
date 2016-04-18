function varargout = biharm(varargin)
%BIHARM   Biharmonic operator of a SPHEREFUN.
%   B = BIHARM(F) returns a SPHEREFUN representing the biharmonic operator 
%   applied to F.
%
%   This is shorthand for BIHARMONIC(F).
%
% See also SPHEREFUN/BIHARMONIC.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = biharm@separableApprox(varargin{:});
end
