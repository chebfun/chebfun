function varargout = biharmonic(varargin)
%BIHARMONIC   Biharmonic operator of a SPHEREFUN.
%   B = BIHARMONIC(F) returns a SPHEREFUN representing the biharmonic 
%   operator applied to F.
%
% See also SPHEREFUN/BIHARM.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = biharmonic@separableApprox(varargin{:});
end
