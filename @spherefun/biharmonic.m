function varargout = biharmonic(varargin)
%BIHARMONIC   Biharmonic operator of a SPHEREFUN2.
%   B = BIHARMONIC(F) returns a SPHEREFUN2 representing the biharmonic 
%   operator applied to F.
%
% See also SPHEREFUN2/BIHARM.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = biharmonic@separableApprox(varargin{:});
end
