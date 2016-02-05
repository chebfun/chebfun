function varargout = biharmonic(varargin)
%BIHARMONIC   Biharmonic operator of a CHEBFUN2.
%   B = BIHARMONIC(F) returns a CHEBFUN2 representing the biharmonic operator 
%   applied to F.
%
% See also CHEBFUN2/BIHARMONIC.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call SEPARABLEAPPROX/BIHARMONIC:
[varargout{1:nargout}] = biharmonic@separableApprox(varargin{:});

end
