function varargout = normalizePivots(varargin)
%NORMALIZEPIVOTS   Scale rows and cols of a CHEBFUN2 so that all pivots are 1.
%
% Additionally, the norm of the kth row and column will be the same.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = normalizePivots@separableApprox(varargin{:});

end
