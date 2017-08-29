function varargout = pivots(varargin)
%PIVOTS   Pivot values of a SPHEREFUN.
% 
%   PIVOTS(F) returns the pivot values used during the construction by the GE
%   algorithm.
%
%   PIVOTS(F, 'normalize'), returns the normalised pivot values. These number
%   are scaled so that the columns and rows have unit 2-norm.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = pivots@separableApprox(varargin{:});

end
