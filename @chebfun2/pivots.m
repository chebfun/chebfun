function varargout = pivots(varargin)
%PIVOTS   Pivot values of a CHEBFUN2.
%
%   PIVOTS(F) returns the pivot values taken during in the constructor by the GE
%   algorithm.
%
%   PIVOTS(F, 'normalize'), returns the normalised pivot values. These numbers
%   are scaled so that the columns and rows have unit 2-norm.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = pivots@separableApprox(varargin{:});

end
