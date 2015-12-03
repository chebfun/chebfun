function varargout = normalizeRowsAndCols(varargin)
%NORMALIZEROWSANDCOLS   Normalize the rows and columns of a CHEBFUN2.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = normalizeRowsAndCols@separableApprox(varargin{:});

end
