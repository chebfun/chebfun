function varargout = normalizeRowsAndCols(varargin)
%NORMALIZEROWSANDCOLS   Normalize the rows and columns of a DISKFUN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = normalizeRowsAndCols@separableApprox(varargin{:});

end