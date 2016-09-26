function varargout = chebpolyplot2(f)
%CHEBPOLYPLOT2   Display bivariate Chebyshev coefficients graphically.
%   CHEBPOLYPLOT2(F) is deprecated. Please use PLOTCOEFFS2.
%
% See also PLOTCOEFFS2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

warning('CHEBFUN:CHEBFUN2:chebpolyplot2:deprecated', ...
        'CHEBPOLYPLOT2 is deprecated. Please use PLOTCOEFFS2 instead.');
warning('off', 'CHEBFUN:CHEBFUN2:chebpolyplot2:deprecated');

[varargout{1:nargout}] = plotcoeffs2(varargin{:});

end
