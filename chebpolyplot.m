function varargout = chebpolyplot(varargin)
%CHEBPOLYPLOT   Display Chebyshev coefficients graphically.
%   CHEBPOLYPLOT(F) is deprecated. Please use PLOTCOEFFS().
%
% See also PLOTCOEFFS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

warning('CHEBFUN:chebpolyplot:deprecated', ...
    'CHEBPOLYPLOT is deprecated. Please use PLOTCOEFFS instead.');
warning('off', 'CHEBFUN:chebpolyplot:deprecated');

[varargout{1:nargout}] = plotcoeffs(varargin{:});

end
