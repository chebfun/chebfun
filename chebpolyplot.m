function varargout = chebpolyplot(varargin)
%CHEBPOLYPLOT   Display Chebyshev coefficients graphically.
%   CHEBPOLYPLOT(F) is deprecated. Please use COEFFSPLOT().
%
% See also COEFFSPLOT.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

warning('CHEBFUN:chebpolyplot:deprecated', ...
    'CHEBPOLYPLOT if deprecated. Please use COEFFSPLOT instead.');
warning('off', 'CHEBFUN:chebpolyplot:deprecated');

[varargout{1:nargout}] = coeffsplot(varargin{:});

end