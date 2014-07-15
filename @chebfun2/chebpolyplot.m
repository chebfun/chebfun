function varargout = chebpolyplot( varargin )
%CHEBPOLYPLOT   Display Chebyshev coefficients of slices graphically.
%   CHEBPOLYPLOT(F) is deprecated. Please use PLOTCOEFFS().
%
% See also PLOTCOEFFS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

warning('CHEBFUN2:chebpolyplot:deprecated', ...
    'CHEBPOLYPLOT is deprecated. Please use PLOTCOEFFS instead.');
warning('off', 'CHEBFUN2:chebpolyplot:deprecated');

[varargout{1:nargout}] = plotcoeffs(varargin{:});


end
