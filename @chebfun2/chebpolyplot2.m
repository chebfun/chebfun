function varargout = chebpolyplot2(f)
%CHEBPOLYPLOT2   Display bivariate Chebyshev coefficients graphically.
%   CHEBPOLYPLOT2(F) is deprecated. Please use PLOTCOEFFS().
%
% See also PLOTCOEFFS2.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

warning('CHEBFUN:CHEBFUN2:chebpolyplot2:deprecated', ...
    'CHEBPOLYPLOT2 is deprecated. Please use PLOTCOEFFS2 instead.');
warning('off', 'CHEBFUN:CHEBFUN2:chebpolyplot2:deprecated');

% TODO: Use this once plotcoeffs2 is working:
% [varargout{1:nargout}] = plotcoeffs2(varargin{:});

X = abs( chebcoeffs2( f ) );  % Absolute value of coefficients. 
X = rot90(X, 2);              % Rotate (MATLAB's convention)
 
% Use a stem3 plot changing the axis to log scale. 
[xx, yy] = meshgrid( 1:size(X,1), 1:size(X,2) );
h = stem3(xx, yy, X, 'fill', 'markerfacecolor', 'k', 'markeredgecolor', 'k');
set(gca, 'ZScale', 'log', 'view', [40 20])
box off

% output handle
if ( nargout ~=0 )
    varargout = { h };
end

end
