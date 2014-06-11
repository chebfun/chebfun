function varargout = chebcoeffsplot2( f )
%CHEBCOEFFSPLOT2   Display bivariate Chebyshev coefficients graphically.
%   CHEBCOEFFSPLOT2(F) plots the bivariate Chebyshev coefficients in a stem3 plot
%   with a semilogy scale.
%
%   H = CHEBCOEFFSPLOT2(F) returns a handle H to the figure.
%
% See also CHEBCOEFFSPLOT, CHEBCOEFFS2.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Compute the bivariate coefficients, truncate them when then fall below
% tolerance for better visual, use stem3.

X = abs( chebcoeffs2( f ) );    % Absolute value of coefficients. 
X = rot90(X, 2);                % Rotate (MATLAB's convention)
 
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
