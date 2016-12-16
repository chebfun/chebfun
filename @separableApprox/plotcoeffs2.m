function varargout = plotcoeffs2( f )
%PLOTCOEFFS2   Display bivariate coefficients graphically.
%   PLOTCOEFFS2(F) plots the bivariate coefficients in a stem3 plot with a
%   semilogy scale.
%
%   H = PLOTCOEFFS2(F) returns a handle H to the figure.
%
% See also PLOTCOEFFS, CHEBCOEFFS2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Compute the bivariate coefficients, truncate them when they fall below
% tolerance for better visual, use stem3.
X = abs( coeffs2( f ) ); % Absolute value of coefficients. 
 
% Use a stem3 plot changing the axis to log scale. 
[xx, yy] = meshgrid( 1:size(X,2), 1:size(X,1) );
h = stem3(xx, yy, X, 'fill', 'markerfacecolor', 'k', 'markeredgecolor', 'k');
set(gca, 'ZScale', 'log', 'view', [40 20])
box off

% Add title and labels
title(gca,  'Bivariate coefficients')
xlabel(gca, 'Degree of polynomial (rows)')
ylabel(gca, 'Degree of polynomial (columns)')
zlabel(gca, 'Magnitude of coefficient')

% output handle
if ( nargout ~=0 )
    varargout = { h };
end

end
