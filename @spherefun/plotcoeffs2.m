function varargout = plotcoeffs2( f )
%PLOTCOEFFS2    Display two-dimensional Fourier coefficients graphically.
%   PLOTCOEFFS2(F) plots the two-dimensional coefficients in a stem3 plot
%   with a semilogy scale.
%
%   H = PLOTCOEFFS2(F) returns a handle H to the figure.
%
% See also SPHEREFUN/PLOTCOEFFS, SPHEREFUN/COEFFS2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Compute the bivariate coefficients, truncate them when they fall below
% tolerance for better visualization, use stem3.
X = abs(coeffs2(f)); % Absolute value of coefficients.
[m, n] = size(X);

% Use a stem3 plot changing the axis to log scale. 
wavem = -ceil((m-1)/2):floor((m-1)/2);
waven = -ceil((n-1)/2):floor((n-1)/2);

[xx, yy] = meshgrid(waven, wavem);
h = stem3(xx, yy, X, 'fill', 'markerfacecolor', 'k', ...
    'markeredgecolor', 'k');
set(gca, 'ZScale', 'log', 'view', [40 20])
box off

% Add title and labels
title(gca, '2D Fourier coefficients')
xlabel(gca, 'Wave number (rows)')
ylabel(gca, 'Wave number (columns)')
zlabel(gca, 'Magnitude of coefficient')

% output handle
if ( nargout ~=0 )
    varargout = { h };
end

end