function varargout = chebpoly2( f )
%CHEBPOLY2  bivariate Chebyshev coefficients
%   X = CHEBPOLY2(F) returns the matrix of bivariate coefficients such that F =
%   sum_i ( sum_j Y(i,j) T_i(y) T_j(x) ), where Y = rot90(X, 2). It is MATLAB
%   convention to flip the coefficients in this silly way.
%
%   [A, D, B] = CHEBPOLY2( f ) returns the same coefficients keeping them in low
%   rank form, i.e., X = A * D * B'.
%
% See also CHEBPOLYPLOT2, CHEBPOLYPLOT.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty( f ) )
    varargout = { [ ] }; 
    return
end

% Get the low rank representation for f:
[cols, d, rows] = cdr(f);

% Get the coeffs of the rows and the columns:
cols_coeffs = chebpoly( cols ).';
rows_coeffs = chebpoly( rows ).';

if ( nargout <= 1 )
    % Return the matrix of coefficients
    varargout = { cols_coeffs * d * rows_coeffs.' }; 
elseif ( nargout <= 3 )
    varargout = {cols_coeffs, d, rows_coeffs};
else
    % Two output variables are not allowed.
    error('CHEBFUN:CHEBFUN2:chebpoly2:outputs', 'Incorrect number of outputs.');
end

end
