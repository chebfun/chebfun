function varargout = coeffs2( f, m, n ) 
% COEFFS2   Bivariate Chebyshev expansion coefficients of f. 
%
%   X = COEFFS2(F) returns the 2D Chebyshev modes of the chebfun2. 
% 
%   [C, D, R] = COEFFS2(F) returns a low rank approximation to the matrix
%   of Chebyshev modes.
% 
%   X = COEFFS2(F, M, N) returns bivariate coefficients as an MxN matrix of
%   Chebyshev modes. 
%
% See also PLOTCOEFFS2, CHEBCOEFFS2, CHEBCOEFFS.

if ( isempty(f) )
    varargout = { [ ] }; 
    return
end

if ( iszero(f) ) 
    [m, n] = length( f );
    varargout = { zeros(n, m) } ; 
    return
end

[cols_coeffs, d, rows_coeffs] = chebcoeffs2( f );

if ( nargin == 2 ) 
    n = m; 
end

if ( nargin > 1 ) 
    [mf, ignored] = size(cols_coeffs); 
    [nf, rf] = size(rows_coeffs); 
    if ( mf <= m ) 
        cols_coeffs = [ cols_coeffs ; zeros(m-mf,rf) ]; 
    else
        cols_coeffs = cols_coeffs(1:m,:);
    end
    if ( mf <= m ) 
        rows_coeffs = [ rows_coeffs ; zeros(n-nf,rf) ]; 
    else
        rows_coeffs = rows_coeffs(1:n,:);
    end
end


if ( nargout <= 1 )
    % Return the matrix of coefficients
    varargout = { cols_coeffs * d * rows_coeffs.' }; 
elseif ( nargout <= 3 )
    varargout = { cols_coeffs, d, rows_coeffs };
else
    % Two output variables are not allowed.
    error('CHEBFUN:CHEBFUN2:coeffs2:outputs', ...
        'Incorrect number of outputs.');
end

end


