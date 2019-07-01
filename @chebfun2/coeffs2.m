function varargout = coeffs2( f, m, n ) 
% COEFFS2   Bivariate Chebyshev expansion coefficients of f. 
%
%   X = COEFFS2(F) returns the 2D Chebyshev coefficients of the chebfun2. 
% 
%   [C, D, R] = COEFFS2(F) returns a low rank approximation to the matrix
%   of Chebyshev modes.
% 
%   X = COEFFS2(F, M, N) returns coefficients as an MxN matrix.
%
% See also PLOTCOEFFS2, CHEBCOEFFS2, CHEBCOEFFS.

if ( isempty(f) )
    if ( nargout <= 1 )
        varargout = { [ ] }; 
    elseif ( nargout == 3 )
        varargout = {[], [], []};
    else
        % Two or four+ output variables are not allowed.
        error('CHEBFUN:CHEBFUN2:coeffs2:outputs', ...
        'Incorrect number of outputs.');
    end
    return
end

% If f is the zero function: 
if ( iszero(f) )
    if ( nargin == 2 ) 
        n = m; 
    elseif ( nargin == 1 )
        % Fix convention between degrees begin [degX, degY] = length(f) and
        % the coefficient matrix being returned in meshgrid form. 
        [n, m] = length( f );   % This line should not be [m,n]=length(f)!
    end
    if ( nargout <=1 )
        varargout = { zeros(m, n) }; 
    elseif ( nargout == 3 )
        varargout = { zeros(m,1), diag(1), zeros(n,1) };
    else
        % Two or four+ output variables are not allowed.
        error('CHEBFUN:CHEBFUN2:coeffs2:outputs', ...
        'Incorrect number of outputs.');
    end
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
        cols_coeffs = [ cols_coeffs ; zeros(m-mf, rf) ];
    else
        cols_coeffs = cols_coeffs(1:m,:);
    end
    if ( nf <= n )
        rows_coeffs = [ rows_coeffs ; zeros(n-nf, rf) ];
    else
        rows_coeffs = rows_coeffs(1:n, :);
    end
end

if ( nargout <= 1 )
    % Return the matrix of coefficients
    varargout = { cols_coeffs * d * rows_coeffs.' }; 
elseif ( nargout == 3 )
    varargout = { cols_coeffs, d, rows_coeffs };
else
    % Two or four+ output variables are not allowed.
    error('CHEBFUN:CHEBFUN2:coeffs2:outputs', ...
        'Incorrect number of outputs.');
end

end