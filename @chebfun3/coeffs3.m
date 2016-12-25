function varargout = coeffs3(f, m, n, p) 
% COEFFS3   Trivariate Cheybshev expansion coefficients of F. 
%   This is a wrapper for CHEBCOEFFS3.
%   C = COEFFS3(F) returns the tensor of trivariate coefficients.
%
%   [core, C, R, T] = COEFFS3(f) returns the same coefficients kept in
%   the Tucker form.
%
%   X = COEFFS3(f, M, N, P) returns coefficients as an M x N x P tensor.
%
% See also CHEBCOEFFS3.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    varargout = {[]};
    return
end

if ( iszero(f) )
    if ( nargin == 1 )
        [m, n, p] = length(f);
    elseif ( nargin == 2 ) 
        n = m; 
        p = m;
    end
    varargout = {zeros(m, n, p)};
    return
end

[core, cols_coeffs, rows_coeffs, tubes_coeffs] = chebcoeffs3(f);

if ( nargin == 2 ) 
    n = m;
    p = m;
end

if ( nargin > 1 ) 
    [mf, m2f] = size(cols_coeffs); 
    [nf, n2f] = size(rows_coeffs); 
    [pf, p2f] = size(tubes_coeffs); 
    if ( mf < m ) 
        cols_coeffs = [cols_coeffs ; zeros(m-mf, m2f)];
    else
        cols_coeffs = cols_coeffs(1:m, :);
    end
    if ( nf < n ) 
        rows_coeffs = [rows_coeffs ; zeros(n-nf, n2f)]; 
    else
        rows_coeffs = rows_coeffs(1:n, :);
    end
    if ( pf < p ) 
        tubes_coeffs = [tubes_coeffs ; zeros(p-pf, p2f)]; 
    else
        tubes_coeffs = tubes_coeffs(1:p, :);
    end    
end

if ( nargout <= 1 )
    % Return the tensor of coefficients.
    varargout = { chebfun3.txm(chebfun3.txm(chebfun3.txm(core, ...
        cols_coeffs, 1), rows_coeffs, 2), tubes_coeffs, 3) }; 
elseif ( nargout <= 4 )
    varargout = {core, cols_coeffs, rows_coeffs, tubes_coeffs};
else
    error('CHEBFUN:CHEBFUN3:coeffs3:outputs', ...
        'Incorrect number of outputs.');
end

end