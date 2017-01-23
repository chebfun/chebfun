function varargout = chebcoeffs3(f)
%CHEBCOEFFS3   Trivariate Chebyshev coefficients
%   C = CHEBCOEFFS3(F) returns the tensor of trivariate coefficients such 
%   that
%   C = sum_{i=0}^{m-1}  sum_{j=0}^{n-1} sum_{k=0}^{p-1} 
%                                      C(i+1,j+1,k+1) T_i(x) T_j(y) T_k(z).
%
%   [CORE, C, R, T] = CHEBCOEFFS3(F) returns the same coefficients kept in
%   the Tucker form.
%
% See also CHEBCOEFFS, CHEBCOEFFS2 and COEFFS3

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    varargout = {[]}; 
    return
end

if ( iszero(f) ) 
    varargout = {0}; 
    return
end

% Get low rank representation of f:
[core, cols, rows, tubes] = tucker(f);

% Get the coeffs of the columns, rows and tubes:
cols_coeffs = chebcoeffs(cols);
rows_coeffs = chebcoeffs(rows);
tubes_coeffs = chebcoeffs(tubes);

if ( nargout <= 1 )
    % Return the tensor of coefficients:
    varargout = {chebfun3.txm(chebfun3.txm(chebfun3.txm(core, ...
        cols_coeffs, 1), rows_coeffs, 2), tubes_coeffs, 3)};
elseif ( nargout == 4 )
    varargout = {core, cols_coeffs, rows_coeffs, tubes_coeffs};
else
    % Two or three outputs are not allowed.
    error('CHEBFUN:CHEBFUN3:chebcoeffs3:outputs', ...
        'Incorrect number of outputs.');
end

end