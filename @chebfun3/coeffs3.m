function varargout = coeffs3(f) 
% COEFFS3   Trivariate Cheybshev expansion coefficients of F. 
%   This is a wrapper for CHEBCOEFFS3.
%   C = COEFFS3(F) returns the tensor of trivariate coefficients.
%
%   [core, C, R, T] = COEFFS3(f) returns the same coefficients kept in
%   the Tucker form.
%
% See also CHEBCOEFFS3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    varargout = {[]};
    return
end

if ( iszero(f) ) 
    [m, n, p] = length(f);
    varargout = {zeros(n, m, p)};
    return
end

if ( nargout <= 1 )
    % Return the tensor of coefficients
    coeffs = chebcoeffs3(f);
    varargout = { coeffs };
elseif ( nargout == 4 )
    [core, cols_coeffs, rows_coeffs, tubes_coeffs] = chebcoeffs3(f);
    varargout = {core, cols_coeffs, rows_coeffs, tubes_coeffs};
else
    error('CHEBFUN:CHEBFUN3:coeffs3:outputs', ...
        'Incorrect number of outputs.');
end

end