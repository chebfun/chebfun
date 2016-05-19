function varargout = coeffs3(f) 
% COEFFS3   Trivariate Cheybshev expansion coefficients of f. 
%   This is a wrapper for chebcoeffs3.

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