function varargout = coeffs2( f ) 
% COEFFS2   Bivariate Cheybshev expansion coefficients of f. 
% 
% Same as CHEBCOEFFS2.
%
% See also PLOTCOEFFS2, CHEBCOEFFS2, CHEBCOEFFS.

if ( isempty(f) )
    varargout = { [ ] }; 
    return
end

if ( iszero(f) ) 
    varargout = { 0 } ; 
    return
end

[cols_coeffs, d, rows_coeffs] = chebcoeffs2( f );

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


