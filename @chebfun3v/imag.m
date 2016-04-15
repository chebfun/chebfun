function F = imag( F )
%IMAG   Imaginary part of a CHEBFUN3V.
%   IMAG(F) returns the imaginary part of a CHEBFUN3V.
%
%   See also CONJ, REAL.

% Empty check:
if ( isempty( F ) )
    return
end

% Take imag part of each component:
F.components = cellfun( @imag, F.components, 'UniformOutput', false );

end
