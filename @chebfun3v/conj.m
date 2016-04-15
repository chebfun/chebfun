function F = conj( F )
%CONJ Complex conjugate of a CHEBFUN3V.
%   CONJ(F) returns the complex conjugate of F. For a complex F, CONJ(F) =
%   REAL(F) - i*IMAG(F).
%
% See also REAL, IMAG. 

% Empty check:
if ( isempty( F ) )
    return
end

% Take conj part of each component:
F.components = cellfun( @conj, F.components, 'UniformOutput', false );

end
