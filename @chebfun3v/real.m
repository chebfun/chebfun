function F = real( F )
%REAL  Real part of a CHEBFUN3V.
%   REAL(F) returns the CHEBFUN3V representing the real part.
%   See also CONJ, IMAG.


% Empty check: 
if ( isempty( F ) )
    return
end

% Take real part of each component:
F.components = cellfun( @real, F.components, 'UniformOutput', false );

end