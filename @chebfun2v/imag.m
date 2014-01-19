function F = imag( F )
%IMAG imaginary part of a chebfun2v 
% 
% IMAG(F) returns the imaginary part of a chebfun2v.
%
% See also CONJ, REAL.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check:
if ( isempty( F ) )
    return
end

% Take conj part of each component:
F.components = cellfun( @imag, F.components, 'UniformOutput', false );

end