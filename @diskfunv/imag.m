function F = imag( F )
%IMAG   Imaginary part of a DISKFUNV.
%   IMAG(F) returns the imaginary part of a DISKFUNV.
%
%   Since all DISKFUNV objects are real-valued, the result of this
%   function is always the original DISKFUNV F.
%
% See also DISKFUNV/CONJ, DISKFUNV/REAL.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty( F ) )
    return
end

% Take imaginary part of each component:
F.components = cellfun( @imag, F.components, 'UniformOutput', false );

end