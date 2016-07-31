function F = imag( F )
%IMAG   Imaginary part of a DISKFUNV.
%   IMAG(F) returns the imaginary part of a DISKFUNV.
%
%   Since only real DISKFUNV objects are supported, the result of this
%   function is always just the original function F.
%
%   See also DISKFUNV/CONJ, DISKFUNV/REAL.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty( F ) )
    return
end

% Take imag part of each component:
F.components = cellfun( @imag, F.components, 'UniformOutput', false );

end
