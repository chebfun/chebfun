function L = laplacian(F)
%LAPLACIAN   Vector Laplacian of a CHEBFUN3V object.
%   LAPLACIAN(F) returns a CHEBFUN3V representing the vector Laplacian of 
%   F. The output L is the vector field of the scalar Laplacian applied to
%   the individual components.
%
% See also CHEBFUN3V/LAP.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(F) )
    L = chebfun3v;
    return
end

L = F;  % Initialise the output to get right domain.

% Call scalar Laplacian for each compoment:
for jj = 1:L.nComponents
    L.components{jj} = laplacian(L.components{jj});
end

end