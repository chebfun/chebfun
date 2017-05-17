function out = isreal(f)
%ISREAL   Test if a CHEBFUN3 object F is real-valued.
%   ISREAL(F) returns logical true if F does not have any nonzero imaginary
%   part and false otherwise.
%  
%   (This is slightly different from the Matlab convention, where isreal(x)
%   is false if x is a complex number whose imaginary part is 0.)
%
% See also CHEBFUN/ISREAL and CHEBFUN2/ISREAL.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    out = true;
    return
end

% Get the low rank representation for f. 
[fCore, fCols, fRows, fTubes] = tucker(f);

% Check individual columns and rows. 
out = isreal(fCore) && isreal(fCols) && isreal(fRows) && isreal(fTubes);

end
