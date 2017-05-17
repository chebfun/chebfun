function out = isreal(f)
%ISREAL   Test if a CHEBFUN3V object F is real-valued.
%   ISREAL(F) returns logical true if F does not have any nonzero imaginary
%   part and false otherwise.
%
%   (This is slightly different from the Matlab convention, where isreal(x)
%   is false if x is a complex number whose imaginary part is 0.)
%
% See also CHEBFUN/ISREAL, CHEBFUN2/ISREAL, and CHEBFUN3/ISREAL.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    out = true;
    return
end

for ii = 1:f.nComponents
    out(ii) = isreal(f.components{ii});
end
out = all(out);

end