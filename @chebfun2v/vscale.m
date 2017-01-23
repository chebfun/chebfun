function out = vscale(F)
%VSCALE   Vertical scale of a CHEBFUN2V.
%   VSCL = VSCALE(F) returns the maximal vertical scale of the components of a
%   CHEBFUN2V object F as determined by evaluating F on a coarse tensor-product
%   grid.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

for jj = 1:F.nComponents
    vscl(jj,1) = vscale(F.components{jj});
end
out = max(vscl);

end
