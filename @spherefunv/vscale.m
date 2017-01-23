function out = vscale(F)
%VSCALE   Vertical scale of a SPHEREFUNV.
%   VSCL = VSCALE(F) returns the maximal vertical scale of the components of a
%   SPHEREFUNV object F as determined by evaluating F on a coarse tensor-product
%   grid.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = max([vscale(F.components{1}), vscale(F.components{2}), ...
    vscale(F.components{3})]);

end
