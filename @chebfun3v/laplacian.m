function L = laplacian(F)
%LAPLACIAN   Vector Laplacian of a CHEBFUN3V object.
%   LAPLACIAN(F) returns a CHEBFUN3V representing the vector Laplacian of F.
%
% See also CHEBFUN3V/LAP.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(F) )
    L = chebfun3v;
    return
end

diff1 = diff(F, 2, 1);
diff2 = diff(F, 2, 2); 
diff3 = diff(F, 2, 3);

% laplacian = F_xx + F_yy + F_zz
L = chebfun3v;
vscales1 = [vscale(diff1.components{1}) + vscale(diff2.components{1}), ...
    vscale(diff3.components{1})];
L.components{1} = chebfun3(@(x,y,z) feval(diff1.components{1},x,y,z) ...
    + feval(diff2.components{1},x,y,z) + feval(diff3.components{1},x,y,z), ...
    'vscaleBnd', vscales1, 'fiberDim', 3);

vscales2 = [vscale(diff1.components{2}) + vscale(diff2.components{2}), ...
    vscale(diff3.components{2})];
L.components{2} = chebfun3(@(x,y,z) feval(diff1.components{2},x,y,z) + ...
    feval(diff2.components{2},x,y,z) + feval(diff3.components{2},x,y,z),...
    'vscaleBnd', vscales2, 'fiberDim', 3);

vscales3 = [vscale(diff1.components{3}) + vscale(diff2.components{3}), ...
    vscale(diff3.components{3})];
L.components{3} = chebfun3(@(x,y,z) feval(diff1.components{3},x,y,z) + ...
    feval(diff2.components{3},x,y,z) + feval(diff3.components{3},x,y,z), ...
    'vscaleBnd', vscales3, 'fiberDim', 3);

end