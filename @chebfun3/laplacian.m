function L = laplacian(f)
%LAPLACIAN   Laplacian of a CHEBFUN3.
%   L = LAPLACIAN(F) returns a CHEBFUN3 representing the Laplacian of F.
%
% See also CHEBFUN3/LAP.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

diff1 = diff(f, 2, 1); 
diff2 = diff(f, 2, 2); 
diff3 = diff(f, 2, 3);

vscales = [vscale(diff1) + vscale(diff2), vscale(diff3)];
% Send two vscales to constructor

L = chebfun3(@(x,y,z) feval(diff1,x, y, z) + feval(diff2,x, y, z) + ...
    feval(diff3,x, y, z) , f.domain, 'vscaleBnd', vscales);
% Instead of calling the compressed_plus, this runs the constructor for addition. 

% L = diff1 + diff2 + diff3; % Using compression_plus
end
