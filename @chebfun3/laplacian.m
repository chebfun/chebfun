function L = laplacian(f)
%LAPLACIAN   Laplacian of a CHEBFUN3.
%
%   L = LAPLACIAN(F) returns a CHEBFUN3 representing the Laplacian of F.
%
%   See also CHEBFUN3/LAP.

diff1 = diff(f, 2, 1); 
diff2 = diff(f, 2, 2); 
diff3 = diff(f, 2, 3);

vscales = [vscale(diff1) + vscale(diff2), vscale(diff3)];
% Send two vscales to constructor

L = chebfun3(@(x,y,z) feval(diff1,x, y, z) + feval(diff2,x, y, z) + ...
    feval(diff3,x, y, z) , f.domain, 'vscaleBnd', vscales, 'fiberDim', 3);
% Instead of calling the compressed_plus, this runs the constructor for addition. 

% L = diff1 + diff2 + diff3; % Using compression_plus
end
