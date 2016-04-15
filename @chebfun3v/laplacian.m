function L = laplacian( F )
%LAPLACIAN Vector Laplacian of a CHEBFUN3V.
%   LAPLACIAN(F) returns a CHEBFUN3V representing the vector Laplacian of F.
%
% See also CHEBFUN3V/LAP.

% Empty check: 
if ( isempty( F ) )
    L = chebfun3v;
    return
end

diff1 = diff(F,2,1); diff2 = diff(F,2,2); diff3 = diff(F,2,3);
%G = diff(Fc{1}, 2, 1) + diff(Fc{2}, 2, 2) + diff(Fc{3}, 2,3); 

% laplacian = F_xx + F_yy + F_zz
L = chebfun3v;
[~, f_fiberDim1] = min( rank(diff1.components{1}) );
[~, f_fiberDim2] = min( rank(diff1.components{2}) );
[~, f_fiberDim3] = min( rank(diff1.components{3}) );
L.components{1} = chebfun3(@(x,y,z) feval(diff1.components{1},x,y,z) + feval(diff2.components{1},x,y,z) + feval(diff3.components{1},x,y,z), 'fiberDim',f_fiberDim1(1));
L.components{2} = chebfun3(@(x,y,z) feval(diff1.components{2},x,y,z) + feval(diff2.components{2},x,y,z) + feval(diff3.components{2},x,y,z), 'fiberDim',f_fiberDim2(1));
L.components{3} = chebfun3(@(x,y,z) feval(diff1.components{3},x,y,z) + feval(diff2.components{3},x,y,z) + feval(diff3.components{3},x,y,z), 'fiberDim',f_fiberDim3(1));

end