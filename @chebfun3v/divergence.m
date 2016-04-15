function f = divergence( F )
%DIVERGENCE   Divergence of a CHEBFUN3V.
%   DIVERGENCE(F) returns the divergence of the CHEBFUN3V i.e.,
%       divergence(F) = F_x + F_y + F_z.
%
% See also DIV.

% Empty check: 
if ( isempty( F ) )
    f = chebfun3;
    return
end

Fc = F.components; 
diff1 = diff( Fc{1}, 1, 1 );
diff2 = diff( Fc{2}, 1, 2 );
diff3 = diff( Fc{3}, 1, 3 );

% [~, f_fiberDim] = min( rank(diff1) );
% f = chebfun3(@(x,y,z) feval(diff1,x, y, z) + feval(diff2,x, y, z) + feval(diff3,x, y, z) , 'fiberDim',f_fiberDim(1));
f = diff1 + diff2 + diff3; % Use compressed plus.

end