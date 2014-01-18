function f = fliplr( f )
%FLIPLR  Flip/reverse a chebfun2 in the x-direction.
%
% G = FLIPLR(F) returns a chebfun2 G with the same domain as F but
% reversed; that is, G(x,y)=F(a+b-x,y), where the domain is [a,b,c,d].
%
% See also FLIPUD.

% Empty check: 
if ( isempty( f ) ) 
    return
end

% Flip the row slices: 
f.rows = flipud( f.rows );

end