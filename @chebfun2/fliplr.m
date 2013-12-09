function f = fliplr(f)
%FLIPLR  Flip/reverse a chebfun2 in the x-direction.
%
% G = FLIPLR(F) returns a chebfun2 G with the same domain as F but
% reversed; that is, G(x,y)=F(a+b-x,y), where the domain is [a,b,c,d].
%
% See also FLIPUD.

if ( isempty( f ) ) % check for empty chebfun2.
    return
end

f.rows = flipud( f.rows );       % Flip the row slices.

end