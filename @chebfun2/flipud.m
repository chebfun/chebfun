function f = flipud( f ) 
%FLIPUD  Flip/reverse a chebfun2 in the y-direction.
%
% G = FLIPUD(F) returns a chebfun2 G with the same domain as F but
% reversed; that is, G(x,y) = F(x,c+d-y), where the domain is [a b] x [c d].
%
% See also FLIPLR, FLIPDIM. 

if ( isempty( f ) )   % check for empty chebfun2
    return
end

f.cols = flipud( f.cols );  % flip the column slices.

end