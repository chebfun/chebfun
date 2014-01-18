function f = flipud( f ) 
%FLIPUD  Flip/reverse a chebfun2 in the y-direction.
%
% G = FLIPUD(F) returns a chebfun2 G with the same domain as F but
% reversed; that is, G(x,y) = F(x,c+d-y), where the domain is [a b] x [c d].
%
% See also FLIPLR, FLIPDIM. 

% Empty check:
if ( isempty( f ) )  
    return
end

% Flip the column slices: 
f.cols = flipud( f.cols ); 

end