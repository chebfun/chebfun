function f = flipdim( f, dim )
%FLIPDIM   Flip/reverse a SEPARABLEAPPROX in a chosen direction.
%   G = FLIPDIM(F, DIM) returns a SEPARABLEAPPROX G with the same domain as F but
%   reversed in a direction, i.e., G(x,y)=F(x, c+d-y). If DIM = 2 (default) then
%   G(x,y) = F(x, c+d-y).  Otherwise DIM = 1 and G(x,y) = F(a+b-x, y). The
%   domain of F is [a, b, c, d].
% 
% See also FLIPLR, FLIPUD.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( numel( dim ) ~= 1 ) 
    error('CHEBFUN:SEPARABLEAPPROX:flipdim:badDim1', 'DIM should be either 1 or 2.')
end

% Use the code in flipud and fliplr. 
if ( dim == 1 )
    f = flipud( f ); 
elseif ( dim == 2 )
    f = fliplr( f ); 
else
    error('CHEBFUN:SEPARABLEAPPROX:flipdim:badDim2', 'Dimension not recognised.');
end

end
