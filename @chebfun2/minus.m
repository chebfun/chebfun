function f = minus(f, g)
%-   Subtraction of two CHEBFUN2 objects.
%   F - G subtracts G from F, where F and G are CHEBFUN2 objects or scalars.
%
% See also PLUS, UMINUS.

% f - g = f + (-g)
f = plus( f, uminus(g) ); 

end
