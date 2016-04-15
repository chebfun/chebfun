function f = minus(f, g)
%-   Subtraction of two CHEBFUN3 objects.
%   F - G subtracts G from F, where F and G are CHEBFUN3 objects or scalars.
%
%   See also PLUS, UMINUS.

% f - g = f + (-g)
f = plus(f, uminus(g));

end