function f = minus(f, g)
%-   Subtraction of two SPHEREFUNS objects.
%   F - G subtracts G from F, where F and G are SPHEREFUN objects or scalars.
%
% See also PLUS, UMINUS.

% f - g = f + (-g)
f = plus(f, uminus(g)); 

end