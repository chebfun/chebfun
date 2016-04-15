function f = minus( f, g )
% - MINUS. Minus of two CHEBFUN3V.  
%   F - G subtracts the CHEBFUN3V F from G componentwise.
%
%   minus(F, G) is called for the syntax f - g.

f = plus( f, uminus( g ) );

end
