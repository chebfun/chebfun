function display(L)
[m, n] = size(L);
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
fprintf('\n  %ix%i block chebmatrix of types:\n\n', m, n)
disp( blockClasses(L) )
end
