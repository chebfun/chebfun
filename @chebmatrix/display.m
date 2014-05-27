function display(L)
%DISPLAY  Print summary of CHEBMATRIX contents.
%
% See also CHEBMATRIX.SPY.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[m, n] = size(L);
fprintf('\n  %i x %i chebmatrix of block types:\n\n', m, n)
disp( blockClasses(L) )

end
