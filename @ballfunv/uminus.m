function g = uminus(f)
%-  BALLFUNV unary minus.
%   -F negates the BALLFUNV F.
%
%   G = UMINUS(F) is called for the syntax '-F'.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = f.comp;
g = ballfunv(-F{1},-F{2},-F{3});
end
