function g = uminus(f)
% UMINUS BALLFUNV unary minus
%   UMINUS(f) is negation of the BALLFUNV f

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = f.comp;
g = ballfunv(-F{1},-F{2},-F{3});
end
