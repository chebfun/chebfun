function f = uminus(f)
%-   Unary minus.
%   -F negates the CHEBFUN F.
%
%    G = uminus(A) is called for the syntax '-A'.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Negate the impulses:
f.impulses(1,:) = -f.impulses(1,:);

% Negate each of the FUNs:
for k = 1:numel(f.funs)
    f.funs{k} = uminus(f.funs{k});
end

end