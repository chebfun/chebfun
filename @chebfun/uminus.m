function f = uminus(f)
%-   CHEBFUN unary minus.
%   -F negates the CHEBFUN F.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Negate the impuleses:
f.impulses(1,:) = -f.impulses(1,:);

% Negate each of the funs:
for k = 1:numel(f.funs)
    f.funs{k} = -f.funs{k};
end

end
