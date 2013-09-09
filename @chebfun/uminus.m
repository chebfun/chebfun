function f = uminus(f)
%-   CHEBFUN unary minus.
%   -F negates the CHEBFUN F.
%
%   G = UMINUS(A) is called for the syntax '-A'.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Handle the empty case.
if ( isempty(f) )
    return
end

% Negate the impulses:
f.impulses = -f.impulses;

% Negate each of the FUNs:
for k = 1:numel(f.funs)
    f.funs{k} = uminus(f.funs{k});
end

end
