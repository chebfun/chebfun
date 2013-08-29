function val = feval(f,x)
%FEVAL Evaluate a SINGFUN at the given points X.
%   For evaluation, the underlying smooth fun is evaluated at X first and
%   then the values computed are scaled by the singular factors.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

%%
% Evaluate the smooth part.
val = feval(f.smoothPart,x);

%%
% Multiply now with the singular factors. 
if ( f.exponents(1) )
    % If there is a non-trivial left singular factor
    val = val.*(1+x).^(f.exponents(1));
end

if ( f.exponents(2) )
    % If there is a non-trivial right singular factor
    val = val.*(1-x).^(f.exponents(2));
end

end
