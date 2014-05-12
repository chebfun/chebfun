function f = fracCalc(f, m)
%FRACCALC   Differentiation or integral of a CHEBFUN for a fractional order.
%  FRACCALC(F, M) is called by DIFF(F, M) and CUMSUM(F, M) when M is not an
%  integer and computes the fractional integral of order M of F, as defined
%  by the Riemann–Liouville integral.
%
% See also CUMSUM, DIFF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Trivial case:
if ( isempty(f) )
    return
end

% Deal with the integer part of M:
for k = 1:m
    f = cumsum(f);
end

% The fractional part:
fracM = m - floor(m);

% Grab some fields from f:
funs = f.funs;
numFuns = numel(funs);

% Loop over each FUN:
for k = 1:numFuns
    f.funs{k} = bndfun.fracCalc(funs(1:k), fracM);
end

end