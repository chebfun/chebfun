% Test file for singfun/conj.m

function pass = test_conj(pref)

% Get preferences.
if ( nargin < 1 )
    pref = singfun.pref;
end

% Pre-allocate pass matrix
pass = zeros(1, 4);

%%
% Real smooth.
f = singfun(@(x) exp(x) );
pass(1) = isequal(conj(f), f);

% Real with exponents.
f = singfun( @(x) 1./((1+x).*(1-x)));
pass(2) = isequal(conj(f), f);

% Purely imaginary.
f = 1i*f;
pass(3) = isequal(conj(f), -f);

% Complex smooth part.
f = singfun(@(x) (sin(x)+1i*cos(x))./((1+x).*(1-x)));
g = singfun(@(x) (sin(x)-1i*cos(x))./((1+x).*(1-x)));
h = conj(f)-g;
pass(4) = iszero(h.smoothPart);
end