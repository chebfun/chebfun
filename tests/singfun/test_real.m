% Test file for singfun/real.m

function pass = test_real(pref)

% Get preferences.
if ( nargin < 1 )
    pref = singfun.pref;
end

% Pre-allocate pass matrix
pass = zeros(1, 4);

% Real smooth.
f = singfun(@(x) exp(x) );
pass(1) = isequal(real(f), f);

% Real with exponents.
f = singfun( @(x) 1./((1+x).*(1-x)));
pass(2) = isequal(real(f), f);

% Purely imaginary.
f = 1i*f;
pass(3) = isequal(real(f), singfun.zeroSingFun());


% Complex smooth part.
f = singfun(@(x) (sin(x)+1i*cos(x))./((1+x).*(1-x)));
g = singfun(@(x) sin(x)./((1+x).*(1-x)));
h = real(f)-g;
pass(4) = iszero(h.smoothPart);
end