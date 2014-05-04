% Test file for singfun/real.m

function pass = test_real(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% Check for empty cases:
f = singfun();
pass(1) = isempty(real(f));

% A purely imaginary smooth SINGFUN should not return a SINGFUN
f = singfun(@(x) 1i*sin(x));
pass(2) = ~isa(real(f), 'singfun');

% Real smooth.
f = singfun(@(x) exp(x) );
pass(1) = isequal(real(f), f.smoothPart);

% Real with exponents.
f = singfun( @(x) 1./((1+x).*(1-x)));
pass(2) = isequal(real(f), f);

% Purely imaginary.
f = 1i*f;
pass(3) = iszero(real(f));

% Complex smooth part.
f = singfun(@(x) (sin(x)+1i*cos(x))./((1+x).*(1-x)));
g = singfun(@(x) sin(x)./((1+x).*(1-x)));
h = real(f)-g;
pass(4) = iszero(h);
end
