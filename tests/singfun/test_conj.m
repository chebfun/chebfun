% Test file for singfun/conj.m

function pass = test_conj(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

%%
% Empty
f = singfun();
pass(1) = isempty(conj(f));

% A smooth SINGFUN conjugated should not return a SINGFUN
f = singfun(@(x) sin(x));
pass(2) = ~isa(conj(f), 'singfun');

% Real smooth.
f = singfun(@(x) exp(x) );
pass(3) = isequal(conj(f), f.smoothPart);

% Real with exponents.
f = singfun( @(x) 1./((1+x).*(1-x)));
pass(4) = isequal(conj(f), f);

% Purely imaginary.
f = 1i*f;
pass(5) = isequal(conj(f), -f);

% Complex smooth part.
f = singfun(@(x) (sin(x)+1i*cos(x))./((1+x).*(1-x)));
g = singfun(@(x) (sin(x)-1i*cos(x))./((1+x).*(1-x)));
h = conj(f)-g;
pass(6) = iszero(h);
end
