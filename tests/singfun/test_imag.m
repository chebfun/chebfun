% Test file for singfun/imag.m

function pass = test_imag(pref)

% Get preferences.
if ( nargin < 1 )
    pref = singfun.pref;
end

% Pre-allocate pass matrix
pass = zeros(1, 4);

% Imaginary smooth.
f = singfun(@(x) 1i*exp(x) );
pass(1) = isequal(1i*imag(f), f.smoothPart);

% Imaginary with exponents.
f = singfun( @(x) 1i./((1+x).*(1-x)));
pass(2) = isequal(1i*imag(f), f);

% Purely real.
f = 1i*f;
pass(3) = isequal(imag(f), singfun.zeroSingFun());

% Complex smooth part.
f = singfun(@(x) (sin(x)+1i*cos(x))./((1+x).*(1-x)));
g = singfun(@(x) cos(x)./((1+x).*(1-x)));
h = imag(f)-g;
pass(4) = iszero(h.smoothPart);
end