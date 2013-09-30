% Test file for singfun/imag.m

function pass = test_imag(pref)

% Get preferences.
if ( nargin < 1 )
    pref = singfun.pref;
end

% Pre-allocate pass matrix
pass = zeros(1, 6);

% Check for empty cases:
f = singfun();
pass(1) = isempty(imag(f));

% A real smooth SINGFUN should not return a SINGFUN
f = singfun(@(x) sin(x));
pass(2) = ~isa(imag(f), 'singfun');

% Imaginary smooth.
f = singfun(@(x) 1i*exp(x) );
pass(3) = isequal(1i*imag(f), f.smoothPart);

% Imaginary with exponents.
f = singfun( @(x) 1i./((1+x).*(1-x)));
pass(4) = isequal(1i*imag(f), f);

% Purely real.
f = 1i*f;
pass(5) = iszero(imag(f));

% Complex smooth part.
f = singfun(@(x) (sin(x)+1i*cos(x))./((1+x).*(1-x)));
g = singfun(@(x) cos(x)./((1+x).*(1-x)));
h = imag(f)-g;
pass(6) = iszero(h);
end