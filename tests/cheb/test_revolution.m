function pass = test_revolution(~)

%% Check results for a cone
f = chebfun(@(x) x, [0,1]);
result = cheb.revolution(f);

% Known exact results
exact.surfaceArea = pi*sqrt(2);
exact.volume = pi/3;
exact.centroidZ = 3/4;
exact.momentOfInertia = pi/10;

% Happy?
pass(1) = abs(result.surfaceArea - exact.surfaceArea) < 1e-15;
pass(2) = abs(result.volume - exact.volume) < 1e-15;
pass(3) = abs(result.centroidZ - exact.centroidZ) < 1e-15;
pass(4) = abs(result.momentOfInertia - exact.momentOfInertia) < 1e-15;

%% Check results for a surface of revolution on an unbounded domain
f = chebfun(@(x) exp(-x), [0,Inf]);
warnstate = warning;
warning('off', 'CHEBFUN:UNBNDFUN:sum:slowDecay')
result = cheb.revolution(f);
warning(warnstate);

% Known exact results
exact.surfaceArea = pi*(sqrt(2) + asinh(1));
exact.volume = pi/2;
exact.centroidZ = 1/2;
exact.momentOfInertia = pi/8;

% Happy?
pass(5) = abs(result.surfaceArea - exact.surfaceArea) < 1e-12;
pass(6) = abs(result.volume - exact.volume) < 1e-12;
pass(7) = abs(result.centroidZ - exact.centroidZ) < 1e-7;
pass(8) = abs(result.momentOfInertia - exact.momentOfInertia) < 1e-12;

end
