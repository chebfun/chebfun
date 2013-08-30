function pass = test_bvp5c(pref)

% [TODO]: Make this test more extensive.

if ( nargin == 0 )
    pref = chebfun.pref();
end

% This is modified from Rodrigo Platte's Jan 2009 test, bvp45ctest.m
d = [0, 4];
y0 = chebfun([1, 0], d, pref);
solinit = bvpinit([0, 1, 2, 3, 4], [1, 0]); 
% Test bvp5c using default tolerance (RelTol = 1e-6)
opts = odeset('RelTol', 1e-6);
y = bvp5c(@twoode, @twobc, y0, opts);         % Chebfun solution
sol = bvp5c(@twoode, @twobc, solinit, opts);  % Matlab's solution
pass(1) = max(max(abs(sol.y' - feval(y,sol.x')))) < 1e-4;

end
