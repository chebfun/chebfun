function pass = test_bvp4c(pref)

% [TODO]: Make this test more extensive.

if ( nargin == 0 )
    pref = chebfun.pref();
end

% This is copied from Rodrigo Platte's Jan 2009 test, bvp45ctest.m
d = [0, 4];
y0 = chebfun([1, 0], d, pref);
solinit = bvpinit([0, 1, 2, 3, 4], [1, 0]); 
% Test bvp4c using default tolerance (RelTol = 1e-3)
y = bvp4c(@twoode, @twobc, y0);         % Chebfun solution
sol = bvp4c(@twoode, @twobc, solinit);  % Matlab's solution
pass(1) = max(max(abs(sol.y' - feval(y,sol.x')))) < 2e-2;

end
