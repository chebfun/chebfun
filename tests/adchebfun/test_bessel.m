% Test file for ADCHEBFUN Bessel functions

function pass = test_bessel

% List of airy functions to test, with different number of inputs
besselj1 = @(f) besselj(1, f);
besselj25 = @(f) besselj(2.5, f);

funcList = {besselj1, besselj25};

% Call the ADCHEBFUN testUnary() method to do the tests.
pass = adchebfun.testUnary(funcList);

end
