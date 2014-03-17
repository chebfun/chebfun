% Test file for ADCHEBFUN trigonometric and related functions, part 4/4

function pass = test_trig4

% List of trigonometric functions to test.
funcList = {@csch, @sec, @secd, @sech, @sin, @sinc, @sind, @sinh, @tan, ...
    @tand, @tanh};

% Call the ADCHEBFUN testUnary() method to do the tests.
pass = adchebfun.testUnary(funcList);

end