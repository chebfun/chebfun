% Test file for ADCHEBFUN trigonometric and related functions, part 1/4

function pass = test_trig1

% List of trigonometric functions to test.
funcList = {@acos, @acosd, @acosh, @acot, @acotd, @acoth, @acsc, @acscd};

% Call the ADCHEBFUN testUnary() method to do the tests.
pass = adchebfun.testUnary(funcList);

end