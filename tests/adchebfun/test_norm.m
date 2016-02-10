% Test file for ADCHEBFUN norm related functions

function pass = test_norm

% List of trigonometric functions to test.
normFun1 = @(u) norm(u, 2);
normFun2 = @(u) norm(u, 2).^2;
normFun3 = @(u) 1./norm(u, 2);
normFun4 = @(u) 1./norm(u, 2).^2;
normFun5 = @(u) norm(u, 2) + norm(diff(u));
normFun5 = @(u) norm(u, 2).^2 + norm(diff(u)).^3;

funcList = {normFun1, normFun2, normFun3, normFun4, normFun5};

% Call the ADCHEBFUN testUnary() method to do the tests.
pass = adchebfun.testUnary(funcList);

end
