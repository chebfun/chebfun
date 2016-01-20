function pass = test_deflationFun(~)

% List of trigonometric functions to test.
x = chebfun('x');
Nop = @(u) diff(u,2) + sin(u);
r = sin(x);
defFun1 = @(u) deflationFun(Nop(u), u, r, 2, .5);
defFun2 = @(u) deflationFun(Nop(u), u, [r exp(r)], 2, .5);
defFun3 = @(u) deflationFun(Nop(u), u, [r exp(r) sin(r)], 2, .5);
defFun4 = @(u) deflationFun(Nop(u), u, [r exp(r) sin(r)], 3, 2);
defFun5 = @(u) deflationFun(Nop(u), u, [r exp(r)], 1, 0);
defFun6 = @(u) deflationFun(Nop(u), u, [r exp(r)], 1, 4);

% defFun1 = @(u) deflationFun(Nop(u), u, [r exp(r)], 2, .5);
% defFun1 = @(u) deflationFun(Nop(u), u, [r exp(r)], 2, .5);
% defFun1 = @(u) deflationFun(Nop(u), u, [r exp(r)], 2, .5);
% 
% defFun2 = @(u) deflationFunH1(Nop(u), u, [r exp(r)], 2, .5);
% normFun1 = @(u) norm(u, 2);
% normFun2 = @(u) norm(u, 2).^2;
% normFun3 = @(u) 1./norm(u, 2);
% normFun4 = @(u) 1./norm(u, 2).^2;
% 
% normFun5 = @(u) normsq(u) + normsq(diff(u));
% normFun6 = @(u) normsq(u) + normsq(diff(u));
% funcList = {normFun1, normFun2, normFun3, normFun4, normFun5, normFun6};

funcList = {defFun1, defFun2, defFun3, defFun4, defFun5, defFun6};
% funcList = {defFun3};
% Call the ADCHEBFUN testUnary() method to do the tests.
pass = adchebfun.testUnary(funcList);

end