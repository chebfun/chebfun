function pass = test_deflationFun(~)

% An operator to test and a function to deflate
x = chebfun('x');
Nop = @(u) diff(u,2) + sin(u);
r = sin(x);

% L2 deflation
defFun1 = @(u) deflationFun(Nop(u), u, r, 2, .5, 'L2');
defFun2 = @(u) deflationFun(Nop(u), u, [r exp(r)], 2, .5, 'L2');
defFun3 = @(u) deflationFun(Nop(u), u, [r exp(r) sin(r)], 2, .5, 'L2');
defFun4 = @(u) deflationFun(Nop(u), u, [r exp(r) sin(r)], 3, 2, 'L2');
defFun5 = @(u) deflationFun(Nop(u), u, [r exp(r)], 1, 0, 'L2');
defFun6 = @(u) deflationFun(Nop(u), u, [r exp(r)], 1, 4, 'L2');

% H1 deflation
defFun7  = @(u) deflationFun(Nop(u), u, r, 2, .5, 'H1');
defFun8  = @(u) deflationFun(Nop(u), u, [r exp(r)], 2, .5, 'H1');
defFun9  = @(u) deflationFun(Nop(u), u, [r exp(r) sin(r)], 2, .5, 'H1');
defFun10 = @(u) deflationFun(Nop(u), u, [r exp(r) sin(r)], 3, 2, 'H1');
defFun11 = @(u) deflationFun(Nop(u), u, [r exp(r)], 1, 0, 'H1');
defFun12 = @(u) deflationFun(Nop(u), u, [r exp(r)], 1, 4, 'H1');


funcList = {defFun1, defFun2, defFun3, defFun4, defFun5, defFun6, ...
    defFun7, defFun8, defFun9, defFun10, defFun11, defFun12};

% Call the ADCHEBFUN testUnary() method to do the tests.
pass = adchebfun.testUnary(funcList);

end