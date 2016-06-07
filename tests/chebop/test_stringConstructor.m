function pass = test_stringConstructor(~)
%TEST_STRINGCONSTRUCTOR    Test that we can pass strings to the CHEBOP
%constructor

%% Linear problem, no x dependence
% String syntax
L1 = chebop('u``+u',[0 5]);
L1.lbc = 1;
L1.rbc = 2;
u1 = L1\0;

% Anonymous function syntax
L2 = chebop(@(u) diff(u,2) + u, [0 5]);
L2.lbc = 1;
L2.rbc = 2;
u2 = L2\0;

% These should be identical!
pass(1) = norm(u1-u2) == 0;

%% Linear problem, x dependence
% String syntax
L1 = chebop('u``+x*u',[0 5]);
L1.lbc = 1;
L1.rbc = 2;
u1 = L1\0;

% Anonymous function syntax
L2 = chebop(@(x,u) diff(u,2) + x.*u, [0 5]);
L2.lbc = 1;
L2.rbc = 2;
u2 = L2\0;

% These should be identical!
pass(2) = norm(u1-u2) == 0;


%% Nonlinear problem, no x dependence
% String syntax
tic
N1 = chebop('u``+sin(u)',[0 5]);
N1.lbc = 1;
N1.rbc = 2;
u1 = N1\0;

% Anonymous function syntax
N2 = chebop(@(u) diff(u,2) + sin(u), [0 5]);
N2.lbc = 1;
N2.rbc = 2;
u2 = N2\0;

% These should be identical!
pass(3) = norm(u1-u2) == 0;

%% Linear problem, x dependence
% String syntax
N1 = chebop('u``+x*sin(u)',[0 5]);
N1.lbc = 1;
N1.rbc = 2;
u1 = N1\0;

% Anonymous function syntax
N2 = chebop(@(x,u) diff(u,2) + x.*sin(u), [0 5]);
N2.lbc = 1;
N2.rbc = 2;
u2 = N2\0;

% These should be identical!
pass(4) = norm(u1-u2) == 0;