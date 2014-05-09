function pass = test_feval_lr
% Check that left/right evaluation with linops behaves as expected.
% Nick Hale, May 2011
% TAD, 3 Feb 2014

% A tolerance to check to
tol = 1e-13;

d = [-1,1];
x = chebfun(@(x) x, d);
s = cos(x+pi/4).*sign(x)+.5;

%%
% Create the linops
L = functionalBlock.feval(0,d);
Ll = functionalBlock.feval(0,d,'left');
Lr = functionalBlock.feval(0,d,'right');

%%
% True left and right points, and an averaged centre point.
cl = -cos(pi/4)+0.5;
cr = cos(pi/4)+0.5;
c = (cl+cr)/2;

%%
% Check the forward operators
pass(1) = abs(L*s-c)<tol;
pass(2) = abs(Ll*s-cl)<tol;
pass(3) = abs(Lr*s-cr)<tol;

%%
% Check when discretized.
Bl = chebmatrix(Ll);
Br = chebmatrix(Lr);
Bl.domain = s.domain;  % must add breakpoint
Br.domain = s.domain;
Bs = chebmatrix(s);

%%
% TODO: We should no longer expect this to work:
% n = [31, 20];
% pass(4) = abs(matrix(Bl,n)*matrix(Bs,n)-cl)<tol;
% pass(5) = abs(matrix(Br,n)*matrix(Bs,n)-cr)<tol;

%%
% Check composition with the derivatives
s = cos(x+pi/4).*abs(x)+.5;
D = operatorBlock.diff(s.domain);
Al = Ll*D;
Ar = Lr*D;
cpl = -sqrt(2)/2;
cpr = sqrt(2)/2;
pass(4) = abs(Al*s-cpl) < 1e3*tol;
pass(5) = abs(Ar*s-cpr) < 1e3*tol;

