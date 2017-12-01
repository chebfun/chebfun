% Test file for gpr.m.

function pass = test_gpr( pref )

% set a seed for the RNG:
seedRNG(23411);
% Generate a few random points to use as test values:
xx = linspace(-1,1,10)'+2e-2*randn(10,1);
% The values to match
yy = exp(xx).*sin(2*xx);

[f,~,fsamples] = gpr(xx,yy,'samples',2);

% Test the accuracy of the posterior mean:
err = norm(f(xx) - yy,Inf);
pass(1) = err < 1e-6;

% Test the accuracy of the samples:
err = norm(fsamples(xx,1) - yy, Inf);
pass(2) = err < 5e-4;

err = norm(fsamples(xx,2) - yy, Inf);
pass(3) = err < 5e-4;

% Test adding noise to the data:
yy = yy + .1*randn(size(yy));
f = gpr(xx,yy,'noise',.1);
pass(4) = std(f(xx)-yy) < .25;


% Test horizontal scale (small interval):
xx = 2e-100*xx-1e-100;
yy = randn(10,1);
f = gpr(xx,yy);
err = norm(f(xx) - yy,Inf);
pass(5) = err < 1e-6;

% Test vertical scale (big function) and horizontal scale(small interval):
yy = 1e100*randn(10,1);
f = gpr(xx,yy);
err = norm((f(xx) - yy)./yy,Inf);
pass(6) = err < 1e-6;

% Test vertical scale (big function, fixed hyperparameter(s)):
xx = linspace(-1,1,10)';
yy = 1e100*sin(3*xx).*exp(xx);
f = gpr(xx,yy,'sigma',1e100,'L',.1);
err = norm((f(xx) - yy)./yy,Inf);
pass(7) = err < 1e-10;

f = gpr(xx,yy,'sigma',1e100);
err = norm((f(xx) - yy)./yy,Inf);
pass(8) = err < 1e-10;

f = gpr(xx,yy,'L',.1);
err = norm((f(xx) - yy)./yy,Inf);
pass(9) = err < 1e-10;


f = gpr(xx,yy);
err = norm((f(xx) - yy)./yy,Inf);
pass(10) = err < 1e-10;


% Test vertical scale (big function):
x = (1:5).^2; y = sin(x);
S = 1e50;
f = gpr(x,y);
ybig = S*y;
fbig = gpr(x,ybig);
err = norm(f-fbig/S);
pass(11) = err < 5e-14;

% Test vertical scale (small function):
ysmall = y/S;
fsmall = gpr(x,ysmall);
err = norm(f-S*fsmall);
pass(12) = err < 5e-14;

% Test horizontal scale (big interval):
x = (1:5).^2; y = sin(x);
S = 1e50;
f = gpr(x,y);
xbig = S*x;
fbig = gpr(xbig,y);
xx = linspace(1,25);
err = norm(f(xx)-fbig(S*xx));
pass(13) = err < 5e-14;

% Test horizontal scale (small interval):
xsmall = x/S;
fsmall = gpr(xsmall,y);
xx = linspace(1,25);
err = norm(f(xx)-fsmall(xx/S));
pass(14) = err < 5e-14;


% Test periodic version of the code:
N = 40;
xx = linspace(-1,1,N)';
xx(2:end-1) = xx(2:end-1)+1e-3*randn(N-2,1);
% The values to match

yy = exp(sin(pi*xx));

[f,~,fsamples] = gpr(xx,yy,'domain',[-1,1],'trig','samples',2);

err = norm(f(xx) - yy,Inf);
pass(15) = err < 1e-6;

% Test the accuracy of the samples:
err = norm(fsamples(xx,1) - yy, Inf);
pass(16) = err < 5e-4;

err = norm(fsamples(xx,2) - yy, Inf);
pass(17) = err < 5e-4;

% Test vertical scale (big function):
x = (1:5).^2; y = sin(x);
S = 1e50;
f = gpr(x,y,'trig');
ybig = S*y;
fbig = gpr(x,ybig,'trig');
err = norm(f-fbig/S);
pass(18) = err < 5e-14;

% Test vertical scale (small function):
ysmall = y/S;
fsmall = gpr(x,ysmall,'trig');
err = norm(f-S*fsmall);
pass(19) = err < 5e-14;

% Test horizontal scale (big interval):
x = (1:5).^2; y = sin(x);
S = 1e50;
f = gpr(x,y,'trig');
xbig = S*x;
fbig = gpr(xbig,y,'trig');
xx = linspace(1,25);
err = norm(f(xx)-fbig(S*xx));
pass(20) = err < 5e-14;

% Test horizontal scale (small interval):
xsmall = x/S;
fsmall = gpr(xsmall,y,'trig');
xx = linspace(1,25);
err = norm(f(xx)-fsmall(xx/S));
pass(21) = err < 5e-14;

end
