%% The Bratu problem
%
% For the BVP
%
%   u'' + lambda*exp(u) = 0,    u(0) = u(1) = 0
%
% it is known that there exist:
%   (i)   two solutions for lambda < lambda_crit ~ 3.5153;
%   (ii)  one solution for lambda = lambda_crit
%   (iii) no solutions for lambda > lambda_crit
%
% We now obtain a bifurcation diagram for this problem.

%% Find the initial solution:
d = [0 1];
lam0 = 0.01;
x = chebfun('x', d);
% Chebop for finding the initial solution:
N = chebop(@(x,u) diff(u,2) + lam0.*exp(u),d);
N.bc = @(x, u)[u(d(1)); u(d(2))];   % Boundary conditions
u0 = N\0;
disp(['Initial solution found. Norm of residual: ' num2str(norm(N.op(x,u0)))])
plot(u0) 
title(sprintf('Initial solution for $\\lambda=%2.2f$', lam0), ...
    'interpreter','latex')

%% Now trace out the curve in (u, lambda) space
% Chebop for tracing the curve:
H = chebop(@(x,u,lam) diff(u,2) + lam.*exp(u),d);
H.bc = @(x,u,lam)[u(d(1)); u(d(2))];
% What we draw on the y-axis of the bifurcation diagram:
measure = @(u) u((d(2)+d(1))/2);
% How many steps along the curve we want to take
maxCounter = 14;
% Maximum allowed stepsize
slmax = 1;
% Manually construct the derivative. Always computing the derivative via AD
% would be expensive, but if we could differentiate the computational graph once
% and for all, we'd get away with cheap derivatives!

% Derivative of differential equation part:
D2 = operatorBlock.diff(d, 2);
A = @(u,lam) D2 + operatorBlock.mult(lam.*exp(u));
g = @(u,lam) exp(u);

% Derivative of boundary conditions:
E = functionalBlock.eval(d);
BCmat = [E(d(1)) 0; E(d(2)) 0];

% Call the path-following method:
[uquasi, lamvec, mvec] = followPath(H, A, g, BCmat, u0, lam0, measure, ...
    'plotOn', 1, 'maxCounter',maxCounter,'slmax',slmax);

%% Plot a nice bifurcation diagram
figure

hfig = figure(1);

xwidth = 800;
ywidth = 400;
set(hfig,'PaperPositionMode','auto')
set(hfig, 'Position', [640 558 xwidth ywidth])

hold on
lamm = [lamvec, mvec]';
fnplt(cscvn(lamm),'b',2)
plot(lamvec,mvec,'ro','markersize',10,'linewidth',2)
lc = 3.5153;

xlabel('$\lambda$','interpreter','latex','fontsize',14')
ylabel('$u\left(\frac{1}{2}\right)$','interpreter','latex','fontsize',14')
yl = ylim;
plot([lc,lc],[yl(1), yl(2)],'k--','linewidth',2)
text(lc-.5,yl(2)-2.75,'$\lambda = \lambda_c \rightarrow$', ...
    'interpreter', 'latex', 'fontsize', 14)
title('Bifurcation diagram for $u''''+\lambda e^u = 0, u(0) = u(1) = 0$', ...
    'interpreter','latex','fontsize',16)
box on
hold off