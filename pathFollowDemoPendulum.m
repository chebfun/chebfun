%% Path following test for a BVP, working with superlinops
% July, 2012
%
%% Setup and find initial solution using chebops
plotOn = 1;

d = [0 10]; % Domain of problem
x = chebfun('x',d);
N = chebop(@(x,u) diff(u,2) + sin(u), d); % Diff. eq. part
lbc = 2;
rbc = 2;
N.bc = @(x,u)[u(d(1))-lbc; u(d(2))-rbc];   % Boundary conditions

u0 = N\0;

disp(['Initial solution found. Norm of residual: ' num2str(norm(N.op(x,u0)))])
%% Iterate along path.
% Begin by finding a tangent direction, then set steplength, then compute
% Newton correction, and repeat.
H = chebop(@(x,u,lam) diff(u,2) + sin(u), d); % Diff. eq. part

lam0 = -feval(diff(u0),d(2))
H.bc = @(x,u,lam)[u(d(1))-lbc,feval(diff(u),d(2))+lam];   % Boundary conditions

measure = @(u) u(d(2));

maxCounter = 48;
D = operatorBlock.diff(d);
D2 = operatorBlock.diff(d, 2);
A = @(u,lam) D2 + operatorBlock.mult(cos(u));
g = @(u,lam) 0*u;

% Derivative of boundary conditions:
E = functionalBlock.eval(d);
BCmat = [E(d(1)) 0; E(d(2))*D 1];

[uquasi lamvec mvec] = followPath(H, A, g, BCmat, u0, lam0, measure, ...
    'plotOn',1, 'maxCounter',maxCounter,'slmax',1);

%% Find what solutions are closest to satisfying original BCs
bcResVec = zeros(numel(uquasi),1);
for uCounter = 1:numel(uquasi)
    bcResVec(uCounter) = norm(feval(uquasi(:,uCounter),d(2))-2);
end
goodSols = [1, 25, 29, 35, 42];

%% Use the good solutions as initial guesses for Newton iteration
r = chebfun;
for rCounter = 1:length(goodSols)
    N.init = uquasi(:,goodSols(rCounter));
    r = [r, N\0];
end
plot(r,'linewidth', 2)

%% Plot a nice diagram with all solutions
lamvec = lamvec(1:46);
mvec = mvec(1:46);
close all
figure(1)
hfig = figure(1);

xwidth = 800;
ywidth = 400;
set(hfig,'PaperPositionMode','auto')
set(hfig, 'Position', [640 558 xwidth ywidth])

hold on
plot(r,'linewidth',2)

xlabel('$t$','interpreter','latex','fontsize',14')
ylabel('$\theta(t)$','interpreter','latex','fontsize',14')

hl = legend('$\theta^{(1)}$', '$\theta^{(2)}$', '$\theta^{(3)}$', ...
    '$\theta^{(4)}$', '$\theta^{(5)}$',-1);
set(hl,'interpreter','latex')
set(hl,'fontsize',14)
xl = xlim;
% ylim([-.2 .4]);
yl = ylim;
set(gca,'xtick',xl(1):2:xl(2))
title('Trajectories of pendulum, satisfying $\theta''''+\sin(\theta) = 0, \, \theta(0) = \theta(10) =2$', ...
    'interpreter','latex','fontsize',16)

plot([0 10],[pi pi],'k--')
plot([0 10],-[pi pi],'k--')
box on
ylim([-3.5, 3.5])
set(gca,'ytick',-3:1.5:3)

hold off

%% Plot a nice bifurcation diagram
figure(1)
hfig = figure(1);

xwidth = 800;
ywidth = 400;
set(hfig,'PaperPositionMode','auto')
set(hfig, 'Position', [640 558 xwidth ywidth])

hold on
lamm = [lamvec, mvec]';
fnplt(cscvn(lamm),'b',2)
plot(lamvec,mvec,'ro','markersize',8,'linewidth',2)
lc = 3.5153;

xlabel('$\lambda$','interpreter','latex','fontsize',14')
ylabel('$\theta(10)$','interpreter','latex','fontsize',14')
xl = xlim;
yl = ylim;
plot([xl(1),xl(2)],[2,2],'k--','linewidth',2)
% set(gca,'xtick',xl(1):xl(2))
text(lc-.5,yl(2)-2.75,'$\lambda = \lambda_c \rightarrow$','interpreter','latex','fontsize',14)
title('Bifurcation diagram for $\theta''''+\sin(\theta) = 0,\, u(0) = 0,\, u''(0) = -\lambda$', ...
    'interpreter','latex','fontsize',16)
box on
hold off
