function pass = test_vanDerPol
dom = [0 40];
mu = 1;
vdpFun = @(y) diff(y, 2) - mu.*(1-y.^2).*diff(y) + y;
anonFun = treeVar.toFirstOrder(vdpFun, chebmatrix(0), dom);
tic
[t,y]=ode113(anonFun, dom, [2 0]);
fprintf('Solution time, conversion from 2nd order to 1st order: %4.4fs.\n', toc);
subplot(1,2,1), plot(t,y(:,1)); title('Automatically converted')

tic
[t,y]=ode113(@vdp1,[0 40],[2 0]);
fprintf('Solution time, built-in demo: %4.4fs.\n', toc);
subplot(1,2,2), plot(t,y(:,1)); title('Built in demo')


%% Crank up mu
mu = 40;
vdpFun = @(y) diff(y, 2) - mu.*(1-y.^2).*diff(y) + y;
anonFun = treeVar.toFirstOrder(vdpFun, chebmatrix(0), dom);
tic
[t,y]=ode113(anonFun, dom, [2 0]);
fprintf('Solution time, conversion 2nd->1st, larger mu: %4.4fs.\n', toc);
figure
plot(t,y(:,1)); title('Automatically converted')

pass = 1;
end