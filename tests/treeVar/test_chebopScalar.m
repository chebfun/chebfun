function pass = test_chebopScalar

mu = 10;
dom = [0, 100];
vdpFun = @(y) diff(y, 2) - mu.*(1-y.^2).*diff(y) + y;
N = chebop(vdpFun, dom);
N.lbc = @(u) [u-2; diff(u)];
y = N\0;
figure
plot(y)
% Show where breakpoints got introduced
hold on, plot(y.domain, y(y.domain),'k*')
shg
pass = 1;
end