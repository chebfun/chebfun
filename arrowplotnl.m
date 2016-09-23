% This code has been hacked together, not yet
% properly written for Chebfun, for the sake of the
% draft of the book "Exploring ODEs".
function arrowplotnl(u,v,dx) % plot chebfun with arrow on end
f = u + 1i*v;
fp = diff(f);
L = .5*diff(domain(f));  % half length of domain
dt = dx/abs(fp(end)); dt = min(dt,L);
for j = 1:1, dt = dt*dx/abs(f(end)-f(end-dt)); dt = min(dt,L); end
df = f(end)-f(end-dt);
seg1 = f(end) + df*chebfun([0 exp(.91i*pi)].');
seg2 = f(end) + df*chebfun([exp(.91i*pi) 0].');
seg3 = f(end) + df*chebfun([0 exp(-.91i*pi)].');
seg4 = f(end) + df*chebfun([exp(-.91i*pi) exp(.91i*pi)].');
plot(join(f,seg1,seg2,seg3,seg4),'color',[0 .45 0])
