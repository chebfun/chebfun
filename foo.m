clc

N = 20;
F = @(x) exp(x);
f = chebfun(F, N);
c = cheb2leg(f.coeffs);

x = [.9 ; -.999999];
x0 = sign(x);
rn = sign(x0);
pn = 1;
dn = 0*x;
s = c(1);
for n = 0:N-2
    dn = ((2*n+1)*(x-x0).*pn + n*dn)./((n+1)*rn);
    pn = (pn + dn).*rn;
    s = s + c(n+2)*pn;
end

s - F(x)


%%
clc

t = acos(x)
t0 = acos(x0)

t = [1e-5;  pi*.9999];
x = cos(t)



dx = x0 - x
idx = t > pi/2;
dx(idx) = -2*cos(t(idx)/2).^2;
dx(~idx) = 2*sin(t(~idx)/2).^2;
dx



sign(dx)
x0

    


