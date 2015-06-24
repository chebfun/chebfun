%f = @(x,y) y.^2.*exp(sin(x) +cos(5*y)-1);
%f = @(x,y) 1+x +y + x.^2;
%f = @(x,y) 
%f = @(x,y) sin(x) + cos(5*y); 
%g = @(th,r) 1+r.*cos(th)+r.*sin(th)+(r.*cos(th)).^2;
% g = @(th,r) (r.*sin(th)).^2.*exp(sin(r.*cos(th))+cos(5*r.*sin(th))-1); 
% g = @(th, r) sin(r.*cos(th))+ cos(5*r.*sin(th));
%t = diskfun(f); 
%u = chebfun2(g, [-pi pi 0 1]); 

%rank(t)
%rank(u)


%choice of pivot values
%clf
%subplot(1,2,1)
%surf(u)
%subplot(1,2,2)
%surf(t)
%clf


%look at cdr coeffs
%[ct dt rt] = cdr(t);
%[cu du ru] = cdr(u);


%semilogy(1./sort(abs(max(dt))), 'rx-')
%hold on
%semilogy(1./sort(abs(max(du))), 'ko-')




f = @(x,y,z) exp(cos(2*pi*x)).*exp(sin(2*pi*y)) + exp(sin(4*z)); 
%f = @(x,y,z) x.*y.*z;
%f = @(x,y,z) cos(1 + 2*pi*(x+y) + 5*sin(pi*z));


g = spherefun(f);g.domain
q = spherefun(f, '2by2');

%colatitude coc


l = @(lam, th) exp(cos(2*pi*(cos(lam).*sin(th)))).*exp(sin(2*pi*sin(lam).*sin(th)))+exp(sin(4*cos(th)));
%l = @(lam, th) cos(lam).*sin(th) .*sin(lam).*sin(th) .*cos(th);
%l = @(lam, th) cos(1+2*pi*(cos(lam).*sin(th)+sin(lam).*sin(th))+5*sin(pi*cos(th)));

m = chebfun2(l, [-pi pi 0 pi]);m.domain

rank(m)
rank(g)

clf
subplot(1,2,1)
surf(g)
subplot(1,2,2)
surf(m)
shg

%look at cdr coeffs
[cm dm rm] = cdr(m);
[cg dg rg] = cdr(g);

size(dm)

clf
semilogy(1./sort(abs(max(dm))), 'rx-')
hold on
semilogy(1./sort(abs(max(dg))), 'ko-')





