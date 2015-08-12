f = @(x,y) y.^2.*exp(sin(x) +cos(5*y)-1);
g = @(th,r) (r.*sin(th)).^2.*exp(sin(r.*cos(th))+cos(5*r.*sin(th))-1); 
%f = @(x,y) 1+x +y + x.^2;
%f = @(x,y) 
%f = @(x,y) sin(x) + cos(5*y); 
%g = @(th,r) 1+r.*cos(th)+r.*sin(th)+(r.*cos(th)).^2;
% g = @(th, r) sin(r.*cos(th))+ cos(5*r.*sin(th));
t = diskfun(f,0); 
u = chebfun2(g, [-pi pi 0 1]); 

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
