function pass = test_BMCsvd( ) 
% Test diskfun BMCsvd() command 

tol = 2e3*chebfunpref().cheb2Prefs.chebfun2eps;

f = diskfun(@(x,y) sin(5*pi*x.*y) + sin(pi*x.^2.*cos(2.*(y-.1))));
s = BMCsvd( f ); 

% Check resolved: 
pass(1) = ( s(end) < 1e3*tol );

% Scale invariant: 
g = 100*f; 
t = BMCsvd( g ); 
pass(2) = norm( s - t/100 ) < tol; 

% orthonormal cols
%h = diskfun(@(t,r) r.*cos(t-.2)+2*r.^2.*sin(2*t)); 

[u, s, v] = BMCsvd(g); 
pass(3) = abs(sum(u(:, 2).*u(:, 27))+sum(u(:,12).*u(:,7))+sum(v(:,28).*v(:,21))+...
    sum(v(:,1).*v(:,2))+sum(u(:,19).*u(:,19))+sum(v(:,1).*v(:,1))-2) < tol;

end 