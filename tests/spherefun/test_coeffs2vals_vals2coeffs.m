function pass = test_coeffs2vals( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps; 
j = 1;
%grab some coeffs
f = spherefun( @(x,y,z) exp(-((x-.4).^2+(y-.9).^2+(z-.1).^2)) ); 

%test 1:
exact = coeffs2(f); 
cfs = spherefun.vals2coeffs(spherefun.coeffs2vals(exact));
pass(j) = norm( exact(:) - cfs(:), 'inf') < tol; 
j = j+1;

%test2: sample doubled grid and see if we get coeffs back: 
f = spherefun(@(x,y,z) sin(pi*x.*y) + sin(pi*x.*z));
[c, d, r] = cdr(f); 
g = c*d*r.'; 
lam = trigpts(50, [-pi, pi]);
th = trigpts(60, [-pi, pi]);
[tt,ll] = meshgrid(th, lam); 
F = g(tt,ll); 
cfs = spherefun.vals2coeffs(F); 
pass(j) = norm(cfs-coeffs2(f, 60, 50), 'inf') < tol;
j = j+1; 

%try the other direction: 
vals = spherefun.coeffs2vals(coeffs2(f)); 
[m,n] = size(vals); 
th = trigpts(n, [-pi, pi]); 
lam = trigpts(m, [-pi, pi]); 
[tt,ll] = meshgrid(th,lam); 
pass(j) = norm(g(tt,ll)-vals, 'inf')< 1e2*tol;
j = j+1;


% test low rank stuff: 
f = randnfunsphere(.3); 
[c, d, r] = cdr(f); 
c = c.coeffs; 
r = r.coeffs;
[u,s,v] = spherefun.coeffs2vals(c, d,r); 
vals = spherefun.coeffs2vals(coeffs2(f)); 
pass(j) = norm( vals -u*s*v.', 'inf') < 100*tol;
j = j+1; 


tt = trigpts(101, [-pi, pi]); 
ll = trigpts(101, [-pi, pi]); 
[c, d,r] = cdr(f); 
c = c(ll); 
r = r(tt); 
[u, s, v] = spherefun.vals2coeffs(c, d, r);
F = coeffs2(f, 101,101);
pass(j) = norm(F-u*s*v.', 'inf') < 3*tol; 


if (nargout > 0)
    pass = all(pass(:));
end

end
