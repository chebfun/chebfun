function pass = test_coeffs2vals( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps; 
j = 1;
%grab some coeffs
f = diskfun(@(x,y) sin(4*((3*x-.4).^2 + y)));

%test 1:
exact = coeffs2(f); 
cfs = diskfun.vals2coeffs(diskfun.coeffs2vals(exact));
pass(j) = norm( exact(:) - cfs(:), 'inf') < tol; 
j = j+1;

%test2: sample doubled grid and see if we get coeffs back: 
f = diskfun(@(x,y) exp(-2*x.^3 - (y-.4).^3)); 
[c, d, r] = cdr(f); 
g = c*d*r.'; 
r = chebpts(51);
t = trigpts(106, [-pi, pi]);
[tt,rr] = meshgrid(t,r); 
F = g(tt,rr); 
cfs = diskfun.vals2coeffs(F); 
pass(j) = norm(cfs-coeffs2(f, 106, 51), 'inf') < tol;
j = j+1; 

%try the other direction: 
vals = diskfun.coeffs2vals(coeffs2(f)); 
[m,n] = size(vals); 
t = trigpts(n, [-pi, pi]); 
r = chebpts(m); 
[tt,rr] = meshgrid(t,r); 
pass(j) = norm(g(tt,rr)-vals, 'inf')< 1e2*tol;
j = j+1;


% test low rank stuff: 
f = randnfundisk(.3); 
[c, d, r] = cdr(f); 
c = c.coeffs; 
r = r.coeffs;
[u,s,v] = diskfun.coeffs2vals(c, d,r); 
vals = diskfun.coeffs2vals(coeffs2(f)); 
pass(j) = norm( vals -u*s*v.', 'inf') < 10*tol;
j = j+1; 


tt = trigpts(100, [-pi, pi]); 
rr = chebpts(51); 
[c, d,r] = cdr(f); 
c = c(rr); 
r = r(tt); 
[u, s, v] = diskfun.vals2coeffs(c, d, r);
F = coeffs2(f, 100,51);
pass(j) = norm(F-u*s*v.', 'inf') < 3*tol; 


if (nargout > 0)
    pass = all(pass(:));
end

end
