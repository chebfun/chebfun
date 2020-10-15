function pass = test_coeffs_vals( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e3*pref.techPrefs.chebfuneps; 
j = 1;
%grab some coeffs
u = diskfun(@(x,y) exp(-cos(pi*(x+y)))); 
v = diskfun(@(x,y) x.*sin(x.*y)); 
f = [u;v];

% test coeffs2: 
[x, y] = coeffs2(f); 
pass(j) = norm(diskfun.coeffs2diskfun(x)-u, inf) < tol;
pass(j+1) = norm(diskfun.coeffs2diskfun(y)-v, inf) < tol;
j = j+2; 

%parameters
[x, y] = coeffs2(f, 50, 60); 
pass(j) = norm(x-coeffs2(u, 50, 60), inf) < tol;
pass(j+1) = norm(y-coeffs2(v, 50, 60), inf) < tol;
j = j+2; 

% test coeffs2diskfunv: 
f2 = diskfunv.coeffs2diskfunv(coeffs2(u), coeffs2(v)); 
pass(j) = norm(f - f2) < tol; 
j = j+1;

%test coeffs2vals: 
[x, y] = coeffs2(f2); 
[u, v] = diskfunv.coeffs2vals(x,y);
pass(j) = norm(diskfun.coeffs2vals(x)-u, 'inf') < tol;
pass(j+1) = norm(diskfun.coeffs2vals(y)-v, 'inf')< tol;
j= j+2;    

%testvals2coeffs: 
[a,b] = diskfunv.vals2coeffs(u,v); 
pass(j) = norm(x - a, 'inf'); 
pass(j+1) = norm(y - b, 'inf'); 

if (nargout > 0)
    pass = all(pass(:));
end
end




