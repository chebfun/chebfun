function pass = test_coeffs_vals( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e3*pref.techPrefs.chebfuneps; 
j = 1;
%grab some coeffs
u = randnfunsphere(.3); 
v = randnfunsphere(.4); 
w = randnfunsphere(.3); 
f = [u;v;w];

% test coeffs2: 
[x, y, z] = coeffs2(f); 
pass(j) = norm(spherefun.coeffs2spherefun(x)-u, inf) < tol;
pass(j+1) = norm(spherefun.coeffs2spherefun(y)-v, inf) < tol;
pass(j+2) = norm(spherefun.coeffs2spherefun(z)-w, inf) < 10*tol;
j = j+3; 

%parameters
[x, y, z] = coeffs2(f, 50, 60); 
pass(j) = norm(x-coeffs2(u, 50, 60), inf) < tol;
pass(j+1) = norm(y-coeffs2(v, 50, 60), inf) < tol;
pass(j+2) = norm(z-coeffs2(w, 50, 60), inf) < tol;
j = j+3; 

% test coeffs2diskfunv: 
f2 = spherefunv.coeffs2spherefunv(coeffs2(u), coeffs2(v),...
    coeffs2(w)); 
pass(j) = norm(f - f2) < tol; 
j = j+1;

%test coeffs2vals: 
[x, y, z] = coeffs2(f2); 
[u, v, w] = spherefunv.coeffs2vals(x,y,z);
pass(j) = norm(spherefun.coeffs2vals(x)-u, 'inf') < tol;
pass(j+1) = norm(spherefun.coeffs2vals(y)-v, 'inf')< tol;
pass(j+2) = norm(spherefun.coeffs2vals(z)-w, 'inf')< tol;
j= j+2;    

%testvals2coeffs: 
[a,b,c] = spherefunv.vals2coeffs(u,v, w); 
pass(j) = norm(x - a, 'inf'); 
pass(j+1) = norm(y - b, 'inf'); 
pass(j+2) = norm(z - c, 'inf'); 

if (nargout > 0)
    pass = all(pass(:));
end
end




