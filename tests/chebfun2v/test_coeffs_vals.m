function pass = test_coeffs_vals( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e3*pref.techPrefs.chebfuneps; 
j = 1;
%grab some coeffs
u = randnfun2(.3); 
v = randnfun2(.4); 
f = [u;v];

% test coeffs2: 
[x, y] = coeffs2(f); 
pass(j) = norm(chebfun2(x, 'coeffs')-u, inf) < tol;
pass(j+1) = norm(chebfun2(y, 'coeffs')-v, inf) < tol;
j = j+2; 

%parameters
[x, y] = coeffs2(f, 50, 60); 
pass(j) = norm(x-coeffs2(u, 50, 60), inf) < tol;
pass(j+1) = norm(y-coeffs2(v, 50, 60), inf) < tol;
j = j+2; 

% test coeffs2chebfun2v: 
f2 = chebfun2v.coeffs2chebfun2v(coeffs2(u), coeffs2(v)); 
pass(j) = norm(f - f2) < tol; 
j = j+1;

%test coeffs2vals: 
[x, y] = coeffs2(f2); 
[u, v] = chebfun2v.coeffs2vals(x,y);
pass(j) = norm(chebfun2.coeffs2vals(x)-u, 'inf') < tol;
pass(j+1) = norm(chebfun2.coeffs2vals(y)-v, 'inf')< tol;
j= j+2;    

%testvals2coeffs: 
[a,b] = chebfun2v.vals2coeffs(u,v); 
pass(j) = norm(x - a, 'inf'); 
pass(j+1) = norm(y - b, 'inf'); 

if (nargout > 0)
    pass = all(pass(:));
end
end




