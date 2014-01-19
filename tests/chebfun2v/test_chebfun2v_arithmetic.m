function pass = chebfun2v_arithmetic
% Check the Chebfun2v constructor for simple arithmetic operations. 
% Alex Townsend, March 2013. 

% These function chosen so that scl does not change. 
f = @(x,y) cos(x); f=chebfun2v(f,f); 
g = @(x,y) sin(y); g=chebfun2v(g,g); 
% exact answers. 
plus_exact = @(x,y) cos(x) + sin(y); plus_exact=chebfun2v(plus_exact,plus_exact); 
minus_exact = @(x,y) cos(x) - sin(y); minus_exact=chebfun2v(minus_exact,minus_exact); 
mult_exact = @(x,y) cos(x).*sin(y); mult_exact=chebfun2v(mult_exact,mult_exact); 
pow_exact = @(x,y) cos(x).^sin(y); pow_exact=chebfun2v(pow_exact,pow_exact);

tol = 1e-14;
try 
    pass(1) = norm(f + g - plus_exact) < tol;
%     pass(2) = norm(f + g - minus_exact) < tol
%     pass(1) = isequal(f + g , plus_exact);
%     pass(2) = isequal(f - g , minus_exact);
%     pass(3) = isequal(f.*g , mult_exact);
%     pass(4) = isequal(f.^g , pow_exact);
    pass = all(pass); 
catch
    pass = 0;
end
end
