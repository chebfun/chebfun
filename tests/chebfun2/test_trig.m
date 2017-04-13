function pass = test_trig( pref ) 
% Test TRIG commands

if ( nargin == 0 ) 
    pref = chebfunpref; 
end

tol = 100*pref.cheb2Prefs.chebfun2eps;

f = {@cos, @sin, @tan, @cosh, @sinh, @tanh, @tand}; 

g = chebfun2(@(x,y) x.*y.^2); 
for jj = 1:numel(f)
    h = f{jj}(g); 
    exact = chebfun2(@(x,y) f{jj}(x.*y.^2)); 
    pass(jj) = ( norm(h - exact) < tol ); 
end

% Make a rank-1 'trig' chebfun2 from coeffs: 
u = rand(10,1); 
v = rand(10,1); 
u = chebfun(u,'trig'); 
v = chebfun(v,'trig'); 
exact = u*v.';
 
coeffs = trigcoeffs(u)*trigcoeffs(v).';
f = chebfun2(coeffs,'coeffs','trig');

pass(jj+1) = norm( f - exact ) < tol; 

end