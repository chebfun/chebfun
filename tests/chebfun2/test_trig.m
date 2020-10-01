function pass = test_trig( pref ) 
% Test TRIG commands

if ( nargin == 0 ) 
    pref = chebfunpref; 
end

tol = 100*pref.cheb2Prefs.chebfun2eps;

f = {@cos, @sin, @tan, @cosh, @sinh, @tanh, @tand}; 

g = chebfun2(@(x,y) x.*y.^2); 
for jj = 1:7
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

pass(8) = norm( f - exact ) < tol; 

% Try a few explicit ones

C0 = zeros(3);

C = C0; C(2,2) = 1;
f = chebfun2(C,'trig','coeffs');
pass(9) = norm( f - 1 ) < tol;

x = chebfun2(@(x,y) x); y = chebfun2(@(x,y) y);
C = C0; C(1,2) = 1;
f = chebfun2(C,'trig','coeffs');
pass(10) = norm( f - exp(-1i*pi*y) ) < tol;

C = C0; C(2,1) = .5; C(2,3) = .5;
f = chebfun2(C,'trig','coeffs');
pass(11) = norm( f - cos(pi*x) ) < tol;

f = chebfun2(@(x,y) sin(2*pi*x),   'trig');
g = chebfun2(@(x,y) cos(2*pi*x)+2, 'trig');
pass(12) = isPeriodicTech( f+1 );
pass(13) = isPeriodicTech( f.^2 );
pass(14) = isPeriodicTech( f./g );

end
