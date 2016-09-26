function pass = test_trigs(pref) 
% Test some trigonometric functions

if ( nargin == 0 ) 
    pref = chebfunpref; 
end
tol = 100*pref.cheb3Prefs.chebfun3eps;

f = {@cos, @sin, @tan, @cosh, @sinh, @tanh, @tand}; 

g = chebfun3(@(x,y,z) x.*y.^2.*z.^3); 
for jj = 1:numel(f)
    h = f{jj}(g); 
    exact = chebfun3(@(x,y,z) f{jj}(x.*y.^2.*z.^3));
    pass(jj) = norm(h - exact) < tol; 
end

end