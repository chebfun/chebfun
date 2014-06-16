function pass = test_trig( pref ) 
% Test TRIG commands

if ( nargin == 0 ) 
    pref = chebfunpref; 
end

tol = 100*pref.eps; 

f = {@cos, @sin, @tan, @cosh, @sinh, @tanh, @tand}; 

g = chebfun2(@(x,y) x.*y.^2); 
for jj = 1:numel(f)
    h = f{jj}(g); 
    exact = chebfun2(@(x,y) f{jj}(x.*y.^2)); 
    pass(jj) = ( norm(h - exact) < tol ); 
end

end