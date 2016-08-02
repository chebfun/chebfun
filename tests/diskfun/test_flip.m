function pass = test_flip( ) 
% Test diskfun flip and rotate functions. 

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

% See the random number generator.

%gaussian shifted right
f = @(x, y) 10*exp(-20*(x-.5).^2-20*(y+.5).^2) ;
h = @(x, y) 10*exp(-20*(-x-.5).^2-20*(y+.5).^2) ;
g = @(x, y) 10*exp(-20*(x-.5).^2-20*(-y+.5).^2) ;
f = diskfun(f); 
h = diskfun(h);
g = diskfun(g);

pass(1) =  ( norm( fliplr(f) - h )  < tol);  
pass(2) = ( norm( flipud(f) - g)  < tol) ; 

%test rotate

th = [0, pi/4, pi/3, 5*pi/7, -pi/6, -pi/3, -pi ]; 
[theta, rad] = meshgrid(trigpts(30, [-pi, pi]), chebpts(30, [0, 1])); 
w = @(t, r) exp(-20*(r.*cos(t)-.5).^2-20.*(r.*sin(t)).^2 ); 
f = diskfun(w, 'polar'); 
for j = 1:length(th) 
    g = rotate(f, th(j));
    pass(j) = ( norm(w(theta-th(j), rad)-g(theta, rad, 'polar'))  < tol);
end



end 

