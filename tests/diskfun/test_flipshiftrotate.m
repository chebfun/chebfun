function pass = test_flipshiftrotate( ) 
% Test diskfun flip, rotate, flipxy and circshift functions. 

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

% See the random number generator.

%gaussian shifted 
f1 = @(x, y) 10*exp(-20*(x-.5).^2-20*(y+.5).^2) ;
h1 = @(x, y) 10*exp(-20*(-x-.5).^2-20*(y+.5).^2) ;
g1 = @(x, y) 10*exp(-20*(x-.5).^2-20*(-y+.5).^2) ;
f = diskfun(f1); 
h = diskfun(h1);
g = diskfun(g1);
w = diskfun(@(x,y)x).*(g +h +f)+diskfun(@(x,y) y);
r1 = @(x,y) x.*(f1(x,y) + h1(x,y) + g1(x,y))+y;
ud = diskfun(@(x,y) r1(x,-y));
lr = diskfun(@(x,y) r1(-x, y)); 
tp = diskfun(@(x,y) r1(y,x)); 
pass(1) =  ( norm( fliplr(f) - h )  < tol);  
pass(2) = ( norm( flipud(f) - g)  < tol) ; 
pass(3) = ( norm( flipud(w) - ud) < tol);
pass(4) = ( norm( fliplr(w) -lr) < tol);
pass(5) = ( norm( flipxy(w)-tp) < tol); 
pass(6) = ( norm( flipdim(w,1)-ud) < tol); 
pass(7) = ( norm( flipdim(w,2)-lr) < tol); 

%test rotate/circshift
th = [0, pi/2, 5*pi/7,  -pi/3, -pi/4, -5*pi ]; 
[theta, rad] = meshgrid(trigpts(30, [-pi, pi]), chebpts(30, [0, 1])); 
w = @(t, r) exp(-20*(r.*cos(t)-.5).^2-20.*(r.*sin(t)).^2 ); 
f = diskfun(w, 'polar'); 
for j = 1:length(th) 
    g = rotate(f, th(j));
    m = circshift(f, th(j)); 
    pass(j+7) = ( norm(w(theta-th(j), rad)-g(theta, rad, 'polar'))  < tol ...
        || norm(g-m) ==0);
end

end 

