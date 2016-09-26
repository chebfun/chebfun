function pass = test_diff( ) 
% Test diff commands

tol = 1e3*chebfunpref().cheb2Prefs.chebfun2eps;
j=1;
% Check empty
u = diskfunv;
f = diffx(u);
pass(j) = isempty(f) & isa(f,'diskfunv'); j=j+1;
f = diffy(u);
pass(j) = isempty(f) & isa(f,'diskfunv'); j=j+1;
f = diff(u, 2, 2); 
pass(j) = isempty(f) & isa(f,'diskfunv'); j=j+1;

%check definitions: 
f = diskfun(@(x,y) cos(3*y.^2).*sin(2*x)); 
g = diskfun(@(x,y) sin(x.^4).*cos(5*y)); 
F = diskfunv(f, g);
G = diskfunv(diff(f, 1,1), diff(g, 1, 1)); 
pass(j) = ( norm(2*diffx(F)-diff(F, 1, 1)-G) < tol); j = j+1;
G = diskfunv(diff(f, 1, 3), diff(g, 1, 3)); 
pass(j) = ( norm(2*diffx(F, 3) - diff(F, 1, 3)-G) < tol); j = j+1;
G = diskfunv(diff(f, 2, 1), diff(g, 2, 1)); 
pass(j) = ( norm(2*diffy(F) - diff(F, 2, 1)-G) < tol); j = j+1; 
G = diskfunv(diff(f, 2, 2), diff(g, 2, 2)); 
pass(j) = ( norm(2*diffy(F, 2)-diff(F, 2, 2)-G)  < tol); j = j+1; 


% Simple examples: 
F = diskfunv(@(x,y) cos(x), @(x,y) sin(y)); 
exact = diskfunv(@(x,y) -sin(x), @(x,y) 0*x); 
pass(j) = (norm( diffx(F) - exact ) < tol); j = j+1;
exact = diskfunv(@(x,y) 0*x, @(x,y) cos(y)); 
pass(j) = ( norm(diffy(F)-exact) < tol ); j = j+1;

% 
f = diskfun(@(x,y) cos(3*y.^2).*sin(2*x)); 
g = diskfun(@(x,y) sin(pi*x.^4).*cos(5*y)); 
F = diskfunv(f, g);
exact = diskfunv( @(x,y) 2*cos(3*y.^2).*cos(2*x), ...
    @(x,y) 4*pi*x.^3.*cos(pi*x.^4).*cos(5*y) ); 
pass(j) = ( norm( diffx(F)-exact )  < tol); j = j+1; 
exact = diskfunv( @(x,y) -6*y.*sin(3*y.^2).*sin(2*x), ...
    @(x,y) -5*sin(5*y).*sin(pi*x.^4)); 
pass(j) = ( norm(diffy(F)-exact) ) < tol;
end