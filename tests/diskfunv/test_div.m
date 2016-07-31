function pass = test_div( ) 
% Test divergence

tol = 1e2*chebfunpref().cheb2Prefs.chebfun2eps;

% Check divergence of an empty diskfunv is an empty diskfun.
u = diskfunv;
f = div(u);
pass(1) = isempty(f) & isa(f,'diskfun');

% Check definition: 
F = diskfunv(@(x,y) cos(x), @(x,y) sin(y)); 
divF = diff(F(1), 1,1) + diff(F(2), 2, 1);   
pass(2) = ( norm(divF - divergence(F) ) < tol );

% Check divergence of the zero field is zero
f = diskfun(@(x,y) 0*x);
u = diskfunv(f,f);
pass(3) = norm(div(u)) < tol; 

% Check divergence of gradient is laplacian: 
f = diskfun(@(x,y) cos(x.*y)); 
pass(4) = ( norm(laplacian( f ) - divergence( gradient( f ) ) ) < tol );

% Check that the div and divergence give the same result.
divergenceu = divergence(u);
divu = div(u); 
pass(5) = isequal(divu,divergenceu);

end