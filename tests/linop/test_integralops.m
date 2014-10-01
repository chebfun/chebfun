function pass = test_integralops
% Test integral operators
%
% Toby Driscoll 28 May 2009
% Nick Hale 6 Jan 2012
% Toby Driscoll 11 March 2014

tol = 1e-10;
method = {@chebcolloc2, @chebcolloc1};
kind = [2 1];

%%
% Fredholm

for k = 1:length(kind)
    d = [0 1];
    x = chebfun(@(x) x, d,'chebkind',kind(k));
    prefs = cheboppref;
    prefs.discretization = method{k};
    F = operatorBlock.fred(@(x,y) sin(2*pi*(x-y)), d);
    A = linop( operatorBlock.eye(d) + F );
    u = chebmatrix( x.*exp(x) );
    f = A*u;
    v = linsolve(A,f,prefs);
    err(k,1) = norm(u{1}-v{1});
end

%%
% Volterra
for k = 1:length(kind)
    
    d = [0,pi];
    x = chebfun(@(x) x, d,'chebkind',kind(k));
    prefs = cheboppref;
    prefs.discretization = method{k};
    V = operatorBlock.volt(@(x,y) x.*y, d);
    f = chebmatrix( x.^2.*cos(x) + (1-x).*sin(x) );
    A = linop( operatorBlock.eye(d) - V );
    u = linsolve(A,f,prefs);
    Au = A*u;
    err(k,2) = norm( u{1} - sin(x) ) ;
    err(k,3) = norm( Au{1} - f{1} ) ;
end

pass = err(:).' < tol;

end
