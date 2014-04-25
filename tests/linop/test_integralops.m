function pass = test_integralops
% Test integral operators
%
% Toby Driscoll 28 May 2009
% Nick Hale 6 Jan 2012
% Toby Driscoll 11 March 2014

tol = 1e-10;
method = {@colloc2, @colloc1};
kind = [2 1];

%%
% Fredholm

for k = 1:length(kind)
    d = [0 1];
    x = chebfun(@(x) x, d,'chebkind',kind(k));
    prefs = cheboppref;
    prefs.discretization = method{k};
    F = operatorBlock.fred(d,@(x,y) sin(2*pi*(x-y)));
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
    V = operatorBlock.volt( d, @(x,y) x.*y );
    f = chebmatrix( x.^2.*cos(x) + (1-x).*sin(x) );
    A = linop( operatorBlock.eye(d) - V );
    u = linsolve(A,f,prefs);
    Au = A*u;
    err(k,2) = norm( u{1} - sin(x) ) ;
    err(k,3) = norm( Au{1} - f{1} ) ;
end

%% TODO: chebop version
%
% % Fredholm
% d = [0,1];
% x = chebfun(@(x) x, d);
% K = @(x,y) sin(2*pi*(x-y));
% A = chebop(@(u) u + fred(K,u), d);
% u = x.*exp(x);
% f = A*u;
% pass(4) = norm(u-A\f) < 1e6*tol;
%
% % Volterra
% d = [0,pi];
% x = chebfun(@(x) x, d);
% K = @(x,y) x.*y;
% A = chebop(@(u) u - volt(K,u), d);
% f = x.^2.*cos(x) + (1-x).*sin(x);
% u = A\f;
%
% pass(5) = norm( u - sin(x) ) < 1e6*tol;
% pass(6) = norm( A*u - f ) < 1e4*tol;

pass = err(:).' < tol;
end