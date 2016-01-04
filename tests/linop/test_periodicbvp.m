function pass = test_periodicbvp
% TAD, 10 Jan 2014

tol = 1e-8; 

%% Building blocks
dom = [-pi pi];
[Z, I, D, C] = linop.primitiveOperators(dom);
[z, E, s] = linop.primitiveFunctionals(dom);
x = chebfun('x', dom);
c = sin(x.^2);
C = operatorBlock.mult(c);   
El = E(dom(1));
Er = E(dom(end));

%% Solve a linear system 
L = [ D^2, -I, sin(x); C, D, chebfun(0,dom); functionalBlock.zero(dom), El, 4 ];
f = [sin(x); chebfun(0,dom); 1 ];
L = addbc(L,'periodic');

%%

type = {@chebcolloc2, @chebcolloc1, @ultraS, @chebcolloc2, @chebcolloc1, @ultraS};
prefs = cheboppref;
prefs.bvpTol = 1e-14;

w = [];
for k = 1:6
    prefs.discretization = type{k};
    w = linsolve(L,f,prefs);

    %%
    % check the ODEs
    w1 = w{1};  w2 = w{2};  w3 = w{3};
    f1 = f{1};  f2 = f{2}; f3 = f{3};
    residual1 = diff(w1,2) - w2 + sin(x).*w3 - f1;
    err(k,1) = norm( residual1 );
    residual2 = c.*w1 + diff(w2) + 0 - f2;
    err(k,2) = norm( residual2 );
    err(k,3) = abs( 0 + feval(w2, dom(1)) + 4*w3 - f3 );

    %%
    % check the BCs
    v = w{2};  u = w{1};
    Du = D*u;  Dv = D*v;
    err(k,4) = abs( u(pi) - u(-pi) );
    err(k,5) = abs( v(pi) - v(-pi) );
    err(k,6) = abs( Du(pi) - Du(-pi) );
    
    %%
    % check continuity 
    err(k,7) = feval(u, pi/2, 'left') - feval(u, pi/2, 'right');
    err(k,8) = feval(v, pi/2, 'left') - feval(v, pi/2, 'right');
    err(k,9) = feval(Du, pi/2, 'left') - feval(Du, pi/2, 'right');
    err(k,10) = feval(u, -pi/2, 'left') - feval(u, -pi/2, 'right');
    err(k,11) = feval(v, -pi/2, 'left') - feval(v, -pi/2, 'right');
    err(k,12) = feval(Du, -pi/2, 'left') - feval(Du, -pi/2, 'right');
    
    if ( k == 2 )
        % introduce breakpoints
        f = [abs(cos(x)); 0*x; 1 ];
    end

end

%% Test for TRIGCOLLOC.

% Domain.
dom = [0 2*pi];

% Differentiation operators.    
D = operatorBlock.diff(dom);
D2 = operatorBlock.diff(dom, 2);

% Multiplication operators.
a1 = chebfun(@(x) 1 + sin(2*x), dom);
A1 = operatorBlock.mult(a1);   
a0 = chebfun(@(x) 1 + cos(x), dom);
A0 = operatorBlock.mult(a0);  

% Linop.
L = linop(D2 + A1*D + A0);

% Rhs.
f = chebfun(@(x) cos(x), dom);

% Solve it with TRIGCOLLOC.
prefs.discretization = @trigcolloc;
u = linsolve(L, f, prefs);
u = u{1};

err(1, 13) = norm(diff(u,2) + a1.*diff(u) + a0.*u - f);
err(2, 13) = abs(u(2*pi) - u(0));
err(3, 13) = abs(feval(diff(u), 2*pi) - feval(diff(u), 0));

%%
pass = err < tol;

end
