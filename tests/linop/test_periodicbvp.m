function pass = test_periodicbvp
% TAD, 10 Jan 2014

tol = 4e-9; 

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

type = {@colloc2, @colloc1, @ultraS, @colloc2, @colloc1, @ultraS};
prefs = cheboppref;
prefs.errTol = 1e-14;

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

pass = err < tol;

end
