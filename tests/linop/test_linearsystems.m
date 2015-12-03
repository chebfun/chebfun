function pass = test_linearsystems
% TAD, 10 Jan 2014

tol = 1e-8; 


%% Building blocks
dom = [-1 1];
I = operatorBlock.eye(dom);
D = operatorBlock.diff(dom);
Z = operatorBlock.zeros(dom);
x = chebfun('x', dom);
c = sin(x.^2);
C = operatorBlock.mult(c);   
E = functionalBlock.eval(dom);
El = E(dom(1));
Er = E(dom(end));

%% Solve a linear system 
L = [ D^2, -I, sin(x); C, D, chebfun(0,dom); functionalBlock.zero(dom), El, 4 ];
f = [(x-1); chebfun(0,dom); 1 ];
B1 = [El, -Er, 0];
B2 = [functionalBlock.sum(dom), El, 0];
B3 = [Er*D, functionalBlock.zero(dom), 0];
B4 = [Er,-El,2];
L = addbc(L,B1,0);
L = addbc(L,B2,1);
L = addbc(L,B3,0);
% L = addbc(L,B4,0);

% %% Solve a linear system 
% L = [ D^2, -I; C D] ;
% f = [(x-1); chebfun(0,dom)];
% B1 = [El, -Er];
% B2 = [functionalBlock.sum(dom), El];
% B3 = [Er*D, functionalBlock.zero(dom)];
% B4 = [Er,-El,2];
% L = addbc(L,B1,0);
% L = addbc(L,B2,1);
% %L = addbc(L,B4,0);

%%

type = {@chebcolloc2, @ultraS, @chebcolloc1, @chebcolloc2, @ultraS, @chebcolloc1};
prefs = cheboppref;
prefs.bvpTol = 1e-13;

for k = 1:6

    prefs.discretization = type{k};
    w = linsolve(L, f, prefs);

    %%
    % check the ODEs
    w1 = w{1};  w2 = w{2};  w3 = w{3};
    f1 = f{1};  f2 = f{2}; f3 = f{3};
    residual1 = diff(w1, 2) - w2 + sin(x).*w3 - f1;
    err(k,1) = norm( residual1 );
    residual2 = c.*w1 + diff(w2) + 0 - f2;
    err(k,2) = norm( residual2 );
    err(k,3) = abs( 0 + feval(w2,dom(1)) + 4*w3 - f3 );

    %%
    % check the BCs
    v = w{2};  u = w{1};
    err(k,4) = abs( u(dom(1)) - v(dom(2)) );
    err(k,5) = abs( sum(u) + v(dom(1)) - 1);
    err(k,6) = abs( feval(diff(u), dom(end)) );
    
    %%
    % check continuity
    Du = D*u;  Dv=D*v;
    err(k,7) = feval(u,1,'left') - feval(u,1,'right');
    err(k,8) = feval(v,1,'left') - feval(v,1,'right');
    err(k,9) = feval(Du,1,'left') - feval(Du,1,'right');
    
    if ( k == 2 )
        f = [abs(x-1); 0*x; 1 ];
    end
end

pass = err < tol;
