close all
clc
clear classes


%% Building blocks
dom = [-2 2];
I = linop.eye(dom);
D = linop.diff(dom);
Z = linop.zeros(dom);
x = chebfun('x', dom);
c = sin(x.^2);
C = linop.mult(c);   
E = linop.eval(dom);
El = E(dom(1));
Er = E(dom(end));


%% Solve a linear system 
L = linop([ D^2, -I, sin(x); C, D, 0*x; linop.zero(dom), El, 4 ]);
f = [abs(x-1); 0*x; 1 ];
B1 = [El, -Er, 0];
B2 = [linop.sum(dom), El, 0];
B3 = [Er*D, linop.zero(dom), 0];
B4 = [Er,-El,2];
L = addbc(L,B1,0);
L = addbc(L,B2,1);
L = addbc(L,B3,0);
%L = addbc(L,B4,0);

%%

type = {@blockColloc2, @blockUS};
w = [];
for k = 1:2
    figure
    wold = w;
    % w = L\f;
    w = linsolve(L, f, type{k});

    %%
    plot(w{1},'b'); hold on
    plot(w{2},'r'); hold off, shg
    w3 = w{3};

    %%
    % check the ODEs
    tol = 1e-10; 
    err(k,1) = norm( diff(w{1},2)-w{2}+sin(x)*w{3} - f{1} );
    pass(k,1) = err(k,1) < tol;
    err(k,2) = norm( c.*w{1} + diff(w{2}) + 0 - f{2} );
    pass(k,2) = err(k,2) < tol;
    err(k,3) = abs( 0 + feval(w{2},dom(1)) + 4*w{3} - f{3} );
    pass(k,3) = err(k,3) < tol;

    %%
    % check the BCs
    v = w{2};  u = w{1};
    err(k,4) = abs( u(-2)-v(2) );
    pass(k,4) = err(k,4) < tol;
    err(k,5) = abs( sum(u)+v(-2) - 1);
    pass(k,5) = err(k,5) < tol;
    err(k,6) = abs( feval(diff(u),dom(end)) );
    pass(k,6) = err(k,6) < tol;
end

pass
err
