function pass = test_linearSystems

tol = 1e-9; 

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
f = [(x-1); 0*x; 1 ];
B1 = [El, -Er, 0];
B2 = [linop.sum(dom), El, 0];
B3 = [Er*D, linop.zero(dom), 0];
B4 = [Er,-El,2];
L = addbc(L,B1,0);
L = addbc(L,B2,1);
L = addbc(L,B3,0);
%L = addbc(L,B4,0);

%%

type = {@colloc2, @ultraS, @colloc2, @ultraS};
w = [];
for k = 1:4
    wold = w;
    w = mldivide(L, f, type{k});

    %%
%     subplot(1, 2, k)
%     plot(w{1},'b'); hold on
%     plot(w{2},'r'); hold off, shg
%     w3 = w{3};

    %%
    % check the ODEs
    err(k,1) = norm( diff(w{1},2)-w{2}+sin(x)*w{3} - f{1} );
    err(k,2) = norm( c.*w{1} + diff(w{2}) + 0 - f{2} );
    err(k,3) = abs( 0 + feval(w{2},dom(1)) + 4*w{3} - f{3} );

    %%
    % check the BCs
    v = w{2};  u = w{1};
    err(k,4) = abs( u(-2)-v(2) );
    err(k,5) = abs( sum(u)+v(-2) - 1);
    err(k,6) = abs( feval(diff(u),dom(end)) );
    
    if ( k == 2 )
        f = [abs(x-1); 0*x; 1 ];
    end
end

err;
pass = err < tol;
