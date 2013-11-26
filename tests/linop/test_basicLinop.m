function pass = test_basicLinop

%% Building blocks
dom = [-2 2];
I = linop.eye(dom);
D = linop.diff(dom);
Z = linop.zeros(dom);
x = chebfun('x', dom);
u = sin(x.^2);
U = linop.mult(u);

%% Solve a linear system
L = linop([ D, -I; I, D ]);
f = [x; 0*x ];
E = linop.eval(dom);
El = E(dom(1));
Er = E(dom(end));
B1 = [El, -Er];
B2 = [linop.sum(dom), El];
L = addbc(L,B1,0);
L = addbc(L,B2,1);


%%

types = {@colloc2,  @ultraS};

for k = 1:2
    
    u = mldivide(L, f, types{k});
    
    %
%     plot(u{1},'b'); hold on
%     plot(u{2},'r'); hold off, shg
    
    %
    % check the ODEs
    err(k, 1) = norm( diff(u{1})-u{2} - f{1} );
    err(k, 2) = norm( u{1} + diff(u{2}) );
    
    %
    % check the BCs
    v = u{2};  u = u{1};
    err(k,3) = abs( u(-2)-v(2) );
    err(k,4) = abs( sum(u)+v(-2) - 1);
    
end

err;
pass = err < 1e-12;
