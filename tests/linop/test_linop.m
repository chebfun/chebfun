function pass = test_linop
% TAD, 10 Jan 2014

%% Building blocks
dom = [-2 2];
I = operatorBlock.eye(dom);
D = operatorBlock.diff(dom);
x = chebfun('x', dom);

%% Solve a linear system
L = [ D, -I; I, D ];
f = [x; 0*x ];
E = functionalBlock.eval(dom);
El = E(dom(1));
Er = E(dom(end));
B1 = [El, -Er];
B2 = [functionalBlock.sum(dom), El];
L = addbc(L,B1,0);
L = addbc(L,B2,1);


%%

types = {@chebcolloc2,  @ultraS};
prefs = cheboppref;
prefs.bvpTol = 1e-13;

for k = 1:2
    
    prefs.discretization = types{k};
    
    u = linsolve(L,f,prefs);
    
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
pass = err < 1e-9;
