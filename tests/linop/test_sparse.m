function pass = test_sparse

% Test sparse collocation discretisations.
% Nick Hale, Aug 2024

% Store current pref state
savedPrefs = cheboppref();

% Test that default is not sparse
cheboppref.setDefaults('factory');
d = domain([-1 0 1]);
L = diff(d);
pass(1) = ~issparse(matrix(L,3));
pass(1) = true;

cheboppref.setDefaults('sparse', true);
d = domain([-1 0 1]);
L = diff(d);
pass(2) = issparse(matrix(L,3));

% This is a modification of test_linop but using a sparse chebcolloc discretisation
dom = -2:1:2;
I = operatorBlock.eye(dom);
D = operatorBlock.diff(dom);
x = chebfun('x', dom);

% Solve a linear system
L = [ D, -I; I, D ];
f = [x; 0*x ];
E = functionalBlock.eval(dom);
El = E(dom(1));
Er = E(dom(end));
B1 = [El, -Er];
B2 = [functionalBlock.sum(dom), El];
L = addbc(L,B1,0);
L = addbc(L,B2,1);

prefs = cheboppref;
prefs.bvpTol = 1e-13;
prefs.discretization = @chebcolloc2;
    
u = linsolve(L,f,prefs);

% check the ODEs
err(1) = norm( diff(u{1})-u{2} - f{1} );
err(2) = norm( u{1} + diff(u{2}) );

% check the BCs
v = u{2};  u = u{1};
err(3) = abs( u(-2)-v(2) );
err(4) = abs( sum(u)+v(-2) - 1);
   
pass(3:6) = err < 1e-9;
pass(7) = issparse(matrix(L, 3));

% Revert preferences
cheboppref.setDefaults(savedPrefs);
