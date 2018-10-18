function pass = test_svds(pref)

% Get preferences:
if ( nargin < 1 )
    pref = cheboppref();
end
tol = 1e1*pref.bvpTol; 

%% Building blocks
dom = [0 pi];
D = operatorBlock.diff(dom);
x = chebfun('x', dom);
c = x.^2;
C = operatorBlock.mult(c);   
z = functionalBlock.zero(dom);
E = functionalBlock.eval(dom);
El = E(dom(1));
Er = E(dom(end));

%% Singular values of differential operator
s_true = (5:-1:0)';
L = linop(D);
[U,S,V] = svds(L, 6, 'bvp', pref);
s = diag(S);
pass(1) = norm(s - s_true, inf) < max(s)*tol;
pass(2) = norm(L*V-U*S) < max(s)*tol;

%% Singular values of self adjoint operator
s_true = (6:-1:1)'.^2;
L = linop(D^2);
L = addbc(L,El,0);
L = addbc(L,Er,0);
[U,S,V] = svds(L, 6, 'bvp', pref);
s = diag(S);
pass(3) = norm(s - s_true, inf) < max(s)*tol;
pass(4) = norm(L*V-U*S) < max(s)*tol;

end
