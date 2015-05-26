function pass = test_eigs_schrodinger(pref)

if ( nargin == 0 )
    pref = cheboppref();
end

% Double-well Schrodinger eigenstates
% Nick Trefethen, November 2010

% A well-known problem in quantum mechanics is the calculation of
% eigenstates of a potential with the shape of a 'double well'.
% Specifically, consider a potential function V(x) defined on [-1,1] by
%
%    V(x) = 1.5 for x in [-.2,.3],  0 otherwise.
%
% We seek eigenmodes of the steady-state ShrÃ¶dinger equation associated
% with this potential, specifically, functions u(x) satisfying
%
%     -0.007u"(x) + V(x)*u(x) = lam*u(x),    u(-1) = u(1) = 0.
%
% for some constant lam. 

x = chebfun(@(x) x);
V = 1.5*(abs(x-0.05) < 0.25);
L = chebop(@(x,u) -0.007*diff(u,2) + V.*u, [-1, 1], 0);
e = eigs(L, 12, pref);

% Results from computation in Chebfun V4:
V4 = [  0.091480998228413
        0.116757122005510
        0.363909308597386
        0.463167687390724
        0.808941736699187
        1.021145960786674
        1.390812031498242
        1.652575851343383
        1.871230031210140
        2.174488704534987
        2.533176594994942
        2.924094539795185];
    
err = norm(V4 - e, inf);
pass = err < 1e-10;

end
