function pass = test_eigs_orrsom(pref)

% Orr-Sommerfeld eigenvalues
% Toby Driscoll and Nick Trefethen, October 2010

if ( nargin == 0 )
    pref = cheboppref();
end

Re = 5772.22;               % Reynolds number
alph = 1;                   % longitudinal Fourier parameter

A = chebop(@(x,u) (diff(u,4) - 2*alph^2*diff(u,2)+alph^4*u)/Re - ...
    2i*alph*u - 1i*alph*(1-x.^2).*(diff(u,2)-alph^2*u), [-1, 1]);
B = chebop(@(x,u) diff(u, 2) - u, [-1 1]);
A.lbc = @(u) [u ; diff(u)];
A.rbc = @(u) [u ; diff(u)];

[V, D] = eigs(A, B, 50, 'LR', pref);
e = diag(D);
e(abs(e) > 1e5) = [];
[ignored, idx] = max(real(e));
e_crit = e(idx);

e_crit_v4 = -0.000078029804093 - 0.261565915010080i;
err = abs(e_crit - e_crit_v4);

tol = 1e-6;
pass(1) = err < tol;

% If We had to remove some entries, then there were spurious eigenvalues..
pass(2) = numel(e) == 50;

end
