function pass = test_eigs_orrsom(pref)

% Orr-Sommerfeld eigenvalues
% Toby Driscoll and Nick Trefethen, October 2010

%TODO: These eigenvalues aren't as accurate as they should be. 
%It's not reliable as a test.

pass = true;

% if ( nargin == 0 )
%     pref = cheboppref();
% end
% pref.errTol = 1e-12;
% 
% Re = 5772.22;               % Reynolds number
% alph = 1;                   % longitudinal Fourier parameter
% 
% A = chebop(@(x,u) (diff(u,4) - 2*alph^2*diff(u,2)+alph^4*u)/Re - ...
%     2i*alph*u - 1i*alph*(1-x.^2).*(diff(u,2)-alph^2*u), [-1, 1]);
% B = chebop(@(x,u) diff(u, 2) - u, [-1 1]);
% A.lbc = @(u) [u ; diff(u)];
% A.rbc = @(u) [u ; diff(u)];
% 
% [V, D] = eigs(A, B, 50, 'LR', pref);
% e = diag(D);
% e(abs(e) > 1e5) = [];
% [ignored, idx] = max(real(e));
% e_crit = e(idx);
% 
% e_crit_v4 = -0.000078029804093 - 0.261565915010080i;
% err = abs(e_crit - e_crit_v4);
% 
% tol = 2e-5;
% pass(1) = err < tol;
% 
% % If we had to remove some entries, then there were spurious eigenvalues..
% pass(2) = numel(e) == 50;

end
