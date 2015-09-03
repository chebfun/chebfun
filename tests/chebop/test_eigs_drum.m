function pass = test_eigs_drum(pref)

% Frequencies of a drum
% Toby Driscoll, November 2010

if ( nargin == 0 )
    pref = cheboppref();
end

% The axisymmetric harmonic vibrations of a circular drum can be described by
% the ODE
%
%   u"(r) + (1/r) u'(r) = -omega^2 u(r),   u'(0)=1, u(1)=0,
%
% where r is the radial coordinate and omega is the frequency of
% vibration. Only discrete positive values of omega are possible, 
% corresponding to the eigenvalues of the differential equation. 

% We multiply the ODE through by r to avoid a potential division by zero.
% This creates a generalized problem in the form A*u = lambda*B*u.
r = chebfun('r', [0, 1]);
A = chebop(0, 1);
A.op = @(r,u) r.*diff(u,2) + diff(u);
A.lbc = 'neumann'; 
A.rbc = 'dirichlet';
B = chebop(0, 1);
B.op = @(r,u) r.*u;

% Then we find the eigenvalues with eigs. It turns out that the omega
% values are also zeros of the Bessel function J0, which gives a way to
% validate the results.
[V, D] = eigs(A, B, pref);
omega = sort(sqrt(-diag(D)));
err(1) = norm(omega - roots( besselj(0, chebfun('r',[0 20]))));
err(2) = 0;
for vCounter = 1:size(V, 2)
    err(2) = err(2) + norm(A(V(:,vCounter)) - ...
        B(V(:,vCounter))*D(vCounter, vCounter)).^2;
end
err(2) = sqrt(err(2));
%%

tol = [1e-9 1e-6];
pass = err < tol;

end
