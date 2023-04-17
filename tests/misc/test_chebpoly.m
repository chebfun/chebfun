% Test file for CHEBPOLY.

function pass = test_chebpoly(pref)

% if ( nargin == 0 )
%     pref = chebfunpref();
% end

tol = 1e-14;

n = [0 1 2 3 5 8];          % Test degrees.
x = linspace(-1, 1, 10)';   % Test points.
x = x([2,end-1]);           % Ignore +/-1 for convenience.

% First kind;
T1 = chebpoly(n, 1);
T1x = feval(T1, x);
T2 = cos(acos(x)*n);
err = T2 - T1x;
pass(1) = norm(err(:), inf) < tol;

% Second kind;
U1 = chebpoly(n, 2);
U1x = feval(U1, x);
U2 = sin(acos(x)*(n+1))./sin(acos(x)*(1+0*n));
err = U2 - U1x;
pass(2) = norm(err(:), inf) < tol;

% Check a couple of mins and maxes
T = chebpoly(27);
pass(3) = abs(max(T)-1) < tol;

% Check a couple of mins and maxes
T = chebpoly(44);
pass(4) = abs(min(T)+1) < tol;

% And a little quasimatrix
T = chebpoly([7:8]);
v = T*[1;1];
pass(5) = abs(v(-1)) < tol;

% Make sure big degrees work
T = chebpoly(1e5);
pass(6) = abs(T(1)-1) < tol;

end
