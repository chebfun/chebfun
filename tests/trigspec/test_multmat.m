function pass = test_multmat(pref)

if ( nargin == 0 ) 
    pref = trigtech.techPref();
end

tol = 10*pref.chebfuneps;

%% CHEBFUN INPUT:

% INPUT: CHEBFUN OF ODD LENGTH M.
% OUTPUT: MATRIX OF ODD SIZE N.
m = 11;
f = chebfun(@(x) sin(pi*x), m, 'trig');
n = 31;
C = zeros(n, 1);
R = zeros(n, 1);
C(2) = -1i/2;
R(2) = 1i/2;
Mexact = toeplitz(C, R);
M = trigspec.multmat(n, f);
pass(1) = norm(full(M) - Mexact) < tol;

% INPUT: CHEBFUN OF ODD LENGTH M.
% OUTPUT: MATRIX OF EVEN SIZE N.
m = 11;
f = chebfun(@(x) sin(pi*x), m, 'trig');
n = 32;
C = zeros(n, 1);
R = zeros(n, 1);
C(2) = -1i/2;
R(2) = 1i/2;
Mexact = toeplitz(C, R);
M = trigspec.multmat(n, f);
pass(2) = norm(full(M) - Mexact) < tol;

% INPUT: CHEBFUN OF EVEN LENGTH M.
% OUTPUT: MATRIX OF ODD SIZE N.
m = 12;
f = chebfun(@(x) sin(pi*x), m, 'trig');
n = 15;
C = zeros(n, 1);
R = zeros(n, 1);
C(2) = -1i/2;
R(2) = 1i/2;
Mexact = toeplitz(C, R);
M = trigspec.multmat(n, f);
pass(3) = norm(full(M) - Mexact) < tol;

% INPUT: CHEBFUN OF EVEN LENGTH M.
% OUTPUT: MATRIX OF EVEN SIZE N.
m = 12;
f = chebfun(@(x) sin(pi*x), m, 'trig');
n = 16;
C = zeros(n, 1);
R = zeros(n, 1);
C(2) = -1i/2;
R(2) = 1i/2;
Mexact = toeplitz(C, R);
M = trigspec.multmat(n, f);
pass(4) = norm(full(M) - Mexact) < tol;

end
