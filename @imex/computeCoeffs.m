function schemeCoeffs = computeCoeffs(~, dt, Lap, ~, S)
%COMPUTECOEFFS   Compute coefficients of an IMEX scheme.
%   SCHEMECOEFFS = COMPUTECOEFFS(K, DT, LAP, ~, S) computes the coefficients
%   needed by the IMEX K from the time-step DT, the laplacian LAP, and the 
%   SPINOPERATOR S.
%
% See also IMEX.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note: For the moment, we only support the LIRK4 scheme.

% Set-up:
nVars = S.numVars;           % number of unknown functions
N = sqrt(size(Lap, 1)/nVars);  % grid points

% Multiplication matrix Tsin2:
m = N;
n = N;
P = speye(m+1); 
P = P(:, 1:m); 
P(1,1) = .5; 
P(m+1,1) = .5;
Q = speye(m+1+4); 
Q = Q(3:m+2,:); 
Q(1,3) = 1; 
Q(1,m+3) = 1;
Msin2 = toeplitz([1/2, 0, -1/4, zeros(1, m+2)]);
Msin2 = sparse(Msin2(:, 3:m+3));
Tsin2 = round(Q*Msin2*P, 15);
In = speye(n);

% Compute LU factorizations of LIRK4 matrices:
Tsin2 = kron(In, Tsin2);
[L, U] = lu(Tsin2); 
[La, Ua] = lu(Tsin2 - 1/4*dt*Lap);
LU = cell(2, 2);
LU{1, 1} = L;
LU{2, 1} = U;
LU{1, 2} = La;
LU{2, 2} = Ua;
schemeCoeffs.multmat = Tsin2;
schemeCoeffs.LU = LU;
schemeCoeffs.L = Lap;

end