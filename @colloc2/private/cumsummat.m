function Q = cumsummat(N)
% CUMSUMMAT  Chebyshev integration matrix.
% Q = CUMSUMMAT(N) is the matrix that maps function values at N Chebyshev
% points to values of the integral of the interpolating polynomial at
% those points, with the convention that the first value is zero.

% TODO: This method is also used in the methods JACPTS and LEGPTS. It should
% probably be made a static method of CHEBTECH (or its subclasses). Thus, it has
% not been reviewed as a part of the LINOP code review. AB, 30/01/14.


% TODO: More efficient implementation?

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.


N = N-1;

persistent cache    % stores computed values for fast return
if isempty(cache), cache = {}; end    % first call

if length(cache) >= N && ~isempty(cache{N})
    Q = cache{N};
    return
else
    cache{N} = [];
end

% Matrix mapping coeffs -> values.
T = cp2cdm(N);

% Matrix mapping values -> coeffs.
Tinv = cd2cpm(N);

% Matrix mapping coeffs -> integral coeffs. Note that the highest order
% term is truncated.
k = 1:N;
k2 = 2*(k-1);  k2(1) = 1;  % avoid divide by zero
B = diag(1./(2*k),-1) - diag(1./k2,1);
v = ones(N,1); v(2:2:end) = -1;
B(1,:) = sum( diag(v)*B(2:N+1,:), 1 );
B(:,1) = 2*B(:,1);

Q = T*B*Tinv;
Q(1,:) = 0;  % make exact
cache{N} = Q;

end

function T = cp2cdm(N)
% Values of Cheb. polys at Cheb nodes, x(n)=-cos(pi*n/N).
theta = pi*(N:-1:0)'/N;
T = cos( theta*(0:N) );
end

function C = cd2cpm(N)
% Three steps: Double the data around the circle, apply the DFT matrix,
% and then take half the result with 0.5 factor at the ends.
theta = (pi/N)*(0:2*N-1)';
F = exp( -1i*theta*(0:2*N-1) );  % DFT matrix
rows = 1:N+1;  % output upper half only
% Impose symmetries on data and coeffs.
C = real( [ F(rows,N+1) F(rows,N:-1:2)+F(rows,N+2:2*N) F(rows,1) ] );
C = C/N;  C([1 N+1],:) = 0.5*C([1 N+1],:);
end
