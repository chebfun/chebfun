function C = cumsum(disc, m)
%CUMSUM   Indefinite integration operator for COLLOC1 discretization.
%   C = CUMSUM(DISC) gives the matrix such that if v=C*u, then u=v' and v=0
%   at the left endpoint, as accurately as possible in Chebyshev polynomial
%   discretization.
%
%   CUMSUM(DISC, M) for positive integer M returns C^M.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Store information about domain and dimensions.
d = disc.domain;
n = disc.dimension;

if ( m == 0 )
    % Trivial case
    C = eye(sum(n));
else
    numIntervals = disc.numIntervals;
    
    % Find the diagonal blocks.
    blocks = cell(numIntervals);
    for k = 1:numIntervals
        len = d(k+1) - d(k);
        blocks{k} = cumsummat(n(k)) * (len/2);  % Scaled cumsummats.
    end
    
    % Assemble!
    C = blkdiag(blocks{:});
    
    % Each subinterval also contributes to the integrals in all the
    % subintervals to its right, creating a triangular structure.
    offset = 0;
    for k = 1:numIntervals
        % Grab the weights for the integral using all of this
        % subinterval.
        row = offset + n(k);
        cols = offset + (1:n(k));
        last = C(row, cols);
        % Copy it to add to the ones that follow.
        offset = row;
        C(offset + 1:end, cols) = repmat(last,[sum(n) - offset, 1]);
    end
    C = C^m;
end

end

function Q = cumsummat(N)
%CUMSUMMAT  Chebyshev integration matrix.
%   Q = CUMSUMMAT(N) is the matrix that maps function values at N Chebyshev
%   points to values of the integral of the interpolating polynomial at those
%   points, with the convention that the first value is zero.

% TODO: More efficient implementation?
% TODO: This is duplicated in a number of places.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

N = N-1;

persistent CACHE  % Stores computed values for fast return
if ( isempty(CACHE) ), CACHE = {}; end  % First call

if ( length(CACHE) >= N && ~isempty(CACHE{N}) )
    Q = CACHE{N};
    return
else
    CACHE{N} = [];
end

% Matrix mapping coeffs -> values.
T = chebtech1.coeffs2vals(eye(N+1));

% Matrix mapping values -> coeffs.
Tinv = chebtech1.vals2coeffs(eye(N+1));

% Matrix mapping coeffs -> integral coeffs. Note that the highest order term is
% truncated.
k = 1:N;
k2 = 2*(k-1);  k2(1) = 1;  % Avoid divide by zero
B = diag(1./(2*k),-1) - diag(1./k2,1);
v = ones(N,1); v(2:2:end) = -1;
B(1,:) = sum( diag(v)*B(2:N+1,:), 1 );
B(:,1) = 2*B(:,1);

Q = T*B(end:-1:1,end:-1:1)*Tinv;
% Make exact:
Q(1,:) = 0;

%Store:
CACHE{N} = Q;

end