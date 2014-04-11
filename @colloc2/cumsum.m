function C = cumsum(disc, m)
%CUMSUM    Indefinite integration operator for COLLOC2 discretization.
%   C = CUMSUM(DISC) gives the matrix such that if v=C*u, then u=v' and v=0
%   at the left endpoint, as accurately as possible in Chebyshev polynomial
%   discretization.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Store information about domain and dimensions.
d = disc.domain;
n = disc.dimension;

if m == 0
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
% CUMSUMMAT  Chebyshev integration matrix.
% Q = CUMSUMMAT(N) is the matrix that maps function values at N Chebyshev
% points to values of the integral of the interpolating polynomial at
% those points, with the convention that the first value is zero.

% TODO: This method is also used in the methods JACPTS and LEGPTS. It should
% probably be made a static method of CHEBTECH (or its subclasses). Thus, it has
% not been reviewed as a part of the LINOP code review. AB, 30/01/14.


% TODO: More efficient implementation?

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
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
T = chebtech2.coeffs2vals(eye(N+1));

% Matrix mapping values -> coeffs.
Tinv = chebtech2.vals2coeffs(eye(N+1));

% Matrix mapping coeffs -> integral coeffs. Note that the highest order
% term is truncated.
k = 1:N;
k2 = 2*(k-1);  k2(1) = 1;  % avoid divide by zero
B = diag(1./(2*k),-1) - diag(1./k2,1);
v = ones(N,1); v(2:2:end) = -1;
B(1,:) = sum( diag(v)*B(2:N+1,:), 1 );
B(:,1) = 2*B(:,1);

Q = T*B(end:-1:1,end:-1:1)*Tinv;
Q(1,:) = 0;  % make exact
cache{N} = Q;

end
