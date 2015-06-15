function B = barymat(y, x, w, s, r, doFlip)
%BARYMAT  Barycentric Interpolation Matrix.
%   BARYMAT(Y, X, W), where Y is a column vector of length M and X and W are
%   column and row vectors of length N, respectively, returns the M*N matrix
%   which interpolates data from the grid X to the grid Y using the 2nd-kind
%   barycentric interpolation formula with barycentric weights W. If W is not
%   supplied it is assumed to be the weights for polynomial interpolation at a
%   2nd-kind Chebyshev grid: W(j) = (-1)^j, W([1, N]) = 0.5*W([1, N]).
%
%   BARYMAT(Y, X, W, S, R) is the same, where S = acos(Y) and R = acos(X). The
%   purpose of this is that Y(j) - X(k) can be more accurately comuted in this
%   'theta space'. This is sometimes referred to as the 'trig trick' in spectral
%   collocation. BARYMAT(Y, X, W, S, R, 1) also performs the 'flipping trick', 
%   which takes advantage of the fact that the smaller entries in R and S can be
%   computed more accurately. Note that X and Y should be symmetric about zero
%   for this work, and it is assumed that S and R are sorted in descending
%   order.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(y) || isempty(x) )
    % Nothing to do here!
    B = []; 
    return
end

if ( size(x, 2) > 1 || size(y, 2) > 1 )
    error('CHEBFUN:barymat:dims', 'Inputs should be column vectors.');
end

% Obtain lengths of inputs:
N = length(x);
M = length(y);

% Nothing to do here!
if ( M == N && all(x == y) )     
    B = eye(N); 
    return
end
    
% Default to the Chebyshev barycentric weights:
if ( nargin < 3 || isempty(w) ) 
    w = ones(1, N); 
    w(2:2:end) = -1; 
    w([1, N]) = 0.5*w([1, N]);
else
    % Ensure w is a row vector:
    w = reshape(w, 1, N);
end

% Repmat(Y-X'):
if ( nargin < 5 )
    B = bsxfun(@minus, y, x.');  
else
    % Use the 'trig trick' that y-x = cos(s)-cos(r) = 2*sin((s+r)/2)*sin((r-s)/2).
    B = 2*bsxfun(@(r, s) sin((s+r)/2).* sin((r-s)/2), r.', s);
end

% Construct the matrix:
if ( M >= 500 && N >= 1000 )            % <-- Experimentally determined.
    % Testing shows BSXFUN is faster in this regime
    B = bsxfun(@rdivide, w, B);         % w(k)/(y(j)-x(k))
    B = bsxfun(@rdivide, B, sum(B, 2)); % Normalisation.
else
    % Else use FOR loops
    for k = 1:N
        B(:,k) = w(k)./B(:,k);          % w(k)/(y(j)-x(k))
    end
    c = 1./sum(B, 2);                   % Normalisation.
    for j = 1:M
        B(j,:) = B(j,:)*c(j);
    end
end

% Where points coincide there will be division by zeros (as with bary.m).
% Replace these entries with the identity:
B(isnan(B)) = 1;

% Flipping trick:
if ( nargin > 5 && doFlip )
    ii = logical(rot90(tril(ones(M, N)), 2)); 
    ii = fliplr(ii);
    rot90D = rot90(B, 2);
    B(ii) = rot90D(ii);
    B(isnan(B)) = 1;
end

end
