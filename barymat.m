function B = barymat(y, x, w)
%BARYMAT  Barycentric Interpolation Matrix.
%   BARYMAT(Y, X, W), where Y is a column vector of length M and X and W are
%   column and row vectors of length N respectively, returns the M*N matrix
%   which interpolates data from the grid X to the grid Y using the 2nd-kind
%   barycentric interpolation formula with barycentric weights W. If W is not
%   supplied it is assumed to be the weights for polynomial interpolation at a
%   2nd-kind Chebyshev grid: W(j) = (-1)^j, W([1, N]) = 0.5*W([1, N]).

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isempty(y) || isempty(x) )
    % Nothing to do here!
    B = []; 
    return
end

if ( min(size(x)) > 1 || min(size(y)) > 1 )
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
if ( nargin < 3 ) 
    w = ones(N, 1); 
    w(2:2:end) = -1; 
    w([1, N]) = 0.5*w([1, N]);
end
% Make sure everything is orientated correctly:
w = w(:).';

% Construct the matrix
if ( M >= 500 && N >= 1000 )     % <-- Experimentally determined.
    
% Testing shows BSXFUN is faster in this regime    
    B = bsxfun(@minus, y, x');   % Repmat(Y-X')
    B = bsxfun(@rdivide, w', B); % w(k)/(y(j)-x(k))
    c = 1./sum(B, 2);            % Normalisation ('denom' in bary-speak).
    B = bsxfun(@times, B, c);
    
else    
    
% Else use FOR loops
    B = bsxfun(@minus, y, x');   % Repmat(Y-X')
    for k = 1:N
        B(:,k) = w(k)./B(:,k);   % w(k)/(y(j)-x(k))
    end
    c = 1./sum(B,2);             % Normalisation ('denom' in bary-speak).
    for j = 1:M
        B(j,:) = B(j,:)*c(j);
    end
    
end

% Where points coincide there will be division by zeros (as with bary.m). 
% Replace these entries with the identity:
B(isnan(B)) = 1;








