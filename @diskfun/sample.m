function varargout = sample(f, varargin)
%SAMPLE   Samples a DISKFUN object on a tensor product grid.
%   X = SAMPLE(F) returns the matrix of values of F(theta, r) on a 
%   Fourier-Chebyshev tensor product grid. (theta, r) are polar coordinates, 
%   with -pi <= theta <= pi and 0 <= r <= 1. 
%
%   [U, D, V] = SAMPLE(F) returns the low rank representation of the
%   values of F on a tensor product grid where X = U * D * V'.
%
%   [U, D, V] = SAMPLE(F,M,N) returns the values of F on an M-by-N
%   tensor product grid.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check. 
if ( isempty(f) )
    varargout = { [] };
    return
end

if ( nargin == 1 ) 
    % Get degrees:
    [m, n] = length(f);
    
elseif ( nargin == 2 ) 
    error('CHEBFUN:DISKFUN:sample:inputs', 'Dimension not specified.'); 
    
else
    m = varargin{ 1 }; 
    n = varargin{ 2 }; 
    if ( (m <= 0) || (n <= 0) )
        error('CHEBFUN:DISKFUN:sample:inputs', ['Number of sample ' ...
             'points must be positive.']);
    end
end

% Get the low rank representation for f. 
[cols, d, rows] = cdr(f);
 
C = sample(cols, max(2*n-1, 1));
C = C(n:end, :);

R = real( sample(rows, m)); 

% Evaluate: 
if ( nargout <= 1 )
    varargout = {C * d * R.'}; 
else
    varargout = {C , d, R}; 
end
    
end