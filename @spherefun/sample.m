function varargout = sample( f, varargin )
%SAMPLE      Samples f on a tensor product grid.
%   X = SAMPLE(F) returns the matrix of values of F on an n-by-m tensor
%   product grid, where n is length of columns of F and m is length of the 
%   rows.
%
%   [U, D, V] = SAMPLE(F) returns the low rank representation of the
%   values of F on an n-by-m tensor product grid. X = U * D * V', where
%   n is the length of the columns of F and m is the length of the rows.
%
%   X = SAMPLE(F, M, N) returns the matrix of values of F on an n-by-m 
%   tensor product grid.
%
%   [U, D, V] = SAMPLE(F, M, N) returns the low rank representation of F 
%   on an n-by-m tensor product grid.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% Empty check. 
if ( isempty(f) )
    varargout = {[]}; 
    return
end

if ( nargin == 1 ) 
    % Get degrees:
    [m, n] = length(f);
elseif ( nargin == 2 ) 
    error('CHEBFUN:SPHEREFUN:sample:inputs', 'Dimension not specified.'); 
else
    m = varargin{1}; 
    n = varargin{2}; 
    if ( (m <= 0) || (n <= 0) )
        error('CHEBFUN:SPHEREFUN:sample:inputs', ['Number of sample ' ...
             'points must be positive.']);
    end
end

% Get the low rank representation for f. 
[cols, d, rows] = cdr(f);

% Use CDR decomposition so we can keep it in low rank form: 
C = real( sample(cols, max(2*n-2,1)) );
C = C([n:2*n-2 1], :);  % Remove doubled up points.
R = real( sample(rows, m) );

% Evaluate: 
if ( nargout <= 1 )
    varargout = {C * d * R.'}; 
else
    varargout = {C, d, R}; 
end
    
end