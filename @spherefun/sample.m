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

% I really wish techs had a "sample" function too that would return m values
% of the tech at it's natural grid.  This has been on the chebfun tracker
% for sometime now.

% Ugly!
if ( n > 1 )
    C = trigtech.coeffs2vals(trigtech.alias(cols.funs{:}.onefun.coeffs, ...
        max(2*n-2, 1)));
else
    C = cols.funs{:}.onefun.values;
end
C = C([n:2*n-2 1], :);  % Remove doubled up points.
R = trigtech.coeffs2vals(trigtech.alias(rows.funs{:}.onefun.coeffs, m)); 

% More ugliness
if ( all(cols.funs{:}.onefun.isReal) && all(rows.funs{:}.onefun.isReal) )
    C = real(C);
    R = real(R);
end

% Evaluate: 
if ( nargout <= 1 )
    varargout = {C * d * R.'}; 
else
    varargout = {C, d, R}; 
end
    
end