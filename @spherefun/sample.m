function varargout = sample( f, varargin )
%SAMPLE      Samples f on a tensor product grid.
%   X = SAMPLE(F) returns the matrix of values of F on a tensor
%   product grid.
%
%   [U, D, V] = SAMPLE(F) returns the low rank representation of the
%   values of F on a tensor product grid. X = U * D * V'.
%
%   [U, D, V] = SAMPLE(F,M,N) returns the values of F on a M-by-N
%   tensor product grid.
%

% Empty check. 
if ( isempty( f ) )
    varargout = { [] }; 
    return
end

if ( nargin == 1 ) 
    % Get degrees:
    [m, n] = length( f );  
elseif ( nargin == 2 ) 
    error('SPHEREFUN:sample:inputs', 'Dimension not specified.'); 
else
    m = varargin{ 1 }; 
    n = varargin{ 2 }; 
end

% Get the low rank representation for f. 
[cols, d, rows] = cdr(f);

% I really wish techs had a "sample" function too that would return m values
% of the tech at it's natural grid.  This has been on the chebfun tracker
% for sometime now.

% Ugly!
C = trigtech.coeffs2vals(trigtech.alias( cols.funs{:}.onefun.coeffs, 2*n ));
C = C([n+1:2*n 1],:);  % Remove doubled up points.
R = trigtech.coeffs2vals(trigtech.alias( rows.funs{:}.onefun.coeffs, m )); 

% More ugliness
if all(cols.funs{:}.onefun.isReal) && all(rows.funs{:}.onefun.isReal)
    C = real(C);
    R = real(R);
end

% Evaluate: 
if ( nargout <= 1 )
    varargout = {C * d * R.'}; 
else
    varargout = {C , d, R}; 
end
    
end
