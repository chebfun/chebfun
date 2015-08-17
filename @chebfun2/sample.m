function varargout = sample( f, varargin )
%SAMPLE      Values of f on a tensor product grid.
%   X = SAMPLE(F) returns the matrix of values of F on a tensor
%   product grid.
%
%   [U, D, V] = SAMPLE(F) returns the low rank representation of the
%   values of F on a tensor product grid. X = U * D * V'.
%
%   [U, D, V] = SAMPLE(F,M,N) returns the values of F on a M-by-N
%   tensor product grid.
%
% See also CHEBCOEFFS2, PLOTCOEFFS2. 

% Empty check. 
if ( isempty( f ) )
    varargout = { [] }; 
    return
end

% For now just call chebpolyval2.  In the future, this function should call
% an appropriate sample function on the underlying tech type.  This would
% require adding "sample" functions at the tech level.
[C,d,R] = chebpolyval2( f, varargin{:} );

% Evaluate: 
if ( nargout <= 1 )
    varargout = {C * d * R.'}; 
else
    varargout = {C , d, R}; 
end

end