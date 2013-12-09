function r = rank(f, tol)
% RANK   Rank of a chebfun2.
%
% RANK(F) produces an estimate of the rank of the approximant F. Note that
% RANK(F)<=LENGTH(F) since 
%
% RANK(F,TOL) is the number of singular values of F greater than TOL/N, 
% where N is the first singular value of F.
%
% See also LENGTH.

if ( isempty( f ) )    % check for an empty chebfun2
    r = [];
    return;
end

if ( nargin == 1 )
    tol = 0;
else
    tol = varargin{1}; 
end

% compute singular values of f. 
s = svd( f ); 

if ( max( s ) == 0  ) %check for zero function. 
    r = 0; 
    return
end

r = find( s/s(1) > tol, 1, 'last'); 

end