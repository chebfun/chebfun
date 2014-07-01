function normA = norms(A, n)
%NORMS   Norms of elements in a CHEBMATRIX object.
%   NORMS(A) computes the Frobenius norm of each element in the CHEBMATRIX
%   object A and returns a matrix the same size as A.
%
%   NORMS(A, 2) or NORMS(A, 'fro') is the same as above.
%
%   NORMS(A, INF) computes the infinity norm of each element in the CHEBMATRIX A.
%
% See also CHEBMATRIX/NORM, CHEBFUN/NORM.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]: Add support for norms of operators (inf x inf blocks).

% Empty CHEBMATRIX has norm 0.
if ( isempty(A) )
    normA = 0;
    return
end

if ( nargin == 1 )
    n = 'fro'; 	% Frobenius norm is the default.
end

% The norm of a CHEBMATRIX with inf x inf block(s) is not supported.
s = cellfun(@(b) min(size(b)), A.blocks);
if ( ~all(isfinite(s(:))) )
    error('CHEBFUN:CHEBMATRIX:norms:notSupported', ...
    'Norms of a chebmatrix with inf x inf block(s) is not supported.')
end

normA = cellfun(@(v) norm(v, n), A.blocks);

end
