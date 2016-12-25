function normA = norm(A, n)
%NORM   Norm of a CHEBMATRIX object.
%   NORM(A) computes the Frobenius norm of the CHEBMATRIX object A, defined as
%   the sum of the squares of the 2-norms of each of the blocks.
%
%   NORM(A, 2) or NORM(A, 'fro') is the same as above.
%
%   NORM(A, INF) return the infinity norm of the CHEBMATRIX A, defined as the
%   maximum infinity norm of each of the blocks.
%
%   NORM(A, -INF) returns the negative infinity norm of the CHEBMATRIX A,
%   defined as the minimum of the negative infinity norm of each of the blocks.
%
% See also CHEBMATRIX, CHEBFUN/NORM.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]: Add support for norms of operators.

% Empty CHEBMATRIX has norm 0.
if ( isempty(A) )
    normA = 0;
    return
end

if ( nargin == 1 )
    n = 'fro'; 	% Frobenius norm is the default.
end

% The norm of a CHEBMATRIX with inf x inf block(s) that are neither 
% CHEBFUN2 nor CHEBFUN3 is not supported:
s = cellfun(@(b) min(size(b)), A.blocks);
t = cellfun(@(b) isa(b, 'chebfun2') || isa(b, 'chebfun3'), A.blocks);
if ( ~all(isfinite(s(:))) )
    if ( ~any(t) )
        error('CHEBFUN:CHEBMATRIX:norm:notSupported', ...
            'Norm of a chebmatrix with operator block(s) is not supported.')
    end
end

% Deal with different cases.
switch n
    
    case {'fro', 2}
        normA = cellfun(@(v) norm(v, 2), A.blocks);
        normA = sqrt(sum(normA(:).^2));
        
    case {inf, 'inf'}
        normA = cellfun(@(v) norm(v, inf), A.blocks);
        normA = max(normA(:));
        
    case {-inf, '-inf'}
        normA = cellfun(@(v) norm(v, -inf), A.blocks);
        normA = min(normA(:));
        
    otherwise
        if ( ~ischar(n) )
            n = num2str(n);
        end
        error('CHEBFUN:CHEBMATRIX:norm:unknown', ...
            'unsupported norm type ''%s''', n)
        
end

end
