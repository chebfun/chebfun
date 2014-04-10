function varargout = matrix(disc, dimension, domain, space)
%MATRIX    Convert operator to matrix using COLLOC discretization.
%   MATRIX(DISC) uses the parameters in DISC to discretize DISC.source as a
%   matrix using COLLOC. 
%
%   MATRIX(DISC, DIM, DOMAIN) overrides the native 'dimension' and 'domain'
%   properties in DISC.
%
%   [PA, P, B, A] = MATRIX(...) returns the additional component matrices
%   resulting from boundary condition manipulations, as described in
%   COLLOC.APPLYCONSTRAINTS.
%
%   See also: COLLOC/INSTANTIATE

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Parse inputs
if ( nargin > 1 )
    disc.dimension = dimension;
    if ( nargin > 2 )
        disc.domain = domain;
    end
end

% Check subinterval compatibility of domain and dimension.
if ( (length(disc.domain) - 1) ~= length(disc.dimension) )
    error('Must specify one dimension value for each subinterval.')
end

L = disc.source;
if ( isa(L, 'chebmatrix') )
    
    if ( nargin < 4 )
        space = max(L.diffOrder, [], 1);
        space = max(space, 0);
    end
    A = instantiate(disc, L.blocks, space);
    
    % We want output on different format depending on whether the source L is a
    % LINOP or another object (most likely a CHEBMATRIX).
    if ( isa(L, 'linop') )
        [out{1:4}] = applyConstraints(disc, A, {space});
    else
        [rows, P] = disc.reduce(A);
        out{1} = cell2mat(rows);
        out{2} = blkdiag(P{:});
        out{3} = [];
        out{4} = A;
%         out{1} = A;
    end
    
    % Parse outputs
    m = max(1, nargout);
    varargout(1:m) = out(1:m);
else
    % The source must be a chebfun, or...?
    % FIXME: Figure out if this is ever called.
    [varargout{1:nargout}] = instantiate(disc,L);
end

end
