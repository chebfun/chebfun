function [M, P, B, A, PS] = matrix(disc, dim, domain)
%MATRIX    Convert operator to matrix using COLLOC discretization.
%   MATRIX(DISC) uses the parameters in DISC to discretize DISC.source as a
%   matrix using COLLOC. 
%
%   MATRIX(DISC, DIM, DOMAIN) overrides the native 'dimension' and 'domain'
%   properties in DISC.
%
%   [PA, P, B, A] = MATRIX(...) returns the additional component matrices
%   resulting from boundary condition manipulations, as described in
%   APPLYCONSTRAINTS.
%
%   See also: INSTANTIATE

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Parse inputs
if ( nargin > 1 )
    disc.dimension = dim;
end
if ( nargin > 2 )
    disc.domain = domain;
end

% Check subinterval compatibility of domain and dimension.
if ( (length(disc.domain) - 1) ~= length(disc.dimension) )
    error('CHEBFUN:chebDiscretisation:matrix:subIntDim', ...
        'Must specify one dimension value for each subinterval.')
end

if ( isa(disc.source, 'chebmatrix') )
    
    % Construct a square representation of each block individually and
    % store in a cell array.
    [A, S] = instantiate(disc);

    % We want output on different format depending on whether the source L is a
    % LINOP or another object (most likely a CHEBMATRIX):
    if ( isa(disc.source, 'linop') )
        % Project rows down, and record the projection matrix as well.
        [PA, P, PS] = disc.reduce(A, S);
                
        % Get constraints:
        B = getConstraints(disc);

        % This should restore squareness to the final matrix.
        M = [ B ; PA ];
        
    else
        % Everything should be of the same dimension.
        M = cell2mat(A);
        
        if ( nargout > 1 )
            error('matrix of a chebmatrix can only return one output.')
        end

        
    end

else
    
    % The source must be a chebfun, or...?
    % Note, this is called by ultraS for functionalBlocks
    M = instantiate(disc);
    
    % Additional outputs (not typically useful)
    if ( nargout > 1 )
        if ( issparse(M) )
            P = speye(size(M));
        else
            P = eye(size(M));
        end
    end
    if ( nargout > 2 )
        B = [];
    end
    A = disc.source;
    
end

end