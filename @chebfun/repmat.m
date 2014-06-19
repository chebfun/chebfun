function G = repmat(f, M, N)
%REPMAT   Replicate and tile a CHEBFUN.
%   REPMAT(F, M, N) or REPMAT(F, [M, N]) creates an array-valued CHEBFUN by
%   tiling copies of F. If F is a column CHEBFUN, then REPMAT(F, 1, N) returns
%   an array-valued CHEBFUN with N*SIZE(F, 2) CHEBFUN columns. If F is a row
%   CHEBFUN, REPMAT(F, M, 1) returns an array-valued CHEBFUN with M*size(F, 1).
%
% See also HORZCAT, VERTCAT, CAT.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs:
if ( nargin == 2 )
    if ( length(M) ~= 2 )
        error('CHEBFUN:CHEBFUN:repmat:notEnoughInputs', ...
            'Requires REPMAT(F, M, N) or REPMAT(F, [M, N]).')
    end
    N = M(2);  
    M = M(1);
end

if ( f(1).isTransposed )
    % REPMAT a row CHEBFUN:
    if ( N ~= 1 )
        error('CHEBFUN:CHEBFUN:repmat:row',...
            'Use REPMAT(F, M, 1) to replicate and tile row CHEBFUN objects.')
    else
        G = repmat({f}, M, 1);
        G = vertcat(G{:});
    end
else
    % REPMAT a column CHEBFUN:
    if ( M ~= 1 )
        error('CHEBFUN:CHEBFUN:repmat:col',...
            'Use REPMAT(F, 1, N) to replicate and tile column CHEBFUN objects.')
    else
        G = repmat({f}, 1, N);
        G = horzcat(G{:});
    end
end
