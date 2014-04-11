function E = expm(disc, t)
%EXPM      Operator exponential for COLLOC discretization.
%   This EXPM is called by LINOP.EXPM to perform propagation of a discrete
%   initial condition via matrix exponential. The returned matrix is the
%   propagator for a discretization of the problem.
%
%   See also LINOP.EXPM.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Dimension information
n = disc.dimension;

if ( t == 0 )
    % Trivial case at t = 0
    E = eye(sum(n));
else
    % We need a copy of the matrix discretized without any side conditions.
    [M, P, B, A] = matrix( disc );
    
    % Reduced and original degrees of freedom
    [mRed, mOrig] = size(P);  
    
    % This step implicitly uses the side conditions in order to lift a reduced
    % discretization to full size.
    Q = [B;P] \ [ zeros(mOrig - mRed, mRed); eye(mRed) ];

    % Propagator of the "reduced" variables: Lift to full size, apply original
    % operator, project down to reduced size, then exponentiate. 
    E = expm(t*P*A*Q);
    
    % Propagator for the original variables: Reduce, propagate, lift. 
    E = Q*E*P;

%     E = expm( [B;P] \ [zeros(size(B));P*A]*t );
end

end