function E = expm(disc,t)

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
d = disc.domain;
n = disc.dimension;

if ( t == 0 )
    E = eye(sum(n));
else
    % We need a copy of the matrix discretized without any side conditions.
    [~,P,B,A] = matrix( disc );
    
    % Assemble the "prolongation" operator.
    [mRed,mOrig] = size(P);
    Q = [B;P] \ [ zeros(mOrig-mRed,mRed); eye(mRed) ];

    % Propagator of the "reduced" variables.
    E = expm(t*P*A*Q);
    
    % Propagator for the original variables. 
    E = Q*E*P;
end

end