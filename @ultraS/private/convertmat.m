function S = convertmat( n, K1, K2 )
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
%CONVERTMAT(A, K1, K2), convert C^(K1) to C^(K2)
S = speye(n);
for s = K1:K2
    S = spconvert(n, s) * S;
end
end
