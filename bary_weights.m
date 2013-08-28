function w = bary_weights(x)
%BARY_WEIGHTS   Barycentric weights
%   W = BARY_WEIGHTS(X) returns scaled barycentric weights for the points X. The
%   weights are scaled such that norm(W, inf) == 1.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

n = length(x);
if ( isreal(x) )
    C = 4/(max(x) - min(x)); % Capacity of interval
else
    C = 1;                   % Scaling by capacity doesn't apply for complex nodes
end

% [TODO]: Why is the top loop not used? (IF 0)

if ( n < 2001 && 0 )         % For small n using matrices is faster.
   V = C*bsxfun(@minus, x, x.');
   V(logical(eye(n))) = 1;
   VV = exp(sum(log(abs(V))));
   w = 1./(prod(sign(V)).*VV).';
   
else                         % For large n use a loop
   w = ones(n,1);
   for j = 1:n
       v = C*(x(j)-x); v(j) = 1;
       vv = exp(sum(log(abs(v))));
       w(j) = 1./(prod(sign(v))*vv);
   end
end

% Scaling
w = w./max(abs(w));

end
