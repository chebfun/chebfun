function w = trigBaryWeights(x)
%TRIGBARYWEIGHTS   Barycentric weights for trigonometric interpolation.
%   W = TRIGBARYWTS(X) returns scaled barycentric weights for the points 
%   in the columm vector X. The weights are scaled such that norm(W, inf) == 1.
% 
% REFERENCES:
%   [1] Berrut, Jean-Paul. "Baryzentrische Formeln zur trigonometrischen 
%   Interpolation (I)." Zeitschrift für angewandte Mathematik und Physik 
%   ZAMP 35.1 (1984): 91-105.
% 
%   [2] Henrici, Peter. "Barycentric formulas for interpolating trigonometric 
%   polynomials and their conjugates." Numerische Mathematik 33.2 
%   (1979): 225-234.
%
% See also TRIGBARY()

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.


% Check input dimension:
[n, m] = size(x);
if ( m > 1 )
    if ( n > 1 )
        error('CHEBFUN:trigBaryWts:matrix', 'Input must be a vector.')
    else
        % Allow a row vector:
        n = m;
        x = x.';
    end
end

% Detect equispaced points:
if ( all(abs(diff(x, 2)) < max(abs(x))*eps) )
    w = ones(size(x));
    w(2:2:end) = -1;
    return
end
    
    

% Compute the weights:
if ( (n < 2001) )              % For small n using matrices is faster.
   V = sin(1/2*bsxfun(@minus, x, x.'));
   V(1:n+1:end) = 1;
   VV = exp(sum(log(abs(V))));
   w = 1./(prod(sign(V)).*VV).';
   
else                           % For large n use a loop.
   w = ones(n,1);
   for j = 1:n
       v = sin(1/2*(x(j) - x)); v(j) = 1;
       vv = exp(sum(log(abs(v))));
       w(j) = 1./(prod(sign(v))*vv);
   end
end

% Scaling:
w = w./norm(w, inf);

end