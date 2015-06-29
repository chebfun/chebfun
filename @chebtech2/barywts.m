function v = barywts(n)
%BARYWTS   Barycentric weights for Chebyshev points of 2nd kind.
%   BARYWTS(N) returns the N barycentric weights for polynomial interpolation on
%   a Chebyshev grid of the 2nd kind. The weights are normalised so that they
%   have infinity norm equal to 1 and the final entry is positive.
%
% See also BARY, CHEBPTS.   

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See Thm. 5.2 of Trefethen, Approximation Theory and Approximation Practice, 
% SIAM, 2013 for more information.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( n == 0 )                      % Special case (no points)
    v = [];
elseif ( n == 1 )                  % Special case (single point)
    v = 1;
else                               % General case
    v = [ones(n-1,1) ; .5];        % Note v(end) is positive.
    v(end-1:-2:1) = -1; 
    v(1) = .5*v(1);
end

end
