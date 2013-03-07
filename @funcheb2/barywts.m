function v = barywts(n)
%BARYWTS % 2nd-kind Chebyshev barycentric weights.
%   BARYWTS(N) returns the N barycentric weights for polynomial interpolation on
%   a Chebyshev grid of the 2nd kind.
%
% See also FUNCHEB2.BARY.m, FUN2.CHEBPTS.   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See Thm. 5.2 of Trefethen: Approximation Theory and Approximation Practice, 
% SIAM, (2013) for more information.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
if ( n == 0 )                      % Special case (no points!)
    v = [];
elseif ( n == 1 )                  % Special case (single point)
    v = 1;
else                               % General case
%     v = [.5 ; ones(n-1,1)];
%     v(2:2:end) = -1; 
%     v(end) = .5*v(end);
    
    v = ones(n,1);
    v(end-1:-2:1) = -1;
    v(end) = 0.5;
    v(1) = v(1)*0.5;
end

end


