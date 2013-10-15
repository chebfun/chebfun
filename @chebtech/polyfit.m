function f = polyfit(y, n)
%POLYFIT   Fit polynomial to a CHEBTECH.
%   F = POLYFIT(Y, N) returns a CHEBTECH F corresponding to the polynomial of
%   degree N that fits the CHEBTECH Y in the least-squares sense.
%
%   Note CHEBTECH/POLYFIT does not not support more than one output argument in
%   the way that MATLAB/POLYFIT does.
%
%   If y is a global polynomial of degree n then this code has an O(n (log n)^2)
%   complexity. If y is piecewise polynomial then it has an O(n^2) complexity.
%
% See also LEG2CHEB, CHEB2LEG.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = y;

if ( n <= length(y) )
    
    % Extract Chebyshev coeficients of y:
    c_cheb = y.coeffs;
    % Convert to Legendre coefficients:
    c_leg = chebtech2.cheb2leg(c_cheb);
    % Truncate to degree N:s
    c_leg = c_leg((end-n+1):end);
    % Convert to Chebyshev coefficients:
    c_cheb = chebtech2.leg2cheb(c_leg);
    % Compute corresponding values on a Chebyshev grid:
    v_cheb = y.coeffs2vals(c_cheb);
    % Update the values, coefficients, and vscale of the CHEBTECH:
    f.coeffs = c_cheb;
    f.values = v_cheb;
    f.vscale = max(abs(v_cheb), [], 1);
    
end
