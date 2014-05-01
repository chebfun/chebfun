function f = polyfit(f, n)
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

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]: 
%   Array-valued CHEBFUN objects.
if ( size(f, 2) > 1 )
    error('CHEBFUN:chebtech:polyfit', ...
        'Array-valued CHEBTECH objects are ot supported by POLYFIT().');
end

if ( n <= length(f) )
    
    % Extract Chebyshev coeficients of y:
    c_cheb = f.coeffs;
    % Convert to Legendre coefficients:
    c_leg = chebtech.cheb2leg(c_cheb);
    % Truncate to degree N:s
    c_leg = c_leg((end-n+1):end);
    % Convert to Chebyshev coefficients:
    c_cheb = chebtech.leg2cheb(c_leg);
    f.coeffs = c_cheb;
    % Update vscale of the CHEBTECH:
    f.vscale = getvscl(f);
    
end
