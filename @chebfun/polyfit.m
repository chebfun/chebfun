function f = polyfit(y, n)
%POLYFIT   Fit polynomial to a CHEBFUN.
%   F = POLYFIT(Y, N) returns a CHEBFUN F corresponding to the polynomial of
%   degree N that fits the CHEBFUN Y in the least-squares sense.
%
%   If Y is a global polynomial of degree n then this code has an O(n (log n)^2)
%   complexity. If Y is piecewise polynomial then it has an O(n^2) complexity.
%
%   F = POLYFIT(X, Y, N, D), where D is a DOMAIN object, returns a CHEBFUN F on
%   the domain D which corresponds to the polynomial of degree N that fits the
%   data (X, Y) in the least-squares sense. X should be a real-valued column
%   vector and Y should be a matrix with size(Y,1) = size(X,1).
%   
%   F = POLYFIT(Y, N) where Y is represented as a periodic TRIGFUN object
%   returns the degree N trigonometric polynomial fit of length 2N+1.
%
%   Note CHEBFUN/POLYFIT does not not support more than one output argument in
%   the way that MATLAB/POLYFIT does.
%
% See also INTERP1.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~isscalar(n) || round(n) ~= n )
    error('CHEBFUN:CHEBFUN:polyfit:input2', 'N must be scalar integer.')
end
    
if ( any(isinf(y.domain)) )
    error('CHEBFUN:CHEBFUN:polyfit:unbounded', ...
        'Unbounded domains are not supported.');
end

if ( n > length(y) && numel(y.funs) == 1 && isa(y.funs{1}.onefun, 'chebtech') )
    % Nothing to do here!
    f = y;
    return
end

if ( isPeriodicTech(y) )    
    % Get the coefficients of the least square approximation:
    c = truncCoeffs(y.funs{1}, n);
    
    % Construct the fit:
    f = chebfun(c, y.domain, 'coeffs', 'periodic');
    return;
    
end

% Compute first n+1 Legendre coeffs:
cleg = legcoeffs(y, n + 1);

% Convert to Chebyshev coeffs:
c = zeros(size(cleg));
for k = 1:size(c, 2)
    c(:,k) = leg2cheb(cleg(:,k));   
end

% Make a CHEBFUN:
f = chebfun(c, domain(y, 'ends'), 'coeffs');    

% Deal with row CHEBFUNs.
if ( y.isTransposed )    
    f = f.';
end

end
