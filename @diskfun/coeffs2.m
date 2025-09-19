function varargout = coeffs2( f, m, n ) 
%COEFFS2   Fourier--Chebyshev coefficients of a DISKFUN. 
% 
%   X = COEFFS2( F ) returns the coefficients of the DISKFUN F in the
%   Fourier--Chebyshev bases, where the columns are the Chebyshev
%   coefficients and the rows are Fourier coefficients.
% 
%   [C, D, R] = COEFFS2( F ) returns a low rank approximation to the
%   coefficients. 
%
%   X = COEFFS2(F, M, N) returns bivariate coefficients with N Chebyshev 
%   coefficients in the radial direction and M Fourier coefficients in
%   angular direction.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.



%empty check
if isempty(f)
    varargout = {}; 
    return
end

% Calculate the CDR decomposition: 
[C, D, R] = cdr( f ); 

if nargin == 1
    % Find the coefficients of each slice: 
    U = C.coeffs; 
    R = R.coeffs;
else
    if nargin == 2
        n = m;
    end
    % Find the coefficients of each slice: 
    U = chebtech2.alias(C.coeffs,n); 
    R = trigtech.alias(R.coeffs,m);
end

% Prepare the output. Keep in low rank form if nargin > 1.
if ( nargout <= 1 ) 
    varargout = {U*D*R.'};
elseif ( nargout <= 3 )
    varargout = {U, D, R};
else
    error('DISKFUN:COEFFS:NARGOUT',...
            'Too many output arguments')
end

end