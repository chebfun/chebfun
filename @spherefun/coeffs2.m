function varargout = coeffs2(f, m, n)
%COEFFS2   Double Fourier coefficients of a SPHEREFUN. 
% 
%   X = COEFFS2(F) returns the 2D Fourier modes of the spherefun, viewed
%   as a doubly periodic function. 
% 
%   [C, D, R] = COEFFS2(F) returns a low rank approximation to the 2D
%   Fourier modes.
% 
%   X = COEFFS2(F, M, N) returns bivariate coefficients with N Fourier 
%   modes in latitude and M Fourier modes in longitude. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% Calculate the CDR decomposition: 
[C, D, R] = cdr(f); 

if ( nargin == 1 )
    % Find the  coefficients of each slice: 
    U = C.coeffs; 
    R = R.coeffs;
else
    if ( nargin == 2 )
        n = m;
    end
    % Find the coefficients of each slice:
    U = trigtech.alias(C.coeffs, n); 
    R = trigtech.alias(R.coeffs, m);
end

% Prepare the output. Keep in low rank form if nargin > 1.
if ( nargout <= 1 ) 
    varargout = { U*D*R.' };
elseif ( nargout <= 3 )
    varargout = { U, D, R };
else
    error('SPHEREFUN:COEFFS:NARGOUT',...
            'Too many output arguments')
end

end