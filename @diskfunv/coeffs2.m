function varargout = coeffs2( f, varargin) 
%COEFFS2   componentwise Fourier--Chebyshev coefficients of a DISKFUNV. 
% 
%   [X, Y] = COEFFS2( F ) returns the coefficients for each component of 
%   the DISKFUNV F in the Fourier--Chebyshev bases, where the columns 
%   are the Chebyshev coefficients and the rows are Fourier coefficients.
%  
%
%   [X, Y] = COEFFS2(F, M, N) returns bivariate coefficients with N Chebyshev 
%   coefficients in the radial direction and M Fourier coefficients in
%   angular direction.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


%empty check
if isempty(f)
    varargout = {}; 
    return
end

F1 = f.components{1};
F2 = f.components{2}; 
if isempty(varargin)
    varargout{1} = coeffs2(F1); 
    varargout{2} = coeffs2(F2); 
else
    varargout{1} = coeffs2(F1, varargin{:}); 
    varargout{2} = coeffs2(F2, varargin{:}); 
end
    

end