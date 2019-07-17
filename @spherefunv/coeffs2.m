function varargout = coeffs2( f, varargin ) 
%COEFFS2   componentwise double Fourier coefficients of a SPHEREFUNV. 
% 
%   [X, Y, Z] = COEFFS2( F ) returns the 2D Fourier modes for each component 
%   of the SPHEREFUNV F, where each component is a doubly-periodic
%   function on the double Fourier sphere.
%  
%   [X, Y] = COEFFS2( F, M, N ) returns bivariate coefficients with N Fourier 
%   modes in the latitude direction and M Fourier modes in the
%   longitude direction.
%
%   See also: SPHEREFUN/COEFFS2

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.



%empty check
if isempty(f)
    varargout = {}; 
    return
end

F1 = f.components{1};
F2 = f.components{2}; 
F3 = f.components{3}; 
if isempty(varargin)
    varargout{1} = coeffs2(F1); 
    varargout{2} = coeffs2(F2);
    varargout{3} = coeffs2(F3); 
else
    varargout{1} = coeffs2(F1, varargin{:}); 
    varargout{2} = coeffs2(F2, varargin{:}); 
    varargout{3} = coeffs2(F3, varargin{:});
end
    

end