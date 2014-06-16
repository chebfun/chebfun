function H = discriminant(f, x, y, varargin)
%DISCRIMINANT the determinant of Hessian of a CHEBFUN2 at (x,y) 
%   H = DISCRIMINANT(F,x,y) returns the determinant of the Hessian of F at
%   (x,y). The gradient of F should be zero at (x,y).
% 
%   H = DISCRIMINANT(F,G,x,y) returnes the determinant of the 'border' Hessian
%   of F at (x,y).
%
%   Note that we cannot represent the Hessian matrix because we do not allow
%   horizontal concatenation of CHEBFUN2 objects.
%
% See also JACOBIAN. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty( f ) )
    H = [];
    return
end

% Mixed second partial derivatives:
fxx = diff(f, 2, 2); 
fyy = diff(f, 2, 1); 
fxy = diff(f, [1 1]);
    
if ( nargin == 3 )      % Standard Hessian
    % Evaluate at (x,y):
    fxx = feval(fxx, x, y);
    fxy = feval(fxy, x, y);
    fyy = feval(fyy, x, y);
    H = fxx.*fyy - fxy.^2;
    
elseif ( nargin == 4 )  % Bordered Hessian
    % Parse user inputs:
    g = x; 
    x = y; 
    y = varargin{1}; 
    
    % Evaluate at (x,y):
    fxx = feval(fxx, x, y);
    fyy = feval(fxy, x, y);
    fxy = feval(fyy, x, y);
    
    % Partial diff and evaluate:
    gx = diff(g, 1, 2); 
    gy = diff(g, 1, 1);
    gx = feval(gx, x, y); 
    gy = feval(gy, x, y);
    
    % Determinant of bordered Hessian.
    H = -gx.*(gx.*fyy - gy.*fxy) + gy.*(gx.*fxy - gy.*fxx);
    
else
    error('CHEBFUN:CHEBFUN2:discriminant:badInput', 'Invalid input arguments.');
end

end
    
