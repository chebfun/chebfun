function h = max( f, g, dim )
%MAX   Maximum value of a SEPARABLEAPPROX in one direction.
%   MAX(f) returns a chebfun representing the maximum of the SEPARABLEAPPROX along the
%   y direction, i.e, MAX(f) = @(x) max( f ( x, : ) )
%
%   MAX(f, [], dim) returns a CHEBFUN representing the maximum of f along the
%   DIM direction. If DIM = 1 is along the y-direction and DIM = 2 is along the
%   x-direction.
%
%   WARNING: This function is not always accurate to full machine precision. 
% 
%   For the global maximum use MAX2.
%
% See also MIN, MAX2, MIN2, MINANDMAX2.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( f ) ) 
    error('CHEBFUN:SEPARABLEAPPROX:max:input', 'SEPARABLEAPPROX is empty');
end

% Default to max of one separableApprox:
if ( nargin < 2 )
    g = []; 
end

% Default to maximum along the y direction: 
if ( nargin < 3 )
    dim = 1;
end

% Do not allow max(F, G): 
if ( nargin > 1 && ~isempty( g ) )
    error('CHEBFUN:SEPARABLEAPPROX:max:twoSeparableApproxInputs', ...
        'Unable to maximize two SEPARABLEAPPROX objects.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We have no idea how to achieve this in an efficient way. This
% is an attempt to return a result, but typically it won't be accurate to 
% more than 4-5 digits. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dom = f.domain;
n = 2049; 
if ( dim == 1 )
    vals = sample(f, n, n); 
    h = chebfun( max( vals ).', dom(1:2), 'splitting', 'on' );
    h = simplify( h.', [], 'globaltol' ); 
elseif ( dim == 2 )
    vals = sample(f, n, n);  
    h = chebfun( max( vals, [], 2 ), dom(3:4), 'splitting', 'on' );
    h = simplify( h, [], 'globaltol' );
elseif ( dim == 0 ) 
    error('CHEBFUN:SEPARABLEAPPROX:max:dim', ...
        'Dimension argument must be a positive integer scalar within indexing range.')
else
   % return the SEPARABLEAPPROX. This is analogous to that MAX() command in
   % MATLAB.
   h = f;  
end

end
