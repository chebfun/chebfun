function g = max( f, g, dim )
%MAX   Maximum value of a CHEBFUN in one direction.
%   MAX(f) returns a chebfun representing the maximum of the CHEBFUN2 along the
%   y direction, i.e, MAX(f) = @(x) max( f ( x, : ) )
%
%   MAX(f, [], dim) returns a CHEFBUN representing the maximum of f along the
%   DIM direction. If DIM = 1 is along the y-direction and DIM = 2 is along the
%   x-direction.
%
%   This function is not considered highly accurate. Expect no more than five
%   or six digits of precision. For the global maximum use MAX2.
%
% See also MIN, MAX2, MIN2, MINANDMAX2.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check: 
if ( isempty(f) ) 
    error('CHEBFUN2:MAX:INPUT', 'CHEBFUN2 is empty');
end

% Do not allow max(F, G): 
if ( nargin > 1 && ~isempty(g) )
    error('CHEBFUN2:max', 'Unable to maximise two CHEBFUN2 objects.');
end

% Default to maximum along the y direction: 
if ( nargin < 1 )  
    dim = 1;
end
dom = f.domain;

% We have no idea how to achieve this in a really fast and efficient way. This
% is an attempt to at least get a result, but there are no guarantees that this
% will be 16 digit accurate.

sample = 2049; 
if ( dim == 1 )
    vals = chebpolyval2(f, sample, sample); 
    g = chebfun( max( vals ).', dom(1:2), 'splitting', 'on' );
    g = simplify( g.' ); 
else
    vals = chebpolyval2(f, sample, sample);  
    g = chebfun( max( vals, [], 2 ), dom(3:4), 'splitting', 'on' );
    g = simplify( g );
end

end
