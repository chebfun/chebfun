function g = max( f, varargin )
%MAX  Maximum value of a chebfun2 in one direction
%
% MAX(f) returns a chebfun representing the maximum of the chebfun2 along
% the y direction, i.e
%
%   MAX(f) = @(x) max( f ( x, : ) )
%
% MAX(f,[],dim) returns a chebfun representing the maximum of f along the
% DIM direction. If DIM = 1 is along the y-direction and DIM = 2 is along
% the x-direction. 
%
% This function is not considered highly accurate.  Expect no more than 
% five or six digits of precision. For the global maximum use MAX2. 
%
% See also MIN, MAX2, MIN2, MINANDMAX2.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check: 
if ( isempty(f) ) 
    error('CHEBFUN2:MAX:INPUT','Chebfun2 is empty');
end

% Default to maximum along the y direction: 
if ( nargin == 1 )  
    dim = 1;
end

% Do not allow max(F,G): 
if ( nargin > 1 )   
    if ( ~isempty(varargin{1}) )
        error('CHEBFUN2:max','Unable to maximise two chebfun2 objects.');
    end
end

if ( nargin == 2 )  % default to maximum along the y direction. 
    dim = 1;
end

% What variable to maximise over? 
if ( nargin == 3 )  
    dim = varargin{ 2 };
end


%%
% We have no idea how to achieve this in a really fast and efficient way.
% This is an attempt to at least get a result, but there are no guarantees
% that this will be 16 digit accurate.  

dom = f.domain;
sample = 2049; 
if ( dim == 1 )
    vals = chebpolyval2(f, sample, sample); 
    g = chebfun( max( vals )', dom(1:2), 'splitting', 'on' );
    g = g.';
    g = simplify( g ); 
else
    vals = chebpolyval2(f, sample, sample);  
    g = chebfun( max( vals, [], 2 ), dom(3:4), 'splitting', 'on' );
    g = simplify( g );
end

end
