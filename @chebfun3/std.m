function g = std(f, varargin)
%STD   Standard deviation of a CHEBFUN3 along one variable.
%   G = STD(F) returns the standard deviation of F in the x-variable 
%   (default). If F is defined on the cuboid [a, b] x [c, d] x [e, g] then
%
%                      b 
%                     /
%     G.^2 = 1/(b-a)  | (F(x,y,z) - mean(F, 1))^2 dx
%                     /
%                     a
%
%   The output G is a CHEBFUN2 object over the rectangle [c, d] x [e, g]. 
%
%   G = STD(F, FLAG, DIM) takes the standard deviation along the x, y, or
%   z directions if DIM = 1, 2, or 3, respectively.  FLAG is ignored and
%   kept in this function so the syntax agrees with the Matlab STD command.
%
% See also CHEBFUN/STD, CHEBFUN2/STD, CHEBFUN2/STD2, CHEBFUN3/STD2 and 
% CHEBFUN3/STD3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(f) )
    g = chebfun2();
    return
end

dom = f.domain; 
if ( nargin < 3 )
    dim = 1;   % default to std over the x-variable.
elseif ( nargin == 3 )
    dim = varargin{ 2 }; 
else 
    error( 'CHEBFUN:CHEBFUN3:std:nargin', 'Too many input arguments.' ); 
end

%   std(X) = sqrt( mean( (X - mean(X) )^2 ) ).
if ( dim == 1 )          % x-variable.
    mx = chebfun3(@(x,y,z) feval(mean(f, 1), y, z), dom);
    g = sqrt(1/(diff(dom(1:2))) * sum((f - mx).^2, 1)) ;
    
elseif ( dim == 2 )      % y-variable.
    my = chebfun3(@(x,y,z) feval(mean(f, 2), x, z), dom);
    g = sqrt(1/(diff(dom(3:4))) * sum((f - my).^2, 2));
    
elseif ( dim == 3 )      % z-variable.
    mz = chebfun3(@(x,y,z) feval(mean(f, 3), x, y), dom);
    g = sqrt(1/(diff(dom(5:6))) * sum((f - mz).^2, 3));
    
else
    error('CHEBFUN:CHEBFUN3:std:dim', ...
        'Third argument should have value 1, 2 or 3.');
end

end
