function g = std( f, varargin )
%STD   Standard deviation of a SEPARABLEAPPROX along one variable.
%   G = STD(F) returns the standard deviation of F in the y-variable (default).
%   That is, if F is defined on the rectangle [a,b] x [c,d] then
%
%                         d 
%                        /
%     std(F)^2 = 1/(d-c) | ( F(x,y) - mean(F,1) )^2 dy
%                        /
%                        c
%
%   G = STD(F, FLAG, DIM) takes the standard deviation along the y-variable if
%   DIM = 1 and along the x-variable if DIM = 2. The FLAG is ignored and kept in
%   this function so the syntax agrees with the Matlab STD command.
%
% See also CHEBFUN/STD, SEPARABLEAPPROX/MEAN.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( f ) ) 
    g = chebfun();
    return
end

dom = f.domain; 
if ( nargin < 3 )
    dim = 1;   % default to std over the y-variable. 
elseif ( nargin == 3 )
    dim = varargin{ 2 }; 
else 
    error( 'CHEBFUN:SEPARABLEAPPROX:std:nargin', 'Too many input arguments.' ); 
end

if ( dim == 1 )          % y-variable.
    mx = chebfun2( @(x,y) feval( mean(f, 2), x ), dom );
    g = sqrt( 1/( diff( dom(3:4) ) ) * sum( ( f - mx ).^2, 1 ) ) ;
elseif ( dim == 2 )      %  x-variable.
    my = chebfun2( @(x,y) feval( mean(f, 1), y), dom );
    g = sqrt( 1/( diff( rect(1:2) ) ) * sum( ( f - my ).^2, 2 ) );
else
    error('CHEBFUN:SEPARABLEAPPROX:std:dim', ...
        'Third argument should have value 1 or 2.');
end

end
