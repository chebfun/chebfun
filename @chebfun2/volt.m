function f = volt( K, v )
%VOLT  Volterra integral operator.
%
% V = VOLT(K,f) returns a row chebfun resulting from the integral
%     
%      f(x) = (K*v)(x) = int( K(x,y) v(y), y=a..x )
%   
% The kernel function K(x,y) must be a smooth chebfun2.
%
% Example:
% 
% f = VOLT(chebfun2(@(x,y) exp(x-y)),chebfun('x'));  
% 
% See also FRED. 

% Copyright 2013 by The University of Oxford and The Chebfun2 Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun2 information.

if ( ~isa( K, 'chebfun2' ) )
    error('CHEBFUN2:volt:input','First argument must be a chebfun2');
end

% Get the low rank representation for f. 
cols = K.cols; 
rows = K.rows; 
piv = K.pivotValues; 
d = 1./piv; 
d(d==inf) = 0;  % set infinite values to zero. 
dom = K.domain; 

if isa( v, 'function_handle' )
    v = chebfun( v, dom(3:4) );  % convert to a chebfun on the right interval. 
end

% Domain compatibility: 
if ( ~domainCheck(cols, v) )
    error('CHEBFUN2:FRED:CHEBDOMAIN','Domain of chebfun and chebfun2 kernel do not match');
end

RR = diag( d ) * rows;

% Cumsum with cols and v:  (This can be sped up.)
f = chebfun( 0, dom(1:2) );
for jj = length(K) : -1 : 1
    CC = cumsum( v .* cols(:,jj) );
    f = f + CC .* RR(:,jj);
end

% Transpose to preserve linear algebra: 
f = f.'; 

end