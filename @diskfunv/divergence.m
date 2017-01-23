function div = divergence( F ) 
%DIVERGENCE  Numerical  divergence of a DISKFUNV. 
%   D = DIVERGENCE( F ) returns the numerical divergence of the
%   DISKFUNV. 
%
% See also DIV, GRAD, CURL.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

% Empty check: 
if ( isempty( F ) )
    div = diskfun();
    return
end

% Extract components: 
Fc = F.components; 

% Calculate the divergence: 
div = diff(Fc{1}, 1) + diff(Fc{2}, 2);

end
