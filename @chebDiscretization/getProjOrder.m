function projOrder = getProjOrder(L)

% TODO: Document

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~isa(L, 'linop') )
    
    % No adjustment
    projOrder = 0;
    
else
    
    
    projOrder = sizeReduction(L);
    
end
    
end