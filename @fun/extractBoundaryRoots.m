function f = extractBoundaryRoots(f)
%EXTRACTBOUNDARYROOTS   Extract boundary roots and represent them by an
%   appropriate ONEFUN.
%
%   G = EXTRACTBOUNDARYROOTS(F) returns a FUN G the SMOOTHPART of whose ONEFUN 
%   is free of roots at the boundary points -1 and 1 and these roots are 
%   represented by a proper EXPONENTS.
%
% See also ROOTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( ~isnumeric(f) )
    
    % Get the boundary values:
    endVals = [get(f, 'lval'); get(f, 'rval')];
    tol = 2e1*eps*get(f, 'vscale');
    
    if ( any(any(abs(endVals) < repmat(tol, 2, 1))) )
        
        if ( issmooth(f) ) 
            
            % If F is a SMOOTHFUN, call EXTRACTBOUNDARYROOTS@SMOOTHFUN:            
            [f.onefun, rootsLeft, rootsRight] = extractBoundaryRoots(f.onefun);

            data.exponents = [rootsLeft rootsRight];
            h = singfun(f.onefun, data);
            f.onefun = h;
            
        elseif ( issing(f) ) 
            
            % If F is a SINGFUN, call EXTRACTBOUNDARYROOTS@SINGFUN:            
            f.onefun = extractBoundaryRoots(f.onefun);

        end
        
    end

end
