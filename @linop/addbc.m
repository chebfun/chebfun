function L = addbc(L, varargin)
%ADDBC  Append to linop constraints (keep existing).

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

if isequal( varargin{1}, 'periodic' )
    if ~isempty( L.constraint )
        warning('Clearing out existing constraints to replace with periodicity.')
    end
    
    L = deriveContinuity(L,true);  % modifies continuity property
    
    % We're going to move the periodic continuity to the constraints, so that
    % we're not fooled into thinking that the interior breakpoints have been
    % done.
    L.constraint = L.continuity;
    L.continuity = linopConstraint();
    
else
    L.constraint = append(L.constraint, varargin{:});
end

end
