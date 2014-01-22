function L = addbc(L, varargin)
%ADDBC     Append to linop constraints.
%   L = ADDBC(L,FUN,VAL) adds a new constraint on the linop L. The
%   functional FUN when applied to a function will be required to be VAL.
%
%   L = ADDBC(L,'periodic') replaces all side conditions with continuity
%   meant to ensure that the function is periodic. 
%
%   See also LINOPCONSTRAINT. 

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

if isequal( varargin{1}, 'periodic' )
    if ~isempty( L.constraint )
        warning('Clearing existing constraints to replace with periodicity.')
    end
    
    L = deriveContinuity(L, true);  % modifies continuity property
    
    % We're going to move the periodic continuity to the constraints, so that
    % we're not fooled into thinking that the interior breakpoints have been
    % done.
    L.constraint = L.continuity;
    L.continuity = linopConstraint();
    
else
    % Append the input constraint to the LINOPCONSTRAINT currently stored in the
    % input linop.
    L.constraint = append(L.constraint, varargin{:});
end

end
