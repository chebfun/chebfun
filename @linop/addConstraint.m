function L = addConstraint(L, varargin)
%ADDCONSTRAINT Append to linop constraints.
%   L = ADDCONSTRAINT(L, FUN, VAL) adds a new constraint on the linop L. The
%   functional FUN when applied to a function will be required to be VAL.
%
%   L = ADDCONSTRAINT(L, 'periodic') replaces all side conditions with continuity
%   meant to ensure that the function is periodic.
%
%   Example:
%     [Z, I, D, C] = linop.primitiveOperators([-1 1]);
%     [z, E, s] = linop.primitiveFunctionals([-1 1]);
%     A = [ D^2+I, D;  Z, D^2-I ];    % 2-by-2 chebmatrix
%     op1 = [ E(-1), z ];   
%     A = addConstraint(A, op1, 1);   % impose u{1}(-1) + 0*u{2} = 1
%     op2 = [ E(1), -E(1) ]; 
%     A = addConstraint(A, op2, 0);   % impose u{1}(1) - u{2}(1) = 1
%     op3 = [ z, E(1)*D ]; 
%     A = addConstraint(A, op3, -1);  % impose 0*u{1} + (u{2}')(1) = -1
%     op4 = [ s, z ];
%     A = addConstraint(A, op4, 1);   % impose sum(u{1}) + 0*u{2} = 1
%
% See also LINOP.ADDCONTINUITY, LINOPCONSTRAINT. 

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

if ( isequal(varargin{1}, 'periodic') )
    if ( ~isempty(L.constraint) )
        warning('CHEBFUN:LINOP:addConstraint:overwrite', ...
            'Clearing existing constraints to replace with periodicity.')
    end
    
    L = deriveContinuity(L, L.domain, true);  % modifies continuity property
    
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
