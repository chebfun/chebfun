function [N, L] = clearPeriodicBCs(N, L)
%CLEARPERIODICBCS    Clear periodic boundary conditions.
%   [N, L] = CLEARPERIODICBCS(N, L, PREF) clears N.BC, L.CONSTRAINT, and
%   L.CONTINUITY. This method only gets called if the current discretization in
%   use is TRIGCOLLOC, and is required since TRIGCOLLOC by construction will
%   only return periodic function as the solution, so there is no need to impose
%   periodic conditions on the discretized linear systems that arise.
%
% See also CHEBOP/DETERMINEDISCRETIZATION, CHEBOP/SOLVEBVP, CHEBOP/EIGS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(N.bc) )
    % N.BC is empty already, so do nothing!
elseif ( isa(N.bc, 'char') && strcmpi(N.bc, 'periodic') )
    % Clear the conditions:
    N.bc = [];
    L.constraint = [];
    L.continuity = [];
else
    error('CHEBFUN:CHEBOP:clearPeriodicBCs:nonperiodic', ...
        'TRIGCOLLOC only works with periodic boundary conditions.');
end

end
