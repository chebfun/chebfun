function [N, L, pref] = clearPeriodicBC(N, L, pref)
%CLEARPERIODICBC    Clear periodic bounadry conditions.
%   [N, L, PREF] = CLEARPERIODICBC(N, L, PREF) clears N.BC, L.CONSTRAINT,
%   and L.CONTINUITY, if PREF.DISCRETIZATION is FOURCOLLOC.
%
% See also CHEBOP/DETERMINEPREF, CHEBOP/SOLVEBVP, CHEBOP/EIGS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check boundary conditions if using FOURCOLLOC.
if ( isequal(pref.discretization, @fourcolloc) )
    if ( isempty(N.bc) )
        % No need to clear the BCs, do nothing!
    elseif ( isa(N.bc, 'char') && strcmpi(N.bc, 'periodic') )
        % FOURCOLLOC uses periodic functions, so there is no need to specify
        % boundary conditions. We clear them out of the chebop object to avoid
        % problems later in the code.
        N.bc = [];
        L.constraint = [];
        L.continuity = [];
    else
        error('CHEBFUN:CHEBOP:solvebvp:bc', ...
            'FOURCOLLOC only works with periodic boundary conditions.');
    end
end

end
