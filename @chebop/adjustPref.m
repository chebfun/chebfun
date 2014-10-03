function [N, L, pref] = adjustPref(N, L, isPrefGiven, pref)
%ADJUSTPREF    Adjust preferences for a CHEBOP object.
%   [N, L, PREF] = ADJUSTPREF(N, L, ISPREFGIVEN, PREF) choses the right
%   discretization PREF.DISCRETIZATION to use to solve an ODE problem or an 
%   eigenvalue problem with periodic boundary conditions, modeled by a
%   CHEBOP N and a LINOP L.
%
% See also CHEBOP/SOLBEBVP, CHEBOP/EIGS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% The user passed a CHEBOPPREF object PREF with PREF.DISCRETIZATION set to 
% FOURCOLLOC, or FOURCOLLOC is the default discretization.
% Since FOURCOLLOC does not support breakpoints, we need to throw an error if
% breakpoints are present. Note that here and below, we look at L.domain rather
% than N.domain, as the domain of the LINOP L will also include any breakpoints
% arising from discontinuous coefficients of N (which we only become aware of
% when we do the linearization).
if ( isequal(pref.discretization, @fourcolloc) && length(L.domain) > 2 )
    error('CHEBFUN:CHEBOP:solvebvp:breakpointsInDomain', ...
        ['Problems with periodic boundary conditions where breakpoints \n', ...
        'are present cannot be solved using the FOURCOLLOC class.\n' ...
        'Please change the discretization to CHEBCOLLOC1/2 or ULTRAS.'])
end

% If the boundary conditions are periodic, change the discretization to 
% FOURCOLLOC, unless one (or more) of the following applies:
% - The user passed a PREF object.
% - The default discretization is ULTRAS.
% - Breakpoint(s) are present.
if ( isa(N.bc, 'char') && strcmpi(N.bc, 'periodic') && ~isPrefGiven ...
        && ~isequal(pref.discretization, @ultraS) && length(L.domain) < 3 )
    pref.discretization = @fourcolloc;
end

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
