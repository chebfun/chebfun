function pref = determineDiscretization(N, L, isPrefGiven, pref)
%DETERMINEDISCRETIZATION    Determine discretization for a CHEBOP object.
%
%   PREFOUT = DETERMINEDISCRETIZATION(N, L, ISPREFGIVEN, PREFIN) choses the
%   correct discretization PREFOUT.DISCRETIZATION to be used when solving
%   problems (BVP/EIGS/EXPM), modeled by a CHEBOP N and a LINOP L.
%
% See also CHEBOP/CLEARPERIODICBC, CHEBOP/SOLVEBVP, CHEBOP/EIGS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Case 1. Determine the discretization if the user wants a discretization using 
% values:
if ( isequal(pref.discretization, @values) )
    % The default for the periodic case if TRIGCOLLOC. But, since TRIGCOLLOC 
    % does not support breakpoints, it will be used only if there are no 
    % breakpoints. Note that here and below, we look at L.domain rather
    % than N.domain, as the domain of the LINOP L will also include any 
    % breakpoints arising from discontinuous coefficients of N (which we only 
    % become aware of when we do the linearization):
    if ( isa(N.bc, 'char') && strcmpi(N.bc, 'periodic') && ...
            length(L.domain) < 3 )
        pref.discretization = @trigcolloc;
    % Otherwise (i.e. periodic + breakpoints or other boundary conditions), 
    % use CHEBCOLLOC2:
    else
        pref.discretization = @chebcolloc2;
    end
    
% Case 2. Determine the discretization if the user wants a discretization using 
% coefficients:
elseif ( isequal(pref.discretization, @coeffs) )
    if ( isa(N.bc, 'char') && strcmpi(N.bc, 'periodic') && ...
            length(L.domain) < 3 )
        pref.discretization = @ultraS;
    else
        pref.discretization = @ultraS;
    end

% Case 3. A function handle referring to a class was passed:
else
    % If TRIGCOLLOC was chosen and there are brekpoints, we need to throw an 
    % error:
    if ( isequal(pref.discretization, @trigcolloc) && length(L.domain) > 2 )
        error('CHEBFUN:CHEBOP:solvebvp:breakpointsInDomain', ...
        ['Problems with periodic boundary conditions where breakpoints \n', ...
        'are present cannot be solved using the TRIGCOLLOC class.\n' ...
        'Please change the discretization to CHEBCOLLOC1/2 or ULTRAS.'])
    end
    % If the boundary conditions are periodic, change the discretization to
    % TRIGCOLLOC, unless one (or more) of the following applies:
    % - The user passed a PREF object.
    % - Breakpoint(s) are present.
    if ( isa(N.bc, 'char') && strcmpi(N.bc, 'periodic') && ~isPrefGiven ...
        && length(L.domain) < 3 )
        pref.discretization = @trigcolloc;
    end
end

end
