function pref = determineDiscretization(N, lengthDom, pref)
%DETERMINEDISCRETIZATION    Determine discretization for a CHEBOP object.
%
%   PREFOUT = DETERMINEDISCRETIZATION(N, LENGTHDOM, PREFIN) choses the
%   correct discretization PREFOUT.DISCRETIZATION to be used when solving
%   problems (BVP/EIGS/EXPM), modeled by a CHEBOP N.
%
% See also CHEBOP/CLEARPERIODICBC, CHEBOP/SOLVEBVP, CHEBOP/EIGS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Case 1. A string was passed:
if ( isa(pref.discretization, 'char') )
    
    % Case 1.1. Determine the discretization if the user wants a discretization 
    % using values:
    if ( strcmpi(pref.discretization, 'values') )
        % The default for the periodic case if TRIGCOLLOC. But, since TRIGCOLLOC
        % does not support breakpoints, it will be used only if there are no
        % breakpoints.
        if ( isa(N.bc, 'char') && strcmpi(N.bc, 'periodic') && lengthDom < 3 )
            pref.discretization = @trigcolloc;
        % Otherwise (i.e. periodic + breakpoints or other boundary 
        % conditions), use CHEBCOLLOC2:
        else
            pref.discretization = @chebcolloc2;
        end
        
    % Case 1.2. Determine the discretization if the user wants a discretization 
    % using coefficients:
    elseif ( strcmpi(pref.discretization, 'coeffs') )
        % Same here with TRIGSPEC and breakpoints:
        if ( isa(N.bc, 'char') && strcmpi(N.bc, 'periodic') && lengthDom < 3 )
            pref.discretization = @trigspec;
        else
            pref.discretization = @ultraS;
        end
    % Error:    
    else
        error('CHEBFUN:CHEBOP:solvebvp:determineDiscretization', ...
        'PREF.DISCRETIZATION should be VALUES or COEFFS.\n')
    end

% Case 2. A function handle referring to a class was passed:
else
    % If TRIGCOLLOC or TRIGSPEC were chosen and there are breakpoints, we need
    % to throw an error:
    if ( ( isequal(pref.discretization, @trigcolloc) || ...
            isequal(pref.discretization, @trigspec) ) && lengthDom > 2 )
        error('CHEBFUN:CHEBOP:solvebvp:breakpointsInDomain', ...
        ['Problems with periodic boundary conditions where breakpoints \n', ...
        'are present cannot be solved with the TRIGCOLLOC/TRIGSPEC class.\n' ...
        'Please change the discretization to CHEBCOLLOC1/2 or ULTRAS.'])
    end
end

end
