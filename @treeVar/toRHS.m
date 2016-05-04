function anonFun = toRHS(systemInfix, varArrays, coeffs,  indexStart, totalDiffOrders)
%TORHS   Convert infix expressions to anonymous function suited for ODE solvers.
%   Calling sequence:
%      ANONFUN = TORHS(SYSTEMINFIX, VARARRAYS, COEFFS, INDEXSTART, ...
%       TOTALDIFFORDERS)
%   where the inputs are
%      SYSTEMFINFIX:    A cell array of strings of infix expressions, each of
%                       which describes the right-hand side of a differential
%                       equation suitable for MATLAB's ODE solvers.
%      VARARRAYS:       A cell array of cell arrays, so that the ith element
%                       of VARARRAYS contains the values of the variables
%                       (scalars and CHEBFUNS) that appear in the ith equation
%                       of SYSTEMINFIX.
%      COEFFS:          A cell array of the coefficients that multiply the
%                       left hand side of SYSTEMINFIX. For example, for the
%                       expression 5*diff(u) + u, COEFFS{1} = 5, while for
%                       sin(x)*diff(u) + u, COEFFS{1} will be the CHEBFUN
%                       sin(x).
%      INDEXSTART:      A vector that denotes at which index we should start
%                       indexing each variable from.
%      TOTALDIFFORDERS: A vector which contains the maximum diffOrder of each
%                       variable that appears in the problem.
%   and the output is
%      ANONFUN:    An anonymous function on a form suitable to be passed to
%                  the MATLAB ODE solvers for describing the first order ODE
%                  system.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Load the variables from VARARRAYS into the current workspace.
loadVariables(varArrays)

% Initialise an empty string, to which we'll add the infix expressions.
infixForm = '';

% Loop through the infix expressions from SYSTEMINFIX to build one big
% expression which we convert to an anonymous function.
for sysCounter = 1:length(systemInfix)
    
    % We need to start by doing the first order reformulation. For example, if
    % we are dealing with diff(u,3) + sin(u), whose first order reformulation is
    %   u'_1 = u_2
    %   u'_2 = u_3
    %   u'_3 = -sin(u_1)
    % we generate the corresponding strings for the RHS of the coupled system,
    % in this case
    %   'u(2); u(3);'
    
    % Initialise an empty string for the first order reformulation:
    varString = '';
    for orderCounter = 1:totalDiffOrders(sysCounter)-1
        varString = [varString, ...
            sprintf('u(%i); ', orderCounter + indexStart(sysCounter))];
    end
    
    % The differential equations can have coefficients not equal to 1 on the
    % left-hand side, e.g. for 5*diff(u) = sin(u). These are stored in the
    % COEFFS cell array, in this example, COEFFS{1} = 5. Since the MATLAB
    % solvers require the ODEs to be specified on the format y' = f(t,y), we
    % need to divide through by the COEFFS.
    
    if ( isnumeric(coeffs{sysCounter}) )
        % If the COEFF operating the LHS of the current differential equation is
        % a scalar, we create a string that picks it out of the COEFFS cell.
        coeffStr = sprintf('coeffs{%i}', sysCounter);

    else
        % If COEFF is not numeric, it must be a CHEBFUN. But that requires us to
        % evaluate it at every point T when we evaluate the ODE fun, so we add
        % the '(t)' part as well, so that it'll be evaluated at the correct
        % coordinate throughout the ODE solver.
        coeffStr = sprintf('coeffs{%i}(t)', sysCounter);

    end
    
    % Add the current equation INFIXFORM, which contains all the equations so
    % far.
    infixForm = [ infixForm, varString , systemInfix{sysCounter}, ...
        './' coeffStr, '; ' ];
end

% Add the @(t,u) to the front of the infix expression, and surround it by [].
infixForm = [ '@(t,u)[', infixForm , ']' ];

% Eval the INFIXFORM string to get the anonymous function.
anonFun = eval(infixForm);

end

function loadVariables(varArrays)
%LOADVARIABLES   Load variables into the memory scope of the caller function.

for arrCounter = 1:length(varArrays)
    varArray = varArrays{arrCounter};
    
    % Loop through the array of variables, and assign into memory of the
    % caller scope:
    for i = 1:size(varArray, 1)
        assignin('caller', varArray{i, 1}, varArray{i, 2})
    end
end
end
