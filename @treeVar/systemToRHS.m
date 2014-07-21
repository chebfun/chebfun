function anonFun = systemToRHS(systemInfix, varArrays, coeffs,  indexStart, totalDiffOrders)
% FEVAL Evaluates an anon with an input argument, similar to f(u) where f
% is an anonymous function and u is the argument.

% Copyright 2011 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Load these variables into workspace
loadVariables(varArrays)

% Add the @(u) part
infixForm = '';
for sysCounter = 1:length(systemInfix)
    varString = '';
    for orderCounter = 1:totalDiffOrders(sysCounter)-1
        varString = [varString, ...
            sprintf('u(%i); ', orderCounter + indexStart(sysCounter))];
    end
    if ( isnumeric(coeffs{sysCounter}) )
        coeffStr = sprintf('coeffs{%i}', sysCounter);

    else
        % If COEFF is not numeric, it must be a CHEBFUN. But that requires us to
        % evaluate it at every point T when we evaluate the ODE fun.
        coeffStr = sprintf('coeffs{%i}(t)', sysCounter);

    end
    infixForm = [infixForm, varString , systemInfix{sysCounter}, ...
        './' coeffStr, '; '];
end
% Get rid of the last ; 
infixForm(end-1:end) = [];
infixForm = ['@(t,u)[', infixForm , ']'];
anonFun = eval(infixForm);

end

function loadVariables(varArrays)

for arrCounter = 1:length(varArrays)
    varArray = varArrays{arrCounter};
    
    for i=1:size(varArray, 1)
        assignin('caller',varArray{i, 1}, varArray{i, 2})
    end
end
end