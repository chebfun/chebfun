function anonFun = toRHS(infix, varArray, coeff)
% FEVAL Evaluates an anon with an input argument, similar to f(u) where f
% is an anonymous function and u is the argument.

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Load these variables into workspace
loadVariables(varArray)

% Add the @(u) part

if ( isnumeric(coeff) )
    infix = ['@(t, u) [u(2); ', infix, './coeff]'];
else
    % If COEFF is not numeric, it must be a CHEBFUN. But that requires us to
    % evaluate it at every point T when we evaluate the ODE fun.
    infix = ['@(t, u) [u(2); ', infix, './coeff(t)]'];
end
anonFun = eval(infix);

end

function loadVariables(varArray)

for i=1:size(varArray, 1)
    assignin('caller',varArray{i, 1}, varArray{i, 2})
end

end