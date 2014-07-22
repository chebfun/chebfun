function anonFun = toAnon(infix, varArray)
% FEVAL Evaluates an anon with an input argument, similar to f(u) where f
% is an anonymous function and u is the argument.

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Load these variables into workspace
loadVariables(varArray)

% Add the @(u) part
infix = ['@(t, u) ', infix];

anonFun = eval(infix);

end

function loadVariables(varArray)

for i=1:size(varArray, 1)
    assignin('caller',varArray{i, 1}, varArray{i, 2})
end

end