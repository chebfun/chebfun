function anonFun = toAnon(infix, varArray)
%TOANON   Convert the infix form of an expression to an anonymous function
%   Calling sequence:
%      ANONFUN = TOANON(INFIX, VARARRAY)
%   where the inputs are:
%      INFIX:      A string, describing a mathematical expression on infix form,
%                  written in MATLAB syntax.
%      VARARRAY:   A MATLAB cell, which contains the name and values of all
%                  variables that appear in INFIX.
%   and the output is:
%      ANONFUN:    An anonymous function, corresponding to the function 
%                  represented in INFIX, with the variables in VARARRAY loaded 
%                  in its workspace.
%
%   This method is used when converting expressions to first order format to
%   obtain an anonymous function that when evaluated, returns the coefficient
%   multiplying the highest order derivative appearing in the expression.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Load these variables into workspace:
loadVariables(varArray)

% Add the @(u) part:
infix = ['@(t,u) ', infix];

% Call eval() on the string to obtain an anonymous function:
anonFun = eval(infix);

end

function loadVariables(varArray)
%LOADVARIABLES   Load variables into the memory scope of the caller function.

% Loop through the array of variables, and assign into memory of the caller
% scope:
for i = 1:size(varArray, 1)
    assignin('caller', varArray{i, 1}, varArray{i, 2})
end

end
