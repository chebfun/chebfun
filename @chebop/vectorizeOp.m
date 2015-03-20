function funOut = vectorizeOp(funIn)
%VECTORIZEOP    Vectorize anonymous functions used for CHEBOP op fields
%
% FUNOUT = VECTORIZEOP(FUNIN), where FUNIN is an anonymous function, returns
% FUNOUT, another anonymous function that corresponds to FUNIN on vectorized
% form.
%
% This method differs from the built-in MATLAB VECTORIZE method in that
% parameters that appear in the original function do get carried over to the 
% vectorized function.
%
% Example:
%   b = 3;
%   fun = @(x) b*x^3;
%   x = 0:0.1:0;
%   % fun(x) % This returns an error, suggesting using .^ instead
%   funVec = chebop.vectorizeOp(fun);
%   funVec(x)
%   funVecBuiltIn = str2func(vectorize(fun))
%   funVecBuiltIn(x) % This returns an error, b is undefined
%
% See also: vectorize


% The input function on vectorized form
fVec = vectorize(funIn);

% Begin by obtaining information about potential variables in funIn
fFunc = functions(funIn);
fWork = fFunc.workspace;

% Load the variables in the input function to the current workspace
loadVariables(fWork)


% Ideally, here, we'd call func2str. However, MATLAB will not recognize the
% variables in the workspace if we call it like that. Hence, we need to call
% EVAL, which will include the variables we've loaded into the workspace above
funOut = eval(fVec);

end

function loadVariables(work)

for wCounter = 1:length(work)
    ws = work{wCounter};
    wsFields = fieldnames(ws);
    for varCounter = 1:length(wsFields)
        varName = wsFields{varCounter};
        assignin('caller', varName, ws.(varName))
    end
end
end