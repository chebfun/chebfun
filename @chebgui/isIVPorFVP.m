function isIorF = isIVPorFVP(guifile, bcInput, allVarString)
%%
dom = guifile.domain;
spaceLoc = strfind(dom, ' ');
firstSpace = min(spaceLoc);
lastSpace = max(spaceLoc);

leftDomStr = dom(2:firstSpace-1);
rightDomStr = dom(lastSpace+1:end-1);

numStr = '[+-]?[0-9]*\.?[0-9]*(e[0-9]+)*';
% Find parts of the boundary conditions input where condition are imposed on the
% left or right ends of the interval.
leftConditionLocs = ...
    regexp(bcInput,[allVarString, '''*\(', leftDomStr, '\)']);
rightConditionLocs = ...
    regexp(bcInput,[allVarString, '''*\(', rightDomStr, '\)']);

% Find any parts of the string where condition of the form u(0) are imposed.
pointEvalConditionLocs = ...
    regexp(bcInput,[allVarString, '''*\(', numStr, '\)']);

% Find where the variables appear at places other then where point evaluation
% takes place. This happens e.g. in conditions such as 'sum(u) = 1':
globalConditionLocs = ...
    regexp(bcInput, [allVarString, '[^''*\(', numStr, '\)]']);

% Below, we like to work with cells which arise in the case where more than one
% condition is imposed. So ensure LEFTCONDITIONLOCS, RIGHTCONDITIONLOCS and
% POINTEVALCONDITIONLOCS are all cells:
leftConditionLocs = vec2cell(leftConditionLocs);
rightConditionLocs = vec2cell(rightConditionLocs);
pointEvalConditionLocs = vec2cell(pointEvalConditionLocs);
globalConditionLocs = vec2cell(globalConditionLocs);

% Find locations in the inputs where either left or right conditions are
% imposed:
leftRigthConditionLocs = ...
    cellfun(@union, leftConditionLocs, rightConditionLocs, ...
    'uniformOutput', false);
interiorConditionLocs = ...
    cellfun(@setdiff, pointEvalConditionLocs, leftRigthConditionLocs, ...
    'uniformOutput', false);
leftConditionsImposed = any(cellfun(@any, leftConditionLocs));
rightConditionsImposed = any(cellfun(@any, rightConditionLocs));
interiorConditionsImposed = ~all(cellfun(@isempty,interiorConditionLocs));
globalConditionsImposed = ~all(cellfun(@isempty, globalConditionLocs));

% If any conditions are imposed at points other than the endpoints, we certainly
% don't have an IVP/FVP. Also, if both left and right endpoint conditions are
% imposed, we're dealing with a BVP, not an IVP/FVP.
if ( interiorConditionsImposed || globalConditionsImposed || ...
        (leftConditionsImposed && rightConditionsImposed) )
    isIorF = 0;
elseif ( leftConditionsImposed )
    % We only have left condition(s) imposed, so we're working with an IVP:
    isIorF = 1;
else
    % We only have right condition(s) imposed, so we're working with an FVP:
    isIorF = 2;
end
%%
end

function out = vec2cell(vec)
if ( ~iscell(vec) )
    out = {vec};
else
    out = vec;
end
end