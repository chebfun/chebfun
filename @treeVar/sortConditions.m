function idx = sortConditions(funIn, domain)
%SORTBCS    Return a vector with indices on how to sort the results of
%           evaluating N.LBC/RBC.

numArgs = nargin(funIn);
args = cell(numArgs, 1);
argsVec = zeros(1, numArgs);

% Populate the args cell
for argCount = 1:numArgs
    argsVec(argCount) = 1;
    args{argCount} = treeVar(argsVec, domain);
    % Reset the index vector
    argsVec = 0*argsVec;
end

% Evaluate FUNIN with the TREEVAR arguments:
bcResults = funIn(args{:});

% Look at the results of evaluating the boundary conditions, find what
% constraint operated on what variable, and what it's diffOrder was:
varList = cell(numArgs, 1);
diffOrderList = varList;
for tCounter = 1:length(bcResults)
    % Current tree we're looking at:
    tempTree = bcResults(tCounter).tree;
    
    % Check whether more than one variable appear in the condition
    if ( sum(tempTree.ID) > 1 )
        error('CHEBFUN:TREEVAR:sortConditions:nonSeparated', ...
            ['For initial value problems, only separated ', ...
            'conditions are supported.']);
    end
    
    % What's the active variable in the current tree (i.e. what variable did the
    % constraint apply to)?
    activeVar = find(tempTree.ID == 1);
    % What's the diffOrder of the current constraint?
    activeDiffOrder = tempTree.diffOrder;
    activeDiffOrder = activeDiffOrder(activeVar);
    
    % Store in a list what variable the current constraint applies
    % to, and what the current diffOrder is:
    varList{activeVar} = [varList{activeVar}, tCounter];
    diffOrderList{activeVar} = [diffOrderList{activeVar}, ...
        activeDiffOrder];
end

% Initalise an index vector to be returned
idx = [];

% Go through the list of what variables appeared in what constraints, and sort
% them based on diffOrders:
for varCounter = 1:numArgs
    [dummy, diffOrderIndex] = sort(diffOrderList{varCounter});
    tempIndex = varList{varCounter}(diffOrderIndex);
    idx = [idx, tempIndex];
end

end