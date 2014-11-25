function idx = sortConditions(funIn, domain)
%SORTCONDITIONS   Returns how the results of evaluating BCs should be sorted.
%   Calling sequence:
%      IDX = SORTCONDITIONS(FUNIN, DOMAIN)
%   where the inputs are:
%      FUNIN:  An anonymous function, that usually specifies N.LBC or N.RBC of a
%              CHEBOP.
%      DOMAIN: The domain of the problem.
%   and the output is
%      IDX:    A vector with indices on how to sort the results of evaluating
%              N.LBC/RBC, so that they match the order required by the MATLAB ODE
%              solvers.
%
%   Example:
%      Assume that for a coupled system, we want to specify the ICs
%         u(0) = 1, u'(0) = 2, v(0) = 0, v'(0) = 3.
%      One anonymous function we could create to specify those conditions and
%      append to a CHEBOP is
%         fun = @(u,v) [u-1; diff(v)-3; v; diff(u)-2]; 
%      However, when we evaluate this function to pick out the initial conditions, 
%      we get the vector
%         [1; 3; 0; 2]
%      which is in the incorrect order to be passed to the MATLAB solvers, which
%      require the conditions to be of the order u, u', v, v'. Calling
%      sortConditions(),
%         idx = treeVar.sortConditions(fun, dom) 
%         idx =
%             1     4     3     2
%      gives the correct order in which the above vector has to be sorted so that
%      the values are in the correct order for MATLAB.

% Check how many unknowns appear in FUNIN.
numArgs = nargin(funIn);
args = cell(numArgs, 1);

% The ID vector to be passed to the TREEVAR constructor.
argsVec = zeros(1, numArgs);

% Populate the args cell with TREEVAR objects.
for argCount = 1:numArgs
    % Set the ID of the current variable to 1:
    argsVec(argCount) = 1;
    % Construct the TREEVAR:
    args{argCount} = treeVar(argsVec, domain);
    % Reset the index vector:
    argsVec = 0*argsVec;
end

% Evaluate FUNIN with the TREEVAR arguments:
bcResults = funIn(args{:});

% Look at the results of evaluating the boundary conditions, find what
% constraint operated on what variable, and what its diffOrder was:
varList = cell(numArgs, 1);
diffOrderList = varList;

% Loop through the resulting syntax trees.
for tCounter = 1:length(bcResults)
    % Current tree we're looking at:
    tempTree = bcResults(tCounter).tree;
    
    % Check whether more than one variable appear in the condition, which we
    % don't support:
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