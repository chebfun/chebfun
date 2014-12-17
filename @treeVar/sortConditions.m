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
%              N.LBC/RBC, so that they match the order required by the MATLAB 
%              ODE solvers.
%
%   Example:
%      Assume that for a coupled system, we want to specify the ICs
%         u(0) = 1, u'(0) = 2, v(0) = 0, v'(0) = 3.
%      One anonymous function we could create to specify those conditions and
%      append to a CHEBOP is
%         fun = @(u,v) [u-1; diff(v)-3; v; diff(u)-2]; 
%      However, when we evaluate this function to pick out the initial 
%      conditions, we get the vector
%         [1; 3; 0; 2]
%      which is in the incorrect order to be passed to the MATLAB solvers, which
%      require the conditions to be of the order u, u', v, v'. Calling
%      sortConditions(),
%         idx = treeVar.sortConditions(fun, dom) 
%         idx =
%             1     4     3     2
%      gives the correct order in which the above vector has to be sorted so 
%      that the values are in the correct order for MATLAB.

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
    
    % Ensure that only initial/final conditions on the form we support appear:
    assert(acceptedCondition(tempTree), ...
        'CHEBFUN:TREEVAR:sortConditions:unsupportedCondition', [ ...
        'Initial/final condition not supported. Please ensure that the\n' ...
        'unknown function(s) appear only once in each condition,\n' ...
        'e.g. not as u + u'' - 1, and that the coefficient of the unknown\n'...
        'function is 1, e.g. not 5*u-1.']);
    
    % What's the active variable in the current tree (i.e. what variable did the
    % constraint apply to)?
    activeVar = find(tempTree.ID == 1);
    
    % What's the diffOrder of the current constraint?
    activeDiffOrder = tempTree.diffOrder;
    activeDiffOrder = activeDiffOrder(activeVar);
    
    % Store in a list what variable the current constraint applies
    % to, and what the current diffOrder is:
    varList{activeVar} = [ varList{activeVar}, tCounter ];
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
    idx = [ idx, tempIndex ];
end

end

function accepted = acceptedCondition(tree)
%ACCEPTEDCONDITION   Check whether we have encountered an unsupported IC
%   
% ACCEPTEDCONDITION(TREE) returns FALSE if we are working with functions such as
%   @(u) 5*u - 1
%   @(u) u + diff(u),
% which are not supported by the automatic first order reformulation. Otherwise,
% it returns TRUE.

% If TREE is not a struct, or it has height 0, it can't cause any harm:
if ( ~isstruct(tree) || ( tree.height == 0 ) )
    accepted = true;

% Deal with the case where the operator is + or -
elseif ( any(strcmp(tree.method, {'plus', 'minus'})) )

    % If TREE is of height 1, and it's operator is + or -, we're OK
    if ( tree.height == 1 )
        accepted = true;
        return
    end
    
    % Otherwise, we have to be more careful. If either argument is not a struct,
    % we only need to check the other. If both arguments are trees, we need
    % ensure that both don't have any IDs, and in the case that happens, check
    % the validity of each subtree.
    if ( ~isstruct(tree.left) )
        accepted = acceptedCondition(tree.right);
    elseif ( ~isstruct(tree.right) )
        accepted = acceptedCondition(tree.left);
    elseif ( any(tree.left.ID) && any(tree.right.ID) )
        accepted = false;
    else
        accepted = acceptedCondition(tree.right) && ...
            acceptedCondition(tree.right);
    end
    
    
% We're happy to support a tree with the diff method if it's height is 1:
elseif ( strcmp(tree.method, 'diff') && ( tree.height == 1 ) )
    accepted = true;
    
% Encountering any other operator implies that we had a condition on the form
% @(u) 5*u-1, which is not supported.
else
    accepted = false;
end

end