function [newTree, derTree] = splitTree(tree, maxOrder)
% Isolate where the highest order derivative appears.

% If the highest order derivative is lower than the maximum order of
% the overall expression, we don't need to investigate it further.
if ( ~isstruct(tree) || tree.diffOrder < maxOrder )
    newTree = tree;
    derTree = [];
else
    switch tree.numArgs
        case 1
            [newTree, derTree] = splitTree(tree, maxOrder);
            
        case 2
            if ( strcmp(tree.method, 'diff') )
                newTree = [];
                derTree = tree;
            else
                [newTreeLeft, derTreeLeft] = ...
                    treeVar.splitTree(tree.left, maxOrder);
                [newTreeRight, derTreeRight] = ...
                    treeVar.splitTree(tree.right, maxOrder);
                
                if ( isempty(newTreeLeft) )
                    newTree = oneTree(newTreeRight, tree.method);
                elseif (isempty(newTreeRight) )
                    newTree = oneTree(newTreeLeft, tree.method);
                else
                    if ( ~isstruct(newTreeLeft) )
                        newDiffOrder = newTreeRight.diffOrder;
                    elseif ( ~isstruct(newTreeRight) )
                        newDiffOrder = newTreeLeft.diffOrder;
                    else
                        newDiffOrder = max(newTreeLeft.diffOrder, ...
                            newTreeRight.diffOrder);
                    end
                    
                    newTree = struct('method', tree.method, ...
                        'numArgs', tree.numArgs, ...
                        'left', newTreeLeft, 'right', newTreeRight, ...
                        'diffOrder', newDiffOrder);
                end
                
                if ( isempty(derTreeLeft) )
                    derTree = oneTree(derTreeRight, tree.method);
                elseif (isempty(derTreeRight) )
                    derTree = oneTree(derTreeLeft, tree.method);
                else
                    derTree = struct('method', tree.method, ...
                        'numArgs', tree.numArgs, ...
                        'left', derTreeLeft, 'right', derTreeRight, ...
                        'diffOrder', max(derTreeLeft.diffOrder, ...
                        derTreeRight.diffOrder));
                end
                
            end
    end
end
end

function ot = oneTree(tree, operator)
if ( strcmp(operator,'minus') )
    ot = struct('method', 'uminus', 'numArgs', 1, 'center', tree, ...
        'diffOrder', tree.diffOrder);
else
    ot = tree;
end
end