function [newTree, derTree] = splitTree(tree, maxOrder)
%SPLITTREE     Split syntax trees into derivative part and non-derivative part

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Isolate where the highest order derivative appears.

% If the highest order derivative is lower than the maximum order of
% the overall expression, we don't need to investigate it further.

% Find what variables we have actually computed derivatives of.
diffVar = maxOrder > 0;

% If our tree don't have any derivatives of highest order that appear in the
% problem, we can return the subtree. But only need to check for those variables
% that we actually differentiate with respect to, e.g. if maxOrder = 0, we can
% return the subtree as well.
if ( ~isstruct(tree) )
    newTree = tree;
    derTree = [];
    return
end

treeDiffOrder = tree.diffOrder;
if ( all(treeDiffOrder(diffVar) < maxOrder(diffVar)) )
    newTree = tree;
    derTree = [];
    return
end

switch tree.numArgs
    case 1
        % We're dealing with a unary operator,
        [newTree, derTree] = splitTree(tree, maxOrder);
        
    case 2
        if ( any(strcmp(tree.method, {'diff','times'})) )
            newTree = [];
            derTree = tree;
        else
            [newTreeLeft, derTreeLeft] = ...
                treeVar.splitTree(tree.left, maxOrder);
            [newTreeRight, derTreeRight] = ...
                treeVar.splitTree(tree.right, maxOrder);
            
            if ( isempty(newTreeLeft) )
                newTree = oneTreeFromRight(newTreeRight, tree.method);
            elseif (isempty(newTreeRight) )
                newTree = newTreeLeft;
            else
                if ( ~isstruct(newTreeLeft) )
                    newDiffOrder = newTreeRight.diffOrder;
                    newHeight = newTreeRight.height;
                elseif ( ~isstruct(newTreeRight) )
                    newDiffOrder = newTreeLeft.diffOrder;
                    newHeight = newTreeLeft.height;
                else
                    newDiffOrder = max(newTreeLeft.diffOrder, ...
                        newTreeRight.diffOrder);
                    newHeight = max(newTreeLeft.height, newTreeRight.height);
                end
                
                newTree = struct('method', tree.method, ...
                    'numArgs', tree.numArgs, ...
                    'left', newTreeLeft, 'right', newTreeRight, ...
                    'diffOrder', newDiffOrder, ...
                    'height', newHeight);
            end
            
            if ( isempty(derTreeLeft) )
                derTree = oneTreeFromRight(derTreeRight, tree.method);
            elseif (isempty(derTreeRight) )
                derTree = derTreeLeft;
            else
                derTree = struct('method', tree.method, ...
                    'numArgs', tree.numArgs, ...
                    'left', derTreeLeft, 'right', derTreeRight, ...
                    'diffOrder', max(derTreeLeft.diffOrder, ...
                    derTreeRight.diffOrder), ...
                    'height', max(derTreeLeft.height, derTreeRight.height));
            end
            
        end
end
end

function ot = oneTreeFromRight(tree, operator)
if ( strcmp(operator,'minus') )
    if ( isstruct(tree) )
        ot = struct('method', 'uminus', 'numArgs', 1, 'center', tree, ...
            'diffOrder', tree.diffOrder, 'ID', tree.ID, ...
            'height', tree.height + 1);
    else
        ot = -tree;
    end
else
    ot = tree;
end
end