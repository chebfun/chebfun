function newTree = expandTree(tree, maxOrder)
% EXPANDTREE. Convert expressions like 5*(diff(u) + u) to 5*diff(u) + 5*u.

% If the input is not a tree, or its diffOrder is less than the
% maxOrder we consider, don't need to do anything.
if ( ~isstruct(tree) || tree.height <= 1 || tree.diffOrder < maxOrder )
    newTree = tree;
elseif ( tree.numArgs == 1 )
    newTree = expandTree(tree.center, maxOrder);
else
    if ( ( strcmp(tree.method, 'diff') && tree.height <= 1 ) )
        newTree = tree;
    elseif ( any(strcmp(tree.method, {'plus','minus'})) )
        % Expand the tree recursively:
        tree.left = treeVar.expandTree(tree.left, maxOrder);
        tree.right = treeVar.expandTree(tree.right, maxOrder);
        newTree = tree;
    elseif ( strcmp(tree.method, 'times') )
        if ( ~isstruct(tree.left) && tree.right.height <= 1 )
            newTree = tree;
            return;
        elseif ( ~isstruct(tree.right) && tree.left.height <= 1 )
            newTree = tree;
            return;
        end
        
        if ( ~isstruct(tree.left) )
            leftNumArgs = 0;
        else
            leftNumArgs = tree.left.numArgs;
        end
        
        if ( ~isstruct(tree.right) )
            rightDiffOrder = 0;
            rightNumArgs = 0;
        else
            rightNumArgs = tree.right.numArgs;
        end
        
        if ( leftNumArgs == 0 || ...
                (strcmp(tree.left.method,'diff') && tree.left.height <= 1 )  )
            leftTree = tree.left;
            leftDiffOrder = 0;
            splitLeft = false;
        elseif ( leftNumArgs == 1 )
            leftTree = expandTree(tree.left, maxOrder);
            leftDiffOrder = leftTree.diffOrder;
            splitLeft = false;
        else
            leftLeft = treeVar.expandTree(tree.left.left, maxOrder);
            leftRight = treeVar.expandTree(tree.left.right, maxOrder);
            splitLeft = true;
        end
        
        if ( rightNumArgs == 0 || ...
                (strcmp(tree.right.method,'diff') && tree.right.height <= 1 )  )
            rightTree = tree.right;
            rightDiffOrder = 0;
            splitRight = false;
        elseif ( rightNumArgs == 1 )
            rightTree = expandTree(tree.right, maxOrder);
            rightDiffOrder = rightTree.diffOrder;
            splitRight = false;
        else
            rightLeft  = treeVar.expandTree(tree.right.left, maxOrder);
            rightRight = treeVar.expandTree(tree.right.right, maxOrder);
            splitRight = true;
        end
        
        if ( splitLeft && ~splitRight )
            newLeft = struct('method','times', 'numArgs', 2, ...
                'left', leftLeft, ...
                'right', rightTree,...
                'diffOrder', max(leftLeft.diffOrder, rightDiffOrder));
            newRight = struct('method','times', 'numArgs', 2, ...
                'left', leftRight, ...
                'right', rightTree, ...
                'diffOrder', max(leftRight.diffOrder, rightDiffOrder));
        elseif ( leftNumArgs < 2 && rightNumArgs == 2 )
            newLeft = struct('method','times', 'numArgs', 2, ...
                'left', leftTree, ...
                'right', rightLeft,...
                'diffOrder', max(leftDiffOrder, rightLeft.diffOrder));
            newRight = struct('method','times', 'numArgs', 2, ...
                'left', leftTree, ...
                'right', rightRight, ...
                'diffOrder', max(leftDiffOrder, rightRight.diffOrder));
        end
        
%         if ( ( isstruct(tree.rightleft) ) && ( tree.left.numArgs == 2 ) )
%             leftLeft = treeVar.expandTree(tree.left.left, maxOrder);
%             leftRight = treeVar.expandTree(tree.left.right, maxOrder);
%             leftTwoArgs = 1;
%         end
%         

        
        newTree = struct('method', 'plus', 'numArgs', 2, ...
            'left', newLeft, 'right', newRight, ...
            'diffOrder', max(newLeft.diffOrder, newRight.diffOrder));
    end
end
end