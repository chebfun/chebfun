function newTree = expandTree(tree, maxOrder)
%EXPANDTREE    Convert expressions like 5*(diff(u) + u) to 5*diff(u) + 5*u.
%
% The role of the EXPANDTREE method is to split up syntax trees, so that
% expressions involving the highest derivative appearing in the equation stand
% alone, that is, not inside parenthesis. For example, EXPANDTREE will convert
% expressions like 5*(diff(u) + u) to 5*diff(u) + 5*u. This is necessary in
% order to be able to put the correct expression on the right-hand side when
% calling MATLAB's built in solvers.
%
% Calling sequence:
%   NEWTREE = EXPANDTREE(TREE, MAXORDER)
% where the inputs are:
%   TREE:       The syntax tree we want to split.
%   MAXORDER:   A vector containing the maximum differential order of each
%               variable that appear in the problem under consideration.
% and the output is:
%   NEWTREE:    A syntax tree where the highest order deriative has been taken
%               out of any parenthesis where other variables with lower
%               differential order appear.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( ~isstruct(tree) || tree.height <= 1 || all(tree.diffOrder < maxOrder) )
    % If the input is not a tree, it is a very short tree, or its diffOrder is
    % less than the maxOrder we consider, don't need to do anything.
    newTree = tree;
    
elseif ( tree.numArgs == 1 )
    % If we have a unary operator, we expand the center node recursively.
    newTree = expandTree(tree.center, maxOrder);
    
else
    % We are dealing with a binary operator.
    if ( ( strcmp(tree.method, 'diff') && tree.height <= 1 ) )
        % We've hit diff(u, k) -- no need to expand the tree further.
        newTree = tree;
        
    elseif ( any(strcmp(tree.method, {'plus', 'minus'})) )
        % We're at a + or a -, so expand the syntax tree recursively:
        tree.left = treeVar.expandTree(tree.left, maxOrder);
        tree.right = treeVar.expandTree(tree.right, maxOrder);
        newTree = tree;
        
    elseif ( strcmp(tree.method, 'times') )
        if ( ~isstruct(tree.left) && tree.right.height <= 1 )
            % We're at 5*u or x.*u.
            newTree = tree;
            return;
        elseif ( ~isstruct(tree.right) && tree.left.height <= 1 )
            % We're at u*5 or u.*x.
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
            leftHeight = 0;
            splitLeft = false;
        elseif ( leftNumArgs == 1 )
            leftTree = treeVar.expandTree(tree.left, maxOrder);
            leftDiffOrder = leftTree.diffOrder;
            leftHeight = leftTree.height;
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
            rightHeight = 0;
            splitRight = false;
        elseif ( rightNumArgs == 1 )
            rightTree = treeVar.expandTree(tree.right, maxOrder);
            rightDiffOrder = rightTree.diffOrder;
            rightHeight = 0;
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
                'diffOrder', max(leftLeft.diffOrder, rightDiffOrder), ...
                'height', max(leftLeft.height, rightHeight) + 1);
            newRight = struct('method','times', 'numArgs', 2, ...
                'left', leftRight, ...
                'right', rightTree, ...
                'diffOrder', max(leftRight.diffOrder, rightDiffOrder) + 1, ...
                'height', max(leftRight.height, rightHeight) + 1);
        elseif ( ~splitLeft && splitRight )
            newLeft = struct('method','times', 'numArgs', 2, ...
                'left', leftTree, ...
                'right', rightLeft,...
                'diffOrder', max(leftDiffOrder, rightLeft.diffOrder), ...
                'height', max(leftHeight, rightLeft.height) + 1);
            newRight = struct('method','times', 'numArgs', 2, ...
                'left', leftTree, ...
                'right', rightRight,...
                'diffOrder', max(leftDiffOrder, rightRight.diffOrder), ...
                'height', max(leftHeight, rightRight.height) + 1);
        elseif ( splitLeft && splitRight )
            disp('Split both?')
        elseif ( ~splitLeft && ~splitRight )
            newLeft = leftTree;
            newRight = rightTree;
        elseif ( leftNumArgs < 2 && rightNumArgs == 2 )
            newLeft = struct('method','times', 'numArgs', 2, ...
                'left', leftTree, ...
                'right', rightLeft,...
                'diffOrder', max(leftDiffOrder, rightLeft.diffOrder), ...
                'height', max(leftHeight, rightLeft.height) + 1);
            newRight = struct('method','times', 'numArgs', 2, ...
                'left', leftTree, ...
                'right', rightRight, ...
                'diffOrder', max(leftDiffOrder, rightRight.diffOrder), ...
                'height', max(leftHeight, rightRight.height) + 1);
        end
         
%         if ( ( isstruct(tree.rightleft) ) && ( tree.left.numArgs == 2 ) )
%             leftLeft = treeVar.expandTree(tree.left.left, maxOrder);
%             leftRight = treeVar.expandTree(tree.left.right, maxOrder);
%             leftTwoArgs = 1;
%         end
%                 
        newTree = struct('method', 'plus', 'numArgs', 2, ...
            'left', newLeft, 'right', newRight, ...
            'diffOrder', max(newLeft.diffOrder, newRight.diffOrder), ...
            'height',max(newLeft.height, newRight.height) + 1);
    end
end
end