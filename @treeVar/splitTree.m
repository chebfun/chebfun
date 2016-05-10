function [newTree, derTree] = splitTree(tree, maxOrder)
%SPLITTREE   Split syntax trees into derivative part and non-derivative part
%   Calling sequence:
%      [NEWTREE, DERTREE] = SPLITTREE(TREE, MAXORDER)
%   where the inputs are:
%      TREE:       The syntax tree to be split.
%      MAXORDER:   A vector that contains the the maximum differential order of
%                  each variable that appears in the problem.
%   and the outputs are
%      NEWTREE:    A syntax tree which describes the factor in which the highest
%                  order derivative appears. E.g. if we split the expression
%                  5*diff(u) + sin(u), NEWTREE is the syntax tree corresponding
%                  to 5*diff(u).
%      DERTREE:    A syntax tree which describes the remaining factors not 
%                  included in NEWTREE. In the example above, DERTREE is the 
%                  syntax tree corresponding to sin(u).

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Find what variables we have actually computed derivatives of.
diffVar = maxOrder > 0;

% If our input is not a syntax tree, we can be certain we don't need to split it
% (this can e.g. happen if TREE is a CHEBFUN or a scalar).
if ( ~isstruct(tree) )
    newTree = tree;
    derTree = [];
    return
end

% Find the diffOrders in the input tree.
treeDiffOrder = tree.diffOrder;

% If our tree doesn't have any derivatives of highest order that appear in the
% problem, we can return the subtree. But only need to check for those variables
% that we actually differentiate with respect to, e.g. if maxOrder = 0, we can
% return the subtree as well.
if ( all(treeDiffOrder(diffVar) < maxOrder(diffVar)) )
    newTree = tree;
    derTree = [];
    return
end

% Split the tree recursively, depending on how many inputs the operator takes.
switch tree.numArgs
    case 1
        % We're dealing with a unary operator.
        [newTree, derTree] = treeVar.splitTree(tree.center, maxOrder);
        
    case 2
        if ( any(strcmp(tree.method, {'diff', 'times', 'rdivide'})) )
            % We have already been through the expandTree() method, which
            % guarantees that there are no terms left in the tree which include
            % highest order derivatives along with other expressions in which
            % the unknown variable appears, e.g. 5*(diff(u) + u). Since we don't
            % reach this point in the code unless the input TREE contains the
            % maximum diffOrder of the problem, we can safely assume that the
            % input TREE is the derivative tree, and return it.
            newTree = [];
            derTree = tree;
        else
            % We're at + or -, split the left and right children trees
            % recursively.
            [newTreeLeft, derTreeLeft] = ...
                treeVar.splitTree(tree.left, maxOrder);
            [newTreeRight, derTreeRight] = ...
                treeVar.splitTree(tree.right, maxOrder);
            
            if ( isempty(newTreeLeft) )
                % We only had a derivative part on the left. Thus, the
                % non-derivative parts will only consist of the non-derivative
                % part of the right tree. However, if the method of our current
                % tree is MINUS(), we must add a UMINUS() in front of
                % NEWTREERIGHT before we can return it. That will be done in the
                % ONETREEFROMRIGHT() method. Simply put, we're converting a
                % binary minus between an empty tree on left and a non-empty
                % tree on the right to a UMINUS on the non-empty tree.
                newTree = oneTreeFromRight(newTreeRight, tree.method);
            elseif (isempty(newTreeRight) )
                % We only had a derivative part on the right, so the
                % non-derivative part only consists of the left part.
                newTree = newTreeLeft;
            elseif ( ~isstruct(newTreeLeft) && ~isstruct(newTreeRight) )
                % Both left and right trees were actually a CHEBFUN or scalar,
                % combine them in a single CHEBFUN/scalar and return as the
                % NEWTREE.
                newTree = eval([tree.method, '(newTreeLeft, newTreeRight)']);
            else
                % Had non derivative parts on both left and right, potentially a
                % combination of CHEBFUN/scalars and syntax trees. Need to
                % combine them.
                if ( ~isstruct(newTreeLeft) )
                    % Left tree was actually a CHEBFUN or scalar, so the new
                    % diffOrders and height only depend on the right tree.
                    newDiffOrder = newTreeRight.diffOrder;
                    newHeight = newTreeRight.height;
                elseif ( ~isstruct(newTreeRight) )
                    % Right tree was actually a CHEBFUN or scalar, so the new
                    % diffOrders and height only depend on the right tree.
                    newDiffOrder = newTreeLeft.diffOrder;
                    newHeight = newTreeLeft.height;
                else
                    % Both left and right tree were actually syntax trees, so
                    % find the maximum heights and diffOrders.
                    newDiffOrder = max(newTreeLeft.diffOrder, ...
                        newTreeRight.diffOrder);
                    newHeight = max(newTreeLeft.height, newTreeRight.height);
                end
                
                % Construct a new syntax tree from the new left and right child
                % trees (the non-derivative parts).
                newTree = struct('method', tree.method, ...
                    'numArgs', tree.numArgs, ...
                    'left', newTreeLeft, 'right', newTreeRight, ...
                    'diffOrder', newDiffOrder, ...
                    'height', newHeight);
            end
            
            if ( isempty(derTreeLeft) )
                % Left child tree only consisted of non-derivative part. Like
                % above, return the right tree, but add a UMINUS() on top of it
                % if needed.
                derTree = oneTreeFromRight(derTreeRight, tree.method);
            elseif (isempty(derTreeRight) )
                % Right child tree only consisted of non-derivative part.
                derTree = derTreeLeft;
            else
                % Combine the left and right derivative subtrees.
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
%ONETREEFROMRIGHT   Convert binary minus to uminus when left tree is empty.
%   Add a UMINUS in front of a syntax tree in the case where the left child tree
%   of MINUS() only consisted of the derivative part. Simply put, we're 
%   converting a binary minus between an empty tree on left and a non-empty tree 
%   on the right to a UMINUS on the non-empty tree.
if ( strcmp(operator, 'minus') )
    % Only need to worry if the original operator was a -.
    if ( isstruct(tree) )
        % If the input is a syntax tree, add a UMINUS on the top of it.
        ot = struct('method', 'uminus', 'numArgs', 1, 'center', tree, ...
            'diffOrder', tree.diffOrder, 'ID', tree.ID, ...
            'height', tree.height + 1);
    else
        % Input was a CHEBFUN or a scalar, can simply negate it.
        ot = -tree;
    end
else
    % Didn't have a MINUS(), so can simply return the input.
    ot = tree;
end
end
