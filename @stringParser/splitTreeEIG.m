function [newTree, lambdaTree, lambdaSign] = splitTreeEIG(treeIn)
%SPLITTREEEIG    Split a syntax tree, specific to EIG problems.
%   [NEWTREE, LAMBDATREE, LAMBDASIGN] = SPLITTREEEIG(TREEIN) splits the syntax
%   tree TREEIN into two trees, NEWTREE and LAMBDATREE. LAMBDATREE contains the
%   syntax tree that the eigenvalue parameter LAMBDA appears in, NEWTREE contains
%   the other part of TREEIN. The value of LAMBDASIGN corresponds to the sign in
%   front of LAMBDA in the original syntax tree TREEIN.
%
% See also: STRONGPARSER/SPLITTREE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Begin by finding the subtree which contains the eigenvalue parameter LAMBDA.
[newTree, lambdaTree, lambdaSign] = findLambda(treeIn, 1);

% Do the basic splitting (converting = into -) in newTree
newTree = stringParser.splitTree(newTree);

end

function [newTree, lambdaTree, lambdaSign] = findLambda(treeIn, lambdaSign)
%FINDLAMBDA   Find where the eigenvalue parameter LAMBDA appears in TREEIN.

% Initialization.
newTree = treeIn;

% Start with empty trees.
leftEmpty = 1;
rightEmpty = 1;
lambdaTree = [];
lambdaTreeLeft = [];
lambdaTreeRight = [];
treeCenter = treeIn.center;

% Check whether we have a syntax tree to the left. If so, look recursively for
% LAMBDA in there.
if ( isfield(treeIn, 'left') )
    [newLeft, lambdaTreeLeft, lambdaSign] = findLambda(treeIn.left, lambdaSign);
    newTree.left = newLeft;
    leftEmpty = 0;
end

% Check whether we have a syntax tree to the right. If so, look recursively for
% LAMBDA in there.
if ( isfield(treeIn, 'right') )
    [newRight, lambdaTreeRight, lambdaSign] = findLambda(treeIn.right, ...
        lambdaSign);
    newTree.right = newRight;
    rightEmpty = 0;
end

% Return a new lambdaTree. If the operator in the center of treeIn is a *, we
% want to return the whole treeIn (e.g. when we see lambda*u). If not, we return
% the latest lambdaTree (e.g. when we see lambda*u+1, and the + is the operator
% of the current tree).

% Start with the left tree if it is nonempty:
if ( ~isempty(lambdaTreeLeft) )
    if ( strcmp(treeCenter{2}, 'OP*') )
        lambdaTree = treeIn;
    % If we have a =, and lambda is on the left, we need to add a unary -
    % switch signs on the lambdaTree.
    elseif ( strcmp(treeCenter{2}, 'OP=') )
        lambdaSign = -1*lambdaSign;
        lambdaTree = lambdaTreeLeft;
    else
        lambdaTree = lambdaTreeLeft;
    end
end

% Then go through the right tree:
if ( ~isempty(lambdaTreeRight) )
    if ( strcmp(treeCenter{2},'OP*') )
        lambdaTree = treeIn;
    elseif ( strcmp(treeCenter{2}, 'UN-') )
        lambdaTree= struct('center', {{'-','UN-'}}, 'right', lambdaTreeRight);
    % If we have a binary -, and lambda is on the right, we need to switch
    % signs on the lambdaTree. Add a unary minus at the top of the tree.
    elseif ( strcmp(treeCenter{2}, 'OP-') )
        lambdaTree= struct('center', {{'-','UN-'}}, 'right', lambdaTreeRight);
    else
        lambdaTree = lambdaTreeRight;
    end
end

% Both left and right trees are empty. We must be at a leaf!
if ( leftEmpty && rightEmpty )
   % We encounter a lambda variable. Replace it by a zero.
   if ( strcmp(treeCenter{2}, 'LAMBDA') )
       lambdaTree = newTree;
       newTree = struct('center', {{'0', 'NUM'}});
   end 
end

end
