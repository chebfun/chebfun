function [newTree, pdeSign] = splitTreePDE(treeIn)
%SPLITTREEPDE   Split a syntax tree, specific to PDE problems.
%   [NEWTREE, PDESIGN] = SPLITTREPDE(TREEIN) goes through the syntax tree
%   TREEIN, and isolates the part where the PDE variable, e.g. u_t, appears.
%   splits the syntax tree TREEIN into two trees, NEWTREE and LAMBDATREE.
%   LAMBDATREE contains the syntax tree that the eigenvalue parameter LAMBDA
%   appears in, NEWTREE contains the other part of TREEIN. The value of
%   LAMBDASIGN corresponds to the sign in front of LAMBDA in the original syntax
%   tree TREEIN.
%
% See also: STRINGPARSER/SPLITTREE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Begin by replacing the subtree which contains the pde_variable with a 0
[newTree, ignored, pdeSign] = findPDE(treeIn, 1);

% Do the basic splitting (converting = into -) in newTree
newTree = stringParser.splitTree(newTree);

end

function [newTree, pdeTree, pdeSign] = findPDE(treeIn, pdeSign)
%FINDPDE    Find where the PDE variable (e.g. u_t) appears in TREEIN.

% Initialization.
newTree = treeIn;
% Create some empty trees.
leftEmpty = 1;
rightEmpty = 1;
pdeTree = [];
pdeTreeLeft = [];
pdeTreeRight = [];
treeCenter = treeIn.center;

% Go recursively through the tree on the left.
if ( isfield(treeIn, 'left') )
    [newLeft, pdeTreeLeft, pdeSign] = findPDE(treeIn.left, pdeSign);
    newTree.left = newLeft;
    leftEmpty = 0;
end

% Go recursively through the tree on the right.
if ( isfield(treeIn, 'right') )
    [newRight, pdeTreeRight, pdeSign] = findPDE(treeIn.right, pdeSign);
    newTree.right = newRight;
    rightEmpty = 0;
end

% Return a new pdeTree. If the operator in the center of treeIn is a *,
% we want to return the whole treeIn (e.g. when we see 1*u_t). If not,
% we return the latest pdeTree (e.g. when we see u_t+1).
%
% Start by looking through the left tree.
if ( ~isempty(pdeTreeLeft) ) 
    if ( strcmp(treeCenter{2}, 'OP*') )
        pdeTree = treeIn;
    else
        pdeTree = pdeTreeLeft;
    end
end

% Now look through the right tree.
if ( ~isempty(pdeTreeRight) )
    if ( strcmp(treeCenter, 'OP*') )
        pdeTree = treeIn;
    % If we have a =, and u_t is on the right, we need to switch signs on
    % the pdeTree.
    elseif ( strcmp(treeCenter{2}, 'OP=') )
        disp('PDE on right')
        pdeSign = -1*pdeSign;
        pdeTree = pdeTreeRight;    
    % If we have a -, and we have u_t on the right we need to switch signs
    % on the pdeTree.
    elseif ( strcmp(treeCenter{2}, 'OP-') || strcmp(treeCenter{2}, 'UN-') )
        pdeSign = -1*pdeSign;
        pdeTree = pdeTreeRight;
    else
        pdeTree = pdeTreeRight;
    end
end

% Both left and right trees are empty. We must be at a leaf! Replace the PDEVAR
% with a 0.
if ( leftEmpty && rightEmpty )
    % We encounter a PDE variable. Replace it by a zero.
    if ( strcmp(treeCenter{2}, 'PDEVAR') )
        pdeTree = newTree;
        newTree = struct('center', {{'0', 'NUM'}});
    end
end

end
