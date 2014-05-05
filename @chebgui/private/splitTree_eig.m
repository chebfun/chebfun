function [newTree lambdaTree lambdaSign] = splitTree_eig(guifile,treeIn)
% SPLITTREE_BVP Split a syntax tree (replace = with -) for a EIG problem

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Begin by finding the subtree which contains lambda

[newTree lambdaTree lambdaSign] = findLambda(treeIn,1);

% Do the basic splitting (converting = into -) in newTree
newTree = splitTree_commas_equalSigns(guifile,newTree);

end

function [newTree lambdaTree lambdaSign] = findLambda(treeIn,lambdaSign)

newTree = treeIn;
leftEmpty = 1; rightEmpty = 1;
lambdaTree = []; lambdaTreeLeft = []; lambdaTreeRight = [];
treeCenter = treeIn.center;

if isfield(treeIn,'left')
    [newLeft lambdaTreeLeft lambdaSign] = findLambda(treeIn.left,lambdaSign);
    newTree.left = newLeft;
    leftEmpty = 0;
end

if isfield(treeIn,'right')
    [newRight lambdaTreeRight lambdaSign] = findLambda(treeIn.right,lambdaSign);
    newTree.right = newRight;
    rightEmpty = 0;
end

% Return a new lambdaTree. If the operator in the center of treeIn is a *,
% we want to return the whole treeIn (e.g. when we see lambda*u). If not,
% we return the latest lambdaTree (e.g. when we see lambda*u+1, and the +
% is the operator of the current tree).
if ~isempty(lambdaTreeLeft)
    if strcmp(treeCenter{2},'OP*')
        lambdaTree = treeIn;
    % If we have a =, and lambda is on the left, we need to add a unary -
    % switch signs on the lambdaTree.
    elseif strcmp(treeCenter{2},'OP=') 
        lambdaSign = -1*lambdaSign;
        lambdaTree = lambdaTreeLeft;
    else
        lambdaTree = lambdaTreeLeft;
    end
end
if ~isempty(lambdaTreeRight)
    if strcmp(treeCenter{2},'OP*')
        lambdaTree = treeIn;
    elseif strcmp(treeCenter{2},'UN-')
        lambdaTree= struct('center',{{'-','UN-'}},'right', lambdaTreeRight);
    % If we have a binary -, and lambda is on the right, we need to switch
    % signs on the lambdaTree. Add a unary minus at the top of the tree.
    elseif strcmp(treeCenter{2},'OP-')
        lambdaTree= struct('center',{{'-','UN-'}},'right', lambdaTreeRight);
    else
        lambdaTree = lambdaTreeRight;
    end
end

if leftEmpty && rightEmpty % Must be at a leaf
    % We encounter a lambda variable. Replace it by a zero.
   if strcmp(treeCenter{2},'LAMBDA')
       lambdaTree = newTree;
       newTree = struct('center',{{'0', 'NUM'}});
   end 
end

end


