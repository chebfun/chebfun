function [newTree pdeSign] = splitTree_pde(guifile,treeIn)
% SPLITTREE_PDE Split a syntax tree (replace = with -) for a PDE

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Begin by replacing the subtree which contains the pde_variable with a 0

[newTree pdeTree pdeSign] = findPDE(treeIn,1);

% Do the basic splitting (converting = into -) in newTree
newTree = splitTree_commas_equalSigns(guifile,newTree);

end

function [newTree pdeTree pdeSign] = findPDE(treeIn, pdeSign)

newTree = treeIn;
leftEmpty = 1; rightEmpty = 1;
pdeTree = []; pdeTreeLeft = []; pdeTreeRight = [];
treeCenter = treeIn.center;

if isfield(treeIn,'left')
    [newLeft pdeTreeLeft pdeSign] = findPDE(treeIn.left,pdeSign);
    newTree.left = newLeft;
    leftEmpty = 0;
end

if isfield(treeIn,'right')
    [newRight pdeTreeRight pdeSign] = findPDE(treeIn.right, pdeSign);
    newTree.right = newRight;
    rightEmpty = 0;
end


% Return a new pdeTree. If the operator in the center of treeIn is a *,
% we want to return the whole treeIn (e.g. when we see 1*u_t). If not,
% we return the latest pdeTree (e.g. when we see u_t+1).
if ~isempty(pdeTreeLeft)
    if strcmp(treeCenter{2},'OP*')
        pdeTree = treeIn;
    else
        pdeTree = pdeTreeLeft;
    end
end
if ~isempty(pdeTreeRight)
    if strcmp(treeCenter,'OP*')
        pdeTree = treeIn;
    % If we have a =, and u_t is on the right, we need to switch signs on
    % the pdeTree.
    elseif strcmp(treeCenter{2},'OP=')
        disp('PDE on right')
        pdeSign = -1*pdeSign;
        pdeTree = pdeTreeRight;    
    % If we have a -, and we have u_t on the right we need to switch signs
    % on the pdeTree.
    elseif strcmp(treeCenter{2},'OP-') || strcmp(treeCenter{2},'UN-')
        pdeSign = -1*pdeSign;
        pdeTree = pdeTreeRight;
    else
        pdeTree = pdeTreeRight;
    end
end


if leftEmpty && rightEmpty % Must be at a leaf.
    % We encounter a PDE variable. Replace it by a zero.
    if strcmp(treeCenter{2},'PDEVAR')
        pdeTree = newTree;
        newTree = struct('center',{{'0', 'NUM'}});
    end
end
end
