function treeOut = bivariate(leftTree, rightTree, method, type)
%BIVARIATE   Construct a syntax tree for bivariate functions.
%
% The BIVARIATE() method is called by most bivariate methods of the TREEVAR
% class, such as minus, plus and rtimes. It constructs the corresponding syntax
% tree, where the appropriate method is at the top, and the syntax trees that
% the method operates on become the left and right nodes.
%
% Calling sequence:
%   TREEOUT = BIVARIATE(LEFTTREE, RIGHTTREE, METHOD, TYPE)
%
% The inputs to this method are:
%   LEFTTREE:   A MATLAB struct for representing the syntax tree of the 
%               mathematical expressions that is the left argument that METHOD
%               operates on.
%   RIGHTTREE:  A MATLAB struct for representing the syntax tree of the 
%               mathematical expressions that is the right argument that METHOD
%               operates on.
%   METHOD:     The name of the method that invokes the call to BIVARIATE().
%   TYPE:       An integer that indicated which of the arguments to METHOD were
%               TREEVAR objects (as opposed to scalars and CHEBFUN objects). The
%               possible values are:
%                   0: Only the left argument was a TREEVAR,
%                   1: Only the right argument was a TREEVAR,
%                   2: Both left and right arguments were TREEVAR objects.   
%
% The output of this method is:
%   TREEOUT: A MATLAB struct that represents the syntax tree of the mathematical
%            expression, obtained once METHOD has operated on LEFTTREE and
%            RIGHTTREE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Construct a new MATLAB struct to be returned.

isPM = any( strcmp(method, {'plus', 'minus'}) );

if ( type == 2 )
    % Both LEFTTREE and RIGHTTREE were TREEVAR objects.
    treeOut = struct('method', method, 'numArgs', 2, ...
        'left', leftTree, 'right', rightTree, ...
        'diffOrder', max(leftTree.diffOrder, rightTree.diffOrder), ...
        'ID', leftTree.ID | rightTree.ID, ...
        'height', max(leftTree.height, rightTree.height) + 1, ...
        'hasTerms', isPM || leftTree.hasTerms || rightTree.hasTerms);
elseif ( type == 1 )
    % Only RIGHTTREE was a TREEVAR.
    treeOut = struct('method', method, 'numArgs', 2, ...
        'left', leftTree, 'right', rightTree, ...
        'diffOrder', rightTree.diffOrder, ...
        'height', rightTree.height + 1, ...
        'ID', rightTree.ID, ...
        'hasTerms', isPM || rightTree.hasTerms);
else
    % Only LEFTTREE was a TREEVAR.
    treeOut = struct('method', method, 'numArgs', 2, ...
        'left', leftTree, 'right', rightTree, ...
        'diffOrder', leftTree.diffOrder, ...
        'height', leftTree.height + 1, ...
        'ID', leftTree.ID, ...
        'hasTerms', isPM || leftTree.hasTerms);
end

end
