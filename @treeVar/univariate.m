function treeOut = univariate(treeIn, method)
%UNIVARIATE   Construct a syntax tree for univariate functions.
%   The UNIVARIATE method is called by most univariate methods of the TREEVAR
%   class, such as COS, EXP and SIN. It constructs the corresponding syntax
%   tree, where the appropriate method is at the top and the syntax tree up to
%   that point becomes the center node.
%
%   Calling sequence:
%      TREEOUT = UNIVARIATE(TREEIN, METHOD)
%
%   The inputs to this method are:
%      TREEIN: A MATLAB struct for representing the syntax tree of a 
%              mathematical expressions that the METHOD operates on.
%      METHOD: The name of the method that invokes the call to UNIVARIATE.
%
%   The output of this method is:
%      TREEOUT: A MATLAB struct that represents the syntax tree of the 
%                mathematical expression, obtained once METHOD has operated on 
%                TREEIN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Construct a new MATLAB struct to be returned.
treeOut = struct('method', method, ...
    'numArgs', 1, 'center', treeIn, ...
    'diffOrder', treeIn.diffOrder, ...
    'height', treeIn.height + 1, ...
    'ID', treeIn.ID, ...
    'hasTerms', treeIn.hasTerms);
end
