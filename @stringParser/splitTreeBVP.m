function treeOut = splitTreeBVP(treeIn)
% SPLITTREE_BVP Split a syntax tree (replace = with -) for a BVP 

% TODO:  Documentation.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/chebfun/ for Chebfun information.

% Only need to do the basic splitting in case of BVPs
treeOut = stringParser.splitTree_commas_equalSigns(treeIn);