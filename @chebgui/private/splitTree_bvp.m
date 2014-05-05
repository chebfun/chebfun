function [treeOut] = splitTree_bvp(guifile,treeIn)
% SPLITTREE_BVP Split a syntax tree (replace = with -) for a BVP 

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Only need to do the basic splitting in case of BVPs
treeOut = splitTree_commas_equalSigns(guifile,treeIn);

% % If the top-center entry is a =, we need to convert that into a - in case
% % we're working with a ODE. Otherwise, do nothing.
% for k = 1:numel(treeOut)
%     treeCenter = treeIn(k).center;
%     if strcmp(treeCenter{2},'OP=')
%         treeOut(k).center = {'-', 'OP-'};
%     end
% end