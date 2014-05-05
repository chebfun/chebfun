function [treeOut] = splitTree_commas_equalSigns(guifile,treeIn)
% SPLITTREE_COMMAS_EQUALSIGNS This is the basic method for splitting syntax
% trees. It looks for individual expressions separated by commas (e.g. u(1)
% = 0, u'(3) = 4), and in each expression, converts the = into a - so that
% we end up with expressions corresponding to u(1)-0, u'(3)-4

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

treeOut = treeIn;

% If the top-center entry is a comma, we have separate expressions. Go
% recursively through them to convert = into -. Since commas have the
% lowest operator precedence of all operators, we know that they can only
% appear at the top of syntax trees. 
%
% Here we also check for =. If we reach the elseif condition, we know we
% haven't got a match for a comma, so that we're in the separate expression
% case. Since = has the second lowest operator precedence after commas,
% they can as well only appear at the top of each individual expression
% syntax trees.
for k = 1:numel(treeOut)
    treeCenter = treeIn(k).center;
    % Check for commas and =. Do nothing otherwise.
    if strcmp(treeCenter{2},'COMMA')
        treeOut(k).left = splitTree_commas_equalSigns(guifile,treeOut.left);
        treeOut(k).right = splitTree_commas_equalSigns(guifile,treeOut.right);
    elseif strcmp(treeCenter{2},'OP=')
        treeOut(k).center = {'-', 'OP-'};
    end
end