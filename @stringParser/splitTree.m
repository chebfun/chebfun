function treeOut = splitTree(treeIn)
%SPLITTREE   Split a syntax tree.
%   TREEOUT = SPLITTREE(TREEIN) returns a new syntax TREEOUT, obtained by
%   replacing a potential comma or an equal sign at the top of TREEIN. This
%   method looks for individual expressions separated by commas (e.g. u(1) = 0,
%   u'(3) = 4), and in each expression, converts the = into a - so that we end
%   up with expressions corresponding to u(1)-0, u'(3)-4.
%
% See also: STRINGPARSER/SPLITTREEEIG, STRINGPARSER/SPLITTREEPDE.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Initialize.
treeOut = treeIn;

% If the top-center entry is a comma, we have separate expressions. Go
% recursively through them to convert = into -. Since commas have the lowest
% operator precedence of all operators, we know that they can only appear at the
% top of syntax trees.
%
% Here we also check for =. If we reach the elseif condition, we know we haven't
% got a match for a comma, so that we're in the separate expression case. Since
% = has the second lowest operator precedence after commas, they can as well
% only appear at the top of each individual expression syntax trees.

for k = 1:numel(treeOut)
    % The center node.
    treeCenter = treeIn(k).center;
    
    % Check for commas and =. Do nothing otherwise.
    if ( strcmp(treeCenter{2}, 'COMMA') )
        % If we have a COMMA, which separates expressions, we call the method
        % recursively on the left and the right syntax trees.
        treeOut(k).left = stringParser.splitTree(treeOut.left);
        treeOut(k).right = stringParser.splitTree(treeOut.right);
        
    elseif ( strcmp(treeCenter{2}, 'OP=') )
        % If we have an = sign, replace it with a -, so that an expression such
        % as 'u(3) = 4' gets converted to 'u(3) - 4'.
        treeOut(k).center = {'-', 'OP-'};
    end
end

end
