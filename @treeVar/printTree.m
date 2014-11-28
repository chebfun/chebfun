function s = printTree(tree, indentStr)
%PRINTTREE   Print a syntax tree.
%   Calling sequence:
%      S = TREEVAR.PRINTTREE(TREE, INDSTR)
%   where the inputs are:
%      TREE:   A MATLAB struct, specifying a syntax tree (usually stored in a
%              TREEVAR)
%      INDSTR: A MATLAB string, that governs the indentation when printing the
%              current node in the syntax tree.
%   and the output is:
%      S:      A MATLAB string, which describes the syntax tree.
%
%   Usually this method is called from within the TREEVAR/PRINT method.
%
%   Example:
%      % First, define a TREEVAR and carry out some operations:
%      u = treeVar();
%      v = cos(u);
%      w = sin(diff(u));
%      t = v + w;
%      % The following are equivalent:
%      print(t)
%      treeVar.printTree(t.tree)
%
% See also TREEVAR.PRINT.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(tree) )
    % Print an empty tree:
    s = sprintf('(empty tree)\n');
    return
end

% Initalize a string to be returned:
s = '';

varString = tree.method;
if ( strcmp(varString, 'constr') )
    % If we're at a constructor leaf, change the text slightly:
    if ( length(tree.ID) == 1)
        % Scalar case.
        varString = 'u';    
    else
        % Systems case, print u and the variable number:
        varString = sprintf('u%i', find(tree.ID == 1));
    end
end

if ( nargin < 2 )
    % We're at the start of the recursion (otherwise, we'd have passed in an
    % INDENTSTR as well).
    
    % At first, we don't do any indentation:
    indentStr = '';
    
    % Print the first method:
    s = [ s, sprintf('%s%s\tdiffOrder: [%s]\n', indentStr, varString, ...
        num2str(tree.diffOrder)) ];
else
    % We're inside the recursion.
    
    if ( indentStr(end) ~= '|')
        % If the current indentation string doesn't end with |, we replace the
        % last whitespace with | when printing the current method. This will be
        % required when we're still setting up the spacing, as we don't want to
        % print a | every time we indent.
        
        % Temporary indentation string for the current method:
        indStrTemp = indentStr;
        indStrTemp(end) = '|';
    else
        indStrTemp = indentStr;
    end
    
    % Print the current method:
    s = [ s, sprintf('%s--%s\tdiffOrder: [%s]\n', indStrTemp, varString, ...
        num2str(tree.diffOrder)) ];
end

% Add whitespace to the indentation string, to be used further in the recusion:
indentStr = [ indentStr, '  ' ];

% Print the remaining syntax tree recursively:
switch tree.numArgs
    case 1
        % Printing univariate methods.
        s = [ s, treeVar.printTree(tree.center, [indentStr, ' ']) ];
    
    case 2
        % Printing bivariate methods.
        
        if ( isstruct(tree.left) )
            % If the left child is a tree, print it recursively:
            s = [ s, treeVar.printTree(tree.left, [indentStr, '|']) ];
        elseif ( isnumeric(tree.left) )
            % Left node was a scalar:
            s = [ s, sprintf('%s|--numerical \tValue: %s\n', ...
                indentStr, num2str(tree.left, '%2.2f')) ];
        else
            % Left node was a CHEBFUN:
            s = [ s, fprintf('%s|--chebfun\n', indentStr) ];
        end
        
        if ( isstruct(tree.right) )
            % If the right child is a tree, print it recursively:
            s = [ s, treeVar.printTree(tree.right, [indentStr, ' ']) ];
        elseif ( isnumeric(tree.right) )
            % Right node is a scalar:
            s = [ s, sprintf('%s|--numerical \tValue: %s\n', ...
                indentStr, num2str(tree.right, '%2.2f')) ];
        else
            % Right node is a CHEBFUN:
            s = [ s, fprintf('%s|--chebfun\n', indentStr) ];
        end       
end

end