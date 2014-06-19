function prefixOut = tree2prefix(syntaxTree)
%TREE2PREFIX   Convert a syntax tree to prefixForm
%   PREFIXOUT = TREE2PREFIX(SYNTAXTREE) returns the cell-array of strings
%   PREFIXOUT that corresponds to the prefix form of SYNTAXTREE. Generally,
%   SYNTAXTREE will be the output of the STRINGPARSER/PARSER method.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:
%   When this method is called, the input has already been through the parser,
%   so we don't have to worry about syntax errors.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Deal with multi-tree case
prefixOut = {};
if ( numel(syntaxTree) > 1 )
    for k = 1:numel(syntaxTree)
        prefixOut{k} = stringParser.tree2prefix(syntaxTree{k}); %#ok<AGROW>
    end
end

% Begin by finding the center node of the (partial) syntax tree
nextCenter = syntaxTree.center;

% Find the type of the center node and the token itself
nextType = char(nextCenter(2));
nextSymbol = char(nextCenter(1));

% Go recursively through the SYNTAXTREE tree, and return the prefix form. We
% make a switch on nextType and determine the action
switch ( nextType )
    
    % The operators +-*/^= take in two arguments, one to the left and one to
    % the right. Furthermore, COMMA will only appear here as a divider of
    % equations (since we get rid of them as arguments in feval, i.e. u(1, left)
    % in the parser, so we know COMMA always has two arguments as well
    case {'OP+', 'OP-', 'OP*', 'OP/', 'OP^', 'OP=', 'OP>', 'OP>=', 'OP<', 'OP<=', 'COMMA'}
        prefixOut = [{nextSymbol, nextType};
            stringParser.tree2prefix(syntaxTree.left);
            stringParser.tree2prefix(syntaxTree.right)];
    % Unary operators only have one argument which is stored to the right
    case 'UN+'
        prefixOut = [{'+', 'UN+'};
            stringParser.tree2prefix(syntaxTree.right)];
    case 'UN-'
        prefixOut = [{'-', 'UN-'};
            stringParser.tree2prefix(syntaxTree.right)];
    % If we get a number or a variable we simply return that. Those types
    % are leaves so the don't have any children.
    case {'NUM',  'VAR',  'INDVAR', 'PDEVAR', 'LAMBDA', 'STR'}
        prefixOut = {nextSymbol, nextType};
    % A function only has one argument which is to the right
    case 'FUNC1'
        prefixOut = [{nextSymbol, 'FUNC1'};
            stringParser.tree2prefix(syntaxTree.right)];
    case 'FUNC2'
        prefixOut = [{nextSymbol, 'FUNC2'};
            stringParser.tree2prefix(syntaxTree.left);
            stringParser.tree2prefix(syntaxTree.right)];
    case 'FUNC3'
        prefixOut = [{nextSymbol, 'FUNC3'};
            stringParser.tree2prefix(syntaxTree.left);
            stringParser.tree2prefix(syntaxTree.right);
            stringParser.tree2prefix(syntaxTree.arg)];        
    otherwise % Only possible token left is a derivative
        prefixOut = [{nextSymbol, nextType};
            stringParser.tree2prefix(syntaxTree.right)];
        
end

end
