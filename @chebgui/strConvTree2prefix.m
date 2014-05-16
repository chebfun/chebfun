function prefixOut = strConvTree2prefix(syntaxTree)
% TREE2PREFIX function takes in output from the parser and returns a array of
% strings containing the expression on prefix form. When we get into this
% function we have put the output through the parser so we don't have to
% worry about syntax errors.

% TODO:  Documentation.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/chebfun/ for Chebfun information.

% Deal with multi-tree case
prefixOut = {};
if ( numel(syntaxTree) > 1 )
    for k = 1:numel(syntaxTree)
        prefixOut{k} = chebgui.strConvTree2prefix(syntaxTree{k});
    end
end

% Begin by finding the center node of the (partial) syntax tree
nextCenter = syntaxTree.center;

% Find the type of the center node and the token itself
nextType = char(nextCenter(2));
nextSymbol = char(nextCenter(1));

% Make a switch on nextType and determine the action
switch ( nextType )
    % The operators +-*/^= take in two arguments, one to the left and one to
    % the right. Furthermore, COMMA will only appear here as a divider of
    % equations (since we get rid of them as arguments in feval, i.e. u(1,left)
    % in the parser, so we know COMMA always has two arguments as well
    case {'OP+', 'OP-', 'OP*', 'OP/', 'OP^', 'OP=', 'OP>', 'OP>=', 'OP<', 'OP<=', 'COMMA'}
        prefixOut = [{nextSymbol, nextType};
            chebgui.strConvTree2prefix(syntaxTree.left);
            chebgui.strConvTree2prefix(syntaxTree.right)];
    % Unary operators only have one argument which is stored to the right
    case 'UN+'
        prefixOut = [{'+', 'UN+'};
            chebgui.strConvTree2prefix(syntaxTree.right)];
    case 'UN-'
        prefixOut = [{'-', 'UN-'};
            chebgui.strConvTree2prefix(syntaxTree.right)];
    % If we get a number or a variable we simply return that. Those types
    % are leaves so the don't have any children.
    case {'NUM',  'VAR',  'INDVAR', 'PDEVAR', 'LAMBDA', 'STR'}
        prefixOut = {nextSymbol, nextType};
    % A function only has one argument which is to the right
    case 'FUNC1'
        prefixOut = [{nextSymbol, 'FUNC1'};
            chebgui.strConvTree2prefix(syntaxTree.right)];
    case 'FUNC2'
        prefixOut = [{nextSymbol, 'FUNC2'};
            chebgui.strConvTree2prefix(syntaxTree.left);
            chebgui.strConvTree2prefix(syntaxTree.right)];
    case 'FUNC3'
        prefixOut = [{nextSymbol, 'FUNC3'};
            chebgui.strConvTree2prefix(syntaxTree.left);
            chebgui.strConvTree2prefix(syntaxTree.right);
            chebgui.strConvTree2prefix(syntaxTree.arg)];        
    otherwise % Only possible token left is a derivative
        prefixOut = [{nextSymbol, nextType};
            chebgui.strConvTree2prefix(syntaxTree.right)];
end
