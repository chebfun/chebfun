function prefixOut = tree2prefix(guifiles,syntaxTree)
% TREE2PREFIX function takes in output from the parser and returns a array of
% strings containing the expression on prefix form. When we get into this
% function we have put the output through the parser so we don't have to
% worry about syntax errors.

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Deal with multi-tree case
prefixOut = {};
if numel(syntaxTree) > 1,
    for k = 1:numel(syntaxTree)
        prefixOut{k} = pree2prefix(guifiles,syntaxTree{k});
    end
end

% Begin by finding the center node of the (partial) syntax tree
nextCenter = syntaxTree.center;

% Find the type of the center node and the token itself
nextType = char(nextCenter(2));
nextSymbol = char(nextCenter(1));

% Make a switch on nextType and determine the action
switch nextType
    % The operators +-*/^= take in two arguments, one to the left and one to
    % the right. Furthermore, COMMA will only appear here as a divider of
    % equations (since we get rid of them as arguments in feval, i.e. u(1,left) in the
    % parser, so we know COMMA always has two arguments as well
    case {'OP+','OP-','OP*','OP/','OP^','OP=','OP>','OP>=','OP<','OP<=','COMMA'}
        prefixOut = [{nextSymbol, nextType}; tree2prefix(guifiles,syntaxTree.left); tree2prefix(guifiles,syntaxTree.right)];
    % Unary operators only have one argument which is stored to the right
    case 'UN+'
        prefixOut = [{'+', 'UN+'}; tree2prefix(guifiles,syntaxTree.right)];
    case 'UN-'
        prefixOut = [{'-', 'UN-'}; tree2prefix(guifiles,syntaxTree.right)];
    % If we get a number or a variable we simply return that. Those types
    % are leafs so the don't have any childs
    case {'NUM', 'VAR', 'INDVAR','PDEVAR','LAMBDA','STR'}
        prefixOut = {nextSymbol, nextType};
    % A function only has one argument which is to the right
    case 'FUNC1'
        prefixOut = [{nextSymbol,'FUNC1'}; tree2prefix(guifiles,syntaxTree.right)];
    case 'FUNC2'
        prefixOut = [{nextSymbol,'FUNC2'}; tree2prefix(guifiles,syntaxTree.left); tree2prefix(guifiles,syntaxTree.right)];
    case 'FUNC3'
        prefixOut = [{nextSymbol,'FUNC3'}; tree2prefix(guifiles,syntaxTree.left); tree2prefix(guifiles,syntaxTree.right); tree2prefix(guifiles,syntaxTree.arg)];        
    otherwise % Only possible token left is a derivative
        prefixOut = [{nextSymbol,nextType}; tree2prefix(guifiles,syntaxTree.right)];
end
