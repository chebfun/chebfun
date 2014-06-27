function [out, varCounter, varArray] = tree2infix(tree, varCounter, varArray)
if ( nargin < 2 )
    varCounter = 1;
    varArray = [];
end

if isempty(tree)
    out = '';
    return
end

switch tree.numArgs
    case 0
        out = 'u(1)';
    case 1
        [tempOut, varCounter, varArray] = ...
            treeVar.tree2infix(tree.center, varCounter, varArray);
        out = sprintf('%s(%s)', tree.method, tempOut);
    case 2
        if ( strcmp(tree.method, 'diff') )
            out = sprintf('u(%i)', tree.right + 1);
            return
        end
        
        if ( isstruct(tree.left) )
            [leftInfix, varCounter, varArray] = ...
                treeVar.tree2infix(tree.left, varCounter, varArray);
        else
            varName = sprintf('var%i', varCounter);
            leftInfix = varName;
            varArray = [varArray; {varName, tree.left}];
            varCounter = varCounter + 1;
        end
        
        if ( isstruct(tree.right) )
            [rightInfix, varCounter, varArray] = ...
                treeVar.tree2infix(tree.right, varCounter, varArray);
        else
            varName = sprintf('var%i', varCounter);
            rightInfix = varName;
            varArray = [varArray; {varName, tree.right}];
            varCounter = varCounter + 1;
        end
        
        out = sprintf('%s(%s, %s)', tree.method, ...
            leftInfix, rightInfix);
end

end