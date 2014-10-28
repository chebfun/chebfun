function s = printTree(tree, ind, indStr)
if ( isempty(tree) )
    s = sprintf('(empty tree)\n');
    return
end

s = '';

if ( nargin < 2 )
    ind = 0;
    indStr = '';
end

if ( ind == 0)
    s = [s, sprintf('%s%s\tdiffOrder: %i\n', indStr, tree.method, ...
        tree.diffOrder)];
else
    if ( indStr(end) ~= '|')
        indStrTemp = indStr;
        indStrTemp(end) = '|';
    else
        indStrTemp = indStr;
    end
    s = [s, sprintf('%s--%s\tdiffOrder: %i\n', indStrTemp, tree.method, ...
        tree.diffOrder)];
end
indStr = [indStr, '  '];
switch tree.numArgs
    case 1
        s = [s, treeVar.printTree(tree.center, ind + 1, [indStr, '|'])];
    case 2
        
        if ( isstruct(tree.left) )
            s = [s, treeVar.printTree(tree.left, ind+1, [indStr, '|'])];
        elseif ( isnumeric(tree.left) )
            s = [s, sprintf('%s|--numerical \tValue: %4.2f\n', indStr, tree.left)];
        else
            s = [s, fprintf('%s|--chebfun\n', indStr)];
        end
        
        if ( isstruct(tree.right) )
            s = [s, treeVar.printTree(tree.right, ind+1, [indStr, ' '])];
        elseif ( isnumeric(tree.right) )
            s = [s, sprintf('%s|--numerical \tValue: %4.2f\n', indStr, tree.right)];
        else
            s = [s, fprintf('%s|--chebfun\n', indStr)];
        end       
end
end