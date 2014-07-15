function printTree(tree, ind, indStr)
if ( isempty(tree) )
    fprintf('(empty tree)\n');
    return
end


if ( nargin < 2 )
    ind = 0;
    indStr = '';
end

if ( ind == 0)
    fprintf('%s%s\tdiffOrder: %i\n', indStr, tree.method, ...
    tree.diffOrder)
else
    if ( indStr(end) ~= '|')
        indStrTemp = indStr;
        indStrTemp(end) = '|';
    else
        indStrTemp = indStr;
    end
    fprintf('%s--%s\tdiffOrder: %i\n', indStrTemp, tree.method, ...
    tree.diffOrder)
end
indStr = [indStr, '  '];
switch tree.numArgs
    case 1
        treeVar.printTree(tree.center, ind + 1, [indStr, '|']);
    case 2

        if ( isstruct(tree.left) )
            treeVar.printTree(tree.left, ind+1, [indStr, '|']);
        elseif ( isnumeric(tree.left) )
            fprintf('%s|--numerical \tValue: %4.2f\n', indStr, tree.left)
        else
            fprintf('%s|--chebfun\n', indStr)            
        end
        
        if ( isstruct(tree.right) )
            treeVar.printTree(tree.right, ind+1, [indStr, ' ']);
        elseif ( isnumeric(tree.right) )
            fprintf('%s|--numerical \tValue: %4.2f\n', indStr, tree.right)
        else
            fprintf('%s|--chebfun\n', indStr)            
        end
        
end
end