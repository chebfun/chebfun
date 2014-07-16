function [out, varCounter, varArray] = tree2infix(tree, eqno, indexStart)
varCounter = 1;
varArray = [];

if isempty(tree)
    out = '';
    return
end

[out, varCounter, varArray] = ...
    toInfix(tree, eqno, indexStart, varCounter, varArray);

end

function [out, varCounter, varArray] = toInfix(tree, eqno, indexStart, varCounter, varArray)

switch tree.numArgs
    case 0
        out = sprintf('u(%i)', indexStart(tree.ID));
    case 1
        [tempOut, varCounter, varArray] = ...
            toInfix(tree.center, eqno, indexStart, varCounter, varArray);
        out = sprintf('%s(%s)', tree.method, tempOut);
    case 2
        if ( strcmp(tree.method, 'diff') )
            out = sprintf('u(%i)', tree.right + indexStart(tree.ID));
            return
        end
        
        if ( isstruct(tree.left) )
            [leftInfix, varCounter, varArray] = ...
                toInfix(tree.left, eqno, indexStart, varCounter, varArray);
        else
            varName = sprintf('eq%i_var%i', eqno, varCounter);
            if ( isnumeric(tree.left) )
                % The left tree is a scalar.
                leftInfix = varName;
            else
                % The left tree must be a chebfun -- need to be able to evaluate
                % it on gridpoints
                leftInfix = [varName, '(t)'];
            end
            varArray = [varArray; {varName, tree.left}];
            varCounter = varCounter + 1;
        end
        
        if ( isstruct(tree.right) )
            [rightInfix, varCounter, varArray] = ...
                toInfix(tree.right, eqno, indexStart, varCounter, varArray);
        else
            varName = sprintf('eq%i_var%i', eqno, varCounter);
            if ( isnumeric(tree.right) )
                % The left tree is a scalar.
                rightInfix = varName;
            else
                % The left tree must be a chebfun -- need to be able to evaluate
                % it on gridpoints
                rightInfix = [varName, '(t)'];    
            end
            varArray = [varArray; {varName, tree.right}];
            varCounter = varCounter + 1;
        end
        
        out = sprintf('%s(%s, %s)', tree.method, ...
            leftInfix, rightInfix);
end

end