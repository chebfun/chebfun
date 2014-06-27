classdef  (InferiorClasses = {?chebfun}) treeVar
    
    properties
        tree
    end
    
    methods
        
        function obj = treeVar(varargin)
            obj.tree = struct('method', 'constr', 'numArgs', 0, 'diffOrder', 0);
        end
        
        function f = cos(f)
            f.tree = struct('method', 'cos', 'numArgs', 1, 'center', f.tree, ...
                'diffOrder', f.tree.diffOrder);
        end
        
        function f = diff(f, k)
            % By default, compute first derivative:
            if ( nargin < 2 )
                k = 1; 
            end
                        f.tree = struct('method', 'diff', 'numArgs', 2, ...
                'left', f.tree, 'right', k, 'diffOrder', f.tree.diffOrder + k);
        end
        
        function f = exp(f)
            f.tree = struct('method', 'exp', 'numArgs', 1, 'center', f.tree, ...
                'diffOrder', f.tree.diffOrder);
        end
        
        function h = minus(f, g)
            h = treeVar();
            if ( ~isa(f, 'treeVar') )
                h.tree = struct('method', 'minus', 'numArgs', 2, ...
                    'left', f, 'right', g.tree, 'diffOrder', g.tree.diffOrder);
            elseif ( ~isa(g, 'treeVar') )
                h.tree = struct('method', 'minus', 'numArgs', 2, ...
                    'left', f.tree, 'right', g, 'diffOrder', f.tree.diffOrder);
            else
                h.tree = struct('method', 'minus', 'numArgs', 2, ...
                    'left', f.tree, 'right', g.tree, ...
                    'diffOrder', max(f.tree.diffOrder, g.tree.diffOrder));
            end
        end
        
        function h = plus(f, g)
            h = treeVar();
            if ( ~isa(f, 'treeVar') )
                h.tree = struct('method', 'plus', 'numArgs', 2, ...
                    'left', f, 'right', g.tree, 'diffOrder', g.tree.diffOrder);           
            elseif ( ~isa(g, 'treeVar') )
                h.tree = struct('method', 'plus', 'numArgs', 2, ...
                    'left', f.tree, 'right', g, 'diffOrder', f.tree.diffOrder);
            else
                h.tree = struct('method', 'plus', 'numArgs', 2, ...
                    'left', f.tree, 'right', g.tree, ...
                    'diffOrder', max(f.tree.diffOrder, g.tree.diffOrder));
            end
        end
        
        function h = power(f, g)
            h = treeVar();
            if ( ~isa(f, 'treeVar') )
                h.tree = struct('method', 'power', 'numArgs', 2, ...
                    'left', f, 'right', g.tree, 'diffOrder', g.tree.diffOrder);
            elseif ( ~isa(g, 'treeVar') )
                h.tree = struct('method', 'power', 'numArgs', 2, ...
                    'left', f.tree, 'right', g, 'diffOrder', f.tree.diffOrder);
            else
                h.tree = struct('method', 'power', 'numArgs', 2, ...
                    'left', f.tree, 'right', g.tree, ...
                    'diffOrder', max(f.tree.diffOrder, g.tree.diffOrder));
            end
        end
        
        function f = sin(f)
            f.tree = struct('method', 'sin', 'numArgs', 1, 'center', f.tree, ...
                'diffOrder', f.tree.diffOrder);
        end
        
        function h = times(f, g)
            h = treeVar();
            if ( ~isa(f, 'treeVar') )
                h.tree = struct('method', 'times', 'numArgs', 2, ...
                    'left', f, 'right', g.tree, 'diffOrder', g.tree.diffOrder);           
            elseif ( ~isa(g, 'treeVar') )
                h.tree = struct('method', 'times', 'numArgs', 2, ...
                    'left', f.tree, 'right', g, 'diffOrder', f.tree.diffOrder);
            else
                h.tree = struct('method', 'times', 'numArgs', 2, ...
                    'left', f.tree, 'right', g.tree, ...
                    'diffOrder', max(f.tree.diffOrder, g.tree.diffOrder));
            end
        end
        
    end
    
    methods ( Static = true )
        [newTree, derTree] = splitTree(tree, maxOrder)
        
        [infix, varCounter, varArray] = tree2infix(tree, varCounter, varArray)
        
        anonFun = toAnon(infix, varArray)
        
        function printTree(tree, ind)
            
            if ( nargin < 2 )
                ind = 0;
            end
            
            spaceStr = repmat(' ', 1, 2*ind);
            fprintf('%s%s.\tdiffOrder: %i\n', spaceStr, tree.method, ...
                tree.diffOrder)
            
            switch tree.numArgs
                case 1
                    treeVar.printTree(tree.center, ind + 1);
                case 2
                    if ( isstruct(tree.left) )
                        treeVar.printTree(tree.left, ind + 1);
                    else
                        fprintf('%s  numerical\n', spaceStr)
                    end
                    
                    if ( isstruct(tree.right) )
                        treeVar.printTree(tree.right, ind + 1);
                    else
                        fprintf('%snumerical\n', spaceStr)
                    end
            end
        end
        
        function out = tree2prefix(tree)
            
            switch tree.numArgs
                case 0
                    out = 'u';
                    
                case 1
                    out = [{tree.method}; ...
                        treeVar.tree2prefix(tree.center)];
                    
                case 2
                    out = [{tree.method}; ...
                        treeVar.tree2prefix(tree.left); ...
                        treeVar.tree2prefix(tree.right)];
            end
        end                
        
    end
    
    
end