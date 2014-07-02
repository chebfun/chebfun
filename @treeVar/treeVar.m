classdef  (InferiorClasses = {?chebfun}) treeVar
    
    properties
        tree
    end
    
    methods
        
        function obj = treeVar(varargin)
            %             obj.tree.method = 'constr';
            %             obj.tree.numArgs = 0;
            %             obj.tree.diffOrder = 0;
            % %             tempTree.method = 'constr';
            % %             tempTree.numArgs = 0;
            % %             tempTree.diffOrder = 0;
            % %             obj.tree = tempTree;
            obj.tree  = struct('method', 'constr', 'numArgs', 0, ...
                'diffOrder', 0, 'height', 0, 'multCoeff', 1);
        end
        
        function f = cos(f)
            %             f.tree.center = f.tree;
            %             f.tree.method = 'cos';
            %             f.tree.numArgs = 1;
            f.tree = struct('method', 'cos', 'numArgs', 1, 'center', f.tree, ...
                'diffOrder', f.tree.diffOrder, 'height', f.tree.height + 1, ...
                'multCoeff', 1);
        end
        
        function f = diff(f, k)
            % By default, compute first derivative:
            if ( nargin < 2 )
                k = 1;
            end
            f.tree = struct('method', 'diff', 'numArgs', 2, ...
                'left', f.tree, 'right', k, 'diffOrder', f.tree.diffOrder + k, ...
                'height', f.tree.height + 1, 'multCoeff', f.tree.multCoeff);
        end
        
        function f = exp(f)
            f.tree = struct('method', 'exp', 'numArgs', 1, 'center', f.tree, ...
                'diffOrder', f.tree.diffOrder, 'height', f.tree.height + 1, ...
                'multCoeff', 1);
        end
        
        function h = minus(f, g)
            h = treeVar();
            if ( ~isa(f, 'treeVar') )
                h.tree = struct('method', 'minus', 'numArgs', 2, ...
                    'left', f, 'right', g.tree, 'diffOrder', g.tree.diffOrder,  ...
                    'height', g.tree.height + 1, 'multCoeff', 1);
            elseif ( ~isa(g, 'treeVar') )
                h.tree = struct('method', 'minus', 'numArgs', 2, ...
                    'left', f.tree, 'right', g, 'diffOrder', f.tree.diffOrder,  ...
                    'height', f.tree.height + 1, 'multCoeff', 1);
            else
                h.tree = struct('method', 'minus', 'numArgs', 2, ...
                    'left', f.tree, 'right', g.tree, ...
                    'diffOrder', max(f.tree.diffOrder, g.tree.diffOrder), ...
                    'height', max(f.tree.height, g.tree.height) + 1, 'multCoeff', 1);
            end
        end
        
        function h = mtimes(f, g)
            if ( isnumeric(f) || isnumeric(g) )
                h = times(f, g);
            else
                error('Dimension mismatch');
            end
            
        end
        
        function h = plus(f, g)
            h = treeVar();
            if ( ~isa(f, 'treeVar') )
                h.tree = struct('method', 'plus', 'numArgs', 2, ...
                    'left', f, 'right', g.tree, 'diffOrder', g.tree.diffOrder, ...
                    'height', g.tree.height + 1, 'multCoeff', 1);
            elseif ( ~isa(g, 'treeVar') )
                h.tree = struct('method', 'plus', 'numArgs', 2, ...
                    'left', f.tree, 'right', g, 'diffOrder', f.tree.diffOrder, ...
                    'height', f.tree.height + 1, 'multCoeff', 1);
            else
                h.tree = struct('method', 'plus', 'numArgs', 2, ...
                    'left', f.tree, 'right', g.tree, ...
                    'diffOrder', max(f.tree.diffOrder, g.tree.diffOrder), ...
                    'height', max(f.tree.height, g.tree.height) + 1, 'multCoeff', 1);
            end
        end
        
        function h = power(f, g)
            h = treeVar();
            if ( ~isa(f, 'treeVar') )
                h.tree = struct('method', 'power', 'numArgs', 2, ...
                    'left', f, 'right', g.tree, 'diffOrder', g.tree.diffOrder, ...
                    'height', g.tree.height + 1);
            elseif ( ~isa(g, 'treeVar') )
                h.tree = struct('method', 'power', 'numArgs', 2, ...
                    'left', f.tree, 'right', g, 'diffOrder', f.tree.diffOrder, ...
                    'height', f.tree.height + 1);
            else
                h.tree = struct('method', 'power', 'numArgs', 2, ...
                    'left', f.tree, 'right', g.tree, ...
                    'diffOrder', max(f.tree.diffOrder, g.tree.diffOrder), ...
                    'height', max(f.tree.height, g.tree.height) + 1);
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
                    'left', f, 'right', g.tree, 'diffOrder', g.tree.diffOrder, ...
                    'height', g.tree.height + 1, 'multCoeff', g.tree.multCoeff.*f);
            elseif ( ~isa(g, 'treeVar') )
                h.tree = struct('method', 'times', 'numArgs', 2, ...
                    'left', f.tree, 'right', g, 'diffOrder', f.tree.diffOrder, ...
                    'height', f.tree.height + 1 , 'multCoeff', f.tree.multCoeff.*g);
            else
                h.tree = struct('method', 'times', 'numArgs', 2, ...
                    'left', f.tree, 'right', g.tree, ...
                    'diffOrder', max(f.tree.diffOrder, g.tree.diffOrder), ...
                    'height', max(f.tree.height, g.tree.height) + 1 , ...
                    'multCoeff', f.tree.multCoeff.*g.tree.multCoeff);
            end
        end
        
    end
    
    methods ( Static = true )
        [newTree, derTree] = splitTree(tree, maxOrder)
        
        [infix, varCounter, varArray] = tree2infix(tree, varCounter, varArray)
        
        anonFun = toAnon(infix, varArray)
        
        anonFun = toRHS(infix, varArray, coeff)
        
        coeff = getCoeffs(infix, varArray)
        
        printTree(tree, ind, indStr)
        
        newTree = expandTree(tree, maxOrder)
        
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
        
        function funOut = toFirstOrder(funIn)
            problemFun = funIn(treeVar());
            
            expTree = treeVar.expandTree(problemFun.tree, ...
                problemFun.tree.diffOrder);
            
            [newTree, derTree] = treeVar.splitTree(expTree, ...
                problemFun.tree.diffOrder);
            
            [infixDer, dummy, varArrayDer] = treeVar.tree2infix(derTree);
            coeffFun = treeVar.toAnon(infixDer, varArrayDer);
            coeffArg = [zeros(1, expTree.diffOrder), 1];
            coeff = coeffFun(coeffArg)
            
            newTree = struct('method', 'uminus', 'numArgs', 1, 'center', newTree);
            [infix, varCounter, varArray] = treeVar.tree2infix(newTree);
            funOut = treeVar.toRHS(infix, varArray, coeff);
        end
        
        
        
    end

end