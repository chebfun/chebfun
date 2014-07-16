classdef  (InferiorClasses = {?chebfun}) treeVar
    
    properties
        tree
    end
    
    methods
        
        function obj = treeVar(varargin)
            if ( nargin > 0 )
                IDvec = varargin{1};
            else
                IDvec = 1;
            end
            obj.tree  = struct('method', 'constr', 'numArgs', 0, ...
                'diffOrder', 0*IDvec, 'height', 0, 'multCoeff', 1, ...
                'ID', logical(IDvec));
        end
        
        function f = cos(f)
            f.tree = f.univariate(f.tree, 'cos');
        end
        
        function f = diff(f, k)
            % By default, compute first derivative:
            if ( nargin < 2 )
                k = 1;
            end
            f.tree = struct('method', 'diff', 'numArgs', 2, ...
                'left', f.tree, 'right', k, ...
                'diffOrder', f.tree.diffOrder + k*f.tree.ID, ...
                'height', f.tree.height + 1, ...
                'ID', f.tree.ID, ...
                'multCoeff', f.tree.multCoeff);
        end
        
        function display(u)
            if ( length(u) == 1 )
                disp('treeVar with tree:')
                disp(u.tree);
            else
                disp('Array-valued treeVar, with trees:');
                for ( treeCounter = 1:length(u) )
                    fprintf('tree %i\n', treeCounter)
                    disp(u(treeCounter).tree);
                end
            end
        end
        
        function f = exp(f)
            f.tree = f.univariate(f.tree, 'exp');
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
                h.tree = treeVar.bivariate(f, g.tree, 'plus', 1);
            elseif ( ~isa(g, 'treeVar') )
                h.tree = treeVar.bivariate(f.tree, g, 'plus', 0);
%                 h.tree = struct('method', 'plus', 'numArgs', 2, ...
%                     'left', f.tree, 'right', g, ...
%                     'diffOrder', f.tree.diffOrder, ...
%                     'ID', f.tree.ID, ...
%                     'height', f.tree.height + 1, ...
%                     'multCoeff', 1);
            else
                h.tree = treeVar.bivariate(f.tree, g.tree, 'plus', 2);
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
            f.tree = f.univariate(f.tree, 'sin');
        end
        
        function h = times(f, g)
            h = treeVar();
            if ( ~isa(f, 'treeVar') )
                h.tree = treeVar.bivariate(f, g.tree, 'times', 1);

%                 h.tree = struct('method', 'times', 'numArgs', 2, ...
%                     'left', f, 'right', g.tree, 'diffOrder', g.tree.diffOrder, ...
%                     'height', g.tree.height + 1);
            elseif ( ~isa(g, 'treeVar') )
                h.tree = treeVar.bivariate(f.tree, g, 'times', 0);
% 
%                 h.tree = struct('method', 'times', 'numArgs', 2, ...
%                     'left', f.tree, 'right', g, 'diffOrder', f.tree.diffOrder, ...
%                     'height', f.tree.height + 1);
            else
                h.tree = treeVar.bivariate(f.tree, g.tree, 'times');
            end
        end
        
    end
    
    methods ( Static = true )
        
        
        function treeOut = univariate(treeIn, method)
            % Construct a syntax tree for univariate functions.
            treeOut = struct('method', method, ...
                'numArgs', 1, 'center', treeIn, ...
                'diffOrder', treeIn.diffOrder, ...
                'height', treeIn.height + 1, ...
                'multCoeff', 1, ...
                'ID', treeIn.ID);
        end
        
        function treeOut = bivariate(leftTree, rightTree, method, type)
            % type
            %   2: Both treeVars
            %   1: Only right treeVar
            %   0: Only left treeVar
            if ( type == 2 )
                treeOut = struct('method', method, 'numArgs', 2, ...
                    'left', leftTree, 'right', rightTree, ...
                    'diffOrder', max(leftTree.diffOrder, rightTree.diffOrder), ...
                    'ID', leftTree.ID | rightTree.ID, ...
                    'height', max(leftTree.height, rightTree.height) + 1);
            elseif ( type == 1 )
                treeOut = struct('method', method, 'numArgs', 2, ...
                    'left', leftTree, 'right', rightTree, ...
                    'diffOrder', rightTree.diffOrder, ...
                    'height', rightTree.height + 1, ...
                    'ID', rightTree.ID);
            else
                treeOut = struct('method', method, 'numArgs', 2, ...
                    'left', leftTree, 'right', rightTree, ...
                    'diffOrder', leftTree.diffOrder, ...
                    'height', leftTree.height + 1, ...
                    'ID', leftTree.ID);
                
                
            end
        end
        
        
        [newTree, derTree] = splitTree(tree, maxOrder)
        
        [infix, varCounter, varArray] = tree2infix(tree, diffOrders, varCounter, varArray)
        
        anonFun = toAnon(infix, varArray)
        
        anonFun = toRHS(infix, varArray, coeff)
        
        coeff = getCoeffs(infix, varArray)
        
        printTree(tree, ind, indStr)
        
        newTree = expandTree(tree, maxOrder)
        
        plotTree(tree, varargin)
        
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
            coeff = coeffFun(coeffArg);
            
            newTree = struct('method', 'uminus', 'numArgs', 1, 'center', newTree);
            [infix, varCounter, varArray] = treeVar.tree2infix(newTree);
            funOut = treeVar.toRHS(infix, varArray, coeff);
        end
        
        
        
    end
    
end