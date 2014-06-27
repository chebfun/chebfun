classdef  (InferiorClasses = {?chebfun}) treeChebfun < adchebfun
    
    properties
        tree
    end
    
    methods
        
        function obj = treeChebfun(varargin)
            obj = obj@adchebfun(varargin{:});
            obj.tree = struct('method', 'constr', 'numArgs', 0, 'diffOrder', 0);
        end
        
        function f = cos(f)
            f = cos@adchebfun(f);
            f.tree = struct('method', 'cos', 'numArgs', 1, 'center', f.tree, ...
                'diffOrder', f.tree.diffOrder);
        end
        
        function f = diff(f, k)
            % By default, compute first derivative:
            if ( nargin < 2 )
                k = 1; 
            end
            
            f = diff@adchebfun(f, k);
            f.tree = struct('method', 'diff', 'numArgs', 2, ...
                'left', f.tree, 'right', k, 'diffOrder', f.tree.diffOrder + k);
        end
        
        function f = exp(f)
            f = exp@adchebfun(f);
            f.tree = struct('method', 'exp', 'numArgs', 1, 'center', f.tree, ...
                'diffOrder', f.tree.diffOrder);
        end
        
        function h = minus(f, g)
            % Call superclass method:
            h = minus@adchebfun(f, g);
            
            if ( ~isa(f, 'treeChebfun') )
                h.tree = struct('method', 'minus', 'numArgs', 2, ...
                    'left', f, 'right', g.tree, 'diffOrder', g.tree.diffOrder);
            elseif ( ~isa(g, 'treeChebfun') )
                h.tree = struct('method', 'minus', 'numArgs', 2, ...
                    'left', f.tree, 'right', g, 'diffOrder', f.tree.diffOrder);
            else
                h.tree = struct('method', 'minus', 'numArgs', 2, ...
                    'left', f.tree, 'right', g.tree, ...
                    'diffOrder', max(f.tree.diffOrder, g.tree.diffOrder));
            end
        end
        
        function h = plus(f, g)
            % Call superclass method:
            h = plus@adchebfun(f, g);
            
            if ( ~isa(f, 'treeChebfun') )
                h.tree = struct('method', 'plus', 'numArgs', 2, ...
                    'left', f, 'right', g.tree, 'diffOrder', g.tree.diffOrder);           
            elseif ( ~isa(g, 'treeChebfun') )
                h.tree = struct('method', 'plus', 'numArgs', 2, ...
                    'left', f.tree, 'right', g, 'diffOrder', f.tree.diffOrder);
            else
                h.tree = struct('method', 'plus', 'numArgs', 2, ...
                    'left', f.tree, 'right', g.tree, ...
                    'diffOrder', max(f.tree.diffOrder, g.tree.diffOrder));
            end
        end
        
        function h = power(f, g)
            % Call superclass method:
            h = power@adchebfun(f, g);
            
            if ( ~isa(f, 'treeChebfun') )
                h.tree = struct('method', 'power', 'numArgs', 2, ...
                    'left', f, 'right', g.tree, 'diffOrder', g.tree.diffOrder);
            elseif ( ~isa(g, 'treeChebfun') )
                h.tree = struct('method', 'power', 'numArgs', 2, ...
                    'left', f.tree, 'right', g, 'diffOrder', f.tree.diffOrder);
            else
                h.tree = struct('method', 'power', 'numArgs', 2, ...
                    'left', f.tree, 'right', g.tree, ...
                    'diffOrder', max(f.tree.diffOrder, g.tree.diffOrder));
            end
        end
        
        function f = sin(f)
            f = sin@adchebfun(f);
            f.tree = struct('method', 'sin', 'numArgs', 1, 'center', f.tree, ...
                'diffOrder', f.tree.diffOrder);
        end
        
        function h = times(f, g)
            % Call superclass method:
            h = times@adchebfun(f, g);
            
            if ( ~isa(f, 'treeChebfun') )
                h.tree = struct('method', 'times', 'numArgs', 2, ...
                    'left', f, 'right', g.tree, 'diffOrder', g.tree.diffOrder);           
            elseif ( ~isa(g, 'treeChebfun') )
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
                    treeChebfun.printTree(tree.center, ind + 1);
                case 2
                    if ( isstruct(tree.left) )
                        treeChebfun.printTree(tree.left, ind + 1);
                    else
                        fprintf('%s  numerical\n', spaceStr)
                    end
                    
                    if ( isstruct(tree.right) )
                        treeChebfun.printTree(tree.right, ind + 1);
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
                        treeChebfun.tree2prefix(tree.center)];
                    
                case 2
                    out = [{tree.method}; ...
                        treeChebfun.tree2prefix(tree.left); ...
                        treeChebfun.tree2prefix(tree.right)];
            end
        end                
        
    end
    
    
end