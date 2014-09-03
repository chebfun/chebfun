classdef  (InferiorClasses = {?chebfun}) treeVar
    
    properties
        tree
        domain
    end
    
    methods
        
        function obj = treeVar(varargin)
            if ( nargin > 0 )
                IDvec = varargin{1};
                obj.domain = varargin{2};
            else
                IDvec = 1;
                obj.domain = [-1 1];
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
                disp('and the domain:')
                disp(u.domain);
            else
                disp('Array-valued treeVar, with trees:');
                for ( treeCounter = 1:length(u) )
                    fprintf('tree %i\n', treeCounter)
                    disp(u(treeCounter).tree);
                    disp('and the domain:')
                    disp(u.domain);
                end
            end
        end
        
        function f = exp(f)
            f.tree = f.univariate(f.tree, 'exp');
        end
        
        function h = minus(f, g)
            h = treeVar();
            if ( ~isa(f, 'treeVar') )
                h.tree = treeVar.bivariate(f, g.tree, 'minus', 1);
            elseif ( ~isa(g, 'treeVar') )
                h.tree = treeVar.bivariate(f.tree, g, 'minus', 0);
            else
                h.tree = treeVar.bivariate(f.tree, g.tree, 'minus', 2);
            end
            h.domain = updateDomain(f, g);
        end
        
        function h = mtimes(f, g)
            if ( isnumeric(f) || isnumeric(g) )
                h = times(f, g);
                h.domain = updateDomain(f, g);
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
            else
                h.tree = treeVar.bivariate(f.tree, g.tree, 'plus', 2);
            end
            h.domain = updateDomain(f, g);
        end
        
        function plot(treeVar)
            treeVar.plotTree(treeVar.tree);
        end
        
        function h = power(f, g)
            h = treeVar();
            if ( ~isa(f, 'treeVar') )
                h.tree = struct('method', 'power', 'numArgs', 2, ...
                    'left', f, 'right', g.tree, 'diffOrder', g.tree.diffOrder, ...
                    'height', g.tree.height + 1, ...
                    'ID', g.tree.ID);
            elseif ( ~isa(g, 'treeVar') )
                h.tree = struct('method', 'power', 'numArgs', 2, ...
                    'left', f.tree, 'right', g, 'diffOrder', f.tree.diffOrder, ...
                    'height', f.tree.height + 1, ...
                    'ID', f.tree.ID);
            else
                h.tree = struct('method', 'power', 'numArgs', 2, ...
                    'left', f.tree, 'right', g.tree, ...
                    'diffOrder', max(f.tree.diffOrder, g.tree.diffOrder), ...
                    'height', max(f.tree.height, g.tree.height) + 1, ...
                    'ID', f.tree.ID | g.tree.ID);
            end
            h.domain = updateDomain(f, g);
        end
        
        function h = rdivide(f, g)
            h = treeVar();
            if ( ~isa(f, 'treeVar') )
                h.tree = treeVar.bivariate(f, g.tree, 'rdivide', 1);
            elseif ( ~isa(g, 'treeVar') )
                h.tree = treeVar.bivariate(f.tree, g, 'rdivide', 0);
            else
                h.tree = treeVar.bivariate(f.tree, g.tree, 'rdivide', 2);
            end
            h.domain = updateDomain(f, g);
        end
        
        function f = sin(f)
            f.tree = f.univariate(f.tree, 'sin');
        end
        
        function h = times(f, g)
            h = treeVar();
            if ( ~isa(f, 'treeVar') )
                h.tree = treeVar.bivariate(f, g.tree, 'times', 1);
            elseif ( ~isa(g, 'treeVar') )
                h.tree = treeVar.bivariate(f.tree, g, 'times', 0);
            else
                h.tree = treeVar.bivariate(f.tree, g.tree, 'times', 2);
            end
            h.domain = updateDomain(f, g);
        end
        
        function dom = updateDomain(f, g)
            if ( isnumeric(f) )
                dom = g.domain;
            elseif ( isnumeric(g) )
                dom = f.domain;
            else
                dom = union(f.domain, g.domain);
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
        
        funOut = toRHS(infix, varArray, coeff, indexStart, totalDiffOrders);
        
        [newTree, derTree] = splitTree(tree, maxOrder)
        
        [infix, varCounter, varArray] = tree2infix(tree, diffOrders, varCounter, varArray)
        
        anonFun = toAnon(infix, varArray)
        
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
        
        [funOut, varIndex, problemDom] = toFirstOrder(funIn, domain)
        
        [funOut, indexStart, problemDom] = toFirstOrderSystem(funIn, domain)
        
        function idx = sortConditions(funIn, domain)
            %SORTBCS    Return a vector with indices on how to sort the results
            %            of evaluating N.LBC/RBC
            
            numArgs = nargin(funIn);
            args = cell(numArgs, 1);
            argsVec = zeros(1, numArgs);
            
            % Populate the args cell
            for argCount = 1:numArgs
                argsVec(argCount) = 1;
                args{argCount} = treeVar(argsVec, domain);
                % Reset the index vector
                argsVec = 0*argsVec;
            end
            
            % Evaluate FUNIN with the TREEVAR arguments:
            bcResults = funIn(args{:});
            
            % Look at the results of evaluating the boundary conditions, find
            % what constraint operated on what variable, and what it's diffOrder
            % was:
            varList = cell(numArgs, 1);
            diffOrderList = varList;
            for tCounter = 1:length(bcResults)
                % Current tree we're looking at:
                tempTree = bcResults(tCounter).tree;
                
                % Check whether more than one variable appear in the condition
                if ( sum(tempTree.ID) > 1 )
                    error('CHEBFUN:TREEVAR:sortConditions:nonSeparated', ...
                        ['For initial value problems, only separated ', ...
                        'conditions are supported.']);
                end
                
                % What's the active variable in the current tree (i.e. what
                % variable did the constraint apply to)?
                activeVar = find(tempTree.ID == 1);
                % What's the diffOrder of the current constraint?
                activeDiffOrder = tempTree.diffOrder;
                activeDiffOrder = activeDiffOrder(activeVar);
                
                % Store in a list what variable the current constraint applies
                % to, and what the current diffOrder is:
                varList{activeVar} = [varList{activeVar}, tCounter];
                diffOrderList{activeVar} = [diffOrderList{activeVar}, ...
                    activeDiffOrder];
            end
            
            % Initalise an index vector to be returned
            idx = [];
            
            % Go through the list of what variables appeared in what
            % constraints, and sort them based on diffOrders:
            for varCounter = 1:numArgs
                [dummy, diffOrderIndex] = sort(diffOrderList{varCounter});
                tempIndex = varList{varCounter}(diffOrderIndex);
                idx = [idx, tempIndex];
            end
            
        end
        
    end
    
end