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
                h.tree = treeVar.bivariate(f, g.tree, 'minus', 1);
            elseif ( ~isa(g, 'treeVar') )
                h.tree = treeVar.bivariate(f.tree, g, 'minus', 0);
            else
                h.tree = treeVar.bivariate(f.tree, g.tree, 'minus', 2);
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
            else
                h.tree = treeVar.bivariate(f.tree, g.tree, 'plus', 2);
            end
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
                h.tree = treeVar.bivariate(f.tree, g.tree, 'times', 2);
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
        
        function funOut = toFirstOrder(funIn, domain)
            % Independent variable on the domain
            t = chebfun(@(t) t, domain);
            arg = treeVar(1);
            
            % If funIn only has one input argument, we just give it a treeVar()
            % argument. Otherwise, the first input will be the independent
            % variable on the domain:
            if ( nargin(funIn) == 1 )
                problemFun = funIn(arg);
            else
                problemFun = funIn(t, arg);
            end
            
            maxDifforder = problemFun.tree.diffOrder;
            
            expTree = treeVar.expandTree(problemFun.tree, maxDifforder);
            
            [newTree, derTree] = treeVar.splitTree(expTree, ...
                problemFun.tree.diffOrder);
            
            [infixDer, dummy, varArrayDer] = treeVar.tree2infix(derTree, 1, 1);
            coeffFun = treeVar.toAnon(infixDer, varArrayDer);
            coeffArg = [zeros(1, expTree.diffOrder), 1];
            
            
            coeff = {coeffFun(t, coeffArg)};
            
            newTree = struct('method', 'uminus', 'numArgs', 1, 'center', newTree);
            [infix, varCounter, varArray] = treeVar.tree2infix(newTree, 1, 1);
            infix = {infix};
            varArray = {varArray};
            funOut = treeVar.toRHS(infix, varArray, coeff, 1, maxDifforder);
        end
        
        function [funOut, indexStart] = toFirstOrderSystem(funIn, domain)
            % Independent variable on the domain
            t = chebfun(@(t) t, domain);
            
            % The first argument to funIn must be the independent time variable.
            % So the number of treeVar arguments needed is one less:
            numArgs = nargin(funIn) - 1;
            args = cell(numArgs, 1);
            argsVec = zeros(1, numArgs);
            
            % Populate the args cell
            for argCount = 1:numArgs
                argsVec(argCount) = 1;
                args{argCount} = treeVar(argsVec);
                % Reset the index vector
                argsVec = 0*argsVec;
            end
            
            % Evaluate FUNIN with the TREEVAR arguments:
            fevalResult = funIn(t, args{:});
            
            % Initialize cells to store the infix forms of the expressions, the
            % coefficients multiplying the highest order derivatives and any
            % variables that appear in the anonymous functions:
            systemInfix = cell(length(fevalResult),1);
            coeffs = systemInfix;
            varArrays = systemInfix;
            
            % First look at all diffOrders to ensure we start with the correct
            % indices. INDEXSTART denotes at which index we should start
            % indexing each variable from. E.g., in the coupled system
            %   [v'' + w; v + w'']
            % we will have v = u(1), v' = u(2), w = u(3), w' = u(4), so
            % INDEXSTART = [1, 3], since v and its derivatives starts getting
            % counted at 1, and w and its derivatives start getting counted at
            % 3.
            %
            % The vector INDEXSTARTDER is similar, but here, we also assign the
            % highest order derivative of each variable it's own index. This is
            % so that later on, we can correctly evaluate the coefficient
            % multiplying the highest order derivative in problem. Thus, in the
            % example above, we'd have
            %   v = u(1), v' = u(2), v'' = u(3), w = u(4), w' = u(5), w'' = u(6)
            % Here, INDEXSTARTDER = [1 4].
            indexStart = zeros(1, numArgs);
            indexStartDer = indexStart;
            totalDiffOrders = indexStart;
            for wCounter = 1:length(fevalResult)
                % We use cumsum() to look at how many derivatives have appeared
                % already in the problem.
                newIndex =  [1 (cumsum(fevalResult(wCounter).tree.diffOrder(1:end-1)) + (1:(numArgs-1)))];
                newIndexDer =  ...
                    [1 (cumsum(fevalResult(wCounter).tree.diffOrder(1:end-1) + 1) + (1:(numArgs-1)))];
                indexStart = max(indexStart, newIndex);
                indexStartDer = max(indexStartDer, newIndexDer);
                totalDiffOrders = max(totalDiffOrders, fevalResult(wCounter).tree.diffOrder);
            end
            
            % COEFFARG will be used to evaluate the functions that gives us
            % information about the coefficients multiplying the highest order
            % derivative in each equation. The vector has to be equally long to
            % the total number of derivatives appearing in the problem; we'll
            % then change one of the entries to 1 at a time to get the
            % coefficient information.
            coeffArg = zeros(1, indexStartDer(end) + totalDiffOrders(end));
            
            % Go through each componenent from the result of evaluating FUNIN,
            % and change it to infix format.
            for wCounter = 1:length(fevalResult)
                
                % The current result we're looking at.
                res = fevalResult(wCounter);
                % Current diffOrders
                diffOrders = res.tree.diffOrder;
                
                % Expand the tree, so that PLUS rather than TIMES is sitting at
                % the top of it.
                expTree = treeVar.expandTree(res.tree, diffOrders);
                
                % Split the tree into derivative part and non-derivative part.
                [newTree, derTree] = treeVar.splitTree(expTree, diffOrders);
                
                % Convert the derivative part to infix form.
                [infixDer, dummy, varArrayDer] = ...
                    treeVar.tree2infix(derTree, wCounter, indexStartDer);
                
                % Find what argument corresponds to the highest derivative one
                % in the current expression we're looking at:
                maxDerLoc = find(expTree.diffOrder == max(diffOrders));
                % Convert the infix form of the expression that gives us the
                % coefficient multiplying the highest order derivative appearing
                % in the expression to an anonymous function we can evaluate:
                coeffFun = treeVar.toAnon(infixDer, varArrayDer);
                
                % Reset coeffArg for next evaluation:
                coeffArg = 0*coeffArg;
                
                % Replace one of the 0s in coeffFun with 1 so that we can
                % evaluate COEFFFUN:
                if ( maxDerLoc == numArgs )
                    % The last variable in the system currently appears in the
                    % highest order derivate.
                    coeffArg(end) = 1;
                else
                    % The variable with index maxDerLoc+1 is the next variable
                    % we need to start numbering at. So subtract 1 for the index
                    % of the highest derivate we're currently interested in.
                    coeffArg(indexStartDer(maxDerLoc+1) - 1) = 1;
                end
                
                % Evaluate the COEFFFUN to the coefficient!
                coeffs{wCounter} = coeffFun(t, coeffArg);
                
                % Now work with the remaining syntax tree of the current
                % expression of interest. We need to negate the syntax tree as
                % we're moving it to the right-hand side. But if it already
                % starts with a unary minus, we can simply remove it rather than
                % doing a double negation:
                % [TODO: Remove double UMINUS]
                newTree = struct('method', 'uminus', ...
                    'numArgs', 1, 'center', newTree);
                % Convert current expression to infix form:
                [infix, varCounter, varArray] = ...
                    treeVar.tree2infix(newTree, wCounter, indexStart);
                % Store the infix form and the variables that appeared in the
                % anonymous function.
                systemInfix{wCounter} = infix;
                varArrays{wCounter} = varArray;
            end
            
            % Convert all the infix expressions, coefficients and variables
            % stored to an anonymous function we can evaluate and use as the RHS
            % of our ODE system:
            funOut = treeVar.toRHS(systemInfix, varArrays, coeffs, ...
                indexStart, totalDiffOrders);

        end
        
    end
    
end