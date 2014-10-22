classdef  (InferiorClasses = {?chebfun}) treeVar
    %TREEVAR    A class for analysing syntax trees of ODEs in Chebfun.
    %
    %   The TREEVAR class allows Chebfun to analyse the syntax trees of ODEs in
    %   CHEBFUN. It's current use is to enable Chebfun to automatically convert
    %   (systems of) higher order ODEs to coupled first order systems. This is
    %   particularly useful for initial-value problems (IVPs), as that allows
    %   Chebfun to call one of the built-in MATLAB solvers for solving IVPs via
    %   time-stepping, rather than globally via spectral methods and Newton's
    %   method in function space.
    %
    %   This class is not intended to be called directly by the end user.
    %
    %   T = TREEVAR(ID, DOMAIN), where ID is a Boolean vector corresponding to the
    %   order of variables in the problem, and DOMAIN is interval that the problem
    %   is specified on, returns the TREEVAR object V, which stores the ID and the
    %   DOMAIN. See example below for how the ID vector is specified.
    %
    %   T = TREEVAR() is the same as above, but with the default ID = 1, and
    %   DOMAIN = [-1, 1]. This is useful for quick testing purposes.
    %
    %   Example 1: Construct TREEVAR object for the scalar IVP
    %       u'' + sin(u) = 0
    %   on the interval [0, 10]:
    %       u = treeVar(1, [0 10]);
    %
    %   Example 2: Construct TREEVAR objects for the coupled IVP
    %       u'' + v = 1, u + v' = x
    %   on the interval [0, 5]:
    %       u = treeVar([1 0], [0 5]);
    %       v = treeVar([0 1], [0 5]);
    %
    % See also CHEBOP, CHEBOP/SOLVEIVP.
    
    % Copyright 2014 by The University of Oxford and The Chebfun Developers.
    % See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TREEVAR class description:
    %
    % The TREEVAR class is used by the CHEBOP class to convert higher order ODEs to
    % coupled systems of first order ODEs, which can then be solved using one of the
    % built-in MATLAB solvers, such as ODE113. This is done by evaluating the
    % (anonymous) functions in the .OP field of the CHEBOP with TREEVAR arguments,
    % which will construct a syntax tree of the mathematical expression describing
    % the operator. By then analysing the syntax tree and restructuring it
    % appropriately, conversion to a first-order system is made possible.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        % The syntax tree of a TREEVAR variable, starting from an initial
        % variable. Each syntax tree is a MATLAB struct, which contains the
        % following fields:
        %   METHOD:    The method leading to the construction of the variable.
        %   NUMARGS:   Number of arguments to the method that constructed the
        %       variable.
        %   DIFFORDER: The differential order of the TREEVAR, which represents
        %       how many times the base variable(s) have been differentiated
        %       when before we arrive at the current TREEVAR. Note that
        %       DIFFORDER is vector valued; for example, the sequence
        %           u = treeVar([1 0 0], [0 1]);
        %           v = treeVar([0 1 0], [0 1]);
        %           w = treeVar([0 0 1], [0 1]);
        %           f = diff(u) + diff(w, 2);
        %       will lead to f.DIFFORDER == [1 0 2].
        %   HEIGHT: The height of the syntax tree, i.e., the number of
        %       operations between the base variables(s) and the current
        %       variable.
        %   MULTCOEFF: The multiplication in front of the variable, which can
        %       either be a CHEBFUN or a scalar. For example, the sequence
        %           u = treeVar();
        %           v = sin(x)*u'l
        %       will have v.multcoeff == sin(x).
        %   ID: A Boolean vector, whose ith element is equal to 1 if the TREEVAR
        %       variable was constructed from the ith base variable, 0
        %       otherwise. For example, the sequence
        %           u = treeVar([1 0 0], [0 1]);
        %           v = treeVar([0 1 0], [0 1]);
        %           w = treeVar([0 0 1], [0 1]);
        %           f = u + 2*w;
        %       will lead to f.ID == [1 0 1].
        tree
        % The domain of the problem we're solving when constructing the TREEVAR
        % objects.
        domain
    end
    
    methods
        
        function obj = treeVar(IDvec, domain)
            if ( nargin > 0 )
                % Store the domain.
                obj.domain = domain;
            else
                % Default ID and domain.
                IDvec = 1;
                obj.domain = [-1 1];
            end
            % Initialise a syntax tree for a base variable.
            obj.tree  = struct('method', 'constr', 'numArgs', 0, ...
                'diffOrder', 0*IDvec, 'height', 0, 'multCoeff', 1, ...
                'ID', logical(IDvec));
        end
        
        function f = abs(f)
            f.tree = f.univariate(f.tree, 'abs');
        end
        
        function f = acos(f)
            f.tree = f.univariate(f.tree, 'acos');
        end
        
        function f = acosd(f)
            f.tree = f.univariate(f.tree, 'acosd');
        end
        
        function f = acot(f)
            f.tree = f.univariate(f.tree, 'acot');
        end
        
        function f = acoth(f)
            f.tree = f.univariate(f.tree, 'acoth');
        end
        
        function f = acsc(f)
            f.tree = f.univariate(f.tree, 'acsc');
        end
        
        function f = acscd(f)
            f.tree = f.univariate(f.tree, 'acscd');
        end
        
        function f = acsch(f)
            f.tree = f.univariate(f.tree, 'acsch');
        end
        
        function f = airy(f)
            f.tree = f.univariate(f.tree, 'airy');
        end
        
        function f = asec(f)
            f.tree = f.univariate(f.tree, 'asec');
        end
        
        function f = asecd(f)
            f.tree = f.univariate(f.tree, 'asecd');
        end
        
        function f = asech(f)
            f.tree = f.univariate(f.tree, 'asech');
        end
        
        function f = asin(f)
            f.tree = f.univariate(f.tree, 'asin');
        end
        
        function f = asind(f)
            f.tree = f.univariate(f.tree, 'asind');
        end
        
        function f = asinh(f)
            f.tree = f.univariate(f.tree, 'asinh');
        end
        
        function f = atan(f)
            f.tree = f.univariate(f.tree, 'atan');
        end
        
        function f = atand(f)
            f.tree = f.univariate(f.tree, 'atand');
        end
        
        function f = atanh(f)
            f.tree = f.univariate(f.tree, 'atanh');
        end
        
        function f = cos(f)
            f.tree = f.univariate(f.tree, 'cos');
        end
        
        function f = cosd(f)
            f.tree = f.univariate(f.tree, 'cosd');
        end
        
        function f = cosh(f)
            f.tree = f.univariate(f.tree, 'cosh');
        end
        
        function f = cot(f)
            f.tree = f.univariate(f.tree, 'cot');
        end
        
        function f = cotd(f)
            f.tree = f.univariate(f.tree, 'cotd');
        end
        
        function f = coth(f)
            f.tree = f.univariate(f.tree, 'coth');
        end
        
        function f = diff(f, k)
            % Derivative of a TREEVAR.
            
            % By default, compute first derivative:
            if ( nargin < 2 )
                k = 1;
            end
            % The derivative syntax tree.
            f.tree = struct('method', 'diff', 'numArgs', 2, ...
                'left', f.tree, 'right', k, ...
                'diffOrder', f.tree.diffOrder + k*f.tree.ID, ...
                'height', f.tree.height + 1, ...
                'ID', f.tree.ID, ...
                'multCoeff', f.tree.multCoeff);
        end
        
        function display(u)
            % Display a TREEVAR.
            
            if ( length(u) == 1 )
                % Scalar case.
                disp('treeVar with tree:')
                disp(u.tree);
                disp('and the domain:')
                disp(u.domain);
            else
                % Systems case.
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
                % (CHEBFUN/SCALAR)-TREEVAR
                h.tree = treeVar.bivariate(f, g.tree, 'minus', 1);
            elseif ( ~isa(g, 'treeVar') )
                % TREEVAR - (CHEBFUN/SCALAR)
                h.tree = treeVar.bivariate(f.tree, g, 'minus', 0);
            else
                % TREEVAR - TREEVAR
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
                % (CHEBFUN/SCALAR)+TREEVAR
                h.tree = treeVar.bivariate(f, g.tree, 'plus', 1);
            elseif ( ~isa(g, 'treeVar') )
                % TREEVAR + (CHEBFUN/SCALAR)
                h.tree = treeVar.bivariate(f.tree, g, 'plus', 0);
            else
                % TREEVAR + TREEVAR
                h.tree = treeVar.bivariate(f.tree, g.tree, 'plus', 2);
            end
            h.domain = updateDomain(f, g);
        end
        
        function plot(treeVar)
            % When we plot a TREEVAR, we plot its syntax tree.
            treeVar.plotTree(treeVar.tree);
        end
        
        function h = power(f, g)
            % POWER of a TREEVAR
            h = treeVar();
            if ( ~isa(f, 'treeVar') )
                % (CHEBFUN/SCALAR).^TREEVAR
                h.tree = struct('method', 'power', 'numArgs', 2, ...
                    'left', f, 'right', g.tree, 'diffOrder', g.tree.diffOrder, ...
                    'height', g.tree.height + 1, ...
                    'ID', g.tree.ID);
            elseif ( ~isa(g, 'treeVar') )
                % TREEVAR.^(CHEBFUN/SCALAR)
                h.tree = struct('method', 'power', 'numArgs', 2, ...
                    'left', f.tree, 'right', g, 'diffOrder', f.tree.diffOrder, ...
                    'height', f.tree.height + 1, ...
                    'ID', f.tree.ID);
            else
                % TREEVAR.^TREEVAR
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
                % (CHEBFUN/SCALAR)./TREEVAR
                h.tree = treeVar.bivariate(f, g.tree, 'rdivide', 1);
            elseif ( ~isa(g, 'treeVar') )
                % TREEVAR./(CHEBFUN/SCALAR)
                h.tree = treeVar.bivariate(f.tree, g, 'rdivide', 0);
            else
                % TREEVAR./TREEVAR
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
                % (CHEBFUN/SCALAR).^*TREEVAR
                h.tree = treeVar.bivariate(f, g.tree, 'times', 1);
            elseif ( ~isa(g, 'treeVar') )
                % TREEVAR.*(CHEBFUN/SCALAR)
                h.tree = treeVar.bivariate(f.tree, g, 'times', 0);
            else
                % TREEVAR.*TREEVAR
                h.tree = treeVar.bivariate(f.tree, g.tree, 'times', 2);
            end
            h.domain = updateDomain(f, g);
        end
        
        function f = uminus(f)
            f.tree = f.univariate(f.tree, 'uminus');
        end
        
        function f = uplus(f)
            f.tree = f.univariate(f.tree, 'uplus');
        end
        
        function dom = updateDomain(f, g)
            % UPDATEDOMAIN    Update domain in case we encounter new breakpoints
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
        
        [funOut, varIndex, problemDom] = toFirstOrder(funIn, rhs, domain)
        
        [funOut, indexStart, problemDom] = toFirstOrderSystem(funIn, rhs, domain)
        
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
