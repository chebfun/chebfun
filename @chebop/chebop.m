classdef (InferiorClasses = {?double}) chebop
%CHEBOP  CHEBOP class for representing operators on functions defined on [a,b].
   
    properties ( GetAccess = 'public', SetAccess = 'public' )
        domain = [];    % Domain of the operator
        op = [];        % The operator
        lbc = [];       % Left boundary condition(s)
        rbc = [];       % Right boundary condition(s)
        bc = [];        % Other/internal/mixed boundary conditions
        init = [];      % Initial guess of a solution
        % Default discretization for linear problems
        discretizationType = @colloc2; 
    end
    
    methods
        
        function N = chebop(op, dom)
            if ( nargin == 0 )
                return
            end
            if ( nargin < 2 )
                % Need to access chebfunpref to create an operator on the
                % default domain if none is passed.
                p = chebpref();
                dom = p.domain;
            end
            
            N.op = op;
            N.domain = dom;
            
        end
        
        function u = mldivide(N, rhs)
            pref = cheboppref;
            u = solvebvp(N, rhs, pref);
        end
        
        function nin = nargin(N)
            nin = nargin(N.op);
        end
        
        function N = set.lbc(N, val)
            nin = nargin(N);
            if isnumeric(val)
                if nin <= 2
                    N.lbc = @(u) u - val;
                else
                    error('CHEBFUN:CHEBOP:SETLBC', ...
                    'Can only assign scalar BCs to scalar problems');
                end
            elseif isa(val,'function_handle')
                if ( ( nin == 1 && nargin(val) == 1) || ( nin == nargin(val) + 1) )
                    N.lbc = val;
                else
                    error('CHEBFUN:CHEBOP:SETLBC', ...
                    'Number of inputs to BCs do not match operator.');
                end
            else
                error('CHEBFUN:CHEBOP:SETLBC', ...
                    'Unsupported format of BCs')
            end
        end
        
        function N = set.rbc(N, val)
            nin = nargin(N);
            if isnumeric(val)
                if nin <= 2
                    N.rbc = @(u) u - val;
                else
                    error('CHEBFUN:CHEBOP:SETRBC', ...
                    'Can only assign scalar BCs to scalar problems');
                end
            elseif isa(val,'function_handle')
                if ( ( nin == 1 && nargin(val) == 1) || ( nin == nargin(val) + 1) )
                    N.rbc = val;
                else
                    error('CHEBFUN:CHEBOP:SETRBC', ...
                    'Number of inputs to BCs do not match operator.');
                end
            else
                error('CHEBFUN:CHEBOP:SETRBC', ...
                    'Unsupported format of BCs')
            end
        end
        
        function N = set.bc(N, val)
            if isnumeric(val)
                error('CHEBFUN:CHEBOP:SETBC', ...
                    'Can not assign numerical BCs to .bc field.');
            elseif isa(val,'function_handle')
                if ( nargin(N) == nargin(val) )
                    N.bc = val;
                else
                    error('CHEBFUN:CHEBOP:SETBC', ...
                    'Number of inputs to BCs must match operator.');
                end
            else
                error('CHEBFUN:CHEBOP:SETRBC', ...
                    'Unsupported format of BCs')
            end
        end
        
        [L, res, isLinear] = linearise(N, x, u, flag);
        
        
    end
        
    methods (Static = true) % These should be private methods as well
        
        [displayFig, displayTimer] = displayInfoInit(u,pref);
        
        displayInfoIter(u, delta, iterNo, normdu, cFactor, errEst, lendu, ...
            lambda, lenu, displayFig, displayTimer, pref);
        
        displayInfoFinal(u, delta, iterNo, errEstDE, errEstBC, displayFig, ...
            displayTimer, pref)
        
        function newRHS = convertToRHS(rhs, residual)
            [numRow, numCol] = size(residual);
            
            
            if ( length(rhs) == 1 )
                % Allow a scalar RHS to be converted to a RHS of correct
                % dimensions
                rhs = repmat(rhs, numRow, numCol);
            elseif ~( all(size(rhs) == [numRow, numCol]) )
                error('CHEBFUN:CHEBOP:CONVERTTORHS', ...
                    'RHS does not match output dimensions of operator.');
            end
            
            rhsBlocks = cell(numRow, numCol);
            resBlocks = residual.blocks;
            
            dom = getDomain(residual);
            
            % Convert numerical values in RHS vector into chebmatrix
            for rhsCounter = 1:numRow
                if isa(resBlocks{rhsCounter}, 'chebfun')
                    % If corresponding block in the residual is a chebfun, the
                    % rhs must also be made to be a chebfun
                    rhsBlocks{rhsCounter} = chebfun(rhs(rhsCounter),dom);
                else
                    % Otherwise, the entry in the chebmatrix will be a scalar
                    rhsBlocks{rhsCounter} = rhs(rhsCounter);
                end
            end
            % Convert the cell-array to a chebmatrix
            newRHS = chebmatrix(rhsBlocks);
        end
    end
    
end

