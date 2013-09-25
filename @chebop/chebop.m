%CHEBOP  Construct an operator on chebfuns.

classdef (InferiorClasses = {?double}) chebop
    
    properties ( GetAccess = 'public', SetAccess = 'public' )
        domain = [];    % Domain of the operator
        op = [];        % The operator
        lbc = [];       % Left boundary condition(2)
        rbc = [];       % Right boundary condition(s)
        bc = [];        % Other/internal/mixed boundary conditions
    end
    
    methods
        
        function N = chebop(op, dom)
            if ( nargin == 0 )
                return
            end
            if ( nargin < 2 )
                dom = chebfun.pref('domain');
            end
            
            N.op = op;
            N.domain = dom;
            
        end
        
        function u = mldivide(N, rhs)
            
            %%
            % Initialise an ADCHEBFUN:
            numVars = nargin(N.op)-1;
            zeroFun = adchebfun(0, N.domain);
            u0 = cell(numVars, 1);
            for k = 1:numVars
                u0{k} = seed(zeroFun, k, numVars);
            end

            %%
            % initialise the dependent variable:
            x = chebfun(@(x) x, N.domain);

            %%
            % Evaluate the operators to get a linearisation:
            w = N.op(x, u0{:});
            L = vertcat(w.jacobian);
            res = vertcat(w.func);

            %%
            % Add BCs
            if ( ~isempty(N.lbc) )
                lbcU0 = feval(N.lbc(u0{:}), N.domain(1));
                L = bc(L, lbcU0.jacobian, -lbcU0.func); %#ok<CPROP>
            end

            if ( ~isempty(N.rbc) )
                rbcU0 = feval(N.rbc(u0{:}), N.domain(end));
                L = bc(L, rbcU0.jacobian, -rbcU0.func); %#ok<CPROP>
            end

            if ( ~isempty(N.bc) )
                bcU0 = N.bc(x, u0{:});
                L = bc(L, bcU0.jacobian, -bcU0.func); %#ok<CPROP>
            end
            
            %%
            % Solve:
            u = L\(rhs-res);
            
        end
        
    end
end

