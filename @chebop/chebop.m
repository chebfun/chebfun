%CHEBOP  Construct an operator on chebfuns.

classdef (InferiorClasses = {?double}) chebop
    
    properties ( GetAccess = 'public', SetAccess = 'public' )
        domain = [];    % Domain of the operator
        op = [];        % The operator
        lbc = [];       % Left boundary condition(s)
        rbc = [];       % Right boundary condition(s)
        bc = [];        % Other/internal/mixed boundary conditions
        init = [];      % Initial guess of a solution
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
        
        u = mldivide(N, rhs);
        
        [L, res, isLinear] = linearise(N, x, u, flag);
        
        
    end
end

