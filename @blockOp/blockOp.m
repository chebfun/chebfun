classdef blockOp 
    
    properties (Access=public)
        % This property is assigned the callable function that does the
        % correct operation when called on chebfuns.
        func = [];
    end
    
    methods

        function A = blockOp(f)
            A.func = f;  
        end
        
        % ** Implementation of abstract methods **
        
        % Required operators.
        
        function D = diff(A,m)
            D = blockOp( @(z) diff(z,m) );
        end
        
        function C = cumsum(A,m)
            C = blockOp( @(u) cumsum(u,m) );
        end
        
        function I = eye(A)
            I = blockOp( @(z) z );
        end
        
        function Z = zeros(A)
            Z = blockOp( @(z) chebfun(0,domain) );
        end
        
        function F = diag(A,f)
            F = blockOp( @(z) times(f,z) );
        end    
        
        % Required functionals.
        
        function S = sum(A)
            S = blockOp( @(z) sum(z) );
        end
        
        function E = feval(A,location,direction)
            if (direction < 0)
                E = blockOp( @(u) feval(u,location,'left') );
            elseif (direction > 0)
                E = blockOp( @(u) feval(u,location,'right') );
            else
                E = blockOp( @(u) feval(u,location) );
            end
        end
        
        function F = inner(A,f)
            F = blockOp( @(z) mtimes(f',z) );
        end

        
        % ** Additional operations **
        
        function C = mtimes(A,B)
            C = blockOp( @(z) A.func(B.func(z)) );
        end
        
        function C = plus(A,B)
            C = blockOp( @(z) A.func(z) + B.func(z) );
        end
        
        function C = uminus(A)
            C = blockOp( @(z) -A.func(z) );
        end
        
    end
    
end
        