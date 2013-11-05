classdef blockMatrix ( Abstract ) < operatorBlockRealization & functionalBlockRealization
    
    properties (Access=public)
        % This property is assigned the callable function that does the
        % correct operation when called on chebfuns.
        matfunc = [];
    end
    
    properties (Dependent)
        result
    end
    
    methods

        function A = blockMatrix(f)
            A.matfunc = f;  
        end
       
        function f = get.result(A)
            f = A.matfunc;
        end
        
       
        function D = diff(A,m)
            D = blockMatrix( @(n) diff(n,m) );
        end
        
        function C = cumsum(A,m)
            C = blockMatrix( @(n) cumsum(n,m) );
        end
        
        function I = eye(A)
            I = blockMatrix( @(n) eye(n) );
        end
        
        function Z = zeros(A)
            Z = blockMatrix( @(z) chebfun(0,domain) );
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
        
        function C = mtimes(A, B)
            C = blockMatrix( @(n) A.matfunc(n) * B.matfunc(n) );
        end
        
        function C = plus(A, B)
            C = blockMatrix( @(n) A.matfunc(n) + B.matfunc(n) );
        end
        
        function B = uminus(A)
            B = blockMatrix( @(n) -A.matfunc(n) );
        end

         
    end
    
end
        