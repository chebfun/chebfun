classdef blockFunction
    
    properties (Access=public)
        % This property is assigned the callable function that does the
        % correct operation when called on chebfuns.
        func = [];
        domain
    end
    
    methods

        function A = blockFunction(varargin)
            if isempty(varargin{1})
                % This call is used to create an empty object of the class,
                % so that its methods will be used to process the stack.
                return
            elseif isa(varargin{1},'linBlock')
                % Convert the given linBlock to its function form by
                % evaluating its stack.
                L = varargin{1};
                dummy = blockFunction([]);
                dummy.domain = L.domain;
                A = L.stack( dummy );
            else
                % Called with data. Create a regular object. 
                A.func = varargin{1};
            end
        end
        
        % ** Implementation of abstract methods **
        
        % Required operators.
        
         
        % ** Additional operations **
        
        function C = mtimes(A,B)
            C = blockFunction( @(z) A.func(B.func(z)) );
        end
        
        function C = plus(A,B)
            C = blockFunction( @(z) A.func(z) + B.func(z) );
        end
        
        function C = uminus(A)
            C = blockFunction( @(z) -A.func(z) );
        end
        
    end
    
    methods
        function D = diff(A,m)
            D = blockFunction( @(z) diff(z,m) );
        end
        
        function C = cumsum(A,m)
            C = blockFunction( @(u) cumsum(u,m) );
        end
        
        function I = eye(A)
            I = blockFunction( @(z) z );
        end
        
        function Z = zeros(A)
            Z = blockFunction( @(z) chebfun(0,A.domain) );
        end
        
        function z = zero(A)
            z = blockFunction( @(u) 0 );
        end

        function F = mult(A,f)
            F = blockFunction( @(z) times(f,z) );
        end    
                
        function S = sum(A)
            S = blockFunction( @(z) sum(z) );
        end
        
        function E = feval(A,location,direction)
            if (direction < 0)
                E = blockFunction( @(u) feval(u,location,'left') );
            elseif (direction > 0)
                E = blockFunction( @(u) feval(u,location,'right') );
            else
                E = blockFunction( @(u) feval(u,location) );
            end
        end
        
        function F = inner(A,f)
            F = blockFunction( @(z) mtimes(f',z) );
        end


    end
    
end
        