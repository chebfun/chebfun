classdef colloc2 < linopDiscretization
    properties 
        size = [];  % arbitrary, but fixed in any one instance
        fundomain = [-1 1];
    end
    
    methods
        function A = colloc2(varargin)
            % Collocation matrix on 2nd kind points.
            
            % COLLOC2(DIM,DOMAIN) returns a dummy object that will propagate the
            % dimension size DIM and function domain DOM throughout the delayed
            % evaluation stack.
            %
            % COLLOC2(A,DIM) realizes the linop A (which knows its domain) at
            % dimension DIM.
            % 
            % COLLOC2([]) returns a dummy object that gives access to static
            % methods. 
            
            if ( nargin > 1 )
                if isa(varargin{1},'linop')
                    L = varargin{1};
                    A.size = varargin{2};
                    A.fundomain = domain(L);
                    A = L.delayFun( A );
                else
                    A.size = varargin{1};
%                     validateattributes(varargin{2},{'numeric'},{'increasing','finite'});
                    A.fundomain = varargin{2};
                end
            end
        end
        
        function d = domain(A)
            d = A.fundomain;
        end
        
        % Building blocks.
        
        D = diff(A,m) 
        C = cumsum(A,m)
        
        function I = eye(A)
            n = dim(A);
            I = eye(sum(n));
        end
        
        function Z = zeros(A)
            n = dim(A);
           Z = zeros(sum(n));
        end
        
        F = diag(A,f)
        
        % Required operators.
        function C = mtimes(A,B)
            C = A*B;
        end
        
        function C = plus(A,B)
            C = A+B;
        end
        
        function B = uminus(A)
            B = -A;
        end
        
        % Required functionals.
        
        S = sum(A)
        
        E = evalAt(A,location,direction)
        
        function F = inner(A,f)
            d = chebmatrix.mergeDomains({A,f});
            [x,w] = points(dim(A),d,2);
            F = w.*f(x);
        end
        
    end
    
    methods (Static)
        % Additional methods
        
        B = resize(A,m,n,domain)
        
        [isDone,epsLevel] = convergeTest(v)
        
        function fx = feval(f,dim,dom)
            if ( nargin < 3 )
                dom = f.domain;
            end
            x = colloc2.points(dim,dom);
            fx = f(x);
        end
        
        [x,w] = points(n,d)
        
    end
    
    methods (Access=private)
        
        D = diffmat(N,k)
       
        Q = cumsummat(N)
                
    end
end