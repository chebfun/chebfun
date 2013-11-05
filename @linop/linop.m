classdef linop
    
    properties
        operator  % chebmatrix
        constraint 
        discretizationType = linBlock.defaultDiscretization;
    end
    
    properties (Dependent)
        domain
        blockDiffOrders
    end
    
    methods
        function L = linop(M, C)
            % TODO: check size, inputs
            if ( isa(M, 'linBlock') )
                M = chebmatrix({M});
            end
            L.operator = M;
            if ( nargin < 2 )
                L.constraint = linopConstraint();
            else
                L.constraint = C;
            end
            
        end
        
        function d = get.domain(L)
            d = L.operator.domain;
        end
        
        function d = get.blockDiffOrders(L)
            [m,n] = size(L);
            d = zeros(m,n);

            for i = 1:m
                for j = 1:n
                    block = L.operator.blocks{i,j};
                    if isa(block,'linBlock')
                        d(i,j) = block.diffOrder;
                    end
                end
            end
        end
        
        function varargout = size(L)
            [varargout{1:nargout}] = size(L.operator);
        end

        
        function L = addbc(L, varargin)
            L.constraint = append(L.constraint, varargin{:});
        end
        
        function L = bc(L, c)
            validateattributes(c, {'linopConstraint'})
            L.constraint = c;
        end
        
        function L = addlbc(L, op, value)      
            if ( nargin < 3 )
                value = 0;
            end
            d = L.operator.domain;
            E = linop.feval(d(1), d);
            L = addbc(L, E*op, value);
        end
        
        function L = addrbc(L, op, value)
            if ( nargin < 3 )
                value = 0;
            end
            d = L.operator.domain;
            E = linop.feval(d(end), d);
            L = addbc(L, E*op, value);          
        end
        
        function u = mldivide(L, f)
            u = linsolve(L, f);
        end
 
    end

    % These are provided as more convenient names than the linBlock equivalents.
    methods (Static)
        function D = diff(varargin)
            D = linBlock.diff(varargin{:});
        end
        
        function C = cumsum(varargin)
            C = linBlock.cumsum(varargin{:});
        end
        
        function I = eye(varargin)
            I = linBlock.eye(varargin{:});
        end
        
        function Z = zeros(varargin)
             Z = linBlock.zeros(varargin{:});
        end
        
        function U = mult(varargin)
            U = linBlock.mult(varargin{:});
        end
        
        function Z = zero(varargin)
            Z = linBlock.zero(varargin{:});
        end
        
        function S = sum(varargin)
            S = linBlock.sum(varargin{:});
        end
        
        function E = feval(varargin)
            E = linBlock.feval(varargin{:});
        end
        
        function E = eval(varargin)
            E = linBlock.eval(varargin{:});
        end
        
        function F = inner(varargin)
            F = linBlock.inner(varargin{:});
        end

        function F = dot(varargin)   % synonym for inner()
            F = linBlock.dot(varargin{:});
        end

    end
    
    methods
        [A, b, dom] = linSystem(L, f, dim, matrixType)
        u = linsolve(L, f, type)
    end
    
    methods (Access = private)
        % Find the differential orders of each equation (row).
        d = getRowDiffOrders(L)
        
        % Find the differential orders of each variable (column).
        d = getColDiffOrders(L)
        
        % Figure out how much to reduce dimension in each equation.
        d = getDownsampling(L)
        
        % Construct operators for generic continuity at each breakpoint.
        C = domainContinuity(L, maxorder)
 
        % Append proper breakpoint continuity conditions to a linear system. 
        L = appendContinuity(L)            
    end
    
end