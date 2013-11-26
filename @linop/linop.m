classdef linop < chebmatrix
    
    properties
        constraint = linopConstraint()
        continuity = linopConstraint()
        discretizer = @colloc2;      % Obtain default from (global?) pref?
    end
    
    properties (Dependent)
        sizeReduction
    end
    
    methods
        function L = linop(M)
            %LINOP   Linop constructor.
            %   Linops are typically onstructed from a LINBLOCK or a CHEBMATRIX.
                
            % TODO: Better argument checking
            L = L@chebmatrix(M);
            
         end
        
        function s = get.sizeReduction(L)
            s = getDownsampling(L);
        end
        
%         function d = get.blockDiffOrders(L)
%             [m, n] = size(L);
%             d = zeros(m,n);
%             for i = 1:m
%                 for j = 1:n
%                     block = L.operator.blocks{i,j};
%                     if isa(block,'linBlock')
%                         d(i,j) = block.diffOrder;
%                     end
%                 end
%             end
%         end
%         
%         function varargout = size(L)
%             [varargout{1:nargout}] = size(L.operator);
%         end
%         
        function L = bc(L, c)
            %BC  Set linop constraints (overwrite existing).
            validateattributes(c, {'linopConstraint'})
            L.constraint = c;
        end
        
        function L = addbc(L, varargin)
            %ADDBC  Append to linop constraints (keep existing).
            L.constraint = append(L.constraint, varargin{:});
        end
        
        function L = addlbc(L, op, varargin)
            %ADDLBC  Append to linop constraints (left BC)
            d = L.domain;
            E = linop.feval(d(1), d);
            L = addbc(L, E*op, varargin{:});
        end
        
        function L = addrbc(L, op, varargin)
            %ADDRBC  Append to linop constraints (right BC)
            d = L.domain;
            E = linop.feval(d(end), d);
            L = addbc(L, E*op, varargin{:});
        end
        
        function u = mldivide(L, f, varargin)
            u = linsolve(L, f, varargin{:});
        end
        
        function A = matrix(L,varargin)
            dsc = L.discretizationType(L,varargin{:});
            A = matrix(dsc);
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
    
end