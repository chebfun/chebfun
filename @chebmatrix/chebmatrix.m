classdef (InferiorClasses = {?chebfun, ?operatorBlock, ?functionalBlock}) chebmatrix
    % No size/compatability checking whatsoever!
    
    properties
        blocks = {};
        
        discretizer = @colloc2;        
    end
    
    properties (Dependent)
        domain
        diffOrder
    end
    
    methods
        
        % Constructor.
        function A = chebmatrix(data)
            if isa(data,'chebmatrix')
                A.blocks = data.blocks;
            elseif isa(data,'chebfun') || isa(data,'linBlock')
                A.blocks = {data};
            elseif iscell(data)
                A.blocks = data;
            end
        end
        
        function d = get.domain(L)
            d = chebmatrix.mergeDomains(L.blocks); 
        end
        
        function d = getDomain(L)
            % DOMAIN(L) returns the domain on which functions are defined for
            % the chebmatrix L.
            d = L.domain;
        end       
        
        function d = get.diffOrder(L)
            d = getDiffOrder(L);
        end
        
        function k = numbc(L)
            % NUMBC(L) returns the number of constraints attached to L.
            k = length(L.constraints);
        end
        
        function t = blockClasses(L)
            t = cellfun(@class, L.blocks, 'uniform', false);
        end
        
        function varargout = size(L, varargin)
            %SIZE Number of blocks within the chebmatrix.
            % S = SIZE(L) returns both dimensions.
            % S = SIZE(L, K) returns Kth dimension (K = 1,2).
            % [M, N] = SIZE(L) returns both as scalars.
            [varargout{1:nargout}] = size(L.blocks, varargin{:});
        end
                
        function t = isempty(L)
            t = isempty(L.blocks);
        end
        
        function display(L)
            [m, n] = size(L);
            fprintf('\n  %ix%i block chebmatrix of types:\n\n', m, n)
            disp( blockClasses(L) )
        end
                
        function C = horzcat(varargin)
            C = cat(2, varargin{:});
        end
        
        function C = vertcat(varargin)
            C = cat(1, varargin{:});
        end
        
        function C = minus(A, B)
            C = plus(A, -B);
        end
        
        function u = mldivide(L, f)
            u = linsolve(linop(L), f);
        end
        
        function out = iszero(f)
            % TODO: Implement this properly (for linearity detection)
            out = false;
        end
        
        function d = getDiffOrder(A)
            d = zeros(size(A));
            for j = 1:numel(A.blocks);
                if ( isa(A.blocks{j}, 'operatorBlock') )
                    d(j) = A.blocks{j}.diffOrder;
                end
            end
        end
        
    end
      
    methods ( Access = private )
        
        % Multiply chebmatrix by scalar.
        C = scalartimes(A, z)
        
    end
    
    methods ( Static )
        
        % Union of all breakpoints, with "fuzzy" equality.
        d = mergeDomains(blocks)       
        
    end
    
end