classdef (InferiorClasses = {?chebfun, ?operatorBlock, ?functionalBlock}) chebmatrix
    % No size/compatability checking whatsoever!
    
    properties
        blocks = {};
        
        % The chebmatrix domain may contain breakpoints that no individual
        % block uses, because it merges the breakpoints.
        domain
        
    end
    
    methods
        
        % Constructor.
        function A = chebmatrix(blocks)
            A.domain = chebmatrix.mergeDomains(blocks);
            A.blocks = blocks;
        end
        
        function d = getDomain(L)
            % DOMAIN(L) returns the domain on which functions are defined for
            % the chebmatrix L.
            d = L.domain;
        end
        
        % Construct a single matrix based on DIM-dimensional blocks.
        function A = discretize(L, varargin)
            A = discretizeBlocks(L, varargin{:});
            A = cell2mat(A);
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
            %
            % S = SIZE(L) returns both dimensions.
            % S = SIZE(L, K) returns Kth dimension (K=1, 2).
            % [M, N] = SIZE(L) returns both as scalars.
            [varargout{1:nargout}] = size(L.blocks, varargin{:});
        end
        
        function varargout = blockSizes(A)
            %BLOCKSIZES Sizes of the blocks within the chebmatrix.
            %
            % BLOCKSIZES(L) returns a cell of 1x2 size vectors.
            % [M, N] = BLOCKSIZES(A) returns two matrices of row/column sizes.
            if ( nargout <= 1 )
                varargout = {cellfun(@size, A.blocks, 'uniform', false)};
            else
                varargout{1} = cellfun(@(x) size(x, 1), A.blocks);
                varargout{2} = cellfun(@(x) size(x, 2), A.blocks);
            end
            
        end
        
        function t = isempty(L)
            t = isempty(L.blocks);
        end
        
        function display(L)
            [m, n] = size(L);
            fprintf('\n  %ix%i block chebmatrix of types:\n\n', m, n)
            disp( blockClasses(L) )
        end
        
        function output = spy(A)
%             data = matrixBlocks(A, 10);
            data = discretizeBlocks(A, 10);
            h = cellplot(data);
            
            % CELLPLOT seems to cover up the text representations of double
            % values. We give a positive z-value so that they sit on top again.
            % And we hide its big ugly box.
            for i = 2:length(h)
                if strcmp(get(h(i), 'type'), 'text')
                    set(h(i), 'position', [0 0 1]+get(h(i), 'position'))
                    set(h(i-1), 'vis', 'off')
                end
            end
            
            if ( nargout > 0 )
                output = h;
            end
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
        
        function varargout = plot(L, varargin)
            
            if ( any(any(cellfun(@(L) isa(L, 'linBlock'), L.blocks))) )
                [varargout{1:nargout}] = spy(L, varargin{:});
            else
                ish = ishold;
                cols = get(gcf, 'DefaultAxesColorOrder');
                h = zeros(size(L.blocks));
                for k = 1:numel(L.blocks)
                    fk = L.blocks{k};
                    if ( ~isa(fk, 'chebfun') )
                        fk = chebfun(fk);
                    end
                    h(k) = plot(fk, varargin{:}); 
                    set(h(k), 'color', cols(k,:));
                    hold on
                end
                if ( ~ish )
                    hold off
                end
                if ( nargout > 0 )
                    varargout{1} = h;
                end
            end
        end
        
    end
    
    methods
        
        %Signatures of externally defined methods.
        
        % Fundamental algebraic operations.
        C = mtimes(A, B)
        C = plus(A, B)
        C = uminus(A)
        
        % Replace each block by its DIM-dimensional discretization.
        A = discretizeBlocks(L, dim, dom, matrixType)
        
        % Concatenation
        C = cat(n, varargin)
        
        % TODO
        B = subsref(A, sr)
        
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