classdef linBlock
    
    properties
        domain = [];
        
        % This stores the delayed evaluation stack. Its argument is an empty
        % instance of the class that determines how the operator is
        % instantiated (collocation matrix, function handle, ...).
        stack = [];
        
        % Track differential order.
        diffOrder = 0;
        
    end
    
    properties (Constant)
        % Used whenever a matrix is required but the type is not specified.
        % It doesn't work as a set/get property, because different linops
        % within one chebmatix can't have different defaults.
        defaultDiscretization = @blockColloc2;
%         defaultDiscretization = @blockUS;
    end
    
    
    methods
        function A = linBlock(domain)
            % A = LINOP()        Domain defaults to [-1, 1].
            % A = LINOP(DOMAIN)
            if ( nargin == 0 )
                domain = [-1 1];
            elseif ( nargin == 1 && isa(domain, 'linBlock') )
                A = domain;
                return
            end
            A.domain = domain;
        end
        
        function d = get.diffOrder(A)
            d = A.diffOrder;
        end
        
        
        function C = minus(A, B)
            % C = A - B
            C = A + (-B);
        end
               
        function C = uminus(A)
            C = A;
            C.stack = @(z) -A.stack(z);
        end
        
        function C = horzcat(varargin)
            kill = cellfun(@isempty, varargin);            
            C = chebmatrix( varargin(~kill) );
        end
        
        function C = vertcat(varargin)
            kill = cellfun(@isempty, varargin);            
            C = chebmatrix( varargin(~kill)' );
        end
        
        
        function L = op(A)
            % OP(A) returns function handle for chebfuns
            L = A.stack( blockOp([]) );
            L = L.func;
        end
        
        function L = coeff(A)
            L = A.stack( blockCoeff([], A.domain) );
            L = [L.coeffs{:}];
        end

        function L = discretize(A, dim, varargin)
            % MATRIX(A, DIM) returns a collocation matrix using A's matrixType property.
            % MATRIX(A, DIM, DOMAIN) overrides the domain stored in A.             
            % MATRIX(A, DIM, DOMAIN, CONSTRUCTOR) overrides the matrixType constructor.            
            p = inputParser;
            addOptional(p, 'domain', A.domain, @isnumeric);
            addOptional(p, 'matrixType', linBlock.defaultDiscretization, @(x) isa(x, 'function_handle'))
%             addOptional(p, 'matrixType', 'blockCoeff', @(x) isa(x, 'function_handle'))
            parse(p, varargin{:})
            
            dom = p.Results.domain;
            matrixType = p.Results.matrixType;
            
            dummy = matrixType([]);
            L = dummy.discretize(A, dim, dom);

        end
        
        function L = matrix(A, varargin)
            L = discretize(A, varargin{:});
        end
                        
    end
        
    methods (Static)
        function D = diff(varargin)
            % LINOP.DIFF  Differentiation operator.
            %
            % LINOP.DIFF returns the first-order differentation operator for
            % functions defined on [-1, 1].
            %
            % LINOP.DIFF(DOMAIN) applies to functions defined on DOMAIN, which
            % may include breakpoints.
            %
            % LINOP.DIFF(DOMAIN, M) is the mth order derivative.
            p = inputParser;
            addOptional(p, 'domain', [-1 1], @isnumeric);
            mcheck = @(m) validateattributes(m, {'numeric'}, {'scalar', 'nonnegative', 'integer'});
            addOptional(p, 'm', 1, mcheck);
            parse(p, varargin{:});
            dom = p.Results.domain;
            m = p.Results.m;
            D = operatorBlock(dom);
            D.stack = @(z) diff(z, m);
            D.diffOrder = m;
        end
        
        function C = cumsum(varargin)
            % LINOP.CUMSUM  Antiderivative operator.
            %
            % LINOP.CUMSUM returns the first-order antiderivative operator for
            % functions defined on [-1, 1]. The result of applying the operator
            % is defined uniquely by having value zero at the left endpoint.
            %
            % LINOP.CUMSUM(DOMAIN) applies to functions defined on DOMAIN, which
            % may include breakpoints.
            %
            % LINOP.CUMSUM(DOMAIN, M) is the mth-repeated antiderivative.
            p = inputParser;
            addOptional(p, 'domain', [-1 1], @isnumeric);
            mcheck = @(m) validateattributes(m, {'numeric'}, {'scalar', 'nonnegative', 'integer'});
            addOptional(p, 'm', 1, mcheck);
            parse(p, varargin{:});
            dom = p.Results.domain;
            m = p.Results.m;
            
            C = operatorBlock(dom);
            C.stack = @(z) cumsum(z, m);
            C.diffOrder = -m;
        end
        
        function I = eye(domain)
            % I = EYE(DOMAIN)   identity operator on the same domain
            if nargin==0, domain = [-1 1]; end
            I = operatorBlock(domain);
            I.stack = @(z) eye(z);
            I.diffOrder = 0;
        end
        
        function Z = zeros(domain)
            % Z = ZEROS(DOMAIN)   zero operator on the same domain
            if nargin==0, domain = [-1 1]; end
            Z = operatorBlock(domain);
            Z.stack = @(z) zeros(z);
            Z.diffOrder = 0;
        end
        
        function Z = zero(domain)
            % Z = ZERO(DOMAIN)   zero functional on the domain
            if nargin==0, domain = [-1 1]; end
            Z = functionalBlock(domain);
            Z.stack = @(z) zero(z);
            Z.diffOrder = 0;
        end

        
        function U = diag(u)
            % D = DIAG(U)  diagonal operator from the chebfun U
            U = operatorBlock(u.domain);
            U.stack = @(z) diag(z, u);
            U.diffOrder = 0;
        end
        
        function S = sum(domain)
            % SUM(DOMAIN)  integration functional on the domain
            if nargin==0, domain = [-1 1]; end                        
            S = functionalBlock(domain);
            S.stack = @(z) sum(z);
            S.diffOrder = -1;
        end
        
        function E = feval(location, varargin)
            p = inputParser;
            addRequired(p, 'location');
            addOptional(p, 'domain', [-1 1]);
            valid = @(x) ischar(x) || isnumeric(x);
            addOptional(p, 'direction', 0, valid);
            parse(p, location, varargin{:});
            location = p.Results.location;
            dom = p.Results.domain;
            direction = p.Results.direction;
            if (location < dom(1)) || (location > dom(end))
                error('Evaluation location is not in the domain.')
            end
            
            % Convert direction argument into a number.
            if ischar(direction)
                if any( strncmpi(direction, {'left', '-'}, 1) )
                    direction = -1;
                elseif any( strncmpi(direction, {'right', '+'}, 1) )
                    direction = +1;
                else
                    error('Direction must be ''left'', ''right'', ''+'', or ''-''.')
                end
            end

            E = functionalBlock(dom);
            E.stack = @(z) feval(z, location, direction);
        end
        
        function E = eval(varargin)
            % EVAL(DOMAIN) returns a function E. The output of E(X) is an
            % evaluator at X in the domain.
            E = @(x) linop.feval(x, varargin{:});
        end
        
        function F = inner(f)
            F = functionalBlock(f.domain);
            F.stack = @(z) inner(z, f);
            F.diffOrder = 0;
        end

        function F = dot(f)   % synonym for inner()
            F = inner(f);
        end

    end
    

    
end