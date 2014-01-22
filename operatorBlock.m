classdef operatorBlock < linBlock
%OPERATORBLOCK  Linear map of function to function.
%   This class is not intended to be called directly by the end user.
%
%   See also LINBLOCK, LINOP, CHEBOP, CHEBOPPREF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer notes
%
% One of the two concrete implementations of the abstract type LINBLOCK.
% Operators can be composed using *, added, exponentiated, and applied to
% CHEBFUN objects.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods
        function A = operatorBlock(domain)
            % OPERATORBLOCK constructor, simply calls the LINBLOCK
            % constructor.
            A = A@linBlock(domain);
        end
        
        function varargout = size(A, dim)
            % SIZE  Size of a OPERATORBLOCK.
            % The commands
            %   S = SIZE(A)
            %   [M, N] = SIZE(A)
            %   P = SIZE(A, K)
            % return the results expected from standard MATLAB syntax.
            
            % An OPERATORBLOCK is always of dimensions Inf x Inf.
            m = [Inf Inf];

            if ( nargin > 1 )
                varargout = {m(dim)};
            elseif ( nargout <= 1 )
                varargout = {m};
            else
                varargout = {m(1) m(2)};
            end
        end

        function C = uminus(A)
            % Unary minus of an OPERATORBLOCK.
            C = operatorBlock(A.domain);
            C.stack = @(z) -A.stack(z);
            C.diffOrder = A.diffOrder; 
        end
        

        function C = mtimes(A, B)
            % *    Operator composition, multiplication or application.
            %
            % C = A*B, where A is a OPERATORBLOCK and B is an OPERATORBLOCK,
            % returns the OPERATORBLOCK C that is the the result of composing
            % the operators A and B.
            %
            % C = A*B, or C = B*A,  where A is a OPERATORBLOCK and B is a
            % scalar, returns the OPERATORBLOCK C that is the the result of
            % multiplying A with B.
            %
            % C = A*F, where A is a OPERATORBLOCK and F is a CHEBFUN, returns
            % the CHEBFUN C which is the the result of applying A to F.
            
            % Which case?
            if ( isa(B, 'chebfun') )
                % Convert C to a callable anonymous function
                C = toFunction(A);
                % And apply it to the CHEBFUN
                C = C(B);
            else
                % A scalar is converted into a constant CHEBFUN, which can then
                % be used to create a multiplication OPERATORBLOCK.
                if ( isnumeric(A) )
                    A = operatorBlock.mult( chebfun(A, B.domain) );
                elseif ( isnumeric(B) )
                    B = operatorBlock.mult( chebfun(B, A.domain) );
                end
                
                % Create an OPERATORBLOCK to be returned.
                C = operatorBlock(A.domain);
                
                % The instantiation class must recognize mtimes as a
                % functional composition.
                C.stack = @(z) A.stack(z) * B.stack(z);

                % Difforder of returned OPERATORBLOCK.
                C.diffOrder = A.diffOrder + B.diffOrder;
            end
        end
        
        function C = plus(A, B)
            % +    Operator addition.
            % 
            % C = A + B, where A and B are both OPERATORBLOCK objects return the
            % OPERATORBLOCK C that is the result of adding the operators A and
            % B.
            %
            % C = A + B, or C = B + A, A is an OPERATORBLOCK and B is a scalar
            % is interpreted as C = A + B*I, where I is the identity operator on
            % the domain A is defined on.

            % Did we get passed a scalar?
            if ( isnumeric(A) )
                A = A*operatorBlock.eye(B.domain);
            elseif ( isnumeric(B) )
                B = B*operatorBlock.eye(A.domain);
            end
            
            % Operator addition.
            dom = chebfun.mergeDomains(A.domain, B.domain);
            C = operatorBlock(dom);
            C.stack = @(z) A.stack(z) + B.stack(z);
            C.diffOrder = max(A.diffOrder, B.diffOrder);
        end
        
        function B = mpower(A, pow)
            % ^    Repeated application of an operator.
            % 
            % C = A^M, where A is an OPERATORBLOCK and M a nonnegative integer,
            % returns the OPERATORBLOCK representing M-fold application of A.
            
            % Check second argument is accepted.
            if ( pow ~= round(pow) || pow < 0 )
                error('Power must be a positive integer.')
            end
            
            % Construct OPERATORBLOCK for repeated application.
            B = operatorBlock.eye(A.domain);
            for i = 1:pow
                B = B*A;
            end
        end
        
        function out = iszero(A)
            % TODO: Implement this.
            out = false;
        end
        
    end
    
    methods (Static = true)
        
        function C = cumsum(varargin)
            % OPERATORBLOCK.CUMSUM  Antiderivative operator.
            %
            % C = OPERATORBLOCK.CUMSUM returns the first-order antiderivative
            % operator C for functions defined on [-1, 1]. The result of
            % applying the operator is defined uniquely by having value zero at
            % the left endpoint.
            %
            % C = OPERATORBLOCK.CUMSUM(DOMAIN) returns the first-order
            % antiderivative operator C which applies to functions defined on
            % DOMAIN, which may include breakpoints.
            %
            % C = OPERATORBLOCK.CUMSUM(DOMAIN, M) is the mth-repeated
            % antiderivative.
            
            % Use inputParser to parse the arguments to the method.
            p = inputParser;
            pref = cheboppref;
            addOptional(p, 'domain', pref.domain, @isnumeric);
            mcheck = @(m) validateattributes(m, {'numeric'}, ...
                {'scalar', 'nonnegative', 'integer'});
            addOptional(p, 'm', 1, mcheck);
            parse(p, varargin{:});
            dom = p.Results.domain;
            m = p.Results.m;
            
            % Create the OPERATORBLOCK with information now available.
            C = operatorBlock(dom);
            C.stack = @(z) cumsum(z, m);
            C.diffOrder = -m;
        end
        
        function D = diff(varargin)
            % OPERATORBLOCK.DIFF  Differentiation operator.
            %
            % D = OPERATORBLOCK.DIFF returns the first-order differentation
            % operator D for functions defined on [-1, 1].
            %
            % D = OPERATORBLOCK.DIFF(DOMAIN) returns the first-order
            % differentation operator D which applies to functions defined on
            % DOMAIN, which may include breakpoints.
            %
            % D = OPERATORBLOCK.DIFF(DOMAIN, M) is the mth order derivative.
            
            % Use inputParser to parse the arguments to the method.
            p = inputParser;
            pref = cheboppref;
            addOptional(p, 'domain', pref.domain, @isnumeric);
            mcheck = @(m) validateattributes(m, ...
                {'numeric'}, {'scalar', 'nonnegative', 'integer'});
            addOptional(p, 'm', 1, mcheck);
            parse(p, varargin{:});
            dom = p.Results.domain;
            m = p.Results.m;
            
            % Create the OPERATORBLOCK with information now available.
            D = operatorBlock(dom);
            D.stack = @(z) diff(z, m);
            D.diffOrder = m;
        end
        
        function I = eye(domain)
            % OPERATORBLOCK.EYE  Identity operator.
            %
            % I = OPERATORBLOCK.EYE(DOMAIN) returns the identity operator for
            % functions on the domain DOMAIN.
            pref = cheboppref;
            if ( nargin == 0 )
                domain = pref.domain; 
            end
            
            % Create the OPERATORBLOCK with information now available.
            I = operatorBlock(domain);
            I.stack = @(z) eye(z);
            I.diffOrder = 0;
        end

        function M = mult(u)
            % OPERATORBLOCK.MULT  Multiplication operator.
            %
            % M = OPERATORBLOCK.MULT(U) returns the multiplication operator from
            % the CHEBFUN U, i.e. the operator that maps a CHEBFUN f(x) to
            % u(x)f(x).
            
            % Create the OPERATORBLOCK with information now available.
            M = operatorBlock(u.domain);
            M.stack = @(z) mult(z, u);
            M.diffOrder = 0;
        end
        
        function Z = zeros(domain)
            % OPERATORBLOCK.ZEROS  Zero operator.
            %
            % Z = OPERATORBLOCK.ZEROS(DOMAIN) returns the zero operator for
            % functions on the domain DOMAIN (i.e., the operator that maps all
            % functions to the zero function on the DOMAIN).
            pref = cheboppref;
            if ( nargin == 0 )
                domain = pref.domain;
            end
            
            % Create the OPERATORBLOCK with information now available.
            Z = operatorBlock(domain);
            Z.stack = @(z) zeros(z);
            Z.diffOrder = 0;
        end
        
        function F = fred(domain, kernel, varargin)
            F = operatorBlock(domain);
            F.stack = @(z) fred(z, kernel, varargin{:});
            F.diffOrder = 0;
        end
        
        function V = volt(domain, kernel, varargin)
            V = operatorBlock(domain);
            V.stack = @(z) volt(z, kernel, varargin{:});
            V.diffOrder = 0;
        end
    end
end
