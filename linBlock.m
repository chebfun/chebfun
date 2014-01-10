classdef linBlock
%LINBLOCK   Linear operator on a single function.
%   This class is not intended to be called directly by the end user.
%
%   See also LINOP, CHEBOP, CHEBOPPREF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer notes
%
% A linBlock is an abstract representation of a linear operator on a single
% function defined on a fixed domain. Its main purpose is to maintain an
% unevaluated stack of algebraic steps operating on predefined building blocks.
% It also keeps track of the differential order of the operator.
%
% The reason for the 'Block' name is that these objects can be a block entry in
% a chebmatrix, for operators that apply to a mixture of multiple functions and
% scalars. The preferred usage externally is to create a 1x1 chebmatrix or linop
% rather than to manipulate the blocks themselves.
%
% The stack is used three different ways within the system: (1) Replace the
% building blocks with appropriate matrices, and then apply the algebra. (I.e.,
% collocation.) (2) Derive the coefficients of each power of the derivative.
% (See toCoeff.) (3) Convert to a callable function that can be applied directly
% to a chebfun. (See toFunction.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    properties
        % The domain of functions that are operated upon.
        domain = [];
        
        % The delayed evaluation stack. Its argument is an empty
        % instance of the class that determines how the operator is
        % discretized.
        stack = [];
                
        % Used to track differential order.
        diffOrder = 0;        
    end
    
    methods
        function A = linBlock(varargin)
            % A = LINBLOCK()        Use the default preference domain. 
            % A = LINBLOCK(B)       Self-return for the same type.
            % A = LINBLOCK(DOMAIN)  Null object on the domain.
            
            if ( nargin == 0 )
                p = chebpopref;
                A.domain = p.domain;
            elseif ( nargin == 1 && isa(varargin{1}, 'linBlock') )
                A = varargin{1};
                return
            else
                A.domain = varargin{1};
            end
        end
        
        function f = toFunction(A)
            % Convert the block to a callable function that can be applied to a
            % chebfun.
            B = blockFunction(A);
            f = B.func;
        end
        
        function c = toCoeff(A)
            % Convert the block to a chebfun of coefficients multiplying the
            % powers of the derivative. 
            B = blockCoeff(A);
            c = B.coeffs;
        end
        
        function C = minus(A, B)
            C = A + (-B);
        end
               
        function C = horzcat(varargin)
            kill = cellfun(@isempty, varargin);   % skip empties         
            C = chebmatrix( varargin(~kill) );
        end
        
        function C = vertcat(varargin)
            kill = cellfun(@isempty, varargin);   % skip empties         
            C = chebmatrix( varargin(~kill)' );
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
            pref = cheboppref;
            addOptional(p, 'domain', pref.domain, @isnumeric);
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
            pref = cheboppref;
            addOptional(p, 'domain', pref.domain, @isnumeric);
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
            pref = cheboppref;
            if nargin==0, domain = pref.domain; end
            I = operatorBlock(domain);
            I.stack = @(z) eye(z);
            I.diffOrder = 0;
        end
        
        function Z = zeros(domain)
            % Z = ZEROS(DOMAIN)   zero operator on the same domain
            pref = cheboppref;
            if nargin==0, domain = pref.domain; end
            Z = operatorBlock(domain);
            Z.stack = @(z) zeros(z);
            Z.diffOrder = 0;
        end
        
        function Z = zero(domain)
            % Z = ZERO(DOMAIN)   zero functional on the domain
            pref = cheboppref;
            if nargin==0, domain = pref.domain; end
            Z = functionalBlock(domain);
            Z.stack = @(z) zero(z);
            Z.diffOrder = 0;
        end
        
        function U = mult(u)
            % D = DIAG(U)  diagonal operator from the chebfun U
            U = operatorBlock(u.domain);
            U.stack = @(z) mult(z, u);
            U.diffOrder = 0;
        end
        
        function S = sum(domain)
            % SUM(DOMAIN)  integration functional on the domain
            pref = cheboppref;
            if nargin==0, domain = pref.domain; end
            S = functionalBlock(domain);
            S.stack = @(z) sum(z);
            S.diffOrder = -1;
        end
        
        function E = feval(location, varargin)
            p = inputParser;
            addRequired(p, 'location');
            pref = cheboppref;
            addOptional(p, 'domain', pref.domain, @isnumeric);
            valid = @(x) ischar(x) || isnumeric(x);
            addOptional(p, 'direction', 0, valid);
            parse(p, location, varargin{:});
            location = p.Results.location;
            domain = p.Results.domain;
            direction = p.Results.direction;
 
            % Sanity check. 
            if (location < domain(1)) || (location > domain(end))
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

            E = functionalBlock(domain);
            E.stack = @(z) feval(z, location, direction);
        end
        
        function E = eval(varargin)
            % EVAL(DOMAIN) returns a function E. The output of E(X) is an
            % evaluator at X in the domain.
            E = @(x) linBlock.feval(x, varargin{:});
        end
        
        function F = inner(f)
            F = functionalBlock(f.domain);
            F.stack = @(z) inner(z, f);
            F.diffOrder = 0;
        end

        function F = dot(f)   % synonym for inner()
            F = inner(f);
        end
        
        function F = fred(domain,kernel,varargin)
            F = operatorBlock(domain);
            F.stack = @(z) fred(z,kernel,varargin{:});
            F.diffOrder = 0;
        end

        function V = volt(domain,kernel,varargin)
            V = operatorBlock(domain);
            V.stack = @(z) volt(z,kernel,varargin{:});
            V.diffOrder = 0;
        end

        
    end
    

    
end