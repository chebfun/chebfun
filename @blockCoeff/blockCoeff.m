classdef blockCoeff
%BLOCKCOEFF   Class to convert linear operator to derivative coefficents.
%   This class is not intended to be called directly by the end user.
%
% See also LINOP, CHEBOP, CHEBOPPREF.
    
% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer notes
%
% This class converts a linBlock object into a list of coefficients for the
% (descending) powers of the derivative in the operator. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        coeffs = [];
        domain
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function A = blockCoeff(varargin)
            % When called with no arguments, the returned object causes the
            % block's stack to be evaluated with these methods to produce
            % coefficients.
            if isempty(varargin{1})
                pref = cheboppref;
                A.domain = pref.domain;
                return
                
            % Calling the constructor with a linBlock argument initiates the
            % process of evaluating the stack with a dummy object of this class.
            elseif isa(varargin{1},'linBlock')
                L = varargin{1};
                dummy = blockCoeff([]);
                dummy.domain = L.domain;
                A = L.stack( dummy );

            % If the constructor is called with data, just make a regular object
            % out of it. 
            else
                f = varargin{1};
                if ( ~iscell(f) )
                    f = {f};
                end
                A.coeffs = f;
                A.domain = varargin{2};
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% OTHER BASIC CONSTRUCTORS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        % These are the basic constructors.
        
        function I = eye(A)
            I = blockCoeff( chebfun(1, A.domain), A.domain );
        end
        
        function I = zeros(A)
            I = blockCoeff( chebfun(0, A.domain), A.domain );
        end
        
        function F = mult(A, f)
            F = blockCoeff( f, A.domain );
        end
        
        function C = cumsum(A, m)
            error('CHEBFUN:BLOCKCOEFF:cumsum:notSupported', ...
                'Conversion of integration to coefficients is not supported.')
        end
        
        function D = diff(A, order)
            
            if ( nargin < 2 )
                order = 1;
            end
            
            if ( ~isempty(A.coeffs) )
                % This syntax is used by the mtimes method to differentiate
                % an existing operator.
                D = A;
            else
                % Initialize with an identity. 
                D = eye(A);
            end
            
            % Differentiate the requested number of times.
            c = D.coeffs;
            for d = 1:order
                m = numel(c);
                c = c([1, 1:m]);
                for k = 2:m
                    c{k} = diff(c{k}) + c{k+1};
                end
                c{m+1} = diff(c{m+1});
            end
            D = blockCoeff(c, A.domain);
            
        end
        
        function S = sum(A)
            error('CHEBFUN:BLOCKCOEFF:sum:notSupported', ...
                'Conversion of integration to coefficients is not supported.')
        end
        
        function E = feval(A, location, direction)
            error('CHEBFUN:BLOCKCOEFF:feval:notSupported', ...
                'Conversion of evaluation to coefficients is not supported.')
        end
        
        function F = inner(f)
            error('CHEBFUN:BLOCKCOEFF:inner:notSupported', ...
                'Conversion of inner product to coefficients is not supported.')
        end
        
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
    
        function C = mtimes(A, B)
            
            if ( isnumeric(A) )
                % Allow multiplying a BLOCKCOEFF with a scalar.
                c = B.coeffs;
                for k = 1:numel(B.coeffs)
                    c{k} = A*c{k};
                end
                
                % Create the result.
                C = blockCoeff(c, B.domain);
                return
            elseif ( isnumeric(B) )
                C = mtimes(B, A);
                return
            end
            
            if ( isempty(A.coeffs) || isempty(B.coeffs) )
                C = blockCoeff([]);
                return
            end
            
            % Initialize constant term of A times the coeffs of B.
            c = B.coeffs;
            for k = 1:numel(B.coeffs)
                c{k} = A.coeffs{end}.*c{k};
            end
            
            % Do a convolution-style multiplication for the rest of A.
            if ( numel(A.coeffs) > 1 )
                z = {0*c{1}};  % an appropriate zero chebfun
            end
            for j = 1:numel(A.coeffs)-1
                B = diff(B);  % differentiate this operator
                c = [z, c];
                for k = 1:numel(B.coeffs)
                    c{k} = c{k} + A.coeffs{end-j}.*B.coeffs{k};
                end
            end
            
            % Create the result.
            C = blockCoeff(c, A.domain);
        end
        
        function C = plus(A, B)
            sA = numel(A.coeffs);
            sB = numel(B.coeffs);
            
            % Empty operand returns empty result. 
            if ( (sA == 0) || (sB == 0) )
                C = blockCoeff([]);
                return
            end
            
            % Get to same length.
            if ( sA < sB )
                z = {0*A.coeffs{1}};
                A.coeffs = [repmat(z, 1, sB-sA), A.coeffs];
                sA = sB;
            elseif ( sB < sA )
                z = {0*B.coeffs{1}};
                B.coeffs = [repmat(z, 1, sA-sB), B.coeffs];
            end
            
            % Do the work.
            c = cell(1, sA);
            for k = 1:sA
                c{k} = A.coeffs{k} + B.coeffs{k};
            end
            C = blockCoeff(c, A.domain);
        end
        
        function A = uminus(A)
            A.coeffs = cellfun( @uminus, A.coeffs, 'uniform', false );
         end
        
        function A = uplus(A)
        end
        
    end
    
end