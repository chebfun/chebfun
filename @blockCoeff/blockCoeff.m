classdef blockCoeff 
    properties (Access=public)
        coeffs = [];
        domain = [-1 1];
    end
    
    methods
        function A = blockCoeff(f,domain)
            if ( ~iscell(f) )
                f = {f};
            end
            A.coeffs = f;
            A.domain = domain;
        end
        
        % Required operators.
        
        function D = diff(A,order)
            
            if ( nargin < 2 )
                order = 1;
            end

            if ( order == 0 )
                D = eye(A);
                return
            end
            
            % If starting from nothing, start with the identity.
            if ( isempty(A.coeffs{1}) )
                A = eye(A);
            end
            
            % Differentiate the correct number of times.
            c = A.coeffs;
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
        
        function I = eye(A)
            I = blockCoeff( chebfun(1,A.domain), A.domain);
        end
        
        function I = zeros(A)
            I = blockCoeff( chebfun(0,A.domain), A.domain);
        end
        
        function F = diag(A, f)
            F = blockCoeff(f,f.domain);
        end
        
        function C = cumsum(A,m)
            error('Not available in coeff form.')
        end
        
        % Additional operations.
        
        function C = mtimes(A, B)
            C = B;
            for k = 1:numel(B.coeffs)
                C.coeffs{k} = A.coeffs{end}.*C.coeffs{k};
            end
            % CHECK: Can the coeffs property ever be empty?
            z = {chebfun(0,A.coeffs{1}.domain)};
            for j = 1:numel(A.coeffs)-1
                B = diff(B);
                C.coeffs = [z, C.coeffs];
                for k = 1:numel(B.coeffs)
                    C.coeffs{k} = C.coeffs{k} + A.coeffs{end-j}.*B.coeffs{k};
                end
            end
        end
        
        function C = plus(A, B)
            sA = numel(A.coeffs);
            sB = numel(B.coeffs);
            % CHECK: Can the coeffs property ever be empty?
            z = {chebfun(0,A.coeffs{1}.domain)};
            if ( sA < sB )
                A.coeffs = [repmat(z, 1, sB-sA), A.coeffs];
                sA = sB;
            elseif ( sB < sA )
                B.coeffs = [repmat(z, 1, sA-sB), B.coeffs];
            end
            c = cell(1, sA);
            for k = 1:sA
                c{k} = A.coeffs{k} + B.coeffs{k};
            end
            C = blockCoeff(c,A.domain);
        end
        
        function A = uminus(A)
            for k = 1:numel(A.coeffs)
                A.coeffs{k} = -A.coeffs{k};
            end
        end

    end
    
end
        
