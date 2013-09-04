classdef coeff 
    properties (Access=public)
        coeffs = [];
    end
    
    methods
        function A = coeff(f)
            if ( ~iscell(f) )
                f = {f};
            end
            A.coeffs = f;
        end
        
        % Required operators.
        
        function D = diff(A,domain)
            m = numel(A.coeffs);
            if ( isempty(A.coeffs{1}) )
                c = {chebfun(1,domain), chebfun(0,domain)};
            else
                A.coeffs
                c = A.coeffs([1, 1:m]);
                for k = 2:m
                    c{k} = diff(c{k}) + c{k+1};
                end
                c{m+1} = diff(c{m+1});
            end
            D = coeff( c );
        end
        
        function I = eye(A,domain)
            I = coeff(chebfun(1,domain));
        end
        
        function I = zeros(A,domain)
            I = coeff(chebfun(0,domain));
        end
        
        function F = diag(A, f)
            F = coeff(f);
        end
        
        function C = cumsum(A,domain)
            error('How will we do this??')
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
            C = coeff(c);
        end
        
        function A = uminus(A)
            for k = 1:numel(A.coeffs)
                A.coeffs{k} = -A.coeffs{k};
            end
        end

    end
    
end
        
