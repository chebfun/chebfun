classdef blockUS
    properties
        size = [];  % arbitrary, but fixed in any one instance
        domain = [-1 1];
    end
    
    methods
        function A = blockUS(varargin)
            % ultraspherical method for ChebT coefficients.
            if ( nargin > 1 )
                if ( isa(varargin{1}, 'linBlock') )
                    L = varargin{1};
                    A.size = varargin{2};
                    A.domain = L.domain;
                    A = L.delayFun( A );  % not sure if this is needed...
                else
                    A.size = varargin{1};
                    A.domain = varargin{2};
                end
            end
        end
        
        function D = diff( A, K )
            %DIFFMAT(N, K, N), computes the kth order US derivative matrix
            
            n = sum(dim(A));
            if ( K > 0 )
                D = spdiags((0:n-1)', 1, n, n);
                for s = 1:K-1
                    D = spdiags(2*s*ones(n, 1), 1, n, n)*D;
                end
            else
                D = speye(n, n);
            end
        end
        
        M = mult( A, f, lambda) 
        
        function S = convert( A, K1, K2 )
            %CONVERTMAT(A, K1, K2), convert C^(K1) to C^(K2)
            
            n = sum(dim(A));
            S = speye(n);
            for s = K1:K2
                S = spconvert(n, s) * S;
            end
        end

        function d = dim(A)
            d = A.size;
        end
    end
    
    methods (Static)
        
        function B = resize(A, m, n, dom, difforder)
            % chop off some rows and columns
            B = A(1:m, :);
            dummy = blockUS(m,[-1,1]);
            for j = 1:difforder-1
                B(:,end) = convert(dummy, j-1, j) * B(:,end);
            end
        end
        
        function [isDone, epsLevel] = testConvergence(v)
            % TODO: (for breakpoints and systems)
            v = full(v);
            f = chebtech2({[], flipud(v)});
            [isDone, epsLevel] = strictCheck(f);
        end
        
        function f_coeffs = discretizeFunction(f, dim, dom)
            if ( nargin < 3 )
                dom = f.domain;
            end
            f_coeffs = flipud(get(f, 'coeffs'));
            n = length(f_coeffs);
            % prolong/truncate.
            if ( n > dim )
                f_coeffs = f_coeffs(1:dim);
            else
                f_coeffs = [f_coeffs ; zeros(dim - n, 1)];
            end
            
            % TODO: Mapping to correct US basis.
        end
        
        L = quasi2USdiffmat(L, dim)
        
        function L = discretize(A, dim, dom, varargin)
            if ( isa(A, 'functionalBlock') )
                L = A.delayFun( blockColloc2(dim, dom) );
                L = flipud(chebtech2.coeffs2vals(L.')).';
            elseif ( isa(A.delayFun, 'blockCoeff') )
                L = blockUS.quasi2USdiffmat(A.delayFun, dim);
            else
                L = A.delayFun( blockCoeff([], dom) );
                L = blockUS.quasi2USdiffmat(L, dim);
            end
        end
        
        function f = makeChebfun(u, dom)
            funs = cell(numel(u), 1);
            for k = 1:numel(u)
                ct = chebtech2({[], flipud(u{k})});
                funs{k} = bndfun(ct, dom);
            end
            f = chebfun(funs);
            u = chebmatrix({f});
        end
    end
end
