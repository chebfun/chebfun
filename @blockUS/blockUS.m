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
                    A = L.stack( A );  % not sure if this is needed...
                else
                    A.size = varargin{1};
                    A.domain = varargin{2};
                end
            end
        end
        
        function D = diff(A, m)
            d = A.domain;
            n = dim(A);
            
            if ( m == 0 )
                D = speye(sum(n));
            else
                numIntervals = length(d)-1;
                
                % Find the diagonal blocks.
                blocks = cell(numIntervals);
                for k = 1:numIntervals
                    len = d(k+1) - d(k);
                    blocks{k} = A.diffmat(n(k), m) * (2/len)^m;
                end
                
                % Assemble.
                D = blkdiag(blocks{:});
            end
        end
        
        M = mult( A, f, lambda)
        
        function S = convert( A, K1, K2 )
            %CONVERT(A, K1, K2), convert C^(K1) to C^(K2)
            d = A.domain;
            n = dim(A);
            numIntervals = length(d) - 1;
            
            % Find the diagonal blocks.
            blocks = cell(numIntervals);
            for k = 1:numIntervals
                blocks{k} = A.convertmat(n(k), K1, K2);
            end
            
            % Assemble.
            S = blkdiag(blocks{:});
        end
        
        function d = dim(A)
            d = A.size;
        end
    end
    
    methods (Static)
        
        function S = convertmat( n, K1, K2 )
            %CONVERTMAT(A, K1, K2), convert C^(K1) to C^(K2)
            
            S = speye(n);
            for s = K1:K2
                S = spconvert(n, s) * S;
            end
        end
        
        
        function D = diffmat( n, m )
            %DIFFMAT(N, K, N), computes the kth order US derivative matrix
            if ( m > 0 )
                D = spdiags((0:n-1)', 1, n, n);
                for s = 1:m-1
                    D = spdiags(2*s*ones(n, 1), 1, n, n)*D;
                end
            else
                D = speye(n);
            end
        end
        
        function B = resize(A, m, n, dom, diffOrder)
            % chop off some rows and columns
            v = [];
            nn = cumsum([0 n]);
            for k = 1:numel(dom)-1
                v = [v m(k) + nn(k) + (1:(n(k)-m(k)))];
            end
            A(v.',:) = [];
            
            % Ensure each block is mapped to the max(difforder) ultra space
            B = A;
            numInts = numel(dom) - 1;
            numBlocks = (size(A, 2) - 1)/sum(n);
            ntmp = [0 n];
            mtmp = [0 m];
            for ii = 1:numBlocks
                doii = diffOrder(ii);
                for jj = 1:numInts
                    indn = (ii-1)*sum(n) + (((jj-1)*ntmp(jj)+1):(jj*n(jj)));
                    indm = ((jj-1)*mtmp(jj)+1):(jj*m(jj));
                    B(indm,indn) = blockUS.convertmat(length(indm), diffOrder(ii), max(diffOrder)-1) * B(indm,indn);
                    B(indm,end) = blockUS.convertmat(length(indm), 0, max(diffOrder)-1) * B(indm,end);
                end
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
            f_coeffs = [];
            f = restrict(f, dom);
            for k = 1:numel(dom)-1
                dimk = dim(k);
                tmp = flipud(get(f.funs{k}, 'coeffs'));
                n = length(tmp);
                % prolong/truncate.
                if ( n > dimk )
                    tmp = tmp(1:dimk);
                else
                    tmp = [tmp ; zeros(dimk - n, 1)];
                end
                f_coeffs = [f_coeffs ; tmp];
            end
        end
        
        L = quasi2USdiffmat(L, dim)
        
        function L = discretize(A, dim, dom, varargin)
            if ( isa(A, 'functionalBlock') )
                % TODO: Just call discretize method?
                L = A.stack( blockColloc2(dim, dom) );
                cumsumdim = [0, cumsum(dim)];
                for k = 1:numel(dom)-1
                    Lk = L(cumsumdim(k) + (1:dim(k)));
                    tmp{k} = flipud(chebtech2.coeffs2vals(Lk.')).';
                end
                L = cell2mat(tmp);
            elseif ( isa(A, 'operatorBlock') )
                A.domain = dom;
                L = blockUS.quasi2USdiffmat(A, dim);
            else
                % TODO: Anything here?
            end
        end
        
        function f = makeChebfun(u, dom)
            funs = cell(numel(u), 1);
            for k = 1:numel(u)
                ct = chebtech2({[], flipud(u{k})});
                funs{k} = bndfun(ct, dom(k:k+1));
            end
            f = chebfun(funs);
            u = chebmatrix({f});
        end
    end
end
