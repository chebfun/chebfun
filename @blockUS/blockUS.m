classdef blockUS
    properties
        size = [];  % arbitrary, but fixed in any one instance
        domain = [-1 1];
    end
    
    methods
        function A = blockUS(varargin)
            % ultraspherical method for ChebT coefficients.
            if ( nargin > 1 )
                if isa(varargin{1},'linBlock')
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
            %DIFFMAT(N,K,N), computes the kth order US derivative matrix
            
            n = sum(dim(A));
            if K > 0
                D = spdiags((0:n-1)',1,n,n);
                for s = 1:K-1
                    D = spdiags(2*s*ones(n,1),1,n,n)*D;
                end
            else
                D = speye(n,n);
            end
            %             if nargin == 3
            %                 dom = varargin{1};
            %                 diffMat = ((2./diff(dom))^k)*diffMat;
            %             end
        end
        
        function S = convert( A, K1, K2 )
            %CONVERTMAT(A,K1,K2), convert C^(K1) to C^(K2)
            
            n = sum(dim(A));
            S = speye(n);
            for s = K1:K2
                S = spconvert(n,s) * S;
            end
        end
        function M = mult(A, f, lambda)
            
            n = sum(dim(A));
            
            % get Chebyshev T coefficients
            a = flipud(get(f,'coeffs'));
            
            if ( numel(a) == 1 ) 
                M = a*speye(n);
                return;
            end
            
            % prolong or truncate coefficients
            if ( numel(a) < n )
                a = [a;zeros(n - numel(a),1)];
            else
                a = a(1:n);  % truncate.
            end
            
            if ( lambda == 0 )
                a = a/2;  % just to make formula easier.
                M = sptoeplitz([2*a(1);a(2:end)],[2*a(1);a(2:end)]);
                H = sphankel(a(2:end));
                sub1 = 2:length(a); sub2 = 1:length(a)-1;
                M(sub1,sub2) = M(sub1,sub2)+ H;
            elseif ( lambda == 1 )
                M = sptoeplitz([2*a(1);a(2:end)],[2*a(1);a(2:end)])/2;
                sub = 1:length(a)-2;
                M(sub,sub) = M(sub,sub) - sphankel(a(3:end)/2);
            else
                % TODO: Add higher-order multiplication matrices.
                error
            end
        end
            function d = dim(A)
                d = A.size;
            end
    end
        
        methods (Static)
            function B = resize(A, m, n, dom)
                % chop off some rows and columns
                B = A(1:m,:);
            end
            function [isDone,epsLevel] = testConvergence(v)
                % TODO: (for breakpoints and systems)
                v = full(v);
                f = chebtech2({[], flipud(v)});
                [isDone, epsLevel] = strictCheck(f);
                %             isDone = 1;
            end
            function f_coeffs = discretizeFunction(f,dim,dom)
                if ( nargin < 3 )
                    dom = f.domain;
                end
                f_coeffs = flipud(get(f,'coeffs'));
                % prolong/truncate.
                if length(f_coeffs) > dim
                    f_coeffs = f_coeffs(1:dim);
                else
                    f_coeffs = [f_coeffs;zeros(dim-length(f_coeffs),1)];
                end
                
                % TODO: Mapping to correct US basis.
            end
            
            L = quasi2USdiffmat(L, dim)
                                
            function L = discretize(A, dim, dom, varargin)
                if ( isa(A, 'functionalBlock') )
                    L = A.delayFun( blockColloc2(dim, dom) );
                    L = flipud(chebtech2.coeffs2vals(L.')).';
                else
                    L = A.delayFun( blockCoeff([], dom) );
                    L = blockUS.quasi2USdiffmat(L, dim);
                end
            end
            
            function f = makeChebfun(u, dom)
                funs = cell(numel(u),1);
                for k = 1:numel(u)
                    ct = chebtech2({[], flipud(u{k})});
                    funs{k} = bndfun(ct, dom);
                end
                f = chebfun(funs);
                u = chebmatrix({f});
            end
        end
    end