classdef colloc2 < linopDiscretization
    properties 
        size = [];  % arbitrary, but fixed in any one instance
        domain = [-1 1];
    end
    
    methods
        function A = colloc2(varargin)
            % Collocation matrix on 2nd kind points.
            
            % COLLOC2(DIM,DOMAIN) returns a dummy object that will propagate the
            % dimension size DIM and function domain DOM throughout the delayed
            % evaluation stack.
            %
            % COLLOC2(A,DIM) realizes the linop A (which knows its domain) at
            % dimension DIM.
            % 
            % COLLOC2([]) returns a dummy object that gives access to static
            % methods. 
            
            if ( nargin > 1 )
                if isa(varargin{1},'linop')
                    L = varargin{1};
                    A.size = varargin{2};
                    A.domain = domain(L);
                    A = L.delayFun( A );
                else
                    A.size = varargin{1};
                    validateattributes(varargin{2},{'numeric'},{'increasing','finite'});
                    A.domain = varargin{2};
                end
            end
        end
        
        % Building blocks.
        
        function D = diff(A,m) 
            d = A.domain;
            if m == 0
                D = linop.eye(d);
            else
                numIntervals = length(d)-1;                
                n = dim(A);
                
                % Find the diagonal blocks.
                blocks = cell(numIntervals);
                for k = 1:numIntervals
                    len = d(k+1) - d(k);
                    blocks{k} = diffmat(n(k),m) * (2/len)^m;
                end
                
                % Assemble.
                D = blkdiag(blocks{:});
            end
        end
        
        function C = cumsum(A,m)
            d = A.domain;
            if m == 0
                C = linop.eye(d);
            else
                numIntervals = length(d)-1;
                n = dim(A);

                % Find the diagonal blocks.
                blocks = cell(numIntervals);
                for k = 1:numIntervals
                    len = d(k+1) - d(k);
                    blocks{k} = cumsummat(n(k)) * (len/2);
                end
                
                % Assemble.
                C = blkdiag(blocks{:});
                
                % Each subinterval also contributes to the integrals in all the
                % subintervals to its right, creating a triangular structure.
                offset = 0;
                for k = 1:numIntervals
                    % Grab the weights for the integral using all of this
                    % subinterval.
                    row = offset + n(k);
                    cols = offset + (1:n(k));
                    last = C(row,cols);   
                    % Copy it to add to the ones that follow.
                    offset = row;
                    C(offset+1:end,cols) = repmat(last,[sum(n)-offset,1]);
                end
                C = C^m;
            end
        end
        
        function I = eye(A)
            n = dim(A);
            I = eye(sum(n));
        end
        
        function Z = zeros(A)
            n = dim(A);
           Z = zeros(sum(n));
        end
        
        function F = diag(A,f)
            x = chebpts(dim(A),f.domain,2);
            F = diag( f(x) );
        end
        
        % Required operators.
        % These are not accessed, because the manipulated items are
        % actually square matrices. The basic blocks above could be changed
        % to return colloc2 objects, in which case these would opearte on
        % some field of the objects.
        
        function C = mtimes(A,B)
            C = A*B;
        end
        
        function C = plus(A,B)
            C = A+B;
        end
        
        function B = uminus(A)
            B = -A;
        end
        
        % Required functionals.
        
        function S = sum(A)
            C = cumsummat(dim(A));
            d = A.domain;
            S = C(end,:) * (d(end)-d(1))/2;
        end
        
        function E = evalAt(A,location,direction)
            n = dim(A);
            % Find the collocation points and create an empty functional. 
            x = chebpts(n,A.domain,2);
            offset = cumsum([0;n(:)]);
            N = offset(end);
            E = zeros(1,N);
           
            % Only one subinterval creates nonzero entries in E.
            intnum = linopDiscretization.whichInterval(location,A.domain,direction);
            active = offset(intnum) + (1:n(intnum));
            E(1,active) = barymat(location,x(active));
 
        end
        
        function F = inner(A,f)
            [x,w] = chebpts(dim(A),domain(f),2);
            F = w.*f(x);
        end
        
    end
    
    methods (Static)
        % Additional methods
        
        function B = resize(A,m,n,domain)
            numint = length(m);
            P = cell(1,numint);
            for k = 1:numint               
                xOut = chebpts(m(k),domain(k:k+1),2);
                xIn = chebpts(n(k),domain(k:k+1),2);
                P{k} = barymat(xOut,xIn);
            end
            B = blkdiag(P{:})*A;
        end
        
        function [isDone,epsLevel] = convergeTest(v)
            
            % The chebfun constructor is happy only when coefficient sizes drop below a
            % level that is tied to machine precision. For solutions of BVPs, this is
            % unrealistic, as the condition number of the problem creates noise in the
            % solution at a higher level. Here we try to detect whether the
            % coefficients have reached a "noise plateau" falling below the given
            % relative threshold. If so, we replace the coefficients on the plateau
            % with zero in order to nudge the constructor to stop.
            
            % Copyright 2011 by The University of Oxford and The Chebfun Developers.
            % See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.
            
            isDone = false;
            epsLevel = eps;
            thresh = 1e-6;  % demand at least this much accuracy
            
            n = length(v);
            if n < 17, return, end
            
            % Convert to Chebyshev coefficients.
            c = chebtech2.vals2coeffs(v);
            %             c = chebpoly( chebfun(v) );
            c = c(end:-1:1);
            
            % Magnitude and rescale.
            ac = abs(c)/min(max(abs(v)),1);
            
            % Smooth using a windowed max to dampen symmetry oscillations.
            maxac = ac;
            for k = 1:8
                maxac = max(maxac(1:end-1),ac(k+1:end));
            end
            
            % If too little accuracy has been achieved, do nothing.
            t = find(maxac<thresh,1);
            if isempty(t) || n-t < 16
                return
            end
            
            % Find where improvement in the windowed max seems to stop, by looking at
            % the derivative of a smoother form of the curve.
            dmax = diff( conv( [1 1 1 1]/4, log(maxac(t:end)) ) );
            mindmax = dmax;
            for k = 1:2
                mindmax = min(mindmax(1:end-1),dmax(k+1:end));
            end
            
            %cut = t+k+8 + find(mindmax < 0.02*min(mindmax), 1, 'last');
            cut = find(mindmax > 0.01*min(mindmax), 3);
            if isempty(cut)
                cut = 1;
            else
                cut = cut(end) + t + k + 3;
            end
            
            % Are we satisfied?
            if cut < length(v)
                isDone = true;
                epsLevel = max( abs(c(cut+1)) );
            end
            
        end
        
    end
    
    methods (Access=private)
        function D = diffmat(N,k)
            % DIFFMAT  Chebyshev differentiation matrix
            % D = DIFFMAT(N) is the matrix that maps function values at N Chebyshev
            % points to values of the derivative of the interpolating polynomial at
            % those points.
            %
            % D = DIFFMAT(N,K) is the same, but for the Kth derivative.
            %
            % The matrices are computed using the 'hybrid' formula of Schneider &
            % Werner [1] and Welfert [2] proposed by Tee [3].
            
            % Copyright 2011 by The University of Oxford and The Chebfun Developers.
            % See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.
            
            % References:
            %  [1] Schneider, C. and Werner, W., "Some new aspects of rational
            %   interpolation", Math. Comp. (47) 285--299, 1986.
            %  [2] Welfert, B. D., "Generation of pseudospectral matrices I", SINUM,
            %   (34) 1640--1657.
            %  [3] Tee, T. W., "An adaptive rational spectral method for differential
            %   equations with rapidly varying solutions", Oxford DPhil Thesis, 2006.
            
            persistent cache    % stores computed values for fast return
            if isempty(cache), cache = {}; end    % first call
            
            if nargin < 2, k = 1; end
            
            if N == 0, D = []; return, end
            if N == 1, D = 0; return, end
            
            if length(cache) >= N && length(cache{N}) >= k && ~isempty(cache{N}{k})
                D = cache{N}{k};
                return
            else
                cache{N}{k} = [];
            end
            
            % construct Chebyshev grid and weights
            x = chebpts(N);
            w = [.5 ; ones(N-1,1)]; w(2:2:end) = -1; w(N) = .5*w(N);
            
            ii = (1:N+1:N^2)';              % indices of diagonal
            Dx = bsxfun(@minus,x,x');       % all pairwise differences
            Dx(ii) = Dx(ii) + 1;            % add identity
            Dxi = 1./Dx;                    % reciprocal
            Dw = bsxfun(@rdivide,w.',w);    % pairwise divisions
            Dw(ii) = Dw(ii) - 1;            % subtract identity
            
            % k = 1
            if ~isempty(cache{N}{1})
                D = cache{N}{1};                            % recover from cache
            else
                D = Dw .* Dxi;
                D(ii) = 0; D(ii) = - sum(D,2);              % negative sum trick
                cache{N}{1} = D;                            % store in cache
            end
            
            if k == 1, return, end
            
            % k = 2
            if k > 1 && ~isempty(cache{N}{2})
                D = cache{N}{2};                            % recover from cache
            elseif k > 1
                D = 2*D .* (repmat(D(ii),1,N) - Dxi);
                D(ii) = 0; D(ii) = - sum(D,2);              % negative sum trick
                cache{N}{2} = D;                            % store in cache
            end
            
            % higher orders
            for n = 3:k
                if ~isempty(cache{N}{n})
                    D = cache{N}{n};
                else
                    D = n*Dxi .* (Dw.*repmat(D(ii),1,N) - D);
                    D(ii) = 0; D(ii) = - sum(D,2);          % negative sum trick
                    cache{N}{n} = D;                        % store in cache
                end
            end
            
            if N < 2^11+2
                siz = whos('cache');
                if siz.bytes > cheboppref('maxstorage')
                    cache = {};
                end
                cache{N}{k} = D;
            end
            
        end
        
        function Q = cumsummat(N)
            % CUMSUMMAT  Chebyshev integration matrix.
            % Q = CUMSUMMAT(N) is the matrix that maps function values at N Chebyshev
            % points to values of the integral of the interpolating polynomial at
            % those points, with the convention that the first value is zero.
            
            % Copyright 2011 by The University of Oxford and The Chebfun Developers.
            % See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.
            
            N = N-1;
            
            persistent cache    % stores computed values for fast return
            if isempty(cache), cache = {}; end    % first call
            
            if length(cache) >= N && ~isempty(cache{N})
                Q = cache{N};
                return
            else
                cache{N} = [];
            end
            
            % Matrix mapping coeffs -> values.
            T = cp2cdm(N);
            
            % Matrix mapping values -> coeffs.
            Tinv = cd2cpm(N);
            
            % Matrix mapping coeffs -> integral coeffs. Note that the highest order
            % term is truncated.
            k = 1:N;
            k2 = 2*(k-1);  k2(1) = 1;  % avoid divide by zero
            B = diag(1./(2*k),-1) - diag(1./k2,1);
            v = ones(N,1); v(2:2:end) = -1;
            B(1,:) = sum( diag(v)*B(2:N+1,:), 1 );
            B(:,1) = 2*B(:,1);
            
            Q = T*B*Tinv;
            Q(1,:) = 0;  % make exact
            cache{N} = Q;
            
        end
        
        function T = cp2cdm(N)
            % Values of Cheb. polys at Cheb nodes, x(n)=-cos(pi*n/N).
            theta = pi*(N:-1:0)'/N;
            T = cos( theta*(0:N) );
        end
        
        function C = cd2cpm(N)
            % Three steps: Double the data around the circle, apply the DFT matrix,
            % and then take half the result with 0.5 factor at the ends.
            theta = (pi/N)*(0:2*N-1)';
            F = exp( -1i*theta*(0:2*N-1) );  % DFT matrix
            rows = 1:N+1;  % output upper half only
            % Impose symmetries on data and coeffs.
            C = real( [ F(rows,N+1) F(rows,N:-1:2)+F(rows,N+2:2*N) F(rows,1) ] );
            C = C/N;  C([1 N+1],:) = 0.5*C([1 N+1],:);
        end
                
    end
end