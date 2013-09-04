classdef colloc2 < linopDiscretization
    properties
        size = [];  % arbitrary, but fixed in any one instance
    end
    
    methods
        function A = colloc2(varargin)
            % Collocation matrix on 2nd kind points.
            
            % COLLOC2(DIM) returns a dummy object that will propagate the
            % dimension size DIM throughout the delayed evaluation stack.
            
            % COLLOC2(A,DIM) realizes the linop A at dimension DIM.
            
            if nargin==1
                A.size = varargin{1};
            
            elseif nargin==2 && isa(varargin{1},'linop')
                L = varargin{1};
                A.size = varargin{2};
                A = L.delayFun( A );
                
            end
        end
        
        % Building blocks.
        
        function D = diff(A,domain)
            D = diffmat(A.size);
            D = D*2/diff(domain);
        end
        
        function C = cumsum(A,domain)
            C = cumsummat(A.size);
            C = C*diff(domain)/2;
        end
        
        function I = eye(A,domain)
            I = eye(A.size);
        end
        
        function Z = zeros(A,domain)
            Z = zeros(A.size);
        end

        function F = diag(A,f)
            x = chebpts(A.size,f.domain,2);
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

        function S = sum(A,domain)
            C = cumsummat(A.size);
            S = C(end,:)*diff(domain)/2;
        end
        
        function E = evalAt(A,domain,loc)
            x = chebpts(A.size,domain,2);
            if strcmp(loc,'left')
                loc = domain(1);
                %E = [1 zeros(1,A.size-1)];
            elseif strcmp(loc,'right')
                loc = domain(2);
                %E = [zeros(1,A.size-1) 1];
            end
            E = barymat(loc,x);
        end
        
        function F = inner(A,f)
            [x,w] = chebpts(A.size,domain(f),2);
            F = w.*f(x);
        end
        
    end
    
    methods (Static)
        % Additional methods
        
        function B = resize(A,m,n)
            B = barymat(m,n)*A;
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
%             c = c(end:-1:1);
            
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
                epsLevel = max( abs(c(cut+1:end)) );
            end
                                    
        end

    end
    
    methods (Access=private)
    end
end