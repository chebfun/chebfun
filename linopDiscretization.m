classdef (Abstract) linopDiscretization 
    
    properties
        domain
        linop
    end
    
    properties (Dependent)
        numIntervals
    end
    
    methods (Abstract)
        A = matrix(disc)
        b = rhs(disc,f)
    end
    
    methods

        function n = get.numIntervals(disc)
            n = length(disc.domain) - 1;
        end        

        
        function [isDone,epsLevel] = testConvergence(disc,values)
            
            isDone = false;
            epsLevel = eps;
            thresh = 1e-6;  % demand at least this much accuracy
            
            n = length(values);
            if n < 17, return, end
            
            % Convert to Chebyshev coefficients.
            f = toFunction(disc,values);
            c = chebpoly(f);
            c = c(end:-1:1);
            
            % Magnitude and rescale.
            ac = abs(c)/min(max(abs(values)),1);
            
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
            if cut < length(values)
                isDone = true;
                epsLevel = max( abs(c(cut+1)) );
            end
            
        end
        
        function [x, disc] = mldivide(disc, A, b)
        	x = A\b;
        end
        
        function disc = deriveContinuity(disc)
            % Find automatic smoothness constraints at domain breakpoints.
            L = disc.linop;
            d = L.diffOrder;
            d = max(d,[],1);
            dom = disc.domain;
            
            cont = linopConstraint();
            
            if ( max(d) > 0 ) && ( length(dom) > 2 )
                
                C = domainContinuity(disc,max(d)-1);
                
                % Each function variable gets a zero functional block; each scalar variable
                % gets a scalar zero.
                z = linBlock.zero(dom);
                Z = {};
                for var = 1:length(d)
                    if isnan(d(var)) || d(var) == 0 % scalar
                        Z = [Z,0];
                    else
                        Z = [Z,z];
                    end
                end
                %    Z = chebmatrix(Z);
                
                for var = 1:length(d)
                    % Skip if this is a scalar variable; it plays no role in continuity.
                    if isnan(d(var)) || d(var) == 0
                        continue
                    end
                    B = Z;
                    for m = 0:d(var)-1
                        for k = 2:length(dom)-1
                            B.blocks{var} = C{m+1,k};
                            cont = cont.append(B,0);
                        end
                    end
                end
            end
            
            disc.linop.continuity = cont;
            
        end        
        
    end
    
    methods ( Access=protected )
        function C = domainContinuity(disc,maxorder)
            % Returns expressions of continuity conditions at
            % the breakpoints of the domain of L.
            %   C{m,k} has the (m-1)th-order derivative at breakpoint k
            
            d = disc.domain;
            A = linop.eye(d);
            D = linop.diff(d,1);
            for m = 0:maxorder
                for k = 2:length(d)-1
                    El = linop.feval(d(k),d,'-');
                    Er = linop.feval(d(k),d,'+');
                    C{m+1,k} = (El-Er)*A;
                end
                A = D*A;
            end
        end

    end
    
end
