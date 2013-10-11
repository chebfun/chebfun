classdef (Abstract) blockDiscretization < operatorBlockRealization & functionalBlockRealization
  
    properties (Abstract)
        size
        domain
    end
    
    methods
        function n = dim(A)
            % This method converts a scalar size into a vector with one entry
            % per subinterval of the domain. For now, we repeat the value of n,
            % but another choice is to equidistribute it.
            %
            % This method should be the only way used to access A.size, in order
            % to enforce consistent behavoir.
            n = A.size;
            numint = length(A.domain)-1;
            if (length(n)==1)
                n = repmat(n,1,numint);
            elseif ( length(n) ~= numint )
                error('Mismatch between provided discretization sizes and the number of subintervals.')
            end
        end
    end
    
    methods (Abstract,Static)
        B = resize(A,m,n,dom)
        [isDone,epsLevel] = testConvergence(v)
        fx = discretizeFunction(f,dim,dom)
    end
    
    methods (Static,Access=protected)
        
        function intnum = whichInterval(location,domain,direction)
            % Which subinterval of the domain is active?
            intnum = find( location >= domain, 1, 'last' );
            
            % linop already screened to make sure the location is in the
            % interval, so we can check for being at the right endpoint.
            if (intnum == length(domain))
                if (direction > 0)
                    error('Evaluation direction is undefined at the location.')
                end
                direction = -1;  % this forces the adjustment below
            end
            
            % Need to adjust if at a breakpoint coming from the left, or if at
            % the right endpoint.
            len = domain(end)-domain(1);
            if (direction < 0) && ( abs( location - domain(intnum) ) < 10*eps*len )
                intnum = intnum - 1;
            end
        end
            

    end
    
end
        