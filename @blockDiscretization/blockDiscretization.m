classdef (Abstract) blockDiscretization < operatorBlockRealization & functionalBlockRealization
    
    properties
        domain = [-1 1]
        dimension = []
        source = []
    end
    
    properties (Dependent)
        numIntervals
    end
    
    methods
        function t = isempty(disc)
            t = isempty(disc.source);
        end
        
        function n = get.numIntervals(disc)
            n = length(disc.domain) - 1;
        end        
    end
    
    methods
        
        function ok = validateParameters(disc)
            ok = true;
            if isempty(disc.dimension)
                error('Dimension of discretization not specified.')
            end
            if ( length(disc.dimension) ~= disc.numIntervals )
                error('Discretization lengths improperly specified for subintervals.')
            end
        end
        
        function intnum = whichInterval(disc,location,direction)
            % Which subinterval of the domain is active?
            intnum = find( location >= disc.domain, 1, 'last' );
            
            % linop already screened to make sure the location is in the
            % interval, so we can check for being at the right endpoint.
            if ( intnum == length(disc.domain) )
                if (direction > 0)
                    error('Evaluation direction is undefined at the location.')
                end
                direction = -1;  % this forces the adjustment below
            end
            
            % Need to adjust if at a breakpoint coming from the left, or if at
            % the right endpoint.
            len = disc.domain(end) - disc.domain(1);
            if ( direction < 0 && abs( location - disc.domain(intnum) ) < 10*eps*len )
                intnum = intnum - 1;
            end
        end
        
        function disc = mergeDomains(disc,varargin)
            if isnumeric(varargin{1})
                d = {varargin{1}};
            else
                % Discard numerical blocks:
                isn = cellfun(@isnumeric,varargin);
                varargin(isn) = [];
                % Find the domains of each block (output is cell):               
                d = cellfun(@(A) A.domain,varargin,'uniform',false);
            end
            
            % Collect the endpoints and take the outer hull.
            leftEnds = cellfun(@(x) x(1),d);
            left = min(leftEnds(:));
            rightEnds = cellfun(@(x) x(end),d);
            right = max(rightEnds(:));
            
            % We want to soften 'equality' relative to the domain length.
            tol = 100*eps*(right-left);
            
            % Check to see if the domain endpoints are genuinely different.
            if ( max( abs(leftEnds(:)-left) ) > tol ) || ( max( abs(rightEnds(:)-right) ) > tol )
                error('Domain endpoints are not compatible.')
            end
            
            % Extract all the interior breakpoints.
            d = cellfun(@(x) x(2:end-1),d,'uniform',false);
            
            % Find the unique ones (sorted).
            breakpoints = cat(2,d{:});
            breakpoints = unique(breakpoints);
            
            if ~isempty(breakpoints)
                % Remove all too close to the left endpoint.
                isClose = ( breakpoints - left < tol );
                breakpoints(isClose) = [];
                
                % Remove all too close to the right endpoint.
                isClose = ( right - breakpoints < tol );
                breakpoints(isClose) = [];
                
                % Remove interior points too close to one another.
                isClose =  find( diff(breakpoints) < tol  );
                breakpoints(isClose) = [];
            end
            
            % Put it all together.
            disc.domain = [left breakpoints right];
            
        end
    end
    
    methods (Abstract)
        %[x,w] = points(disc)   % appropriate?
        values = toValues(disc,f)
        f = toFunction(disc,values)
        A = discretize(disc)
    end
    
end
