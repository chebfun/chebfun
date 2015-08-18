function N = set(N, varargin)
%SET   Set CHEBOP2 properties.
%   N.lbc = F sets left boundary conditions. If is a DOUBLE or a CHEBFUN, 
%   then left Dirichlet conditions with boundary data F are set. 
%   If F = @(x,u) f(x,u), then the conditions f(x,u) = 0 are set along the 
%   left edge.  
% 
%   N.rbc = F, N.ubc = F, N.dbc = F are the same as N.lbc = F but for the
%   right, top, and bottom boundary conditions, respectively.
% 
%   N.domain = F sets the domain of the CHEBOP2. F is expected to be a
%   vector of length 4. 
% 
%   N.xorder = Kx sets the differential order of N in the x-variable to Kx.
%   N.yorder = Ky sets the differential order of N in the y-variable to Ky.
%
%   N.op = fh sets the partial differential operator to fh. 
% 
%   N.U = U, N.S = S, N.V = V, sets the rank approximation to the partial
%   differential operator. 

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Empty check.
if ( isempty(N) ) 
    N = []; 
    return
end

rect = N.domain; 
propertyArgIn = varargin;

while ( length(propertyArgIn) >= 2 )
    prop = propertyArgIn{1};
    val = propertyArgIn{2};
    propertyArgIn = propertyArgIn(3:end);
    
    switch prop
        case {'dom', 'domain'}
            if isa(val, 'double')
                N.domain = val;
            end
        case 'bc'
            N.lbc = chebop2.createBC(val, rect(3:4));
            N.rbc = N.lbc;
            N.ubc = chebop2.createBC(val, rect(1:2));
            N.dbc = N.ubc;
        case 'lbc'
             N.lbc = chebop2.createBC(val, rect(3:4));
        case 'rbc'
             N.rbc = chebop2.createBC(val, rect(3:4));
        case 'ubc'
             N.ubc = chebop2.createBC(val, rect(1:2));
        case 'dbc'
             N.dbc = chebop2.createBC(val, rect(1:2));
        case 'op'
            if isa(val, 'function_handle')
                % Do nothing.
            else
                error('CHEBFUN:CHEBOP2:set:opType', ...
                    'Operator must be a function handle.')
            end
            N.op = val;
        case 'xorder'
            N.xorder = val; 
        case 'yorder'
            N.yorder = val; 
        case 'U'
            N.U = val;
        case 'S'
            N.S = val;
        case 'V'
            N.V = val;
        otherwise
            error('CHEBFUN:CHEBOP2:set:unknownProp', ...
                'Unknown CHEBOP2 property.')
    end
    
end

end
