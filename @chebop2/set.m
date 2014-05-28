function N = set(N,varargin)
%SET   Set chebop properties.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isempty(N) )   % empty check  
    N = []; return; 
end

rect = N.domain; 

propertyArgIn = varargin;




while length(propertyArgIn) >= 2,
    prop = propertyArgIn{1};
    val = propertyArgIn{2};
    propertyArgIn = propertyArgIn(3:end);
    switch prop
        case {'dom','domain'}
            if isa(val,'double')
                N.domain = val;
            end
        case 'bc'
            N.lbc = chebop2.createbc(val, rect(3:4));
            N.rbc = N.lbc;
            N.ubc = chebop2.createbc(val, rect(1:2));
            N.dbc = N.ubc;
        case 'lbc'
             N.lbc = chebop2.createbc(val, rect(3:4));
             N.lbcshow = val;
        case 'rbc'
             N.rbc = chebop2.createbc(val, rect(3:4));
             N.rbcshow = val;
        case 'ubc'
             N.ubc = chebop2.createbc(val, rect(1:2));
             N.ubcshow = val;
        case 'dbc'
             N.dbc = chebop2.createbc(val, rect(1:2));
             N.dbcshow = val;
        case 'op'
            if isa(val,'function_handle')
                % Do nothing
            else
                error('CHEBOP:set:opType','Operator must by a function handle.')
            end
            N.op = val;
            if ~iscell(val)
                N.opshow = {char(val)};
            else
                N.opshow = cellfun(@char,val,'uniform',false);
            end
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
            error('CHEBOP:set:unknownprop','Unknown chebop property')
    end
end
end