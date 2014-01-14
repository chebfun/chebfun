function g = power(f, b)
%.^   CHEBFUN power.
%   F.^G returns a CHEBFUN F to the scalar power G, a scalar F to the CHEBFUN
%   power G, or a CHEBFUN F to the CHEBFUN power G. F and or G may be complex.
%
%   H = POWER(F, G) is called for the syntax 'F .^ G'.
%
% See also SQRT, COMPOSE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty case.
if ( isempty(f) || isempty(b) )
    g = chebfun();
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHEBFUN .^ CHEBFUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
if ( isa(f, 'chebfun') && isa(b, 'chebfun') ) 
    
    % Call COMPOSE(): (Note, COMPOSE() checks that the domains match)
    g = compose(f, @power, b);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHEBFUN .^ constant %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ( isa(f, 'chebfun') )                     % CHEBFUN .^ constant
    
    if ( b == 0 )                       % Trivial case
        % Constant CHEBFUN:
        vals = ones(1, min(size(f)));
        g = chebfun(vals, f.domain([1,end]));
        
    elseif ( b == 1 )                   % Identity
        g = f;
        
    elseif ( b == 2 )                   % Square
        % Call TIMES():
        g = f.*f;
        
    elseif ( ( b > 0 ) && ( round(b) == b ) )   % Positive integer
        
        % If SINGFUN is involved:
        if ( issing(f) )
            % If singfun is involved, treat each piece individually:
            g = f;
            numFuns = numel(f.funs);
            for k = 1:numFuns
                g.funs{k} = power(f.funs{k}, b);
            end
            g.impulses = g.impulses.^b;
        else
            % Result will be smooth. Call COMPOSE():
            g = compose(f, @(x) power(x, b));
        end
        
    elseif ( b == .5 )                  % Sqrt
        % Call SQRT():
        g = sqrt(f);
        
    else                                % General case
        
        % Add breaks at the appropriate roots of f:
        if ( isreal(f) )
            f = addBreaksAtRoots(f);
        else
            % Add breaks at the roots of the imaginary part of F to account for
            % the discontinuity in POWER along the negative real semi-axis due 
            % to the branch cut.
            f = addBreaksAtRoots(f, 'imag');
        end
        
        % Loop over each piece individually:
        numFuns = numel(f.funs);
        g = f;
        for k = 1:numFuns
            g.funs{k} = power(f.funs{k}, b);
        end
        
        % Update the impulses:
        g.impulses = g.impulses.^b;
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%% constant .^ CHEBFUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    
    % Call COMPOSE():
    g = compose(b, @(x) power(f, x));

end
