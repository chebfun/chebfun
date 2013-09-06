function g = power(f, b)
% .^   Chebfun power.
%   F.^G returns a CHEBFUN F to the scalar power G, a scalar F to the CHEBFUN
%   power G, or a CHEBFUN F to the CHEBFUN power G. F and or G may be complex.
%
%   H = POWER(F, G) is called for the syntax 'F .^ G'.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty case.
if ( isempty(f) || isempty(b) )
    g = chebfun();
    return
end

% [TODO]: This might need to be changed to include SINGFUN features.
    
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
        
    elseif ( b == .5 )                  % Sqrt
        % [TODO]: implement this.
        warning('Not yet implemented.')
        g = compose(f, @(x) power(x, b));
%         g = sqrt(f);
        
    elseif ( b == 2 )                   % Square
        % Call TIMES():
        g = f.*f;
        
    elseif ( (b > 0) && (round(b) == b) )   % Positive integer
        % Result will be smooth. Call COMPOSE():
        g = compose(f, @(x) power(x, b));
        
    else                                % General case
        % [TODO]: implement this properly (i.e., using SINGFUN).
        warning('Not yet implemented.')
        g = compose(f, @(x) power(x, b));
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%% constant .^ CHEBFUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    % Call COMPOSE():
    g = compose(b, @(x) power(f, x));

end
