function s = rdivide(f, g)
%./   Divide DELTAFUNS with DELTAFUNS
%   RDIVIDE(F, G) computes the pointwise division F./G. The operation is only
%   defined if G does not have any delta functions and G has no roots.
%
% See also LDIVIDE, TIMES.
%
% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Case of empty arguments:
if ( isempty(f) || isempty(g) )
    % Return an empty DELTAFUN:
    s = deltafun();
    return
end

% Check if inputs are other than DELTAFUNs, CLASSICFUNs or doubles.
if ( (~isa(f, 'deltafun') && ~isa(f, 'classicfun') && ~isa(f, 'double')) || ...
     (~isa(g, 'deltafun') && ~isa(g, 'classicfun') && ~isa(g, 'double')) )   
    error('CHEBFUN:DELTAFUN:rdivide:rdivide', ...
        'Input can only be a DELTAFUN, a CLASSICFUN or a double')
end

%% Reciprocal: ( CLASSICFUN or DOUBLE ). / DELTAFUN
if ( isa(f, 'double') || isa( f, 'classicfun') && isa(g, 'deltafun') )    
    if ( anyDelta(g) )
        error('CHEBFUN:DELTAFUN:rdivide:rdivide',
            'Division by delta functions is not defined.');
    end
    % A smooth function is returned in this case:
    s = f ./ g.funPart;
    return
end

%% DELTAFUN./DOUBLE
if ( isa(f, 'deltafun') && isa(g,'double') )
    % Copy the DELTAFUN in the output:
    s = f;
    % Divide the funPart with the double g:
    s.funPart = 1/g * f.funPart;
    % Update the delta functions:
    s.deltaMag = 1/g * s.deltaMag;
    if ( ~anyDelta(s) )
        s = s.funPart;
    end
    return
end

%% DELTAFUN./CLASSICFUN
if ( isa(f, 'deltafun') && isa(g, 'classicfun') )
    % Take reciprocal and make a DELTAFUN:
    g = deltafun( 1./g, [], [] );
    % Now, multiply:
    s = f .* g;
    return
end

%% DELTAFUN./DELTAFUN
if ( isa(f, 'deltafun') && isa(g, 'deltafun') )
    if ( anyDelta(g) )
        error('CHEBFUN:DELTAFUN:rdivide:rdivide',
            'Division by delta functions is not defined.');
    end
    % Take reciprocal of the funPart only and make a DELTAFUN:
    g = deltafun( 1 ./ g.funPart, [], []);
    % Now, multiply:
    s = f .* g;
    return
end

end
