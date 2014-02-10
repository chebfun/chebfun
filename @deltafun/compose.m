function f = compose(f, op, g, pref)
%COMPOSE   Compose a DELTAFUN with an affine operator or another DELTAFUN.
%   H = COMPOSE(F, OP) returns a BNDFUN representing OP(F) where F is also a
%   BNDFUN object, and OP is a function handle of a linear function of the form
%   op = @(x) a*x + b for some constants a and b. No check is made on the
%   funciton handle and it is assumed that OP is of the correct form.
%
%   H = COMPOSE(F, OP, PREF) where OP is a function handle, passes the options 
%   passed by the preferences structure PREF to the composition of the FUNPART
%   of F.
%
% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 3 )
    pref = chebpref;
end

% TODO: NH: What about COMPOSE(f, op, g)?

if ( isa(op, 'function_handle') )
    % OP is a function handle:
    f.funPart = compose(f.funPart, op, [], pref);
    
    deltaLoc = f.deltaLoc;
    deltaMag = f.deltaMag;
    dom = f.funPart.domain;
    if ( ~isempty(deltaLoc) )     
        % Determine a and b of the linear map:
        b = op(0);
        a = op(1) - b;
        % Apply the shifting:
        deltaLoc = (-b + deltaLoc)/a;
        
        % Discard delta function that are now outside the domain:
        idx = deltaLoc >= dom(1) & deltaLoc <= dom(2);
        deltaLoc = deltaLoc(idx);
        deltaMag = deltaMag(:, idx);
        
        % Apply the scaling:
        deltaMag = 1/abs(a)*deltaMag;
        
        % Assign back:
        if ( isempty(deltaLoc) || isempty(deltaMag) )
            f.deltaLoc = [];
            f.deltaMag = [];
        else
            f.deltaLoc = deltaLoc;
            f.deltaMag = deltaMag;
        end
    end        
else
    error( 'DELTAFUN:COMPOSE', 'OP must be a function handle.');
end

end
