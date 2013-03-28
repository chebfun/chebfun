function vals = jumpvals(funs, ends, op)
% Updates the values at breakpoints, i.e., the first row of imps. If there is a
% singular point, op is evaluated in order to obtain a value at the breakpoint.

% Copyright 2011 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

nfuns = numel(funs);
% vals = 0*ends;

% Trivial case:
if ( nfuns == 1 && isempty(funs) )
    vals = []; 
    return
end

% Endpoint values:
vals(1,:) = get(funs{1}, 'lval');
vals(nfuns+1,:) = get(funs{nfuns}, 'rval');

if ( numel(ends) == 2 )
    vals = vals(:).';
    return
end

if ( nargin > 2 ) 
    % Function handle provided

    if ( iscell(op) )
        for k = 2:nfuns
            % Convert from cell to op:
            opk = op{k}; 
            % Function handle is numeric, and of no use!
        
            if ( isa(opk, 'double') )
                vals(k,:) = (get(funs{k-1}, 'rval') + get(funs{k}, 'lval'))/2;
            else
                vals(k,:) = feval(opk, ends(k));
            end
        end
    else
        vals(2:nfuns,:) = feval(op, ends(2:nfuns).');
    end
        
    
else
    % Function handle is not provided
    
    for k = 2:nfuns-1
        vals(k,:) = (get(funs{k-1}, 'rval') + get(funs{k}, 'lval'))/2;
    end
    
end
vals = vals(:).';
