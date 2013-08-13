function f = restrict(f, s)
%RESTRICT   Restrict a SINGFUN to a subinterval.
%   RESCTRICT(F, S) returns a SINGFUN that is restricted to the subinterval
%   [S(1), S(2)] of [-1, 1]. Note that since SINGFUN objects only live on
%   [-1, 1], a linear change of variables is implicitly applied.
%
%   If length(S) > 2, i.e., S = [S1, S2, S3, ...], then RESCTRICT(F, S) returns
%   an array of CHEBTECH objects, where the entries hold F restricted to each of
%   the subintervals defined by S.
%
%   If F is an array-valued function, say [F1, F2], then the restrict(F, S =
%   [S1, S2, S3]) returns the array-valued CHEBTECH {restrict(F1,S).
%   restrict(F2, S)}.
%
%   Note that restrict does not 'simplify' its output.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    return
end

% Check if s is actually a subinterval:
if ( (s(1) < -1) || (s(end) > 1) || (any(diff(s) <= 0)) )
    error('CHEBFUN:SINGFUN:restrict:badinterval', 'Not a valid interval.')
elseif ( (numel(s) == 2) && all(s == [-1, 1]) )
    % Nothing to do here!
    return
end

% Number of subintervals:
numInts = numel(s) - 1;

% Preallocate a cell
f = cell(1, numInts);

% Set the bound which is deemed 'finite'
fntbnd = realmax;

% function values at the endpoints of the subintervals:
[lval, rval] = op(s);

if ( ( lval < fntbnd ) && ( rval < fntbnd ) )
    
    for j = 1:numInts
        
        % define the new operator which will be evaluated in the subinterval by
        % the singfun ctor.
        op = @(x) feval(f, ((1-x)*s(j)+(1+x)*s(j+1))/2);
        
        % call the singfun constructor
        g = singfun( op, zeros(1,2), {'none', 'none'}, [] );
        
        % put in cell
        f{j} = g;
    end
    
elseif ( ( s(1) == -1 ) && ( rval < fntbnd ) )
    
    % define the new operator
    op = @(x) feval(f, ((1-x)*s(1)+(1+x)*s(2))/2);
    
    % call the singfun constructor
    g = singfun( op, [f.exponents(1) 0], {'sing', 'none'}, [] );
    f{j} = g;
    
    for j = 2:numInts
        
        % define the new operator
        op = @(x) feval(f, ((1-x)*s(j)+(1+x)*s(j+1))/2);
        
        % call the singfun constructor
        g = singfun( op, zeros(1,2), {'none', 'none'}, [] );
        
        % put in cell
        f{j} = g;
    end
    
    
elseif ( ( lval < fntbnd ) && ( s(end) == 1 ) )
    
    if numInts >1
        
        for j = 1:numInts-1
            
            % define the new operator
            op = @(x) feval(f, ((1-x)*s(j)+(1+x)*s(j+1))/2);
            
            % call the singfun constructor
            g = singfun( op, zeros(1,2), {'none', 'none'}, [] );
            
            % put in cell
            f{j} = g;
        end
        
    end
    
    % define the new operator
    op = @(x) feval(f, ((1-x)*s(end-1)+(1+x)*s(end))/2);
    
    % call the singfun constructor
    g = singfun( op, [0 f.exponents(2)], {'none', 'sing'}, [] );
    
    % put in cell
    f{end} = g;
    
else
    
    for j = 1:numInts
        
        % define the new operator
        op = @(x) feval(f, ((1-x)*s(j)+(1+x)*s(j+1))/2);
        
        % call the singfun constructor
        g = singfun( op, [], {'sing', 'sing'}, [] );
        
        % put in cell
        f{j} = g;
    end
end

% If there is only one subinterval, return a singfun, instead of a cell
if ( numInts == 1 )
    f = f{:};
end

end