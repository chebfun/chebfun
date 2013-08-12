function f = defineInterval(f, subInt, g)
%DEFINEINTERVAL   Supply a new definition for a chebfun on a subinterval.
%
%   F = DEFINEINTERVAL(F, S, G) redefines the CHEBFUN F by the CHEBFUN or double
%   G in the interval [S(1), S(end)] in F.domain. If F is array-vaued then G
%   should have the same number of columns, i.e., size(F, 2) = size(G, 2), and
%   if G is a CHEBFUN it must be defined on a domain containing [S(1), S(end)].
%   Any new breakpoints S(2:end-1) are also introduced into F.
%
%   An equivalent syntax is F{S(1), S(2), ..., S(end)} = G.
%
% See also CHEBFUN/SUBSASGN, CHEBFUN/RESTRICT, CHEBFUN/DEFINEPOINT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Not a valid subdomain:
if ( any(diff(subInt) < 0) )
    error('CHEBFUN:defineInterval:invalidDomain', 'Not a valid domain.');
end

% TODO: Ensure f and g have the same number of colums.
numColsF = size(f.funs{1}, 2);

% Define at a single point:
if ( subInt(1) == subInt(end) )
    if ( isa(g, 'chebfun') )
        g = feval(g, subInt(1));
    end
    f = definePoint(f, subInt(1), g);
    return
end

% Convert a scalar or empty input to a chebfun.
if ( isnumeric(g) )
    % Allow scalar expansion:
    if ( size(g, 2) == 1 )
            g = repmat(g, 1, numColsF);
    end
    if ( numel(g) == 0 )
        % An empty CHEBFUN:
        g = chebfun;
    else
        % An constant CHEBFUN:
        g = chebfun(g, subInt);
    end
else
    % Restrict G to SUBINT: (Will throw an error if not a valid subdomain).
    g = restrict(g, subInt);
end

% Trivial case.
if ( isempty(f) )
    f = g;
    return
end

% Make sure dimensions add up:
if ( ~isempty(g) && size(g, 2) ~= numColsF )
    error('CHEBFUN:defineInterval:numCols', ...
        'Dimensions of matrices being concatenated are not consistent.');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%% INSERTION/OVERWRITING %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( ~isempty(g) )                               
    
    if ( subInt(end) < f.domain(1) )             % Extension to the left
        % Append FUNs, domain, and impulses:
        padding = chebfun(0, [g.domain(end), f.domain(1)]);
        f.domain = [g.domain, f.domain];
        f.funs = [g.funs, padding.funs, f.funs];
        f.impulses = [g.impulses ; f.impulses];
        
    elseif ( subInt(1) > f.domain(end) )         % Extension to the right
        % Append FUNs, domain, and impulses:
        padding = chebfun(0, [f.domain(end), g.domain(1)]);
        f.domain = [f.domain, g.domain];
        f.funs = [f.funs, padding.funs, g.funs];
        f.impulses = [f.impulses ; g.impulses];
        
    else                                         % SubInt intersects f.domain
        fLeft = restrict(f, [f.domain(1), subInt(1)]);
        fRight = restrict(f, [subInt(end), f.domain(2)]);
        % Insert FUNs, domain, and impulses:
        f.funs = [fLeft.funs, g.funs, fRight.funs];
        f.domain = [fLeft.domain(1:end-1), g.domain, fRight.domain(2:end)];
        f.impulses = [fLeft.impulses(1:end-1,:,:) ; g.impulses ; ...
            fRight.impulses(2:end,:,:)];
        
    end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELETION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    if ( subInt(2) < f.domain(1) || subInt(1) > f.domain(2) )
        error('CHEBFUN:define:badremoveinterval',...
            'Interval to be removed is outside the domain.')
        
    else
        % Restrict from the left and right:
        fLeft = restrict(f, [f.domain(1), subInt(1)]);
        fRight = restrict(f, [subInt(end), f.domain(end)]);
        if ( isempty(fRight) )
            % Delete part of the domain from the right:
            f = fLeft;
        elseif ( isempty(fLeft) )
            % Delete part of the domain from the left:
            f = fRight;
        else
            % Deletion strictly inside the domain - slide to the left.
            newEnds = fRight.domain - fRight.domain(1) + fLeft.domain(end);
            % Update maps in FUN objects:
            for j = 1:numel(fRight.funs)  
                fRight.funs{j} = changeMap(fRight.funs{j}, newEnds(j:j+1));
            end
            % Insert FUNs, domain, and impulses:
            f.funs = [fLeft.funs, fRight.funs];
            f.domain = [fLeft.domain(1:end-1), newEnds];
            f.impulses = [fLeft.impulses(1:end-1,:,:) ; fRight.impulses];
        end
        
    end
end

end