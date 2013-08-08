function r = roots(f, varargin)
%ROOTS   Roots of a CHEBFUN.
%   ROOTS(F) returns the roots of F in its domain of definition. By default,
%   roots are returned at jumps in F which pass through zero, and if F is
%   identically zero on a part of its domain thena single roots is returned at
%   the midpoint. Each of these behaviours can be modified using the optional
%   inputs described below:
%
%   ROOTS(F, 'nojump') prevents ROOTS() from returning points where F changes
%   sign due to a jump discontinuity, such as roots(chebfun(@sign, 'splitting',
%   'on')).
%
%   ROOTS(F, 'nozerofun') prevents ROOTS() from returning a zero at the midpoint
%   of the domain of F when F if identically zero, such as ROOTS(chebfun(0)).
% 
%   ROOTS(F, 'norecursion') deactivates the recursion procedure used to compute
%   roots (see the Guide 3: Rootfinding and minima and maxima for more
%   information on this recursion procedure).
%
%   ROOTS(F, 'all') returns the roots of all the FUN objects representing the
%   smooth pieces of F. Note that by default this disables recursion, and so is
%   equivalent to ROOTS(F, 'all', 'norecursion').
%
%   ROOTS(F, 'complex') returns the real and complex roots of the FUN objects
%   representing the smooth pieces of F that are determined to be non-spurious.
%   (See CHEBELLIPSEPLOT). This capability may remove some spurious roots that
%   can appear if using ROOTS(F, 'all'). ROOTS(F, 'complex') is equivalent to
%   ROOTS(F, 'complex', 'recursion').
%
%   ROOTS(F, 'all', 'recursion') and ROOTS(F, 'complex','norecursion') can be
%   used to activate and deactivate the recursion procedure respectively, to
%   compute the roots as explained in the 'all' and 'complex' modes.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with the trivial empty case:
if ( isempty(f) )
    return
end

rootsPref = parseInputs(f, varargin{:});

% Set horizontal and vertical scales:
el = epslevel(f);
hs = hscale(f);
vs = vscale(f);
htol = el*hs;
vtol = el*vs;
dom = f.domain;

% Initialise vector to store roots:
r = [];

% Zero impulses are roots.
if ( any(abs(f.impulses(1,:,1))) < vs*htol )
    % Left impulses is zero: (or sufficiently close)
    r = dom(1);
end

funs = f.funs;
nFuns = numel(funs);
for k = 1:nFuns

    % Get the roots of the current fun:
    rk = roots(funs{k}, rootsPref);

    % Trim out roots that are repeated on either side of the breakpoint:
    if ( ~isempty(r) )
        rk(abs(r(end) - rk) < htol) = [];
    end

    % Append new roots to r:
    r = [ r ; rk ]; %#ok<AGROW>

    % Are any roots at the next breakpoint?
    if ( ~isempty(r) && abs(r(end) - dom(k+1)) < htol )
        % The next breakpoint is already in the roots list.

    elseif ( abs(f.impulses(k+1,:,1)) < vtol )
        % The next impulse is zero (or small). Add a root:
        r = [ r ; dom(k+1) ];

    elseif ( rootsPref.jumpRoot && k < nFuns && ...
                get(funs{k}, 'rval') * get(funs{k+1}, 'lval') <= 0 )
        % The solution jumps in sign across the break. Add a root:
        r = [ r ; dom(k+1) ];

    end

end

end

function rootsPref = parseInputs(f, varargin)
% Parse the preferences. See documentation in ROOTS().

% Defaults:
rootsPref = struct('all', 0, 'recurse', 1, 'prune', 0,  'zeroFun', 1, ...
    'jumpRoot', 1);

if ( ~isreal(f) )
    % 'jumpRoots' only makes sense for real-valued functions, so disable it:
    rootsPref.jumpRoot = false;
end

% Parse the inputs:
recurseHasBeenSet = 0;
for k = 1:numel(varargin)
    argin = lower(varargin{k});
    switch ( argin )
        case 'all'
            rootsPref.all = 1;
            rootsPref.prune = 0;
            if ( ~recurseHasBeenSet )
                rootsPref.recurse = 0;
            end
        case 'complex'
            rootsPref.prune = 1;
            rootsPref.all = 1;
        case 'zerofun'
            rootsPref.zeroFun = 1;
        case 'nozerofun'
            rootsPref.zeroFun = 0;
        case 'jump'
            rootsPref.jumpRoot = 1;
        case 'nojump'
            rootsPref.jumpRoot = 0;
        case 'recursion'
            rootsPref.recurse = 1;
            recurseHasBeenSet = 1;
        case 'norecursion'            
            rootsPref.recurse = 0;
        otherwise
            error('CHEBFUN:roots:UnknownOption', 'Unknown option in ROOTS.')
    end
end

end







