function f = compose(f, op, g, pref)
%COMPOSE   Composition of CHEBTECH objects.
%   COMPOSE(F, OP) returns a CHEBTECH representing OP(F) where F is also a
%   CHEBTECH object, and OP is a function handle.
%
%   COMPOSE(F, OP, G) returns a CHEBTECH representing OP(F, G) where F and G
%   are CHEBTECH objects, and OP is a function handle.
%
%   COMPOSE(F, G) returns a CHEBTECH representing G(F), where both F and G are
%   also CHEBTECH objects. If the range of F is not in [-1, 1] then an error is
%   thrown.
%
%   COMPOSE(F, OP, G, PREF) or COMPOSE(F, OP, [], PREF) uses the options passed
%   by the preferences structure PREF to build the returned CHEBTECH.  In
%   particular, one can set PREF.REFINEMENTFUNCTION to be a function which takes
%   advantage of F and possibly OP or G being CHEBTECH objects.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs:
if ( nargin < 4 )
    pref = f.techPref();
else
    pref = f.techPref(pref);
end

if ( (nargin < 3) || isempty(g) )
    nfuns = 1;
else
    nfuns = 2;
end

% Set some preferences:
vscale = f.vscale;
pref.minPoints = max(pref.minPoints, length(f));
pref.sampleTest = false;
% pref.eps = max(pref.eps, f.epslevel);
pref.eps = updateEpslevel(pref.eps, f.epslevel, [], length(f));

if ( nfuns == 2 )
    if ( size(f, 2) ~= size(g, 2) )
        error('CHEBFUN:CHEBTECH:compose:dim', ...
              'Matrix dimensions must agree.')
    end

    % Grab some data from G:
    vscale = max(vscale, g.vscale);
    pref.minPoints = max(pref.minPoints, length(g));
%     pref.eps = max(pref.eps, g.epslevel);
    pref.eps = updateEpslevel(pref.eps, g.epslevel, [], length(g));
        
elseif ( isa(op, 'chebtech') )
    % If OP is a CHEBTECH, we grab some of its data:
    if ( (size(op, 2) > 1) && (size(f, 2) > 1) )
        error('CHEBFUN:CHEBTECH:compose:arrval', ...
              'Cannot compose two array-valued CHEBTECH objects.')
    end
    
    values = f.coeffs2vals(f.coeffs);
    if ( norm(values(:), inf) > 1 )
        error('CHEBFUN:CHEBTECH:compose:range', ...
              'The range of f is not contained in the domain of g.')
    end

    vscale = max(vscale, op.vscale);
    pref.minPoints = max(pref.minPoints, length(op));
%     pref.eps = max(pref.eps, op.epslevel);
    pref.eps = updateEpslevel(pref.eps, op.epslevel, [], length(op));
    
end

% Use a naive evaluation procedure if a custom refinement has not been passed.
if ( ischar(pref.refinementFunction) )
    if ( nfuns == 1 )
        op = @(x) feval(op, feval(f, x));
    else
        op = @(x) feval(op, feval(f, x), feval(g, x));
    end
end

% % Lower epslevel so that we do not get too much degredation on repeated calls
% % to COMPOSE():
% pref.eps = max(eps, pref.eps/100);

% Make CHEBTECH object:
f = f.make(op, vscale, f.hscale, pref);

% % Throw a warning: (Removed by NH Apr 2014. See #282)
% if ( ~f.ishappy )
%     warning('CHEBFUN:CHEBTECH:compose:convfail', ...
%         [ 'Composition with ', func2str(op), ...
%           ' failed to converge with ', int2str(length(f)), ' points.' ]);
% end

end

function el = updateEpslevel(el1, el2, n1, n2)
%UPDATEEPSLEVEL   Update .epslevel property.
%   UPDATEEPSLEVEL(EL1, EL2) returns max(EL1, EL2).
%
%   UPDATEEPSLEVEL(EL1, EL2, N1, N2) returns MAX(EL1/T1^(2/3), EL2/T2^(2/3)),
%   where T1 and T2 are the estimated 'tailLengths' of the corresponding
%   CHEBETCHS F1 and F2 as defined by CHEBTECH.CLASSICCHECK(). This is needed as
%   subsequent calls to CLASSICCHECK() such as in CMOPOSE() result in
%   degredation of accuracy as the tolerance is relaxed as the discretization
%   length is increase.

% TODO: Is this the right thing to do?

if ( nargin < 3 || isempty(n1) )
    n1 = 1;
end
if ( nargin < 4 || isempty(n2) )
    n2 = 1;
end
testLength = min(n1, max(5, round((n1-1)/8)));
testLength2 = min(n2, max(5, round((n2-1)/8)));
% el = max(el1, el2);
el = max(el1/testLength^(2/3), el2/testLength2^(2/3));

end
