function F = restrict(F, newDomain)
%RESTRICT   Restrict a CHEBFUN object to a subinterval.
%   G = RESTRICT(F, [S1, S2]) returns a CHEBFUN G defined on the interval [S1,
%   S2] which agrees with F on that interval. Any interior breakpoints in
%   F.DOMAIN within [S1, S2] are kept in G.DOMAIN.
%
%   G = RESTRICT(F, S), where S is a row vector, will introduce additional
%   interior breakpoints at S(2:end-2).
%
%   In both cases, if S(1) > S(end), S(1) < F.domain(1), or S(end) >
%   F.domain(end), then an error is returned. If S is empty or a scalar, then an
%   empty CHEBFUN G is returned.
%
%   Note that G will not be 'simplified'. If this is required, call G =
%   SIMPLIFY(RESTRICT(F)), or G = F{S}.
%
% See also OVERLAP, SUBSREF, DEFINE, SIMPLIFY.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Tweak the domain of a quasimatrix input:
[F, newDomain] = tweakDomain(F, newDomain);

% Loop over the columns:
for k = 1:numel(F)
    F(k) = columnRestrict(F(k), newDomain);
end

end

function f = columnRestrict(f, newDomain)

% Ignore duplicate entries in newDomain:
newDomain([false, ~diff(newDomain)]) = [];

% Empty case:
if ( isempty(f) )
    return
elseif ( isempty(newDomain) || numel(newDomain) == 1 )
    f = chebfun();
    return
end

% Cast a periodic chebfun to a regualr chebfun:
if ( isPeriodicTech(f) )
    f = chebfun(f);
end


% Grab domain from f:
oldDomain = f.domain;

% Trivial case:
if ( (numel(newDomain) == 2) && domainCheck(f, newDomain([1, end])) )    
    % The domains are the same!
    return
end

% Check the domain is valid:
if ( (newDomain(1) < oldDomain(1)) || (newDomain(end) > oldDomain(end)) || ...
        any(diff(newDomain) < 0) )
    % newDom is not a valid subinterval of oldDom!
    error('CHEBFUN:CHEBFUN:restrict:subdom', 'Not a valid subdomain.');
end

% Obtain FUN cell and pointValues from f:
funs = f.funs;
pointValues = f.pointValues;

% Discard intervals to the left:
discardIntsLeft = oldDomain(2:end) < newDomain(1);
oldDomain([discardIntsLeft, false]) = [];
funs(discardIntsLeft) = [];
pointValues([discardIntsLeft, false],:,:) = [];

% Discard intervals to the right:
discardIntsRright = oldDomain(1:end-1) > newDomain(end);
oldDomain([false, discardIntsRright]) = [];
funs(discardIntsRright) = [];
pointValues([false, discardIntsRright],:,:) = [];

% Take the union of the new and old domains:
if ( ~isempty( oldDomain(2:end-1) ) ) % Required due to Matlab union() behavior.
    % NB: the old endpoints will either be in newDomain or are not required.
    newDomain = union(oldDomain(2:end-1), newDomain);
end
numFuns = numel(funs);

% Initialise storage for new FUN objects and pointValues:
newFuns = cell(1, numel(newDomain)-1);

% Loop through each fun and restrict as required:
l = 0;
for k = 1:numFuns
    % Find the breaks which correspond to the kth fun:
    subsIdx = (newDomain >= oldDomain(k)) & (newDomain <= oldDomain(k+1));
    if ( sum(subsIdx) == 2 )
        % This interval is already a FUN: (i.e., no new breaks to introduce)
        newFuns{l+1} = restrict(funs{k}, newDomain(subsIdx));
        l = l + 1;
    else
        numSubs = sum(subsIdx)-1;
        if ( numSubs > 0 )
            % Restrict the FUN at the corresponding break points:
            newFuns(l+(1:numSubs)) = restrict(funs{k}, newDomain(subsIdx));
            l = l + numSubs;
        end
    end
end

% Update the pointValues:
newPointValues = chebfun.getValuesAtBreakpoints(newFuns);
% Restore existing pointValues:
[mask, locB] = ismember(oldDomain, newDomain);
locB(~logical(locB)) = [];
newPointValues(locB,:) = pointValues(mask,:);

% Attach data to CHEBFUN to return as output:
f.domain = newDomain;
f.funs = newFuns;
f.pointValues = newPointValues;

end
