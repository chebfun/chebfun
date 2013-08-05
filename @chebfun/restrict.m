function f = restrict(f, newDomain)
%RESTRICT    Restrict a chebfun to a subinterval.
%   G = RESTRICT(F, S) returns a chebfun G with domain [S(1), S(end)] which
%   agrees with F on that interval. G will contain breakpoints at the same
%   location as F within [S(1), S(end)]. Additional interior breakpoints can be
%   included by passing S(2), ..., S(end)-1. G = F{S} is an equivalent syntax.
%
%   If [S(1) >= S(end)] then an error is returned.
%
% See also OVERLAP, SUBSREF, DEFINE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty case:
if ( isempty(f) || isempty(newDomain) )
    return
end

% Grab domain from f:
oldDomain = f.domain;

% Trivial case:
if ( (numel(newDomain) == 2) && isequal(newDomain, oldDomain([1 end])) )
    % The domains are the same!
    return
end

% Check the domain is valid:
if ( (newDomain(1) < oldDomain(1)) || (newDomain(end) > oldDomain(end)) || ...
        any(diff(newDomain) < 0) )
    % newDom is not a valid subinterval of oldDom!
    error('CHEBFUN:restrict:subdom', 'Not a valid subdomain.');
end

% Obtain funs and impulses from f:
funs = f.funs;
imps = f.impulses;

% Discard intervals to the left:
discardIntsLeft = oldDomain(2:end) < newDomain(1);
oldDomain([discardIntsLeft false]) = [];
funs(discardIntsLeft) = [];
imps([discardIntsLeft, false],:,:) = [];

% Discard intervals to the right:
discardIntsRright = oldDomain(1:end-1) > newDomain(end);
oldDomain([false discardIntsRright]) = [];
funs(discardIntsRright) = [];
imps([false, discardIntsRright],:,:) = [];

% Take the union of the new and old domains:
if ( ~isempty( oldDomain(2:end-1) ) )
    % NB:  the old endpoints will either be in newDomain or are not required.
    newDomain = union(oldDomain(2:end-1), newDomain);
end
numFuns = numel(funs);

% Initialise storage for new funs and impulses:
newFuns = cell(1, numel(newDomain)-1);
newImps = zeros(numel(newDomain), size(imps, 2), size(imps, 3));

% Loop through each fun and restrict as required:
l = 0;
for k = 1:numFuns
    % Find the breaks which correspond to the kth fun:
    subsIdx = (newDomain >= oldDomain(k)) & (newDomain <= oldDomain(k+1));
    if ( sum(subsIdx) == 2 )
        % This interval is already a FUN: (i.e., no new breaks to introduce)
        newFuns(l+1) = funs(k);
        l = l + 1;
    else
        numSubs = sum(subsIdx)-1;
        if ( numSubs > 0 )
            % Restrict the fun at the corresponding break points:
            newFuns(l+(1:numSubs)) = restrict(funs{k}, newDomain(subsIdx));
            l = l + numSubs;
        end
    end
end

% Update the impulses:
newImps(:,:,1) = chebfun.jumpVals(newFuns);
% Restore existing impulses:
[mask, locB] = ismember(oldDomain, newDomain);
locB(~logical(locB)) = [];
newImps(locB,:,:) = imps(mask,:,:);

% Attach data to chebfun to return as output:
f.domain = newDomain;
f.funs = newFuns;
f.impulses = newImps;

end
