function exponents = findSingExponents(op, isSingEnd, singType, pref)
%FINDEXPONENTS Endpoint singularity detection by sampling values.
%  Private method of SINGFUN.

% TODO Documentation

% Copyright 2013 by The University of Oxford and The Chebfun De

tol = pref.singfun.eps;

%[TODO: At the moment we are assuming the same kind
% of singularity at both end points.
if ( strcmpi( singType{1}, 'pole') )
    exponents = -singfun.findPoleOrder(op, isSingEnd);
else
    if ( ~strcmpi(singType{1}, 'branch') )
        warning('CHEBFUN:singfun:findSingExponents:unknownPref',...
            'Blowup preference "%s" unknown; using default',type)
    end
    exponents = -singfun.findBranchOrder(op, isSingEnd);
end
