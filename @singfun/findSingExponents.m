function exponents = findSingExponents(op, singType)
%FINDEXPONENTS   Endpoint singularity detection by sampling values.
%  Private method of SINGFUN.

% TODO Documentation

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

tol = singfun.pref.singfun.eps;
exponents = zeros(1,2);
% loop through each end
singEnd = {'left', 'right'};
for k = 1:2
    if ( strcmpi( singType{k}, 'pole') )
        exponents(k) = singfun.findPoleOrder(op, singEnd{k});
    else
        if strcmpi( singType{k}, 'sing'  )
            exponents(k) = singfun.findSingOrder(op, singEnd{k});
        else
            if strcmpi( singType{k}, 'branch' )
                exponents(k) = singfun.findBranchOrder(op, singEnd{k});
            else
                if strcmpi( singType{k}, 'none' )
                    exponents(k) = 0;                    
                else
                    error('CHEBFUN:SINGFUN:findSingExponents:unknownPref',...
                        'singType "%s" unknown', singType{k})
                end
            end
        end
    end
end