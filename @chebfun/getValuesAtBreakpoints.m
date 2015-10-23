function vals = getValuesAtBreakpoints(funs, ends, op)
%GETVALUESATBREAKPOINTS   Determine values between neighbouring FUN objects.
%   VALS = GETVALUESATBREAKPOINTS(FUNS, ENDS, OP) returns the values at
%   breakpoints ENDS between FUN objects. If OP can be evaluated at ENDS then
%   VALS = OP(ENDS). Otherwise VALS(j) is the average of the right and left
%   limits of its neighbouring funs for interior breaks and the limits from the
%   left and right for the VALS(1) and VALS(end), respectively.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Determine the number of intervals:
numFuns = numel(funs);

% Trivial empty case:
if ( (numFuns == 1) && isempty(funs{1}) )
    vals = [];
    return
end

% Determine the number of columns:
numCols = size(funs{1}, 2);

% Initialise vals:
vals = zeros(numFuns+1, numCols);

if ( (nargin < 3) || isnumeric(op) || iscell(op) )
    % Function handle not provided.

    vals(1,:) = lval(funs{1});
    % Take the mean of the FUNs on either side of the break:
    for k = 2:numFuns
        vals(k,:) = ( rval(funs{k-1}) + lval(funs{k})) / 2;
    end
    vals(numFuns+1,:) = rval(funs{numFuns});

else
    % Function handle provided.

    % Evaluate the function handle at the breaks:
    vals(1:numFuns+1,:) = feval(op, ends.');
    
    %% 
    % When functions are defined on unbounded domains, indeterminate form may
    % occur and the evaluation of the function handle may return NaNs. For
    % example, NaN is return when we try to evaluate x*exp(x) at -Inf.
    % Therefore, we need to replace all NaNs by the appropriate function values
    % obtained using FUNs, i.e., by interpolation:
    
    mask = isnan(vals(1:numFuns+1,:));
    if ( any(mask(1,:)) )
        lvals = lval(funs{1});
        vals(1,mask(1,:)) = lvals(mask(1,:));
    end
    
    for k = 2:numFuns
        if ( any(mask(k,:)) )
            lrvals = (rval(funs{k-1}) + lval(funs{k}))/2;
            vals(k,mask(k,:)) = lrvals(mask(k,:));
        end
    end
    
    if  ( any(mask(numFuns+1,:)) )
        rvals = rval(funs{numFuns});
        vals(numFuns+1,mask(numFuns+1,:)) = rvals(mask(numFuns+1,:));
    end
    
end

end
