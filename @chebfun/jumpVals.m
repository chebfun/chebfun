function jVals = jumpVals(funs, ends, op)
%JUMPVALS   Determine values between FUN objects of neighbouring domains.
%   JVALS = JUMPVALS(FUNS, ENDS, OP) returns the values at breakpoints ENDS
%   between FUN objects, which will be stored in the first row of the .IMPULSES
%   field of the CHEBFUN with .funs  = FUNS.
%
%   If OP can be evaluated at ENDS then JVALS = OP(ENDS). Otherwise JVALS(j) is
%   the average of the right and left limits of its neighbouring funs for
%   interior breaks and the left and right limits for the JVALS(1) and
%   JVALS(end), respectively.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Determine the number of intervals:
numFuns = numel(funs);

% Trivial empty case:
if ( (numFuns == 1) && isempty(funs) )
    jVals = [];
    return
end

% Determine the number of columns:
numCols = size(funs{1}, 2);

% Initialise jVals:
jVals = zeros(numFuns+1, numCols);

if ( (nargin < 3) || isnumeric(op) || iscell(op) ) % Function handle not
                                                   % provided.

    jVals(1,:) = get(funs{1}, 'lval');
    % Take the mean of the funs on either side of the break:
    for k = 2:numFuns
        jVals(k,:) = (get(funs{k-1}, 'rval') + get(funs{k}, 'lval'))/2;
    end
    jVals(numFuns+1,:) = get(funs{numFuns}, 'rval');

elseif ( (nargin == 3) && ~iscell(op) ) % Function handle provided.

    % Evaluate the function handle at the breaks:
    jVals(1:numFuns+1,:) = feval(op, ends.');

end

end
