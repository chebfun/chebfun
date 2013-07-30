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
if ( numFuns == 1 && isempty(funs) )
    jVals = [];
    return
end

% Initialise jVals:
jVals = zeros(numFuns+1, size(funs{1}, 2));

if ( nargin < 3 || isnumeric(op) || iscell(op) )    % Function handle is not provided.

    jVals(1,:) = get(funs{1}, 'lval');
    % Take the mean of the funs on either side of the break:
    for k = 2:numFuns
        jVals(k,:) = (get(funs{k-1}, 'rval') + get(funs{k}, 'lval'))/2;
    end
    jVals(numFuns+1,:) = get(funs{numFuns}, 'rval');

elseif ( nargin == 3 && ~iscell(op) ) % Function handle provided.

    % Evaluate the function handle at the breaks:
    jVals(1:numFuns+1,:) = feval(op, ends.');

% elseif ( iscell(op) )                 % Cell given. Check each entry.
%     % Loop over each endpoint:
%     for k = 1:numFuns + 1
%         % Convert from cell to op:
%         opk = op{min(k, numFuns)};
%
%         if ( isa(opk, 'double') )
%             if ( k == 1 )
%                 jVals(1,:) = get(funs{1}, 'lval');
%             elseif ( k == numFuns+1 )
%                 jVals(numFuns+1,:) = get(funs{numFuns}, 'rval');
%             else
%                 % Take the mean of the funs on either side of the break:
%                 jVals(k,:) = (get(funs{k-1}, 'rval') + get(funs{k}, 'lval'))/2;
%             end
%         else
%             % Evaluate the function handle at the break:
%             jVals(k,:) = feval(opk, ends(k));
%         end
%     end
%
end

% [TODO]: Document how .impulses are stored.
jVals = reshape(jVals.', 1, numFuns+1, size(jVals, 2));
% jVals = reshape(jVals.', numFuns+1, size(jVals, 2), 1);

end