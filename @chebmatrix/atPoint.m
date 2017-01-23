function out = atPoint(A, x, varargin)
%FEVAL   Left-multiply a CHEBMATRIX by a point evaluation.
%   Each row of a CHEBMATRIX A returns either a function or a scalar value.
%   FEVAL(A, X) essentially pre-multiplies A with a point evaluation functional
%   for each row that corresponds to a function output. Rows with scalar outputs
%   are not affected.
%
%   FEVAL(A, X, DIRECTION) uses the direction implied by the third argument for
%   the evaluation (important at a breakpoint in the domain).
%
% See also LINOP.FEVAL. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = A.blocks;
[m, n] = blockSizes(A);
for i = find( isinf(m(:, 1)') )
    for j = 1:size(A, 2)
        item = out{i,j};
        if ( isa(item, 'chebfun') )
            out{i,j} = item(x);
        elseif ( isa(item, 'linBlock') )
            out{i,j} = linBlock.feval(x, item. domain,varargin{:}) * item;
        end
    end
end

% Return a matrix if it's all scalar values. 
if ( cellfun(@isnumeric, out) )
    out = cell2mat(out);
else
    out = chebmatrix(out);
end

end
