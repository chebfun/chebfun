function out = feval(f, x, varargin)
%FEVAL   Evaluate a CHEBMATRIX of CHEBFUN objects and doubles.
%   FEVAL(F, x) evaluates the CHEBMATRIX F at the entries of the matrix X. Since
%   it doesn't make sense to evaluate an operator at a point, an error if thrown
%   if F contains a LINBLOCK object.

% % Throw an error if there are LINBLOCKs:
% if ( any(cellfun(@(b) isa(b, 'linBlock'), f.blocks)) )
%     error('Cannot FEVAL() a linBlock');
% end

% % Store the evaluation of each block in a cell for now:
% outCell = cell(size(f.blocks));
% 
% % Loop over the blocks:
% for k = 1:numel(f.blocks)
%     fk = f.blocks{k};
%     if ( isnumeric(fk) )
%         % Scalars evaluate to themselves:
%         outCell{k} = repmat(fk, size(x)); % Scalar expansion.
% %         outCell{k} = fk;                  % No scalar expansion.
%     else
%         % Call the FEVAL() method of the block:
%         outCell{k} = feval(fk, x, varargin{:});
%     end
% end
% 
% % Convert the cell to a matrix. Dimensions will be OK.
% out = cell2mat(outCell);

b = blockSizes(f);
if ( ~isscalar(x) && any(any(cellfun(@(b) isinf(b(:,2)), b))) )
    error('Can only evaluate operators at scalar values.')
end

out = f.blocks;
for k = 1:numel(f.blocks)
    fk = f.blocks{k};
    if ( isnumeric(fk) )
        % Scalars evaluate to themselves:
        out{k} = repmat(fk, size(x)); % Scalar expansion.
%         outCell{k} = fk;                  % No scalar expansion.
    elseif ( isa(fk, 'functionalBlock') )
        % Do nothing
    elseif  ( isa(fk, 'operatorBlock') )
        out{k} = linBlock.feval(x, fk.domain)*fk;
    elseif ( ~fk.isTransposed )
        % Call the FEVAL() method of the block:
        out{k} = feval(fk, x, varargin{:});
    end
end

if ( cellfun(@isnumeric, out) )
    out = cell2mat(out);
else
    out = chebmatrix(out);
end


end