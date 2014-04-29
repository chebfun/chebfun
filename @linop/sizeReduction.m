function [reduce, d, dRow, dVar] = sizeReduction(L)
%SIZEREDUCTION   Deduce row down-projection dimensions.
%   Each boundary and continuity constraint in a linop forces a reduction in the
%   total number of rows in the discrete operator, so that the composite is
%   square. The reduction is found by down-projection of the result of applying
%   the operator.
%
%   SIZEREDUCTION(L) returns a vector of dimensions by which each column of the
%   system should be down-projected in order to end with a square system. If the
%   differential orders of the variables give the correct total result, they are
%   used; otherwise the reductions are spread as evenly as possible.
%
%   [REDUCE, D, DROW, DVAR] = SIZEREDUCTION(L) also returns the
%   differential order of each variable, the differential order of each
%   system row (equation), and the differential order of each system
%   variable (column).

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

[dummy, d, dRow, dVar] = sizeReduction@chebmatrix(L);
reduce = dVar;           % Nick H's choice;
reduce = max(reduce, 0); % Nick H's choice;

%%
% TODO: Figure out what is going on here!
% 
% [m, n] = size(L);
% dRow(~dRow) = NaN;
% dVar(~dVar) = NaN;
% 
% % This will determine how much is downsampled. 
% numBC = length(L.constraint);
% 
% % Each variable requires a downsampling contribution equal to its
% % differential order. The total diff. orders of the variables might
% % not equal the total of the rows.
% totalRow = sum( dRow(~isnan(dRow)) );
% totalVar = sum( dVar(~isnan(dVar)) );
% 
% totalVar = numBC;  % Alex's cho ice
% if ( totalRow == totalVar )
%     reduce = dRow;
% else
%     % The only reasonable thing is to spread out the D.O. reductions
%     % as evenly as possible among the rows. 
%     
%     % TODO: Update for Nick's choice.
%     reduce = NaN(1,m);
%     isOp = ~isnan(dRow);
%     numOp = sum(isOp);
%     k = floor( totalVar / numOp );
%     reduce(isOp) = k;  % even distribution of order
%     i = rem(totalVar, numOp);  % leftover to be made up
%     if i > 0
%         % Increment the first i rows
%         incr = find(isOp,i);
%         reduce(incr) = reduce(incr) + 1;
%     end
% end
% 
% reduce(isnan(reduce)) = 0;

end
