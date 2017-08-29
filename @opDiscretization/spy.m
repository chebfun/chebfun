function spy(disc, varargin)
%TODO: Either remove this method, or uncomment lines below and ensure this
%method works. The following is broken:
%   D = chebop(@(x,u) diff(u));
%   spy(chebcolloc2(linop(D)))

error('CHEBFUN:OPDISCRETIZATION:spy:noSupported', ...
    'OPDISCRETIZATION currently does not support spy().');
end
% %SPY    Visualize a LINOP.
% %   SPY(DISC) creates a picture of the nonzero pattern of the given LiNOP
% %   discretization DISC. Block boundaries are indicated by gray lines, and
% %   side condition rows are marked off by dashed lines (boundary and
% %   continuity conditions).
% %
% %   SPY(DISC,DIM) uses the dimension vector DIM to override the dimension
% %   property of DISC.
% %
% %   SPY(DISC,S) or SPY(DISC,DIM,S) allows modification of the plot
% %   attributes, as with the built in method.
% %
% % See also LINOP.SPY.
% 
% % Copyright 2017 by The University of Oxford and The Chebfun Developers.
% % See http://www.chebfun.org for Chebfun information.
% 
% % Obtain domain information.
% dom = disc.domain;
% 
% % Parse the inputs:
% if ( ( nargin > 1 ) && isnumeric(varargin{1}) )
%     dim = varargin{1};
%     if ( numel(dim) == 1 )
%         % Scalar expand the dimension here.
%         dim = repmat(dim, 1, length(dom) - 1);
%     end
% 
%     disc.dimension = dim;
%     varargin = varargin(2:end);
% end
% dim = disc.dimension;
%     
% plotOpt = varargin;
% 
% % Override hold state.
% holdState = ishold;
% 
% % Check whether we need to derive continuity conditions.
% L = disc.source;
% if ( isempty(L.continuity) )
%     L = deriveContinuity(L, dom);
% end
% disc.source = L;
% 
% % Spy the matrix, with a useful label.
% data = matrix(disc);
% spy(data, plotOpt{:}), hold on
% s =  sprintf('%i,', disc.dimension);    % list of sizes
% s = [ 'discretization = [', s(1:end-1), ']' ];
% xlabel(s)
% 
% % Find the number of constraints and continuity conditions
% nbc = length(L.constraint);
% ncon = length(L.continuity);
% numInts = numel(dom) - 1;
% 
% % Find all block sizes, substituting in the discretization size for Inf.
% [m, n] = blockSizes(L);
% m(isinf(m)) = sum(dim);
% n(isinf(n)) = sum(dim) + numInts*disc.dimAdjust(isinf(n));
% 
% % Draw vertical block boundaries.
% cscol = cumsum(n(1,:));
% coldiv = cscol(1:end-1) + 1/2;
% rowmax = size(data, 1);
% plot([coldiv ; coldiv], [ 0; rowmax+1]*ones(size(coldiv)), 'color', [.6 .6 .6])
% 
% % Draw horizontal block boundaries. Account for the down-sampling of each row.
% csrow = cumsum( m(:, 1)');                          % remove down-sampling
% rowdiv = csrow(1:end - 1) + 1/2;                    % boundary after each block
% rowdiv = nbc + ncon + rowdiv;                       % offset from top rows
% colmax = size(data, 2);
% % plot([0; colmax+1]*ones(size(rowdiv)), [rowdiv; rowdiv], 'color', [.6 .6 .6])
% 
% % Draw horizontal BC and continuity boundaries.
% y = nbc + 1/2;
% plot([0; colmax+1], [y; y], '--', 'color', [.6 .6 .6])
% if ( ncon > 0 )
%     plot([0; colmax + 1], [y + ncon; y + ncon], '--', 'color', [.6 .6 .6])
% end
% 
% % Reset hold state.
% if ( ~holdState )
%     hold off
% end
% 
% end
% 
