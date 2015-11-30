function spy(L, varargin)
%SPY    Visualize a LINOP.
%   SPY(A) creates a picture of the nonzero pattern of the default
%   discretization of the LINOP A. Block boundaries are indicated by gray lines,
%   and side condition rows are marked off by dashed lines (boundary and
%   continuity conditions).
%
%   SPY(A, S) allows modification of the SPY plot as with the built in method.
%   
%   SPY(A, 'dimension', DIM, ...) uses the dimension vector DIM and SPY(A,
%   'disc', DISCTYPE, ...) uses the discretization DISCTYPE for the
%   visualization. All optional inputs can be used in tandem.
%
%   SPY(A, PREFS, ...), where PREFS is a CHEBOPPREF object, modifies the default
%   discretization type and dimension.
%
% See also LINOP.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Obtain domain information.
dom = L.domain;

% Look for a preference input:
isPref = cellfun(@(v) isa(v, 'cheboppref'), varargin);
if ( any(isPref) )
    prefIdx = find(isPref, 1);
    pref = varargin{prefIdx};
    varargin(prefIdx) = [];
    dim = pref.dimension(1);
else
    % Take the default preferences:
    pref = cheboppref();
    if ( length(dom) == 2 )
        dim = 10;
    else
        dim = 6;
    end
end

k = 1;
% Parse out other inputs:
while ( k < numel(varargin) )
    if ( strncmpi(varargin{k}, 'dimension', 3) )
        dim = varargin{k+1};
        varargin(k:k+1) = [];
    elseif ( strncmpi(varargin{k}, 'domain', 3) )
        dom = L.mergeDomains(dom, varargin{k+1});
        varargin(k:k+1) = [];
    elseif ( strncmpi(varargin{k}, 'discretization', 4) )
        discType = varargin{k+1};
        if ( ischar(discType) )
            discType = eval(['@' discType]);
        end
        pref.discretization = discType;
        varargin(k:k+1) = [];        
    else
        k = k + 1;
    end
end        

if ( numel(dim) == 1 )
    % Scalar expand the dimension here.
    dim = repmat(dim, 1, length(dom) - 1);
end

% Override hold state.
holdState = ishold;

disc = pref.discretization(L, dim, dom);

% Check whether we need to derive continuity conditions.
if ( isempty(L.continuity) )
    L = deriveContinuity(L, dom);
end
disc.source = L;

% Spy the matrix, with a useful label.
data = matrix(disc);
spy(data, varargin{:}), hold on
s =  sprintf('%i,', dim);    % list of sizes
s = [ 'discretization = [', s(1:end-1), ']' ];
xlabel(s)

% Find the number of constraints and continuity conditions
nbc = length(L.constraint);
ncon = length(L.continuity);
numInts = numel(dom) - 1;

% Find all block sizes, substituting in the discretization size for Inf.
[m, n] = blockSizes(L);
m(isinf(m)) = sum(dim);
n(isinf(n)) = sum(dim) + numInts*disc.dimAdjust(isinf(n));

% Draw vertical block boundaries.
cscol = cumsum(n(1,:));
coldiv = cscol(1:end-1) + 1/2;
rowmax = size(data, 1);
plot([coldiv ; coldiv], [ 0; rowmax+1]*ones(size(coldiv)), 'color', [.6 .6 .6])

% Draw horizontal block boundaries. Account for the down-sampling of each row.
csrow = cumsum( m(:, 1)');                          % remove down-sampling
rowdiv = csrow(1:end - 1) + 1/2;                    % boundary after each block
rowdiv = nbc + ncon + rowdiv;                       % offset from top rows
colmax = size(data, 2);
plot([0; colmax+1]*ones(size(rowdiv)), [rowdiv; rowdiv], 'color', [.6 .6 .6])

% Draw horizontal BC and continuity boundaries.
y = nbc + 1/2;
plot([0; colmax+1], [y; y], '--', 'color', [.6 .6 .6])
if ( ncon > 0 )
    plot([0; colmax + 1], [y + ncon; y + ncon], '--', 'color', [.6 .6 .6])
end

% Reset hold state.
if ( ~holdState )
    hold off
end

end

% TODO: The lines below do filled-block spy-plots for LINOP. Probably want to
% remove this.
% 
% function spy(L, varargin)
% %SPY    Visualize a LINOP.
% %   SPY(A) creates a picture of the sparsity and constraint patterns of the
% %   LINOP A. 
% %
% %   SPY(A, S) allows modification of the plot attributes, as with the built
% %   in method.
% %   
% % See also OPDISCRETIZATION.SPY.
% 
% % Copyright 2015 by The University of Oxford and The Chebfun Developers.
% % See http://www.chebfun.org/ for Chebfun information.
% 
% % Obtain domain information.
% dom = L.domain;
% 
% % Look for a preference input:
% isPref = cellfun(@(v) isa(v, 'cheboppref'), varargin);
% if ( any(isPref) )
%     prefIdx = find(isPref, 1);
%     pref = varargin{prefIdx};
%     varargin(prefIdx) = [];
%     dim = pref.dimension(1);
% else
%     % Take the default preferences:
%     pref = cheboppref();
%     if ( length(dom) == 2 )
%         dim = 10;
%     else
%         dim = 6;
%     end
% end
% 
% j = 1;
% % Parse out other inputs:
% while ( j < numel(varargin) )
%     if ( strncmpi(varargin{j}, 'dimension', 3) )
%         dim = varargin{j+1};
%         varargin(j:j+1) = [];
%     elseif ( strncmpi(varargin{j}, 'domain', 3) )
%         dom = L.mergeDomains(dom, varargin{j+1});
%         varargin(j:j+1) = [];
%     elseif ( strncmpi(varargin{j}, 'discretization', 4) )
%         discType = varargin{j+1};
%         if ( ischar(discType) )
%             discType = eval(['@' discType]);
%         end
%         pref.discretization = discType;
%         varargin(j:j+1) = [];        
%     else
%         j = j + 1;
%     end
% end        
% 
% if ( numel(dim) == 1 )
%     % Scalar expand the dimension here.
%     dim = repmat(dim, 1, length(dom) - 1);
% end
% 
% % Override hold state.
% holdState = ishold;
% 
% disc = pref.discretization(L, dim, dom);
% 
% % Check whether we need to derive continuity conditions.
% if ( isempty(L.continuity) )
%     L = deriveContinuity(L, dom);
% end
% disc.source = L;
% 
% % Get the matrix.
% data = matrix(disc);
% %spy(data, varargin{:}), hold on
% % s =  sprintf('%i,', dim);    % list of sizes
% % s = [ 'discretization = [', s(1:end-1), ']' ];
% % xlabel(s)
% 
% % TODO: Better documentation.
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
% hold on
% 
% % Draw horizontal block boundaries. Account for the down-sampling of each row.
% csrow = cumsum( m(:, 1)');                          % remove down-sampling
% rowdiv = csrow(1:end - 1) + 1/2;                    % boundary after each block
% rowdiv = nbc + ncon + rowdiv;                       % offset from top rows
% colmax = size(data, 2);
% plot([0; colmax+1]*ones(size(rowdiv)), [rowdiv; rowdiv], 'color', [.6 .6 .6])
% 
% % Draw horizontal BC and continuity boundaries.
% y = nbc + 1/2;
% plot([0; colmax+1], [y; y], '--', 'color', [.6 .6 .6])
% if ( ncon > 0 )
%     plot([0; colmax + 1], [y + ncon; y + ncon], '--', 'color', [.6 .6 .6])
% end
% 
% % Draw filled blocks:
% % TODO: Not fill empty BC blocks.
% % TODO: Not fill empty subblocks for piecewise domains
% cscol = [0 cscol];
% csrow = ncon + nbc + [0 csrow];
% a = .65;
% b = .35;
% 
% % BC and continuity blocks. 
% for j = 0:nbc-1
%     for k = 1:(numel(cscol)-1)
%         fill([cscol(k)+a cscol(k)+a cscol(k+1)+b cscol(k+1)+b], ...
%             j+b+[a b b a], 'b', 'facealpha', 0.8, 'edgealpha', 0)
%     end
% end
% for j = 0:ncon-1
%     for k = 1:(numel(cscol)-1)
%         fill([cscol(k)+a cscol(k)+a cscol(k+1)+b cscol(k+1)+b], ...
%             nbc+j+b+[a b b a], 'b', 'facealpha', 0.8, 'edgealpha', 0)
%     end
% end
% 
% % Main operator blocks.
% for k = 1:(numel(cscol)-1)
%     for j = 1:(numel(csrow)-1)
%         if ( myiszero(L.blocks{j,k}) )
%             continue
%         end
%         fill([cscol(k)+a cscol(k)+a cscol(k+1)+b cscol(k+1)+b], ...
%             [csrow(j)+a csrow(j+1)+b csrow(j+1)+b csrow(j)+a], 'b', ...
%             'facealpha', 1, 'edgealpha', 0)
%     end
% end
% 
% axis equal, axis square, axis tight
% set(gca,'xtick',[],'ytick',[])  % TODO: Ticks by blocks
% 
% 
% % Reset hold state.
% if ( ~holdState )
%     hold off
% end
% 
% end
% 
% function out = myiszero(f)
% if ( isnumeric(f) )
%     out = any(f);
% else
%     out = iszero(f);
% end
% end
