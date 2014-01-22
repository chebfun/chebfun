function spy(L,dim)
%SPY    Visualize a linop.
%   SPY(L) creates a picture of the nonzero pattern of the default
%   discretization of L. Block boundaries are indicated by gray lines, and
%   side condition rows are marked off by dashed lines (boudnary and
%   continuity conditions).
%   
%   SPY(L,DIM) uses the dimension vector DIM to create the picture.
%
%   See also LINOP.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

if ( nargin < 2 )
    if ( length(L.domain) == 2 )
        dim = 10;
    else
        dim = repmat(6,1,length(L.domain)-1);
    end
end

disc = L.prefs.discretization(L);
disc.dimension = dim;

if isempty(L.continuity)
    disc.source = deriveContinuity(L);
end

% Spy the matrix, with a useful label.
M = matrix(disc);
spy(M)
s =  sprintf('%i,',dim);    % list of sizes
s = [ 'discretization = [', s(1:end-1), ']' ];
xlabel(s)

nbc = length(L.constraint);
ncon = length(L.continuity);

% Reverse engineer to see where the block boundaries are.
[m,n] = blockSizes(L);
m(isinf(m)) = sum(dim);   % use discretization sizes
n(isinf(n)) = sum(dim);

% Override hold state.
hs = ishold;
hold on

% Find all block sizes, substituting in the discretization size for Inf.
[m,n] = blockSizes(L);
m(isinf(m)) = sum(dim);  
n(isinf(n)) = sum(dim);

% Draw vertical block boundaries.
cscol = cumsum(n(1,:));
coldiv = cscol(1:end-1) + 1/2; 
rowmax = sum(m(:,1));
plot([coldiv;coldiv],[0;rowmax+1],'color',[.6 .6 .6])

% Draw horizontal block boundaries. Account for the down-sampling of each
% row.
csrow = cumsum( m(:,1)' - sizeReduction(L)' );  % remove down-sampling
rowdiv = csrow(1:end-1) + 1/2;   % boundary after each block
rowdiv = nbc + ncon + rowdiv;    % offset from top rows
colmax = sum(n(1,:));
plot([0;colmax+1],[rowdiv;rowdiv],'color',[.6 .6 .6])

% Draw horizontal BC and continuity boundaries.
y = nbc+1/2;
plot([0;colmax+1],[y;y],'--','color',[.6 .6 .6])
if ( ncon > 0 )
    plot([0;colmax+1],[y+ncon; y+ncon],'--','color',[.6 .6 .6])
end

% Reset hold state.
if ( ~hs )
    hold off
end

end