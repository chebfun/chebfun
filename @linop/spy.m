function spy(L,dim,prefs)
%SPY    Visualize a linop.
%   SPY(L) creates a picture of the nonzero pattern of the default
%   discretization of L. Block boundaries are indicated by gray lines, and
%   side condition rows are marked off by dashed lines (boundary and
%   continuity conditions).
%
%   SPY(L,DIM) uses the dimension vector DIM to create the picture.
%
%   SPY(L,DIM,PREFS) uses a preferences structure or object like that created by
%   CHEBOPPREF. This allows you to change the type of discretization used for
%   the visualization.
%
%   See also LINOP, CHEBOPPREF.

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Set defaults as needed.
if ( nargin < 3 )
    prefs = cheboppref;
    if ( nargin < 2 )
        if ( length(L.domain) == 2 )
            dim = 10;
        else
            dim = repmat(6, 1, length(L.domain) - 1);
        end
    end
end

disc = prefs.discretization(L);
disc.dimension = dim;

% Check whether we need to derive continuity conditions.
if isempty(L.continuity)
    disc.source = deriveContinuity(L);
end

% Spy the matrix, with a useful label.
M = matrix(disc);
spy(M)
s =  sprintf('%i,', dim);    % list of sizes
s = [ 'discretization = [', s(1:end-1), ']' ];
xlabel(s)

% Find the number of constraints and continuity conditions
nbc = length(L.constraint);
ncon = length(L.continuity);

% Reverse engineer to see where the block boundaries are.
[m, n] = blockSizes(L);
m(isinf(m)) = sum(dim);   % use discretization sizes
n(isinf(n)) = sum(dim);

% Override hold state.
holdState = ishold;
hold on

% Find all block sizes, substituting in the discretization size for Inf.
[m, n] = blockSizes(L);
m(isinf(m)) = sum(dim);
n(isinf(n)) = sum(dim);

% Draw vertical block boundaries.
cscol = cumsum(n(1, :));
coldiv = cscol(1:end-1) + 1/2;
rowmax = sum(m(:, 1));
plot([coldiv;coldiv], [0; rowmax+1]*ones(size(coldiv)), 'color', [.6 .6 .6])

% Draw horizontal block boundaries. Account for the down-sampling of each row.
sizeRedux = getProjOrder(L);                        % could be a row or column?
csrow = cumsum( m(:, 1)' - sizeRedux(:)' );         % remove down-sampling
rowdiv = csrow(1:end - 1) + 1/2;                    % boundary after each block
rowdiv = nbc + ncon + rowdiv;                       % offset from top rows
colmax = sum(n(1, :));
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
