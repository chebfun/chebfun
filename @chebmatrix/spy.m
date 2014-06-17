function spy(A, varargin)
%SPY    Visualize a CHEBMATRIX.
%   SPY(A) creates a picture of the nonzero pattern of the default
%   discretization of the CHEBMATRIX A. Block boundaries are indicated by gray
%   lines.
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
% Example:
%   f = chebmatrix({diag(x)*operatorBlock.diff cos(x) ; functionalBlock.sum 2})
%   spy(f, 'xr', 15, 'disc', @ultraS, 'dom', [-1 0 1], 'dim', 18)
%
% See also CHEBMATRIX, CHEBMATRIX.MATRIX, CHEBOPPREF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Obtain domain information.
dom = A.domain;

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
        dom = A.mergeDomains(dom, varargin{k+1});
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

% Discretize and do a regular spy plot.
data = matrix(A, dim, dom, pref.discretization);
spy(data, varargin{:}); hold on
s =  sprintf('%i,', dim);    % list of sizes
s = [ 'piecewise dimension = [', s(1:end-1), ']' ];
xlabel(s)

% Find all block sizes, substituting in the discretization size for Inf.
[m, n] = blockSizes(A);
m(isinf(m)) = sum(dim);  
n(isinf(n)) = sum(dim);

% Draw horizontal block boundaries.
if ( size(A.blocks, 1) > 1 )
    csrow = cumsum(m(:,1)');
    rowdiv = csrow(1:end-1) + 1/2;   % insert boundary after each block
    colmax = sum(n(1,:));
    plot([0 ; colmax+1], [rowdiv ; rowdiv], 'color', [.6 .6 .6])
end

% Draw vertical block boundaries.
if ( size(A.blocks, 2) > 1 )
    cscol = cumsum(n(1,:));
    coldiv = cscol(1:end-1) + 1/2;
    rowmax = sum(m(:,1));
    plot([coldiv ; coldiv], [0 ; rowmax+1], 'color', [.6 .6 .6])
end

% Clean up hold state.
if ( ~holdState )
    hold off
end

end
