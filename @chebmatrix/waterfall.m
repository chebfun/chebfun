function varargout = waterfall(varargin)
%WATERFALL   Waterfall plot for CHEBMATRIX object.
%   WATERFALL(U), or WATERFALL(U, T) where LENGTH(T) = MIN(SIZE(U)), plots a
%   "waterfall" plot of the CHEBMATRIX U. If U cannot be converted to a
%   QUASIMATRIX (i.e., if it contains INFxINF blocks), then an error is thrown.
%
%   WATERFALL(U, T, PROP1, VAL1, PROP2, VAL2, ...) allows additional plotting
%   options. For further details see CHEBFUN/WATERFALL.
%
%   WATERFALL(U, ..., 'EdgeColors', COLS) allows specification of the Edge
%   Colors for each of the rows in U. If COLS is a standard color string (e.g.,
%   'b') or a 1x3 vector, this functions the same as 'EdgeColor'. If COLs is a
%   cell array or a matrix with the same number of rows as U, the kth row of U
%   is plotted in the color COLS{k} or COLS(k,:).
%
% See also PLOT, PLOT3.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% % First input might be a figure handle:
% [cax, varargin] = axescheck(varargin{:});
% if ( ~isempty(cax) )
%     axes(cax);
% end

% First input is now the CHEBMATRIX:
u = varargin{1};
varargin(1) = [];

% Check that we have a valid CHEBMATRIX to plot:
s = cellfun(@(b) min(size(b)), u.blocks);
if ( ~all(isfinite(s)) )
    error('CHEBFUN:CHEBMATRIX:waterfall:dim', ...
        ['WATERFALL only supports CHEBMATRICES whose blocks have at least ', ...
         'one finite dimension.']);
end

% Parse inputs:
col = [];
k = 1;
while ( k <= numel(varargin) )
    if ( strcmpi(varargin{k}, 'edgecolors') )
        col = varargin{k+1};
        varargin(k:k+1) = [];
    else
        k = k+1;
    end
end

% Transpose if we are given an column CHEBMATRIX:
numRows = size(u,1);
if ( (numRows > 1) && (size(u,2) == 1) )
    u = u.';
    numRows = size(u,1);
end

% Initialise:
h = cell(numRows, 1);
colorData = {};
holdState = ishold();

% If we were given what looks like a system, create colordata
if ( all(size(u) > 1) )
    % Use default colors
    col = get(0, 'DefaultAxesColorOrder');
end

% Loop over the rows:
for k = 1:numRows
    
    % Assign colors if given:
    if ( ~isempty(col) )
        if ( iscell(col) )
            colorData = {'edgecolor', col{k}};
        elseif ( size(col,1) == 1 )
            colorData = {'edgecolor', col};
        elseif ( isnumeric(col) )
            colorData = {'edgecolor', col(k,:)};
        end
    end
    
    % Convert to a quasimatrix:
    uk = u.blocks(k,:);     % Get the blocks of this row.
    uk = quasimatrix(uk);   % Convert them to a quasimatrix.
    
    % Call CHEBFUN/WATERFALL():
    h{k} = waterfall(uk, varargin{:}, colorData{:});
    hold on
    
end

% Reset HOLD state:
if ( ~holdState )
    hold off
end

% Output handles:
if ( nargout > 1 )
    h = cell2mat(h);
    varargout{1} = h;
end

end
