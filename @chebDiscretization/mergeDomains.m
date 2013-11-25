function disc = mergeDomains(disc,varargin)

if isnumeric(varargin{1})
    d = varargin(1);
else
    % Discard numerical blocks:
    isn = cellfun(@isnumeric,varargin);
    varargin(isn) = [];
    % Find the domains of each block (output is cell):
    d = cellfun(@(A) A.domain,varargin,'uniform',false);
end

% Collect the endpoints and take the outer hull.
leftEnds = cellfun(@(x) x(1),d);
left = min(leftEnds(:));
rightEnds = cellfun(@(x) x(end),d);
right = max(rightEnds(:));

% We want to soften 'equality' relative to the domain length.
tol = 100*eps*(right-left);

% Check to see if the domain endpoints are genuinely different.
if ( max( abs(leftEnds(:)-left) ) > tol ) || ( max( abs(rightEnds(:)-right) ) > tol )
    error('Domain endpoints are not compatible.')
end

% Extract all the interior breakpoints.
d = cellfun(@(x) x(2:end-1),d,'uniform',false);

% Find the unique ones (sorted).
breakpoints = cat(2,d{:});
breakpoints = unique(breakpoints);

if ~isempty(breakpoints)
    % Remove all too close to the left endpoint.
    isClose = ( breakpoints - left < tol );
    breakpoints(isClose) = [];
    
    % Remove all too close to the right endpoint.
    isClose = ( right - breakpoints < tol );
    breakpoints(isClose) = [];
    
    % Remove interior points too close to one another.
    isClose =  find( diff(breakpoints) < tol  );
    breakpoints(isClose) = [];
end

% Put it all together.
disc.domain = [left breakpoints right];

end
