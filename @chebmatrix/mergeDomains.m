function d = mergeDomains(varargin)

% foo
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

if iscell(varargin{1})
    % A cell was passed in rather than multiple arguments.
    blocks = varargin{1};
else
    blocks = varargin;
end

% Discard numerical blocks:
isnum = cellfun(@isnumeric,blocks);
blocks(isnum) = [];

% Find the domains of each remaining block (output is cell).
d = cellfun(@(x) x.domain,blocks,'uniform',false);

d = chebfun.mergeDomains(d{:});

end
