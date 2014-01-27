function out = horzcat(varargin)
%HORZCAT   Horizontal concatenation of CHEBFUN objects.
%   [A B] horizontally concatenates the CHEBFUN objects A and B to form an
%   array-valued CHEBFUN or an array of CHEBFUN objects (depending on whether
%   the interior breakpoints of A and B match or not). [A,B] does the same. Any
%   number of CHEBFUN objects can be concatenated within one pair of brackets.
%
% See also VERTCAT, CAT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.
% error('CHEBFUN:horzcat:noSupport', 'HORZCAT of a CHEBFUN is not yet supported.');

% TODO: Document quasimatrix vs array-valued CHEBFUN

% [TODO]: Vertical concatenation. (Chebmatrix / quasimatrix).

% Remove empties:
empties = cellfun(@isempty, varargin);
if ( all(empties) )
    out = varargin{1};
    return
else
    varargin(empties) = [];
end

% Find the locations of the CHEBFUN objects in the inputs:
if ( numel(varargin) == 1 )
    out = varargin{1};
    return
end

% Promote doubles to CHEBFUN objects:
chebfunLocs = cellfun('isclass', varargin, 'chebfun');
chebfun1 = varargin{find(chebfunLocs, 1, 'first')};

% Horizontal concatenation of row CHEBFUN objects produces a CHEBMATRIX:
if ( chebfun1(1).isTransposed )
    out = chebmatrix(varargin);
    return
end

% Promote doubles to CHEBFUN objects:
domain1 = chebfun1.domain;
doubleLocs = find(~chebfunLocs);
for k = doubleLocs
    varargin{k} = chebfun(varargin{k}, domain1);
end

%%%%%%%%%%%%%%%%%%%%%%%% Deal with qausimatrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numElements = cellfun(@numel, varargin);
if ( any(numElements > 1) )
    args = {};
    for k = 1:numel(varargin)
        args = [args, num2cell(varargin{k})]; %#ok<AGROW>
    end
    out = horzcat(args{:});
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Deal with domains %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grab the domains of each of the inputs:
allDomainsCell = cellfun(@(f) f.domain, varargin, 'UniformOutput', false);

% Ensure that the domains match:
domainEnds = allDomainsCell{1}([1 end]);
if ( any(cellfun(@(d) any(d([1 end]) ~= domainEnds), allDomainsCell)) )
    error('CHEBFUN:horzcat:domains', 'Inconsistent domains.');
end
% [TODO]: checkDomains?

% Check to see if interior breakpoints differ:
differentBreakpoints = false;
if ( any(diff(cellfun(@(d) length(d), allDomainsCell))) )
    differentBreakpoints = true;
else
    tol = max(cellfun(@(f) hscale(f).*epslevel(f), varargin));
    if ( any(cellfun(@(d) any(d - domain1) > tol, allDomainsCell)) )
        differentBreakpoints = true;
    end
end

% TODO: Also check to see if an input is a SINGFUN.
isSingular = any(cellfun(@issing, varargin));

%%%%%%%%%%%%%%%%%%%%%%%%%%% FORM A QUASIMATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( differentBreakpoints || isSingular )  % (form a quasimatrix)
    isArrayCheb = cellfun(@(f) isa(f, 'chebfun') && size(f, 2) > 1, varargin);
    if ( any(isArrayCheb) )
        % Break up array-valued CHEBFUNs into single columns:
        args = {};
        for k = 1:numel(varargin)
            args = [args, num2cell(varargin{k})]; %#ok<AGROW>
        end
    else
        args = varargin;
    end
    numCols = numel(args);
    % Initialise CHEBFUN array:
    clear out
    out(1, numCols) = chebfun();
    % Assign columns/rows:
    for k = 1:numCols
        out(k) = args{k};
    end
    
%%%%%%%%%%%%%%%%%%%%%%% FORM AN ARRAY-VALUED CHEBFUN %%%%%%%%%%%%%%%%%%%%%%%%%%%

else % (form an array-valued CHEBFUN) 
    
    % Concatenate the FUNs:
    out = varargin{1};
    numInts = numel(out.domain) - 1;
    for k = 1:numInts
        funs = cellfun(@(f) f.funs{k}, varargin, 'UniformOutput', false);
        out.funs{k} = horzcat(funs{:});
    end
    % Concatenate impulses:
    imps = cellfun(@(f) f.impulses, varargin, 'UniformOutput', false);
    out.impulses = cell2mat(imps);

end

end
