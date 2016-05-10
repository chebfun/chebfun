function out = horzcat(varargin)
%HORZCAT   Horizontal concatenation of CHEBFUN objects.
%   [A B] horizontally concatenates the column CHEBFUN objects A and B to form
%   an array-valued CHEBFUN or an array of CHEBFUN objects (depending on whether
%   the interior breakpoints of A and B match or not). [A,B] does the same. Any
%   number of CHEBFUN objects can be concatenated within one pair of brackets.
%
%   If one of A or B is a scalar or a row vector, it is cast to a constant
%   CHEBFUN.
%
%   If A and B are row CHEBFUN objects then C = [A B] will be a CHEBMATRIX. See
%   CHEBFUN/VERTCAT for more details.
%
%   [A1, A2, ...] concatenates multiple objects.
%
% See also VERTCAT, CAT, QUASIMATRIX, CHEBMATRIX.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

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
isCheb = cellfun('isclass', varargin, 'chebfun');
chebfun1 = varargin{find(isCheb, 1, 'first')};

% Check transpose state:
if ( ~all(chebfun1(1).isTransposed == ...
        cellfun(@(f) double(f(1).isTransposed), varargin(isCheb)) ) )
    error('CHEBFUN:CHEBFUN:horzcat:transpose', ...
        'Dimensions of matrices being concatenated are not consistent. ');
end

% Horizontal concatenation of row CHEBFUN objects produces a CHEBMATRIX:
if ( chebfun1(1).isTransposed )
    args = cellfun(@transpose, varargin, 'UniformOutput', false);
    out = vertcat(args{:})';
    return
end

% Promote doubles to CHEBFUN objects:
domain1 = chebfun1.domain;
doubleLocs = find(~isCheb);
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
    error('CHEBFUN:CHEBFUN:horzcat:domains', 'Inconsistent domains.');
end
% [TODO]: checkDomains?

% Check to see if interior breakpoints differ:
differentBreakpoints = false;
if ( any(diff(cellfun(@(d) length(d), allDomainsCell))) )
    differentBreakpoints = true;
else
    tol = max(cellfun(@(f) hscale(f)*eps, varargin));
    if ( any(cellfun(@(d) any(d - domain1) > tol, allDomainsCell)) )
        differentBreakpoints = true;
    end
end

% TODO: Also check to see if an input is a SINGFUN or DELTAFUN or TRIGFUN:
isSingular = any(cellfun(@issing, varargin));
isDelta = any(cellfun(@isdelta, varargin));
isAllPeriodic = all(cellfun(@isPeriodicTech, varargin));
isQuasiPeriodic = any(cellfun(@isPeriodicTech, varargin)) && ~isAllPeriodic;

%%%%%%%%%%%%%%%%%%%%%%%%%%% FORM A QUASIMATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( differentBreakpoints || isSingular || isDelta || isQuasiPeriodic )  % (form a quasimatrix)
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
    % Concatenate pointValues:
    imps = cellfun(@(f) f.pointValues, varargin, 'UniformOutput', false);
    out.pointValues = cell2mat(imps);

end

end
