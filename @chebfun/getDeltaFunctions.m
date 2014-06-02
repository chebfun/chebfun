function [deltaMag, deltaLoc] = getDeltaFunctions(f)
%GETDELTAFUNCTIONS   Get the delta functions in a CHEBFUN F.
%   [DELTAMAG, DELTALOC] = GETDELTAFUNCTIONS(F) the returns the delta
%   functions within a CHEBFUN F. DELTAMAG is a matrix containing the (signed)
%   magnitude of the delta functions and their higher derivatives, while
%   DELTALOC is a an array containing the location of delta functions.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Initialize variables:
deltaMag = [];
deltaLoc = [];
maxRows = 0;
numCols = 0;
numFuns = numel(f.funs);

% Loop through the FUNs of F:
for j = 1:numFuns
    % If the current FUN is a deltafun:
    if ( isa(f.funs{j}, 'deltafun') )
        d = f.funs{j};
        % If there are deltafunctions in it:
        if ( ~isempty(d.deltaMag) )
            m = size(d.deltaMag, 1);
            % Get the number of columns and rows:
            numCols = numCols + size(d.deltaMag, 2);
            if ( m > maxRows )
                maxRows = m;
            end                        
        end        
    end
end

if ( maxRows > 0 )
    deltaMag = zeros(maxRows, numCols);
    deltaLoc = zeros(1, numCols);
    k = 1;
    for j = 1:numFuns
        % If the current FUN is a deltafun:
        if ( isa(f.funs{j}, 'deltafun') )
            d = f.funs{j};
            if ( ~isempty(d.deltaMag) )
                % Get the delta functions:
                sz = size(d.deltaMag);
                deltaMag(1:sz(1), k:k+sz(2)-1) = d.deltaMag;
                deltaLoc(1, k:k+sz(2)-1) = d.deltaLoc;
                k = k + sz(2);
            end
        end
    end
    
    % Merge delta functions if they are at the same location and delete trivial
    % columns and rows:
    [deltaMag, deltaLoc] = deltafun.mergeColumns(deltaMag, deltaLoc);
    [deltaMag, deltaLoc] = deltafun.cleanColumns(deltaMag, deltaLoc);
    deltaMag = deltafun.cleanRows(deltaMag);
end

end
