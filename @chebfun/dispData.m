function [name, data] = dispData(f)
%DISPDATA   Useful information for DISPLAY.
%   [NAME, DATA] = DISPDATA(F) collects extra information from the given CHEBFUN
%   F and the output will be used by DISPLAY. 
%
% See also DISPLAY.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Number of pieces (i.e., funs) information:
numFuns = numel(f.funs);

%% Collect information about the exponents:

% Loop over each FUN:
for j = 1:numFuns
    tmp = dispData(f.funs{j});
    if  ( ~isempty(tmp) ), tmp = [tmp{:}]; end
    info{j} = tmp;
end

%% Organize the output which will be passed to DISPLAY:

% Preallocate the extra information:
name = ' ';
data = cell(j,1);

% Loop through each fun, look for duplicate names. Insert space if no duplicates.
for j = 1:numFuns
    if ( isempty(info{j}) ), continue, end
    for k = 1:numel(info{j})
        namejk = info{j}(k).name;
        datajk = info{j}(k).data;
        name = [name, namejk];
        z = repmat(' ', size(datajk))
        for l = 1:(j-1)
            data{l,1} = [data{l,1}, z];
        end
        data{j,1} = [data{j,1}, datajk];
        for l = (j+1):numFuns
            isz = true;
            if ( ~isempty(info{l}) )
                for m = 1:numel(info{l})
                    if ( strcmp(namejk, info{l}(m).name) )
                        data{l,1} = [data{l,1}, info{l}(m).data];
                        isz = false;
                        info{l}(m) = [];
                        break
                    end
                end
            end
            isz
            if ( isz )
                data{l,1} = [data{l,1}, z];
            end
        end
    end
end

data

% More information for F can be appended to DATA:

end
