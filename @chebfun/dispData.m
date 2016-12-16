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

% Preallocation:
exps = zeros(numFuns, 2);

% Loop over each FUN:
for j = 1:numFuns
    infoJ = dispData(f.funs{j});
    if ( ~isempty(infoJ) )
        numInfo = numel(infoJ);
        for k = 1:numInfo
            infoJK = infoJ{k};
            if ( strcmpi('exponents', infoJK.name) )
                exps(j,:) = infoJK.data;
            end
        end
    end
end

%% Organize the output which will be passed to DISPLAY:

% Preallocate the extra information:
name = ' ';
data = cell(j,1);
    
% So far, the only extra information is exponents:
if ( any( exps(:) ) )
    name = '  endpoint exponents';
    for j = 1:numFuns
        data{j} = ['        ' '[' num2str(exps(j,1), '%2.2g') '      ' ...
            num2str(exps(j,2), '%2.2g') ']' '  '];
    end
elseif ( ~isempty(infoJ) )
    name = infoJ{1}.name;
end

% More information for F can be appended to INFO:

end
