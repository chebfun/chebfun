function [name, data] = dispInfo(f)
%DISPINFO   Useful information for DISPLAY.
%   INFO = DISPINFO(F) collects extra information from the given CHEBFUN F and
%   the output INFO will be used by DISPLAY. 
%
% See also DISPLAY.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Number of pieces (i.e., funs) information:
numFuns = numel(f.funs);

%% Collect information about the exponents:

% Preallocation:
exps = zeros(numFuns, 2);

% Loop over each FUN:
for j = 1:numFuns
    infoJ = dispInfo(f.funs{j});
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
if any( exps(:) )
    name = '  exponents';
    for j = 1:numFuns
        data{j} = ['  ' num2str(exps(j,:), '%5.2g') '  '];
    end
end    

% More information for F can be appended to INFO:

end