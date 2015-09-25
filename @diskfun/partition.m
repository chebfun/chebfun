function [fp,fm] = partition(f)
% PARTITION  Partition a diskfun into its even/periodic odd/anti-periodic
% parts.
%
% [fp,fm] = partition(f) partitions f into two diskfuns fp & fm with the
% following properties:
% 
%   fp has a CDR decomposition such that C is even and R is pi periodic
%   fm has a CDR decomposition such that C is odd and R is pi anti-periodic
%
% See also COMBINE

if ~isa(f,'diskfun')
    error('DISKFUN:partition:unknown',['Undefined function ''partition'' for ' ...
        'input argument of type %s.'], class(f));
end

if isempty(f)
    fp = diskfun();
    fm = diskfun();
    return
end

% Do the even-pi-periodic case first
id = f.idxPlus;

if isempty(id)
    fp = diskfun();
else
    fp = f;
    fp.cols = fp.cols(:,id);
    fp.rows = fp.rows(:,id);
    fp.pivotValues = fp.pivotValues(id, :);
    fp.pivotLocations = fp.pivotLocations(id, :); %fix locations and indices (s.b. two columns)
    fp.pivotIndices = fp.pivotIndices(id, :);      % not sure why fixing this causes other issues
    fp.idxPlus = 1:length(id);
    fp.idxMinus = [];
end

% Now odd case
id = f.idxMinus;

if isempty(id)
    fm = diskfun();
else
    fm = f;
    fm.cols = fm.cols(:,id);
    fm.rows = fm.rows(:,id);
    fm.pivotValues = fm.pivotValues(id);
    fm.pivotLocations = fm.pivotLocations(id, :);
    fm.pivotIndices = fm.pivotIndices(id, :);
    fm.idxMinus = 1:length(id);
    fm.idxPlus = [];
    fm.nonZeroPoles = 0;
end

end

