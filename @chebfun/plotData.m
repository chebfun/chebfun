function data = plotData(f)
%PLOTDATA   Useful data values for plotting a CHEBFUN object.
%   OUT = PLOTDATA(F) returns a struct containing data that can be used for
%   plotting F. In particular, DATA.xLine and DATA.fLine are for plotting smooth
%   curves (usually passed to plot with '-'), DATA.xPoints and DATA.yPoints
%   contain the (x, F(x)) data used to represent F, and DATA.xJumps and
%   DATA.fJumps are the linear connections between discontinuous pieces.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the data from the FUN objects:
data = struct('xLine', [], 'fLine', [], 'xPoints', [], 'fPoints', [], ...
    'xJumps', [], 'fJumps', []);
% Loop over each FUN for Line and Points data:
for k = 1:numel(f.funs)
    dataNew = plotData(f.funs{k});
    myNaN = NaN(1, size(dataNew.fLine, 2));
    data.xLine = [data.xLine ; NaN ; dataNew.xLine];
    data.fLine = [data.fLine ; myNaN ; dataNew.fLine];
    data.xPoints = [data.xPoints ; NaN ; dataNew.xPoints];
    data.fPoints = [data.fPoints ; myNaN ; dataNew.fPoints];
end

% Loop over each FUN for Jumps data:
for k = 1:(numel(f.funs) - 1)
    data.xJumps = [data.xJumps, NaN, f.funs{k}.domain(2) f.funs{k+1}.domain(1)];
    data.fJumps = [data.fJumps, myNaN', get(f.funs{k}, 'rval')', ...
        get(f.funs{k+1}, 'lval')'];
end

end
