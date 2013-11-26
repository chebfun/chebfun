function [displayFig, displayTimer] = displayInfo(mode, varargin)
% Utility routine for displaying iteration progress in the solve functions.

switch mode
    case 'init'
        [displayFig, displayTimer] = chebop.displayInfoInit(varargin{:});
    case 'iter'
        chebop.displayInfoIter(varargin{:});
    case 'final'
        chebop.displayInfoFinal(varargin{:});
end

end

