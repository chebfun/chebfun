function varargout = displayInfo(mode, varargin)
%DISPLAYINFO    Utility routine for displaying iteration progress.
% Calling sequence:
%   [DISPLAYFIG, DISPLAYTIME] = DISPLAYINFO(MODE, VARARGIN)
% where
%   MODE:       What phase in Newton iteration are we in?
%   VARARGIN:   Various useful information passed from the Newton solver
%               to this method.
% and
%   DISPLAYFIG:     Handle to the figure plots are drawn.
%   DISPLAYTIME:    Handle to the Matlab timer used for pausing the desired time
%                   between plots (cf. t = tic).
%
% See also: chebop/displayInfoFinal, chebop/displayInfoInit, 
%           chebop/displayInfoIter, chebop/displayInfoLinear.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developers note: This method was introduced to allow the nonlinear solvers to
% use the same syntax for displaying information, whether the system is running
% in GUI mode or command-line mode.
%
% This method will call one of three different methods listed below, depending
% on where in the solution process the iteration currently is (start of
% iteration, during iteration, or after iteration finishes).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Display information depending on what phase in the Newton iteration we're in.
switch ( mode )
    
    % If we get passed an initial guess that satisfies the BCs.
    case 'exactInitial'
        chebop.displayInfoExactInitial(varargin{:});
        
    % Start of iteration
    case 'init'
        [displayFig, displayTimer] = chebop.displayInfoInit(varargin{:});
        varargout{1} = displayFig;
        varargout{2} = displayTimer;
        
    % During Newton iteration
    case 'iter'
        varargout{1} = chebop.displayInfoIter(varargin{:});
        varargout{2} = false;

    % Once iteration is over
    case 'final'
        chebop.displayInfoFinal(varargin{:});
        
    % Display special information in case of linear problems
    case 'linear'
        chebop.displayInfoLinear(varargin{:});
    

end

end

