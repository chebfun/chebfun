function varargout = arrowplot(f,varargin)
%ARROWPLOT   Chebfun plot with arrowhead at end
%   ARROWPLOT(F, G), where F and G are CHEBFUNs with the same domain, plots the
%   curve (F,G) in the plane with an arrowhead.
%
%   ARROWPLOT(F, G, 'multi', n) plots the curve (F,G) in the plane with n
%   arrowheads placed on the curve.
%
%   ARROWPLOT(F), where F is a complex CHEBFUN, plots the curve
%   (real(F),imag(F)) in the plane with an arrowhead.
%
%   Arguments can also be quasimatrices, in which case several curves are
%   plotted.
%
%   Plotting options can be passed in the usual fashion. For example,
%   'markerSize' will change the size of the arrowhead (the default is 6).
%   Furthermore, there is an option 'ystretch' for manually adjusting the
%   slope of the arrowhead, so that the slope is 'ystretch' times greater
%   than the computed slope dictates. This is useful for plots which are 
%   rescaled, or have their axes adjusted after plotting.
%
% Examples:
%
%   t = chebfun('t',[0,6]);
%   f = sin(t); g = cos(t), arrowplot(f,g);
%
%   h = exp((-.2+3i)*t); arrowplot(h,'color','r')
%
%   A = []; for k = 1:3, A = [A exp((-.1*k+1i)*t)]; end, arrowplot(A)
%
%   A = []; for k = 1:3, A = [A exp((-.2*k+1i)*t)]; end, arrowplot(A,'multi',5)

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse arguments. If first input argument is a real CHEBFUN, next parameter
% must a CHEBFUN as well. If first input parameter is a complex CHEBFUN, rest of
% arguments must be options for the plot.
if isreal(f)
    assert(nargin > 1, 'CHEBFUN:CHEBFUN:arrowplot:nargin', ...
        ['Arrowplot expects two real valued CHEBFUNs, or a complex CHEBFUN '...
        'as arguments'])
    g = varargin{1};
    varargin(1) = [];
    f = f + 1i*g;
end

% Did the user pass in ystretch argument?
findStretch = find(strcmp(varargin,'ystretch'));
if ~isempty(findStretch)
    ystretch = varargin{findStretch + 1};
    varargin(findStretch:findStretch+1) = [];
else
    ystretch = 1;
end

% Did the user pass a markersize argument? This needs to be treated slightly
% differently, as it's not officially an option for Matlab annotations, so we
% work around it by grabbing the value to pass to 'HeadLength' and 'HeadWidth'/
findMS = find(strcmpi(varargin,'markersize'));
if ~isempty(findMS)
    markerSize = varargin{findMS + 1};
    varargin(findMS:findMS+1) = [];
else
    markerSize = 6;
end

% Did the user pass a multi-head argument? Find how many arrowheads have been
% requested.
findMulti = find(strcmpi(varargin,'multi'));
if ~isempty(findMulti)
    numArrows = varargin{findMulti + 1};
    varargin(findMulti:findMulti+1) = [];
else
    numArrows = 1;
end

% Evaluate the derivative of input CHEBFUNs to get slope information
fp = diff(f);

% Get domain information, and evaluate functions and derivatives at the correct
% point(s):
fdom = domain(f);

% Points to evaluate functions at for arrow location/slopes
evalPts = linspace(fdom(1), fdom(end), numArrows + 1);
evalPts = evalPts(2:end);

% Eval at the points
fPts = feval(f, evalPts);
fpPts = feval(fp, evalPts);

% Normalise the slopes
for aCounter=1:length(fpPts)
    if norm(fpPts(aCounter),2) ~= 0
        fpPts(aCounter) = 0.001*fpPts(aCounter)/norm(fpPts(aCounter),2);
    end
end

% Plot the complex valued CHEBFUNs. Note that we need treat real valued chebfuns
% in a special way. The only way we end up with a real f is if the g input from
% above was 0 (since then, f = f + 1i*g from above just gives the input f). In
% that case, plot(f) won't call the desired complex chebfun plot. See #2106.
if ( isreal(f) )
    pp = plot(f, 0*f, varargin{:});
else
    pp = plot(f,varargin{:});
end

% Loop through the pieces, and plot arrows at the points requested. Don't plot
% any arrows in the case where the CHEBFUN column is a zero chebfun:
isz = iszero(f);

% Loop through the points
for aCounter=1:length(fpPts)
    h(aCounter) = annotation('arrow','Visible','off');
    if ( ~isz(ceil(aCounter/numArrows)) )
        set(h(aCounter),'parent', gca, ...
            'position', [real(fPts(aCounter)) imag(fPts(aCounter)) ...
            real(fpPts(aCounter)) ystretch*imag(fpPts(aCounter))], ...
            'HeadLength', markerSize, 'HeadWidth', markerSize, ...
            'Color', pp(ceil(aCounter/numArrows)).Color, ...
            'Visible', 'on', varargin{:});
    end
end

% Return a handle to the graphic objects if method called with an output. Note
% that we need to concatenate the lines in PP and arrows in H.
if ( nargout == 1)
    varargout{1} = [pp(:); h(:)];
end

end
