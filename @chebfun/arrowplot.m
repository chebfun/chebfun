function varargout = arrowplot(f,varargin)
%ARROWPLOT   Chebfun plot with arrowhead at end
%   ARROWPLOT(F,G), where F and G are CHEBFUNs with the same domain, plots the
%   curve (F,G) in the plane with an arrowhead.
%
%   ARROWPLOT(F), where F is a complex CHEBFUN, plots the curve
%   (real(F),imag(F)) in the plane with an arrowhead.
%
%   Arguments can also be quasimatrices, in which case several curves are
%   plotted.
%
%   Plotting options can be passed in the usual fashion. For example,
%   'markerSize' will change the size of the arrowhead plotted. Furthermore,
%   there is an option 'ystretch', for manually adjusting the slope of the
%   arrowhead, so that the slope is 'ystretch' times greater than the computed
%   slope dictates. This is useful for plots which are rescaled, or have their
%   axes adjusted after plotting.
%
% Examples:
%
%   t = chebfun('t',[0,6]);
%   f = sin(t); g = cos(t), arrowplot(f,g);
%
%   h = exp((-.2+3i)*t); arrowplot(h,'color','r')
%
%   A = []; for k = 1:3, A = [A exp((-.1*k+1i)*t)]; end, arrowplot(A)

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse arguments. If first input argument is a real CHEBFUN, next parameter
% must a CHEBFUN as well. If first input parameter is a complex CHEBFUN, rest of
% arguments must be options for the plot.
if isreal(f)
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

% Did the user pass in markersize argument? This needs to be treated slightly
% differently, as it's not officially an option for Matlab annotations, so we
% workaround it by grabbing the value to pass to 'HeadLength' and 'HeadWidth'/
findMS = find(strcmpi(varargin,'markersize'));
if ~isempty(findMS)
    markerSize = varargin{findMS + 1};
    varargin(findMS:findMS+1) = [];
else
    markerSize = 7;
end


% Evaluate the derivative of input CHEBFUNs to get slope information
fp = diff(f);

% Get domain information, and evaluate functions and derivatives at the correct
% point:
fdom = domain(f);
fend = feval(f, fdom(end));
fpend = feval(fp, fdom(end));

% Normalise the slopes
for aCounter=1:length(fpend)
    fpend(aCounter) = 0.001*fpend(aCounter)/norm(fpend(aCounter),2);
end

% Plot the complex valued CHEBFUNs
pp = plot(f,varargin{:});

% Loop through the pieces, and plot arrows at the end
for aCounter=1:length(fend)
    % Create a vector of arrow annotations, then set their properties below
    h(aCounter) = annotation('arrow');
    set(h(aCounter),'parent', gca, ...
        'position', [real(fend(aCounter)) imag(fend(aCounter)) ...
            real(fpend(aCounter)) ystretch*imag(fpend(aCounter))], ...
        'HeadLength', markerSize, 'HeadWidth', markerSize, ...
        'Color', pp(aCounter).Color, varargin{:});
end

% Return a handle to the graphic objects if method called with an output. Note
% that we need to concat the lines in PP and arrows in H.
if ( nargout == 1)
    varargout{1} = [pp(:); h(:)];
end

end