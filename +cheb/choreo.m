function choreo
%CHOREO   Compute planar choreographies of the n-body problen.
%    CHOREO computes a planar choreogpraphy using an hand-drawn initial guess. 
%   
% It uses trignometric interpolation and quasi-Newton methods. See [1] for more
% details.
%
% [1] H. Montanelli, N. I. Gushterov, Computing planar and spherical
% choreographies, submitted to Physica D.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

close all
LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';
lw = 2; ms = 35; fs = 18;
format long, format compact
n = input('How many bodies? (e.g. 5) '); 
w = input('Angular velocity? (non-integer, e.g. 1.2) ');
k = 15;
N = k*n;
dom = [0 2*pi];

while 1   

% Hand-drawn initial guess:
  xmax = 2.5; ymax = 2.5;
  clf, hold off, axis equal, axis([-xmax xmax -ymax ymax]), box on
  set(gca,'FontSize',fs), title('Draw a curve!')
  h = imfreehand;
  z = getPosition(h);
  delete(h)
  z = z(:,1) + 1i*z(:,2);
  q0 = chebfun(z,dom,'trig');
  q0 = chebfun(q0,dom,N,'trig');
  c0 = trigcoeffs(q0);
  c0(1+floor(N/2)) = 0;
  q0 = chebfun(c0,dom,'coeffs','trig');
    
% Solve the problem:
  hold on, plot(q0,'.-b',LW,lw), drawnow
  c0 = trigcoeffs(q0);
  c0 = [real(c0);imag(c0)];
  A0 = actiongradeval(c0,n,w);
  fprintf('\nInitial acion: %.6f\n',A0)
  options = optimoptions('fminunc');
  options.Algorithm = 'quasi-newton';
  options.HessUpdate = 'bfgs'; 
  options.GradObj = 'on';
  options.Display = 'off';
  [c,A,~, ~,G] = fminunc(@(x)actiongradeval(x,n,w),c0,options);
  fprintf('Action after optimization: %.6f\n',A)
  fprintf('Norm of the gradient: %.3e\n',norm(G))

% Plot the result:
  c = c(1:N) + 1i*c(N+1:2*N);
  c(1+floor(N/2)) = 0;
  q = chebfun(c,dom,'coeffs','trig');
  hold on, plot(q,'r',LW,lw)
  hold on, plot(q(2*pi*(0:n-1)/n),'.k',MS,ms), title('')
  s = input('Want to see the planets dancing? (Yes=1, No=0) '); 
  if s == 1
    hold off, plotonplane(q,n,w), pause
  end
  
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,G] = actiongradeval(c,n,w)
%ACTIONGRADEVAL   Compute the action and its gradient.

% Set up:
  N = length(c)/2; 
  u = c(1:N); 
  v = c(N+1:end);
  c = u + 1i*v;
  if ( mod(N, 2) == 0 )
    k = (-N/2:N/2-1)';
  else
    k = (-(N-1)/2:(N-1)/2)';
  end
  [t,s] = trigpts(N,[0 2*pi]);
  if nargin < 3
      w = 0;
  end

% Evaluate action A:
  c = bsxfun(@times,exp(1i*k*(0:n-1)*2*pi/n),c(:,1));
  vals = ifft(N*c);
  dc = 1i*k.*c(:,1);
  if ( mod(N, 2) == 0 )
    dc(1,1) = 0; % zero the k=-N/2 coefficient
  end
  dv = ifft(N*dc);
  K = abs(dv + 1i*w*vals(:,1)).^2; 
  U = 0;
  for i = 2:n
    U = U + 1./abs(vals(:,1)-vals(:,i)); 
  end
  A = n/2*s*(K + U); % trapezoidal rule

if nargout > 1
% Evaluate gradient G:
  G = zeros(2*N, 1);
  % Loop over the bodies for the A_U contribution:
  for j = 1:n-1
    a = bsxfun(@times,1-cos(k'*j*2*pi/n),cos(t*k')) + ...
        bsxfun(@times,sin(k'*j*2*pi/n),sin(t*k'));
    b = bsxfun(@times,-1+cos(k'*j*2*pi/n),sin(t*k')) + ...
        bsxfun(@times,sin(k'*j*2*pi/n),cos(t*k'));
    f = sqrt((a*u + b*v).^2 + (-b*u + a*v).^2);
    df = bsxfun(@times,a*u + b*v,a) + bsxfun(@times,b*u - a*v,b);
    G(1:N) = G(1:N) - n/2*bsxfun(@rdivide,df,f.^3)'*s';
    df = bsxfun(@times, a*u + b*v, b) + bsxfun(@times,-b*u + a*v,a);
    G(N+1:2*N) = G(N+1:2*N) - n/2*bsxfun(@rdivide,df,f.^3)'*s';
  end
  % Add the A_K contribution:
  if ( mod(N, 2) == 0 )
    G = G + 2*pi*n*[w;(k(2:end)+w);w;(k(2:end)+w)].^2.*[u;v];
  else
    G = G + 2*pi*n*[(k+w);(k+w)].^2.*[u;v];
  end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotonplane(q,n,w)
%PLOTONPLANE   Plot a planar choreography.

xmax = 2.5; ymax = 2.5;
LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';
lw = 2; ms = 35; fs = 18;
dom = [0 2*pi]; T = 2*pi; dt = .1;
xStars = 2*xmax*rand(250,1)-xmax;
yStars = 2*ymax*rand(250,1)-ymax;
for t = dt:dt:10*T
    clf, hold off
    fill(xmax*[-1 1 1 -1 -1],ymax*[-1 -1 1 1 -1],'k')
    axis equal, axis([-xmax xmax -ymax ymax]), box on
    set(gca,FS,fs)
    hold on, plot(xStars,yStars,'.w',MS,4)
    if w ~= 0
        q = chebfun(@(t)exp(1i*w*dt).*q(t),dom,length(q),'trig');
    end
    if t < T
        hold on, plot(q,'r',LW,lw)
    end
    for j = 1:n
        hold on, plot(q(t+2*pi*j/n),'.',MS,ms)
    end
    drawnow
end

end
