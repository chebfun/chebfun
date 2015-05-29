function choreosphere
%CHOREOSPHERE   Compute spherical choreographies of the n-body problen.
%    CHOREO computes a spherical choreogpraphy using an hand-drawn initial guess. 
%   
% It uses trignometric interpolation, stereographic projection and quasi-Newton 
% methods. See [1] for more details.
%
% [1] H. Montanelli, N. I. Gushterov, Computing planar and spherical
% choreographies, submitted to Physica D.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

close all
LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';
lw = 2; ms = 30; fs = 18;
format long, format compact
n = input('How many bodies? (e.g. 5) '); 
R = input('Radius? (e.g. 2) ');
w = input('Angular velocity? (non-integer, e.g. 1.2) ');
k = 15;
N = k*n;
dom = [0 2*pi];

while 1   

% Hand-drawn initial guess:
  xmax = 2; ymax = xmax;
  clf, hold off, axis equal, axis([-xmax xmax -ymax ymax]), box on
  set(gca,FS,fs), title('Draw a curve!')
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
  A0 = actiongradevalsphere(c0,n,w,R);
  fprintf('\nInitial acion: %.6f\n',A0)
  options = struct();
  options.Method = 'lbfgs';
  options.Display = 'off';
  c = minFunc(@(x)actiongradevalsphere(x,n,w,R),c0,options);
  [A,G] = actiongradeval(c,n,w);
  fprintf('Action after optimization: %.6f\n',A)
  fprintf('Norm of the gradient: %.3e\n',norm(G))
  
% Plot the result:
  c = c(1:N) + 1i*c(N+1:2*N);
  q = chebfun(c,dom,'coeffs','trig');
  hold on, plot(q,'r',LW,lw)
  hold on, plot(q(2*pi*(0:n-1)/n),'.k',MS,ms), title('')
  s = input('Want to see the planets dancing? (Yes=1, No=Enter) '); 
  if s == 1
    hold off, plotonsphere(q,n,w,R), pause
  end
  
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,G] = actiongradevalsphere(c,n,w,R)
%ACTIONGRADEVALSPHERE   Compute the action and its gradient.

% Parse inputs:
  if nargin < 3
    w = 0;
    R = 1;
  end
  if nargin < 4
    R = 1;
  end

% Set up:
  N = length(c)/2; 
  u = c(1:N); 
  v = c(N+1:end);
  c = u + 1i*v;
  if ( mod(N, 2) == 0 )
    k = (-N/2:N/2-1)';
    kz = [0,-N/2+1:N/2-1]';
  else
    k = (-(N-1)/2:(N-1)/2)';
    kz = k;
  end
  [t,s] = trigpts(N,[0 2*pi]);

% Evaluate action:
  c = bsxfun(@times,exp(1i*k*(0:n-1)*2*pi/n),c(:,1));
  vals = ifft(N*c);
  dc = 1i*k.*c(:,1);
  if ( mod(N, 2) == 0 )
      dc(1,1) = 0;
  end
  dvals = ifft(N*dc);
  Kt = (2*R^2*abs(dvals)./(R^2 + abs(vals(:,1)).^2)).^2;
  Kr = 0;
  U = 0;
  for j = 1:n
    vj = abs(vals(:,j));
    Kr = Kr + R^2 - (-R^3+R*vj.^2).^2./(R^2+vj.^2).^2;
    for l = 1:j-1
      vi = abs(vals(:,l));
      d = 2*R*abs(vals(:,l)-vals(:,j))./(sqrt(R^2+vi.^2).*sqrt(R^2+vj.^2));
      U = U - 1/R*(1-d.^2/2)./(d.*sqrt(1-d.^2/4));
    end
  end
  A = s*(n/2*Kt + w^2/2*Kr - U);

if nargout > 1
% Evaluate gradient, A_Kr and A_U contributions:
  G = zeros(2*N, 1);
  for j = 0:n-1
    aj = bsxfun(@times,cos(k'*j*2*pi/n),cos(t*k')) - ...
      bsxfun(@times,sin(k'*j*2*pi/n),sin(t*k'));
    bj = -bsxfun(@times,cos(k'*j*2*pi/n),sin(t*k')) - ...
      bsxfun(@times,sin(k'*j*2*pi/n),cos(t*k'));
    g = (aj*u + bj*v).^2 + (-bj*u + aj*v).^2;
    f = (-R^3+R*g)./(R^2+g);
    dg = 2*(bsxfun(@times,aj*u + bj*v,aj) + bsxfun(@times,-bj*u + aj*v,-bj));
    df = bsxfun(@rdivide,2*dg*R^3,(R^2+g).^2);
    G(1:N) = G(1:N) + w^2/2*bsxfun(@times,-2*df,f)'*s';
    dg = 2*(bsxfun(@times,aj*u + bj*v,bj) + bsxfun(@times,-bj*u + aj*v,aj));
    df = bsxfun(@rdivide,2*dg*R^3,(R^2+g).^2);
    G(N+1:2*N) = G(N+1:2*N) + w^2/2*bsxfun(@times,-2*df,f)'*s';
    for l = 0:j-1
      al = bsxfun(@times,cos(k'*l*2*pi/n),cos(t*k')) - ...
        bsxfun(@times,sin(k'*l*2*pi/n),sin(t*k'));
      bl = -bsxfun(@times,cos(k'*l*2*pi/n),sin(t*k')) - ...
        bsxfun(@times,sin(k'*l*2*pi/n),cos(t*k'));
      a = al - aj;
      b = bl - bj;
      C = sqrt((a*u + b*v).^2 + (-b*u + a*v).^2);
      Bl = sqrt(R^2 + (al*u + bl*v).^2 + (-bl*u + al*v).^2);
      Bj = sqrt(R^2 + (aj*u + bj*v).^2 + (-bj*u + aj*v).^2);
      d = 2*R^2*C./(Bl.*Bj);
      dC = bsxfun(@times,a*u + b*v,a) + bsxfun(@times,-b*u + a*v,-b);
      dC = bsxfun(@rdivide,dC,C);
      dBl = bsxfun(@times,al*u + bl*v,al) + bsxfun(@times,-bl*u + al*v,-bl);
      dBl = bsxfun(@rdivide,dBl,Bl);
      dBj = bsxfun(@times,aj*u + bj*v,aj) + bsxfun(@times,-bj*u + aj*v,-bj);
      dBj = bsxfun(@rdivide,dBj,Bj);
      diffd =  2*R^2*(bsxfun(@rdivide,dC,Bl.*Bj) - ...
        bsxfun(@rdivide,dBl,Bl.^2.*Bj./C) - bsxfun(@rdivide,dBj,Bl.*Bj.^2./C));
      G(1:N) = G(1:N) + ...
          1/R*bsxfun(@rdivide,-8*R^4*diffd,d.^2.*(4*R^2-d.^2).^(3/2))'*s';
      dC = bsxfun(@times,a*u + b*v,b) + bsxfun(@times,-b*u + a*v,a);
      dC = bsxfun(@rdivide,dC,C);
      dBl = bsxfun(@times,al*u + bl*v,bl) + bsxfun(@times,-bl*u + al*v,al);
      dBl = bsxfun(@rdivide,dBl,Bl);
      dBj = bsxfun(@times,aj*u + bj*v,bj) + bsxfun(@times,-bj*u + aj*v,aj);
      dBj = bsxfun(@rdivide,dBj,Bj);
      diffd =  2*R^2*(bsxfun(@rdivide,dC,Bl.*Bj) - ...
        bsxfun(@rdivide,dBl,Bl.^2.*Bj./C) - bsxfun(@rdivide,dBj,Bl.*Bj.^2./C));
      G(N+1:2*N) = G(N+1:2*N) + ...
          1/R*bsxfun(@rdivide,-8*R^4*diffd,d.^2.*(4*R^2-d.^2).^(3/2))'*s';
    end
  end
  
% Add the A_Kt contribution:
  a = bsxfun(@times,kz',cos(t*kz'));
  b = bsxfun(@times,kz',sin(t*kz'));
  f = (a*u - b*v).^2 + (b*u + a*v).^2;
  df = 2*(bsxfun(@times,a*u - b*v,a) + bsxfun(@times,b*u + a*v,b));
  a = cos(t*k');
  b = sin(t*k');
  h = R^2 + (a*u - b*v).^2 + (b*u + a*v).^2;
  dh = 2*(bsxfun(@times,a*u - b*v,a) + bsxfun(@times,b*u + a*v,b));
  G(1:N) = G(1:N) + 4*R^4*(n/2)*(bsxfun(@rdivide,df,h.^2)'*s' - ...
    2*bsxfun(@rdivide,dh,h.^3./f)'*s');
  a = bsxfun(@times,kz',cos(t*kz'));
  b = bsxfun(@times,kz',sin(t*kz'));
  df = 2*(bsxfun(@times,b*u + a*v,a) + bsxfun(@times,a*u - b*v,-b));
  a = cos(t*k');
  b = sin(t*k');
  dh = 2*(bsxfun(@times,b*u + a*v,a) + bsxfun(@times,a*u - b*v,-b));
  G(N+1:2*N) = G(N+1:2*N) + 4*R^4*(n/2)*(bsxfun(@rdivide,df,h.^2)'*s' -...
     2*bsxfun(@rdivide,dh,h.^3./f)'*s');
end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotonsphere(q,n,w,R)
%PLOTONSPHERE   Plot a spherical choreography.

if nargin < 3
    w = 0;
    R = 1;
end
LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';
lw = 2; ms = 30; fs = 18;
XTick = [-R,0,R]; YTick = XTick; ZTick = XTick;
N = 200;
Ns = 100;
ns = 50;
[Xs,Ys,Zs] = sphere(Ns);
Xs = R*Xs; Ys = R*Ys; Zs = R*Zs;
XXs = Xs(1:ns:end,:); YYs = Ys(1:ns:end,:); ZZs = Zs(1:ns:end,:);
XXs = XXs'; YYs = YYs'; ZZs = ZZs';
T = 2*pi;
dt = .1;
tt = trigpts(N, [0 2*pi]);
z = q(tt);
[X,Y,Z] = plane2sphere(z,R);
for t = dt:dt:10*T
  clf, hold off
  surf(Xs,Ys,Zs,'EdgeColor','none'), alpha(.15)
  hold on, plot3(Xs(:,1:ns:end),Ys(:,1:ns:end),Zs(:,1:ns:end),'k',LW,1e-10)
  hold on, plot3(XXs,YYs,ZZs,'k',LW,1e-10), view(-37.5,50)
  if t > 2*T
    view(0,90)
  end
  if t < 3*T
    if w ~= 0
      M = [X';Y';Z'];
      Rot = [cos(w*dt),-sin(w*dt),0;sin(w*dt),cos(w*dt),0;0,0,1];
      M = Rot*M;
      X = M(1,:)'; Y = M(2,:)'; Z = M(3,:)';
    end
  hold on, plot3(X,Y,Z,'r',LW,lw)
  end
  for j=1:n
    z = q(t+2*pi*j/n);
    [XX,YY,ZZ] = plane2sphere(z,R);
    if w ~= 0
      M = [XX';YY';ZZ'];
      Rot = [cos(w*t),-sin(w*t),0;sin(w*t),cos(w*t),0;0,0,1];
      M = Rot*M;
      XX = M(1,:)'; YY = M(2,:)'; ZZ = M(3,:)';
    end
    hold on, plot3(XX,YY,ZZ,'.',MS,ms)
  end
  box on, set(gca,FS,fs), set(gca,'XTick',XTick,'YTick',YTick,'ZTick',ZTick)
  axis([-R R -R R -R R]), axis equal
  drawnow
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x1,x2,x3] = plane2sphere(z,R)
%PLANE2SPHERE   Inverse stereographic projection.

  if nargin < 2
    R = 1;
  end
  x1 = 2*R^2*real(z)./(R^2 + abs(z).^2);  
  x2 = 2*R^2*imag(z)./(R^2 + abs(z).^2);
  x3 = R*(-R^2 + abs(z).^2)./(R^2 + abs(z).^2); 
  
end

