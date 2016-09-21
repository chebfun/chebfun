function choreosphere
%CHOREOSPHERE   Compute spherical choreographies of the curved n-body problem.
%    CHEB.CHOREOSPHERE computes a spherical choreography using hand-drawn
%    initial guesses. At the end of the computation, the user is asked to press
%    <1> to display the motion of the planets. To stop the program, press
%    <CTRL>-<C>.
%
% Spherical choreographies are periodic solutions of the n-body problem on the
% sphere in which the bodies share a single orbit. This orbit can be fixed or
% rotating with some angular velocity relative to an inertial reference frame.
%
% The algorithm uses trigonometric interpolation, stereographic projection and
% quasi-Newton methods. See [1] for details.
%
% Example: at the prompt, specify 5 bodies, angular rotation 0, radius 2.
% Then either draw a curve (if your machine has imfreehand) or click
% in 10 or 15 points (if it doesn't) roughly along a figure-8.  Type
% <Enter>, and in a few seconds you will see a choreograpy.
%
% [1] H. Montanelli and N. I. Gushterov, Computing planar and spherical
% choreographies, SIAM Journal on Applied Dynamical Systems 15 (2016),

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

close all
LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';
lw = 2; ms = 30; fs = 18;
format long, format compact
n = input('How many bodies? (e.g. 5) ');
R = input('Radius? (e.g. 2) ');
w = input('Angular velocity? (e.g. 0 or 1.2) ');
k = 15;
N = k*n;
dom = [0 2*pi];

while 1
    
    % Hand-drawn initial guess:
    xmax = 2; ymax = xmax;
    clf, hold off, axis equal, axis([-xmax xmax -ymax ymax]), box on
    set(gca,FS,fs), title('Draw a curve!')
    if ( exist('imfreehand') == 2 )
        % Try to use IMFREEHAND (need the image processing toolbox):
        h = imfreehand();
        z = getPosition(h);
        delete(h)
        z = z(:,1) + 1i*z(:,2);
        q0 = chebfun(z,dom,'trig');
        q0 = chebfun(q0,dom,N,'trig');
        c0 = trigcoeffs(q0);
        c0(1+floor(N/2)) = 0;
        q0 = chebfun(c0,dom,'coeffs','trig');
    else
        % Otherwise, use GINPUT:
        x = []; y = []; button = 1;
        disp('Input points with mouse, press <enter> for final point.')
        while ( button == 1 )
            [xx,yy,button] = ginput(1);
            x = [x; xx]; y = [y; yy]; %#ok<*AGROW>
            hold on, plot(xx,yy,'xb',MS,ms/2), drawnow
        end
        z = x + 1i*y;
        q0 = chebfun(z,dom,'trig');
        q0 = chebfun(q0,dom,N,'trig');
    end
    
    % Solve the problem:
    hold on, plot(q0,'.-b',LW,lw), drawnow
    c0 = trigcoeffs(q0);
    c0 = [real(c0);imag(c0)];
    [A0,G0] = actiongradevalsphere(c0,n,w,R);
    fprintf('\nInitial action: %.6f\n',A0)
    options = optimoptions('fminunc');
    options.Algorithm = 'quasi-newton';
    options.HessUpdate = 'bfgs';
    options.GradObj = 'on';
    options.Display = 'off';
    [c,A,~, ~,G] = fminunc(@(x)actiongradevalsphere(x,n,w,R),c0,options);
    
    % Plot the result:
    fprintf('Action after optimization: %.6f\n',A)
    fprintf('Norm of the gradient: %.3e\n',norm(G)/norm(G0))
    c = c(1:N) + 1i*c(N+1:2*N);
    q = chebfun(c,dom,'coeffs','trig');
    hold on, plot(q,'r',LW,lw)
    hold on, plot(q(2*pi*(0:n-1)/n),'.k',MS,ms), title('')
    s = input('Want to see the planets dancing? (Yes=<1>, No=<0>) ');
    if ( s == 1 )
        hold off, plotonsphere(q,n,w,R), pause
    end
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A, G] = actiongradevalsphere(c,n,w,R)
%ACTIONGRADEVALSPHERE  Compute the action and its gradient on the sphere.

% Set up:
if ( nargin < 3 )
    w = 0;
    R = 2;
end
if ( nargin < 4 )
    R = 2;
end
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

% Evaluate action A:
c = bsxfun(@times,exp(1i*k*(0:n-1)*2*pi/n),c(:,1));
vals = ifft(N*c);
dc = 1i*k.*c(:,1);
if ( mod(N, 2) == 0 )
    dc(1,1) = 0;
end
dvals = ifft(N*dc);
v1 = abs(vals(:,1));
U = 0;
for j = 2:n
    vj = abs(vals(:,j));
    D = 2*R^2*abs(vals(:,1)-vals(:,j))./(sqrt(R^2+v1.^2).*sqrt(R^2+vj.^2));
    U = U - 1/R*(2*R^2-D.^2)./(D.*sqrt(4*R^2-D.^2));
end
K = (2*R^2*abs(dvals+1i*w*vals(:,1))./(R^2+v1.^2)).^2;
A = s*n/2*(K-U);

if ( nargout > 1 )
    % Initialize gradient G:
    G = zeros(2*N,1);
    
    % Loop over the bodies for the AU contribution:
    a0 = cos(t*k');
    b0 = -sin(t*k');
    r0 = sqrt(R^2 + (a0*u + b0*v).^2 + (-b0*u + a0*v).^2);
    dr0du = bsxfun(@times,a0*u + b0*v,a0) + bsxfun(@times,-b0*u + a0*v,-b0);
    dr0du = bsxfun(@rdivide,dr0du,r0);
    dr0dv = bsxfun(@times,a0*u + b0*v,b0) + bsxfun(@times,-b0*u + a0*v,a0);
    dr0dv = bsxfun(@rdivide,dr0dv,r0);
    for j = 1:n-1
        c = bsxfun(@times,cos(k'*j*2*pi/n),cos(t*k')) - ...
            bsxfun(@times,sin(k'*j*2*pi/n),sin(t*k'));
        d = -bsxfun(@times,cos(k'*j*2*pi/n),sin(t*k')) - ...
            bsxfun(@times,sin(k'*j*2*pi/n),cos(t*k'));
        r = sqrt(R^2 + (c*u + d*v).^2 + (-d*u + d*v).^2);
        a = a0 - c;
        b = b0 - d;
        f = (a*u + b*v).^2 + (-b*u + a*v).^2;
        D = 2*R^2*sqrt(f)./(r0.*r);
        dr = bsxfun(@times,c*u + d*v,c) + bsxfun(@times,-d*u + c*v,-d);
        dr = bsxfun(@rdivide,dr,r);
        df = 2*(bsxfun(@times,a*u + b*v,a) + bsxfun(@times,-b*u + a*v,-b));
        dD = 2*R^2*(bsxfun(@rdivide,df,r0.*r.*(2*sqrt(f))) - ...
            bsxfun(@rdivide,dr0du,r0.^2.*r./sqrt(f)) - ...
            bsxfun(@rdivide,dr,r0.*r.^2./sqrt(f)));
        G(1:N) = G(1:N) + ...
            n/(2*R)*bsxfun(@rdivide,-8*R^4*dD,D.^2.*(4*R^2-D.^2).^(3/2))'*s';
        dr = bsxfun(@times,c*u + d*v,d) + bsxfun(@times,-d*u + c*v,c);
        dr = bsxfun(@rdivide,dr,r);
        df = 2*(bsxfun(@times,a*u + b*v,b) + bsxfun(@times,-b*u + a*v,a));
        dD = 2*R^2*(bsxfun(@rdivide,df,r0.*r.*(2*sqrt(f))) - ...
            bsxfun(@rdivide,dr0dv,r0.^2.*r./sqrt(f)) - ...
            bsxfun(@rdivide,dr,r0.*r.^2./sqrt(f)));
        G(N+1:2*N) = G(N+1:2*N) + ...
            n/(2*R)*bsxfun(@rdivide,-8*R^4*dD,D.^2.*(4*R^2-D.^2).^(3/2))'*s';
    end
    
    % Add the AK contribution:
    g1 = bsxfun(@times,kz'+w,cos(t*k'));
    g2 = bsxfun(@times,kz'+w,sin(t*k'));
    g = (g1*u - g2*v).^2 + (g2*u + g1*v).^2;
    h = (R^2 + (a0*u + b0*v).^2 + (-b0*u + a0*v).^2).^2;
    dg = 2*(bsxfun(@times,g1*u - g2*v,g1) + bsxfun(@times,g2*u + g1*v,g2));
    dh = (bsxfun(@times,a0*u + b0*v,a0) + bsxfun(@times,-b0*u + a0*v,-b0));
    dh = bsxfun(@times,dh,4*sqrt(h));
    G(1:N) = G(1:N) + 2*R^4*n*(bsxfun(@rdivide,dg,h)'*s' - ...
        bsxfun(@rdivide,dh,h.^2./g)'*s');
    dg = 2*(bsxfun(@times,g2*u + g1*v,g1) + bsxfun(@times,g1*u - g2*v,-g2));
    dh = (bsxfun(@times,a0*u + b0*v,b0) + bsxfun(@times,-b0*u + a0*v,a0));
    dh = bsxfun(@times,dh,4*sqrt(h));
    G(N+1:2*N) = G(N+1:2*N) + 2*R^4*n*(bsxfun(@rdivide,dg,h)'*s' - ...
        bsxfun(@rdivide,dh,h.^2./g)'*s');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotonsphere(q,n,w,R)
%PLOTONPLANE   Plot a spherical choreography.

if ( nargin < 3 )
    w = 0;
    R = 1;
end
LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';
lw = 2; ms = 30; fs = 18;
XTick = [-R,0,R]; YTick = XTick; ZTick = XTick;
N = 500;
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
    hold on, plot3(XXs,YYs,ZZs,'k',LW,1e-10), view(-45.5,46)
    if ( t > 5*T )
        view(0,90)
    end
    if ( t < 3*T )
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

if ( nargin < 2 )
    R = 1;
end
x1 = 2*R^2*real(z)./(R^2 + abs(z).^2);
x2 = 2*R^2*imag(z)./(R^2 + abs(z).^2);
x3 = R*(-R^2 + abs(z).^2)./(R^2 + abs(z).^2);

end
