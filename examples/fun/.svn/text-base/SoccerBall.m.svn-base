%% Chebfun Soccerball
% Filomena Di Tommaso,  24th July 2013

%%
% (Chebfun example fun/SoccerBall.m)
% [Tags: #Chebfun2,  #soccer]

%%
% Our aim is to draw a soccerball.
%
% The vertices of the faces of the soccerball are the vertices of a truncated
% icosahedron centered at the origin and they are all the even permutations 
% $$
% \begin{array}{c}
% \left( 0,\pm 1,\pm 3\phi \right) \\ 
% \left( \pm 2,\pm \left( 1+2\phi \right) ,\pm \phi \right) \\ 
% \left( \pm 1,\pm \left( 2+\phi \right) ,\pm 2\phi \right)
% \end{array}
% $$
% where $\phi =\frac{1+\sqrt{5}}{2}$ is the golden mean [1]. Moreover, it is
% easy to determine which vertices of the list identify the exagonal and
% pentagonal faces of the soccerball. (Matlab vertices\_faces.m)

%%
% Once we identify the faces it is necessary to draw the edges of the faces by
% distinguishing three cases:
% 
%  # edges belonging from parallels
%  # edges belonging from meridians
%  # slanted edges.

%%
% Let us observe that two points are on the same parallel if they have the
% same $z$ coordinate and they are on a meridian if they are alligned on the $%
% y $ direction. Finally, we need to consider separately the case of North and
% South pole.
% 
% At that point we have all the faces of our soccerball and we want to color
% in black the pentagonal one. In order to do this we consider a grid of
% points for each pentagonal face (actually we consider half of the pentagonal
% faces since the other half are in antipodal position) and we establish if
% each point of the grid is outside or inside the pentagon. The procedure to
% establish the membership of a point $p$ in a pentagon is the following: we
% subdivide the pentagon in three triangle and we calculate the spherical
% barycentric coordinates [2] of $p$ with respect to each of the triangle, $p$
% lies inside the pentagon if all the spherical barycentric coordinates are
% positive.

LW = 'linewidth'; FS = 'fontsize'; MS = 'markersize';

% Define some constants
alpha = (1+sqrt(5))/2;
d = [-2*pi, 2*pi, -2*pi, 2*pi];
r = sqrt(9*alpha+10);
theta = chebfun2(@(t, p) t, d);
phi = chebfun2(@(t, p) p, d);

% Parametric sphere
X = r*cos(theta).*cos(phi);
Y = r*cos(theta).*sin(phi);
Z = r*sin(theta);

% Plot of the surface
colormap([1 1 1])
surf(X,  Y,  Z)
set(gca, 'color', .8*[1 1 1])
axis equal
hold on

%% Vertices and Faces

% Vertices of the soccerball
V = [0, 1, 3*alpha;
    0, 1, -3*alpha;
    0, -1, 3*alpha;
    0, -1, -3*alpha;
    2, 1+2*alpha, alpha;
    2, 1+2*alpha, -alpha;
    2, -(1+2*alpha), alpha;
    2, -(1+2*alpha), -alpha;
    -2, 1+2*alpha, alpha;
    -2, 1+2*alpha, -alpha;
    -2, -(1+2*alpha), alpha;
    -2, -(1+2*alpha), -alpha;
    1, 2+alpha, 2*alpha;
    1, 2+alpha, -2*alpha;
    1, -(2+alpha), 2*alpha;
    1, -(2+alpha), -2*alpha;
    -1, 2+alpha, 2*alpha;
    -1, 2+alpha, -2*alpha;
    -1, -(2+alpha), 2*alpha;
    -1, -(2+alpha), -2*alpha]; 
V = [V; V(:, 2), V(:, 3), V(:, 1);V(:, 3), V(:, 1), V(:, 2)];

% Hexagonal faces
Fexa = [1, 3, 50, 58, 54, 46;
    1, 3, 49, 57, 53, 45;
    13, 5, 33, 25, 53, 45;
    49, 15, 7, 34, 26, 57;
    46, 17, 9, 35, 27, 54;
    50, 19, 11, 36, 28, 58;
    11, 19, 15, 7, 22, 24;
    22, 8, 16, 20, 12, 24;
    20, 52, 60, 32, 40, 12;
    16, 8, 38, 30, 59, 51;
    9, 17, 13, 5, 21, 23;
    10, 18, 14, 6, 21, 23;
    10, 39, 31, 56, 48, 18;
    55, 47, 14, 6, 37, 29;
    41, 25, 33, 37, 29, 43;
    41, 26, 34, 38, 30, 43;
    42, 27, 35, 39, 31, 44;
    42, 28, 36, 40, 32, 44;
    2, 4, 51, 59, 55, 47;
    2, 4, 52, 60, 56, 48];

% pentagonal faces
Fpenta = [1, 45, 13, 17, 46;
    15, 49, 3, 50, 19;
    7, 22, 8, 38, 34;
    26, 41, 25, 53, 57;
    37, 6, 21, 5, 33;
    54, 27, 42, 28, 58;
    11, 24, 12, 40, 36;
    35, 39, 10, 23, 9;
    16, 51, 4, 52, 20;
    18, 2, 14, 47, 48;
    44, 31, 32, 56, 60;
    43, 29, 30, 55, 59];

%%

for i = 1:20
    for k = 1:6
        P = V(Fexa(i, k), :);
        Q = V(Fexa(i, mod(k, 6)+1), :);
        
        PQ = sqrt((P(1)-Q(1))^2+(P(2)-Q(2))^2+(P(3)-Q(3))^2);
        rc = sqrt(P(1)^2+P(2)^2+P(3)^2);
        
        % coefficient of the plane through P, Q and O
        A = P(2)*Q(3)-Q(2)*P(3);
        B = Q(1)*P(3)-P(1)*Q(3);
        C = P(1)*Q(2)-Q(1)*P(2);
        
        % angles which identify point P
        beta1 = asin(P(3)/rc);
        alpha1 = asin(P(2)/(rc*cos(beta1)));
        
        % angles which identify point Q
        beta2 = asin(Q(3)/rc);
        alpha2 = asin(Q(2)/(rc*cos(beta2)));
        
        
        if abs(sum([rc*cos(beta1).*cos(alpha1), rc*cos(beta1).*sin(alpha1), rc*sin(beta1)]-P))>10^-5
            beta1 = pi-beta1;
            alpha1 = asin(P(2)/(rc*cos(beta1)));
        end
        
        if abs(sum([rc*cos(beta2).*cos(alpha2), rc*cos(beta2).*sin(alpha2), rc*sin(beta2)]-Q))>10^-5
            beta2 = pi-beta2;
            alpha2 = asin(Q(2)/(rc*cos(beta2)));
        end
        
        % Point for which abs(alpha1-alpha2)<10^-5 are point on the same meridian
        if P(3) ~= Q(3) && abs(alpha1-alpha2)>10^-5  
            % Hexagon at the north and south pole
            if (i == 1 && k == 2)||(i == 20 && k == 2)
                alpha1 = -alpha1;
            end
            if (i == 1 && k == 6) ||(i == 20 && k == 6)
                alpha2 = -alpha2;
            end
            
            v = chebfun(@(v) v, [min(alpha1, alpha2), max(alpha1, alpha2)]);
            t = atan((-A*cos(v)-B*sin(v))/C);
            
            if abs(t(alpha1)-beta1)>10^-5
                t = atan((-A*cos(v)-B*sin(v))/C)+pi;
            end
            if abs(t(alpha2)-beta2)>10^-5
                t = atan((-A*cos(v)-B*sin(v))/C)+pi;
            end
            
            xc = rc*cos(t).*cos(v);
            yc = rc*cos(t).*sin(v);
            zc = rc*sin(t)+0*v;
        else
            if  abs(P(3)-Q(3)) < 10^-5 && Fexa(i, k) ~= 1 && Fexa(i, k) ~= 2
                rc = sqrt(P(1)^2+P(2)^2);
                              
                % angle between P and Q determined by the "Law of chords"
                angle = 2*asin(PQ/(2*rc));
                
                beta1 = asin(P(2)/rc);
                
                if abs(rc*cos(beta1)-P(1))>10^-10
                    beta1 = pi-beta1;
                end
                
                beta2 = asin(Q(2)/rc);
                if abs(rc*cos(beta2)-Q(1))>10^-10
                    beta2 = pi-beta2;
                end
                
                alpha0 = beta1;
                
                if abs(sum([rc*cos(alpha0+angle), rc*sin(alpha0+angle)]-Q(1:2)))>10^-10
                    alpha0 = beta2;
                end
                t = chebfun(@(t) t, [alpha0, alpha0+angle]);
                
                xc = rc*cos(t);
                yc = rc*sin(t);
                zc = P(3)+0*t;
            else
                angle = 2*asin(PQ/(2*rc));
                
                t = chebfun(@(t) t, [min(beta1, beta2), min(beta1, beta2)+angle]);
                if i == 19 || i == 20
                    t = chebfun(@(t) t, [beta1-angle, beta1]);
                end
                
                v = alpha1;
                xc = rc*cos(t).*cos(v);
                yc = rc*cos(t).*sin(v);
                zc = rc*sin(t)+0*v;
            end
        end
        plot3(xc, yc, zc, 'k', LW, 3)
    end
end

%% Colouring the faces
% Let us now color in black the pentagonal face

vertex = [1, 2, 3;3, 4, 5;1, 3, 5];

for n = 1:6
    if n == 1 || n == 2 || n == 4 || n == 6
        xp = [V(Fpenta(n, :), 1);V(Fpenta(n, 1), 1)];
        yp = [V(Fpenta(n, :), 2);V(Fpenta(n, 1), 2)];
        
        Rx = linspace(min(xp)-0.2, max(xp)+0.2, 40);
        Ry = linspace(min(yp)-0.2, max(yp)+0.2, 40);
        [rx, ry] = meshgrid(Rx, Ry);
        npr = length(Rx);
        
        % d = [min(xp)-0.2, max(xp)+0.2, min(yp)-0.2, max(yp)+0.2];
        
        rz = sqrt(r^2-rx.^2-ry.^2);
        
        for k = 1:3
            Vx = V(Fpenta(n, vertex(k, :)), 1);
            Vy = V(Fpenta(n, vertex(k, :)), 2);
            Vz = V(Fpenta(n, vertex(k, :)), 3);
            
            % Spherical barycentric coordinates
            detM = Vx(1)*(Vy(2)*Vz(3)-Vy(3)*Vz(2))-Vy(1)*(Vx(2)*Vz(3)-Vx(3)*Vz(2))+Vz(1)*(Vx(2)*Vy(3)-Vx(3)*Vy(2));
            detM1 = rx.*(Vy(2)*Vz(3)-Vy(3)*Vz(2))-ry.*(Vx(2)*Vz(3)-Vx(3)*Vz(2))+rz.*(Vx(2)*Vy(3)-Vx(3)*Vy(2));
            detM2 = -rx.*(Vy(1)*Vz(3)-Vy(3)*Vz(1))+ry.*(Vx(1)*Vz(3)-Vx(3)*Vz(1))-rz.*(Vx(1)*Vy(3)-Vx(3)*Vy(1));
            detM3 = rx.*(Vy(1)*Vz(2)-Vy(2)*Vz(1))-ry.*(Vx(1)*Vz(2)-Vx(2)*Vz(1))+rz.*(Vx(1)*Vy(2)-Vx(2)*Vy(1));
            
            b1 = detM1/detM;
            b2 = detM2/detM;
            b3 = detM3/detM;
            
            black_point = (b1 >= 0).*(b2 >= 0).*(b3 >= 0);
            
            plot3(rx(black_point == 1), ry(black_point == 1), rz(black_point == 1), '.k', MS, 20);
            plot3(rx(black_point == 1), ry(black_point == 1), -rz(black_point == 1), '.k', MS, 20);
        end
    else
        yp = [V(Fpenta(n, :), 2);V(Fpenta(n, 1), 2)];
        zp = [V(Fpenta(n, :), 3);V(Fpenta(n, 1), 3)];
        
        Ry = linspace(min(yp)-0.2, max(yp)+0.2, 40);
        Rz = linspace(min(zp)-0.2, max(zp)+0.2, 40);
        [ry, rz] = meshgrid(Ry, Rz);
        npr = length(Ry);
        rx = sqrt(r^2-ry.^2-rz.^2);
        
        for k = 1:3
            Vx = V(Fpenta(n, vertex(k, :)), 1);
            Vy = V(Fpenta(n, vertex(k, :)), 2);
            Vz = V(Fpenta(n, vertex(k, :)), 3);
            
            
            % Spherical barycentric coordinates
            detM = Vx(1)*(Vy(2)*Vz(3)-Vy(3)*Vz(2))-Vy(1)*(Vx(2)*Vz(3)-Vx(3)*Vz(2))+Vz(1)*(Vx(2)*Vy(3)-Vx(3)*Vy(2));
            
            detM1 = rx.*(Vy(2)*Vz(3)-Vy(3)*Vz(2))-ry.*(Vx(2)*Vz(3)-Vx(3)*Vz(2))+rz.*(Vx(2)*Vy(3)-Vx(3)*Vy(2));
            detM2 = -rx.*(Vy(1)*Vz(3)-Vy(3)*Vz(1))+ry.*(Vx(1)*Vz(3)-Vx(3)*Vz(1))-rz.*(Vx(1)*Vy(3)-Vx(3)*Vy(1));
            detM3 = rx.*(Vy(1)*Vz(2)-Vy(2)*Vz(1))-ry.*(Vx(1)*Vz(2)-Vx(2)*Vz(1))+rz.*(Vx(1)*Vy(2)-Vx(2)*Vy(1));
            
            b1 = detM1/detM;
            b2 = detM2/detM;
            b3 = detM3/detM;
            
            black_point = (b1 >= 0).*(b2 >= 0).*(b3 >= 0);
            
            plot3(rx(black_point == 1), ry(black_point == 1), rz(black_point == 1), '.k', MS, 20);
            plot3(-rx(black_point == 1), ry(black_point == 1), rz(black_point == 1), '.k', MS, 20);
        end
    end
end

%% References:
% [1] http://en.wikipedia.org/wiki/Truncated\_icosahedron
% [2] Christoph F\"{u}nfzig: Spherical techniques and their applications in a
% scene graph system: collision detection and occlusion culling. University of
% Braunschweig - Institute of Technology (2007), ISBN 978-3-86727-111-0, pp.
% 1-143


