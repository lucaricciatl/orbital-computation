%clear the environment
clc
clear

%---constants---%
t0 = 0;
tf = 1000;
dt = 0.1;
%---initialize objects---%

%---set grahichs properties---%
[fig,ax] = set_graphichs_properties();

%---generate objects---%
ss = solar_system.generate_solar_system();

%---initialize bodies classes---%
solar_init(ss);

%---initialize spaceship---%
spaceship = solar_system.spaceship();
spaceship.compute_orbital_elements();

%---plot static elements---%
plot_orbits(ss,ax);

%---trasferimento #1---%
x_sp = spaceship.x0;
y_sp = spaceship.y0;
z_sp = spaceship.z0;
%plot3(x_sp, y_sp,z_sp, 'o','MarkerSize',8);
%initialPos = [x_sp, y_sp , z_sp];
initialVel = [0,0.2,0];
%[distance, vf ] = orbit_transfer.longpatch(10000,ss.earth.mass,initialPos,initialVel);
hold on

%plot3(distance(:,1), distance(:,2),distance(:,3), 'b' ,'LineWidth', 1);
%plot3(distance(end,1), distance(end,2),distance(end,3), 'o','MarkerSize',8);

%dynamic plot
% for t = 1: tf
%     pause(0.0000000000001)
%     plot(distance(t,1), distance(t,2),'.')
% end

%---CICLO DI DISEGNO---%
fn = fieldnames(ss);

linkdata on
brush on
fn = fieldnames(ss);
num_struct = numel(fn);
x = zeros(1,num_struct);
y = zeros(1,num_struct);
z = zeros(1,num_struct);

for k=1:num_struct
    body =ss.(fn{k});%
    x(k) = body.x0;
    y(k) = body.y0;
    z(k) = body.z0;
    
end


o = [0 0 0];  %# Origin
a = [5,0,0];  %# Vector 1
b = [0 5 0];  %# Vector 2
c = [0,0,5];      %# Resultant
starts = zeros(3,3);
ends = [a;b;c];
p = plot3(ax,x,y,z,'o','MarkerSize',4 ,'MarkerFaceColor','white','HandleVisibility','off');
quiver3(starts(:,1), starts(:,2), starts(:,3), ends(:,1), ends(:,2), ends(:,3),'MarkerSize',8 ,'HandleVisibility','off')

p.XDataSource = 'x(:)';
p.YDataSource = 'y(:)';
p.ZDataSource = 'z(:)';

t = 1: 36525 %un secolo

%2451545 ---> primo gennaio 2000
%2451545 + 365.25 * 29 + 281 ---> 7 settembre 2029
%2451545 + 365.25 * 30 + 281 + 8 ---> 15 settembre 2030
%2451545 + 365.25 * 32 + 122 +15 ---> 15 aprile 2032
%2451545 + 365.25 * 34+ 122 + 183 + 15 ---> 15 ottobre 2034

% %%%%CICLO DI DISEGNO 
for Teph = 2451545 : 1 : 2451545 + 365.25 * 29 + 281
    T = (Teph-2451545)/36525;
    for k=1:num_struct
        body = ss.(fn{k});
        pos = body.stepf(T);
        x(k) = pos(1);
        y(k) = pos(2);
        z(k) = pos(3);
    end
    refreshdata
    %drawnow
    pause(0.1)
end
    

% T1 = (2451545 + 365.25 * 29 + 281 -2451545)/36525;
% print("posizione della terra al 7 settembre 2029")
% earthpos1 = ss.earth.stepf(T1)
% 
% T2 = (2451545 + 365.25 * 30 + 281 + 8)/36525;
% print("posizione di venere al 15 settembre 2030")
% venuspos1 = ss.venus.stepf(T2)
% 
% T3 = (2451545 + 365.25 * 32 + 122 +15)/36525;
% print("posizione della terra al 15 aprile 2032")
% earthpos2 = ss.earth.stepf(T3)
% 
% T4 = (2451545 + 365.25 * 32 + 122 +15)/36525;
% print("posizione di giove al 15 ottobre 2034")
% jupiterpos1 = ss.jupiter.stepf(T4)

%2)disegna posizione nel tempo della navicella


% 
% symstep();
%timepos = solve_solar(t0,dt,tf,ss);

%---compute influence sphere---%
% [inf_xe,inf_ye,inf_ze]  = ss.earth.compute_influence_sphere();
% [inf_xm,inf_ym,inf_zm] = ss.mars.compute_influence_sphere();
% [inf_xj,inf_yj,inf_zj] = ss.jupiter.compute_influence_sphere();
% [inf_xme,inf_yme,inf_zme] = ss.mercury.compute_influence_sphere();
% [inf_xv,inf_yv,inf_zv] = ss.venus.compute_influence_sphere();
% [inf_xsat,inf_ysat,inf_zsat] = ss.saturn.compute_influence_sphere();
% [inf_xu,inf_yu,inf_zu] = ss.uranus.compute_influence_sphere();
% [inf_xn,inf_yn,inf_zn] = ss.neptune.compute_influence_sphere();
%---draw influence spheres---%
% s1 = surface(inf_xe,inf_ye,inf_ze,'facealpha',0.5);
% s2 = surface(inf_xm,inf_ym,inf_zm,'facealpha',0.5);
% s3 = surface(inf_xj,inf_yj,inf_zj,'facealpha',0.5);
% s4 = surface(inf_xme,inf_yme,inf_zme,'facealpha',0.5);
% s5 = surface(inf_xv,inf_yv,inf_zv,'facealpha',0.5);
% s6 = surface(inf_xsat,inf_ysat,inf_zsat,'facealpha',0.5);
% s7 = surface(inf_xu,inf_yu,inf_zu,'facealpha',0.5);
% s8 = surface(inf_xn,inf_yn,inf_zn,'facealpha',0.5);
% s1.EdgeColor = "none";
% s1.FaceColor = "white";
% s2.EdgeColor = "none";
% s2.FaceColor = "white";
% s3.EdgeColor = "none";
% s3.FaceColor = "white";
% s4.EdgeColor = "none";
% s4.FaceColor = "white";
% s5.EdgeColor = "none";
% s5.FaceColor = "white";
% s6.EdgeColor = "none";
% s6.FaceColor = "white";
% s7.EdgeColor = "none";
% s7.FaceColor = "white";
% s8.EdgeColor = "none";
% s8.FaceColor = "white";
%---legend---%


%---ciclo di disegno---%
%plotta orbite
% TIMESCALE = 5 / pi*10^7;        %1anno= (pi*10^7s) == 10 s
% plot = plot3(ax,0,0,0,'o','HandleVisibility','off');%plotta posizione
% 
% 
% fn = fieldnames(ss);
% num_struct = numel(fn);
% 
% S = gobjects(1,num_struct);
% for j = 1:ceil(tf-t0)
%     X = zeros(1,num_struct);
%     Y = zeros(1,num_struct);
%     Z = zeros(1,num_struct);
%     for k=1:num_struct
%         body = ss.(fn{k});
%         [x,y,z] = body.step_2body(dt);
% %         x =timepos(1,k,j); 
% %         y =timepos(2,k,j);
% %         z =timepos(3,k,j); 
%         X(k) = x;
%         Y(k) = y;
%         Z(k) = z;
%         disp(X); %plotta nan perchè non è inizializzatat bene la posizione.
%         disp(Y); %plotta nan perchè non è inizializzatat bene la posizione.
%         disp(Z); %plotta nan perchè non è inizializzatat bene la posizione.
% %         disp(body.influence_radius);
% %          [xs,ys,zs] = sphere(100);
% %          xs = 100 * xs + x;
% %          ys = 100 * ys + y;
% %          zs = 100 * zs + z;
% %        S(k) = surface(ax,xs,ys,zs,'facealpha',0.5,'HandleVisibility','off');
%         %plotta sfera influna
%     end
%     plot.XData = X; %sostituisci coordinate x
%     plot.XData = Y; %sostituisci coordinate y
%     plot.XData = Z; %sostituisci coordinate z
%     %plotta astronave
%     %aggiorna astronave
% pause(0.05)     % pause to control animation speed
% end



% [xmat,ymat,zmat] = peaks(100); %qui ci si mettono le posizioni nel tempo da plottare
%  xvec = xmat(:);
%  yvec = ymat(:);
%  zvec = zmat(:);
%  comet3(xvec,yvec,zvec)
