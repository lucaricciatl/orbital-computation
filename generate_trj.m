
clc;
clear;

[fig,ax] = set_graphichs_properties();

sun = solar_system.body();
sun.name = "sun";
sun.mass = 0;

sun.diameter = 0;
sun.x= 0;
sun.y= 0;
sun.z= 0;
sun.vx= 0;
sun.vy=  0;
sun.vz = 0;

traveller1 = solar_system.body() ;                     
traveller1.name = 'traveller1' ;

s = struct();
s.earth = solar_system.earth();
s.venus = solar_system.venus();
s.jupiter = solar_system.jupiter();
plot_orbits(s,ax)


ind = 'pro';
av = 'retro';

G = 6.67*10^-11; %[m^3/(kg*s^2)]
mu_sun = 1.327565122000000e+11;   %[km^3/(s^2)]
AU_KM = 1/149597828; %[km/au]
DR = 2*pi/360; %%deg to radiants conversion units

J_time_terra1 = 2462386.5;
J_time_venus1 = 2462606.5;
J_time_terra2 = 2463337.5;
J_time_jupiter1 = 2464250.5;

earthpos1 = s.earth.stepf(2462386.5)';
venuspos1 = s.venus.stepf(2462606.5)';
earthpos2 =  s.earth.stepf(2463337.5)';
jupiterpos1 = s.jupiter.stepf(2464250.5)';
jupiterpos2 = s.jupiter.stepf(2464210.5)'; 
jupiterspeed = (jupiterpos1-jupiterpos2)/AU_KM/(10*3600*24);
days1=373;
days2=567;
days3= 365*2.5;
daysh = days1+days2+days3;
v_terra1 = 27.5;
v_venere = 42.2;
%orbit_transer.lambert(R1, R2, dt, tolerance)
o = [0 0 0];  %# Origin
a = [5,0,0];  %# Vector 1
b = [0 5 0];  %# Vector 2
c = [0,0,5];      %# Resultant
starts = zeros(3,3);
ends = [a;b;c];
v_terra2 = 38.1;
quiver3(starts(:,1), starts(:,2), starts(:,3),...
        ends(:,1), ends(:,2), ends(:,3), ...
        'MarkerSize',8 , ...
        'HandleVisibility','off')

%---
%plotta le posizioni dei pianeti
plot3(ax,earthpos1(1),earthpos1(2),earthpos1(3), ...
      'o','MarkerSize',8 ,'HandleVisibility','off', ...
      'color','r')
  
plot3(ax,venuspos1(1),venuspos1(2),venuspos1(3),'o','MarkerSize',8 ,'HandleVisibility','off','color','g')
plot3(ax,earthpos2(1),earthpos2(2),earthpos2(3),'o','MarkerSize',8 ,'HandleVisibility','off','color','b')
plot3(ax,jupiterpos1(1),jupiterpos1(2),jupiterpos1(3),'o','MarkerSize',8 ,'HandleVisibility','off','color','m')


%%PRIMO TRASFRIMENTO TERRA1-VENERE
[v1_1, v2_1] = lib.lambert(earthpos1/AU_KM,venuspos1/AU_KM,conv.days_to_second(days1),av,mu_sun);
[orb1, t1] = orbit_transfer.intpl_orbit(days1 , earthpos1/ AU_KM , v1_1 , mu_sun);
plt1 = plot3(ax,orb1(:,1)* AU_KM,orb1(:,2)* AU_KM,orb1(:,3)* AU_KM,'r','linewidth',2,'HandleVisibility','off')

%%SECONDO TRASFRIMENTO VENERE-TERRA2
[v1_2, v2_2] = lib.lambert(venuspos1 / AU_KM,earthpos2 / AU_KM, conv.days_to_second(days2),av,mu_sun);
[orb2, t2] = orbit_transfer.intpl_orbit(days2, venuspos1/ AU_KM , v1_2 ,mu_sun);
plt2 = plot3(ax,orb2(:,1)*AU_KM,orb2(:,2)* AU_KM,orb2(:,3)* AU_KM,'g','linewidth',2,'HandleVisibility','off')

%%TERZO TRASFRIMENTO TERRA2-GIOVE
[v1_3, v2_3] = lib.lambert(earthpos2 / AU_KM,jupiterpos1/ AU_KM, conv.days_to_second(days3),av,mu_sun);
[orb3, t3] = orbit_transfer.intpl_orbit(days3, earthpos2/ AU_KM , v1_3 ,mu_sun);
plt3 = plot3(ax,orb3(:,1)*AU_KM,orb3(:,2)* AU_KM,orb3(:,3)* AU_KM,'b','linewidth',2,'HandleVisibility','off')

%%TRASFRIMENTO HOMANN TERRA1
[v1_h1, v2_h1] = lib.lambert(earthpos1 / AU_KM,jupiterpos1/ AU_KM, conv.days_to_second(daysh),av,mu_sun);
[orb1h, th1] = orbit_transfer.intpl_orbit(daysh, earthpos1/ AU_KM , v1_h1 ,mu_sun);
plt4 =plot3(ax,orb1h(:,1)*AU_KM,orb1h(:,2)* AU_KM,orb1h(:,3)* AU_KM,'w','linewidth',2,'HandleVisibility','off');
%%TRASFRIMENTO HOMANN TERRA2
[v1_h2, v2_h2] = lib.lambert(earthpos2 / AU_KM,jupiterpos1/ AU_KM, conv.days_to_second(daysh),ind,mu_sun);
[orb2h, th2] = orbit_transfer.intpl_orbit(daysh, earthpos2/ AU_KM , v1_h2 ,mu_sun);
plt5 = plot3(ax,orb2h(:,1)*AU_KM,orb2h(:,2)* AU_KM,orb2h(:,3)* AU_KM,'y','linewidth',2,'HandleVisibility','off')
v_terra1_hohmann = 22;
v_terra2_hohmann = 30;

COMPLEXORBIT1 = vertcat([orb1(:,1)* AU_KM,orb1(:,2)* AU_KM,orb1(:,3)* AU_KM], ...
                        [orb2(:,1)*AU_KM,orb2(:,2)* AU_KM,orb2(:,3)* AU_KM] );
COMPLEXORBIT1 = vertcat(COMPLEXORBIT1,[orb3(:,1)*AU_KM,orb3(:,2)* AU_KM,orb3(:,3)* AU_KM]);
            
COMPLEXORBIT2 = [orb1h(:,1)*AU_KM,orb1h(:,2)* AU_KM,orb1h(:,3)* AU_KM];

COMPLEXORBIT3 = [orb2h(:,1)*AU_KM,orb2h(:,2)* AU_KM,orb2h(:,3)* AU_KM];
                

L1 = legend([plt1,plt2,plt3,plt4,plt5], ...
            'first trasansfer earth1-venus1', ...
            'second trasansfer venus1-earth2', ...
            'third trasansfer earth2-jupiter1', ...
            'homann trasansfer earth1-jupiter1', ...
            'homann trasansfer earth2-jupiter1 retrograde' ...
            );


%plotta la prima traiettoria completa nel tempo
sizearr1 = size(orb1(:,1));
sizearr2 = size(orb2(:,1));
sizearr3 = size(orb3(:,1));
sizearr1 = sizearr1(1);
sizearr2 = sizearr2(1);
sizearr3 = sizearr3(1);
sizearr = sizearr1+sizearr2+sizearr3;

x = zeros(1,1);
y = zeros(1,1);
z = zeros(1,1);
plt = plot3(ax,x,y,z,'o','MarkerSize',8,'MarkerFaceColor','r','MarkeredgeColor','black','HandleVisibility','off');

plt.XDataSource = 'x(1,:)';
plt.YDataSource = 'y(2,:)';
plt.ZDataSource = 'z(3,:)';
for k = 2 : 200 : sizearr 
    x(1,k) = COMPLEXORBIT1(k,1);
    y(2,k) = COMPLEXORBIT1(k,2);
    z(3,k) = COMPLEXORBIT1(k,3);
    refreshdata
    %drawnow
    pause(0.01)
    k = k;
end



%plotta la prima traiettoria  homann completa nel tempo
sizearr1 = size(orb1h(:,1));

sizearr = sizearr1;

x = zeros(1,1);
y = zeros(1,1);
z = zeros(1,1);
plt = plot3(ax,x,y,z,'o','MarkerSize',8,'MarkerFaceColor','w','MarkeredgeColor','black','HandleVisibility','off');

plt.XDataSource = 'x(1,:)';
plt.YDataSource = 'y(2,:)';
plt.ZDataSource = 'z(3,:)';
for k = 2 : 200 : sizearr 
    x(1,k) = COMPLEXORBIT2(k,1);
    y(2,k) = COMPLEXORBIT2(k,2);
    z(3,k) = COMPLEXORBIT2(k,3);
    refreshdata
    %drawnow
    pause(0.01)
    k = k;
end

%plotta la seconda traiettoria  homann completa nel tempo
sizearr1 = size(orb2h(:,1));

sizearr = sizearr1;

x = zeros(1,1);
y = zeros(1,1);
z = zeros(1,1);
plt = plot3(ax,x,y,z,'o','MarkerSize',8,'MarkerFaceColor','y','MarkeredgeColor','black','HandleVisibility','off');

plt.XDataSource = 'x(1,:)';
plt.YDataSource = 'y(2,:)';
plt.ZDataSource = 'z(3,:)';
for k = 2 : 200 : sizearr 
    x(1,k) = COMPLEXORBIT3(k,1);
    y(2,k) = COMPLEXORBIT3(k,2);
    z(3,k) = COMPLEXORBIT3(k,3);
    refreshdata
    %drawnow
    pause(0.01)
    k = k;
end




trasferimento = {'terra1-vener1';'venere1-terra2';'terra2-giove1';'terra1-giove1 homann';'terra2-giove1 homann'};
posizione_iniziale = [earthpos1;venuspos1;earthpos2;earthpos1;earthpos1];
posizione_finale = [venuspos1;earthpos2;earthpos1;jupiterpos1;jupiterpos1];
velocita_partenza = [v1_1; v1_2; v1_3; v1_h1;v1_h2];
velocita_arrivo = [v2_1; v2_2; v2_3; v2_h1;v2_h2];
giorni_impiaegati = [days1 ;days2; days3; daysh;daysh];
v_tot_final=97.62;
T1 = table(trasferimento, ...
          posizione_iniziale, ...
          posizione_finale, ...
          velocita_partenza, ...
          velocita_arrivo, ...
          giorni_impiaegati ...
          )

text = {'deltaV(Non coindice con il carburante consumato)'};
trasferimento = {'1º patch Terra1-Venere1 + Δv fuga dalla terra';'Δv 1º flyby su 1º Venere1';'Δv 2º flyby su Terra2';'Δv Trasferimento alla Hohmann Terra1-Giove1 ';'Δv Trasferimento alla Hohmann Terra2-Giove1 (retrogrado) ';'Δv Terra2-Io(compresi Δv interni alla SOI)';'Totale Δv per missione su Io'};
deltaV = [norm(v1_1)-v_terra1 ; norm(v1_2)-v_venere ; norm(v1_3)-v_terra2 ; norm(v1_h1)-v_terra1_hohmann ; norm(v1_h2)+v_terra2_hohmann; 8.18 ;norm(v1_1)+norm(v1_2)+norm(v1_3)-v_tot_final];





T2 = table(trasferimento, ... 
          deltaV)



%% Venere Close-up (317) FLYBY
T0 = J_time_venus1 - 1;
T1 = J_time_venus1 + 1;
mu_venus = G*10^-9 * s.venus.mass;
dfl = 1;
vpl = (orb1(end,1:3)-orb1(end-dfl,1:3))/(t1(end)-t1(end-dfl))  ;
[e,a,v_out1] = flyby_NDM(v2_1 ,vpl ,orb1(end-dfl,1:3)/AU_KM,mu_venus);
deltaV1 = norm(v_out1-v1_1);


fly1 = solar_system.body();
fly1.name = 'fly1';
fly1.a =a * AU_KM;
fly1.e = e;
fly1.i = 0;
fly1.an = 0;
fly1.pa = 0;
fly1.l = 0;
fly1.attractor = s.venus;

sv = struct();
sv.venus = s.venus;
sv.fly1 = fly1;

%set graphics
fig1 = figure;
ax1=axes;
opengl software
hold on
light
axis tight
shading interp
lighting gouraud 
grid on
xlabel('x[au]')
ylabel('y[au]')
zlabel('z[au]')
set(gca,'Color','k')
set(gcf,'Color','black');
set(gca,'xcolor','w')
set(gca,'ycolor','w')
set(gca,'zcolor','w')
view(3)
width=1920/2;
height=1080/2;
set(gcf,'position',[0,0,width,height])
xlim([-0.0000006 0.0000006])
ylim([-0.0000006 0.0000006])
zlim([-0.0000006 0.0000006])
daspect([1 1 1])
set(groot,'DefaultFigureGraphicsSmoothing','on')
    
[x,y,z] = sv.fly1.orbitpath();
plot3(0,0,0,'o')
R_SOI = 660400*AU_KM;
pos_body_sphere = [0,0,0];
% [X,Y,Z] = sphere(100) ;
% surf = surface(ax1,R_SOI* X, R_SOI* Y, R_SOI*Z,'facealpha',0.2,'HandleVisibility','off');
% surf.EdgeColor = "none";
% surf.FaceColor = "yellow";
body_sphere_AU(1,pos_body_sphere)
plot3(x,y,z);
R_PL = s.venus.diameter/2*AU_KM;
[X,Y,Z] = sphere(100) ;
surf = surface(ax1,R_PL* X, R_PL* Y, R_PL*Z,'facealpha',0.9,'HandleVisibility','off');
surf.EdgeColor = "none";
surf.FaceColor = "yellow";
plot3(x,y,z);
hold off

%% Terra Close-up (317) FLYBY
T0 = J_time_terra2 - 1;
T1 = J_time_terra2 + 1;
mu_earth = G*10^-9 * s.earth.mass;
dfl = 10;
vpl = (orb2(end,1:3)-orb2(end-dfl,1:3))/(t2(end)-t2(end-dfl))  ;
[e,a,v_out2] = flyby_NDM(v2_2 ,vpl ,orb2(end-dfl,1:3)/AU_KM,mu_earth);
deltaV1 = norm(v_out2-v1_3);

fly2 = solar_system.body();
fly2.name = 'fly2';
fly2.a =a * AU_KM;
fly2.e = e;
fly2.i = 0;
fly2.an = 0;
fly2.pa = 0;
fly2.l = 0;
fly2.attractor = s.earth;

sv = struct();
sv.earth = s.earth;
sv.fly2 = fly2;
%set graphics
fig1 = figure;
ax2=axes;
opengl software
hold on
light
axis tight
shading interp
lighting gouraud 
grid on
xlabel('x[au]')
ylabel('y[au]')
zlabel('z[au]')
set(gca,'Color','k')
set(gcf,'Color','black');
set(gca,'xcolor','w')
set(gca,'ycolor','w')
set(gca,'zcolor','w')
view(3)
width=1920/2;
height=1080/2;
set(gcf,'position',[0,0,width,height])
xlim([-0.0000003 0.0000003])
ylim([-0.0000003 0.0000003])
zlim([-0.0000003 0.0000003])
daspect([1 1 1])
set(groot,'DefaultFigureGraphicsSmoothing','on')
 pos_body_sphere = [0,0,0];   

[x,y,z] = sv.fly2.orbitpath();
plot3(0,0,0,'o')
R_SOI = 0.01;
% [X,Y,Z] = sphere(100) ;
% surf = surface(ax2,R_SOI* X, R_SOI* Y, R_SOI*Z,'facealpha',0.2,'HandleVisibility','off');
% surf.EdgeColor = "none";
% surf.FaceColor = "yellow";
body_sphere_AU(2,pos_body_sphere)
plot3(x,y,z);
R_PL = s.earth.diameter/2*AU_KM;
[X,Y,Z] = sphere(100) ;
surf = surface(ax2,R_PL* X, R_PL* Y, R_PL*Z,'facealpha',0.9,'HandleVisibility','off');
surf.EdgeColor = "none";
surf.FaceColor = "yellow";
plot3(x,y,z);
hold off


%% plot Cattura su Giove (Marco)
figure(10)
[traj,delta_v,hyp] = capture_hyp2(3,jupiterpos1,jupiterspeed,orb3(end-2:end,1:3),45000,v2_3,s.jupiter.an,s.jupiter.i);
[park_jup,t_V2] = park_in(3, jupiterpos1, 45000, s.jupiter.an,s.jupiter.i,[0,0,0], [2034,10,15,0,0,0], [2034,11,5,0,0,0]);
hold off

%%%si è scelta un'orbita di parcheggio maggiore per poter fare un trasf. alla Hohmann di diminuzione dell'orbita
%%%
figure(11)
%[traj,delta_v,hyp] = capture_hyp2(3,jupiterpos1,jupiterspeed,orb3(end-2:end,1:3),450000,v2_3,s.jupiter.an,s.jupiter.i);

[park_jup,t_V2] = park_in(3, jupiterpos1, 450000, s.jupiter.an,s.jupiter.i,[0,0,0], [2034,10,15,0,0,0], [2034,11,5,0,0,0]);
%                 park_in(obj_id, body_pos, radius, ra,i, p_ref, start, finish)
jupiter_diameter_m = s.jupiter.diameter*149597828*1000;
jupiter_radii_m = jupiter_diameter_m/2;
%io orbita a 350000 dalle nubi di Giove
%FARE ORBITA IO CON CENTRO IN 0,0 E RAGGIO 421800
hold on
%calcolo con lambert la delta_v per la manovra di Hohmann
t_io = 60*24*3600;
G_KM_io = 6.67408e-20;  %km^3/Kg^2
jupiter_mass = 1.8981e+27;   %kg
pl_mu = G_KM_io * jupiter_mass; %km^3/s^2
% [v1_io, v2_io] = lib.lambert([0,450000,0],[0,-421800,0],t_io,'pro',pl_mu)
[orbit4, total_time, last_vel, total_deltaV]= homan(s.jupiter.mass, jupiter_radii_m, 450000000);


%cattura su IO
%arriva lo stesso giorno (homan ci sta 0.4099 giorni) --position:421800---0.0028
hold on
io_position = [0,421800,0];
io_speed1 = [-( G_KM_io*jupiter_mass/421800^3) ,0,0]; %[km/s]
diff_orbit = [0, 0.03*10^5, 0; 0, 0.034*10^5, 0];
io_speed2 = [-( G_KM_io*(0.08931*10^24)/(421800+0.03*10^5)^3) ,0,0];
%[traj,delta_v,hyp] = capture_Io(3,io_position,io_speed,diff_orbit,200,v2_3,0,0);
%                    capture_hyp2(goal_id,pl_r0,v_arr,orbit,park_r, v_in,i, RA)
[park_jup,t_V2] = park_in(5, io_position, 200, 0,0,[0,0,0], [2034,11,15,0,0,0], [2034,11,31,0,0,0]);


%%PLOT FUGA DELLA TERRA
solar_system_and_time

      