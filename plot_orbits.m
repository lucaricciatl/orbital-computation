function plot_orbits(ss,ax)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%for iteration over data structure.
fn = fieldnames(ss);
num_struct = numel(fn);
alpha = 0.6;
%preallocate
L = gobjects(1,num_struct+4) ;  %preallocate plots

%per controllare che le orbite siano allineate ai nodi
% [xsu,  ysu] = meshgrid(-50:0.05:50); % Generate x and y data
% zsu = zeros(size(xsu, 1)); % Generate z data
% surf(xsu, ysu, zsu,'facealpha',0.2,'HandleVisibility','off',"FaceColor","white","EdgeColor", "none") % Plot the surface

for k=1:num_struct
    body = ss.(fn{k});
    %---generate orbits---%
    [xk,yk,zk] = body.orbitpath();  
    %---plot orbits---%
    col = [rand,rand,rand];
    L(k) = plot3(ax,xk,yk,zk,"Color",col,"DisplayName",body.name);   
    L(k).Color(4) = alpha;
    %---compute apoapsis/periapsis/ascending/descending---%
    [xp,yp,zp] = body.compute_periapsis(); 
    disp([xp,yp,zp]);
    [xa,ya,za] = body.compute_apoapsis(); 
    [xan,yan,zan] = body.compute_ascending_node();
    [xdn,ydn,zdn] = body.compute_descending_node();
    %---plot apoapsis/periapsis/ascending/descending---%
    pa = plot3(ax,xp,yp,zp,'v','Color','black','MarkerSize',5,"MarkerFaceColor",col,"DisplayName","periapsis",'HandleVisibility','off');   
    aa = plot3(ax,xa,ya,za,'^','Color','black','MarkerSize',5,"MarkerFaceColor",col,"DisplayName","apoapsis",'HandleVisibility','off');  
    an = plot3(ax,xan,yan,zan,'d','Color','black','MarkerSize',5,"MarkerFaceColor",col,"DisplayName","ascending node",'HandleVisibility','off');
    dn = plot3(ax,xdn,ydn,zdn,'d','Color','black','MarkerSize',5,"MarkerFaceColor",col,"DisplayName","descending node",'HandleVisibility','off');    

end
% L(num_struct)= pa;
% L(num_struct)= aa;
% L(num_struct)= an;
% L(num_struct)= dn;
% legend(L,"TextColor",'#FFFFFF');

%---plot sun---%
% [xs,ys,zs] = ss.sun.bodysphere();
% for i=1:5:50
% s = surface(ax,i*xs,i*ys,i*zs,'facealpha',0.2,'HandleVisibility','off');
% s.EdgeColor = "none";
% s.FaceColor = "yellow";
% end

end

 