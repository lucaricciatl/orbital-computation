function sun = sun()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    solar_system.solar_system_parameters();
    sun = solar_system.body();
    sun.name = "sun";
    sun.mass = sun_mass;
    sun.diameter = sun_diameter;
    sun.x= sun_x0;
    sun.y= sun_y0;
    sun.z= sun_z0;
    sun.vx= sun_vx0;
    sun.vy=  sun_vy0;
    sun.vz = sun_vz0;
    
    sun.a = sun_a;
    sun.i = sun_i;
    sun.e= sun_e;
    sun.an = sun_an;
    sun.pa = sun_pa;
    sun.l = sun_l;
end

