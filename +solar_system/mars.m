function mars = mars()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    solar_system.solar_system_parameters();
    mars = solar_system.body();
    mars.name = "mars";
    mars.mass = mars_mass;
    mars.diameter = mars_diameter;
    mars.a = mars_a;
    mars.i = mars_i;
    mars.e= mars_e;
    mars.an = mars_an;
    mars.pa = mars_pa;
    mars.l = mars_l;
    
    mars.a_dot = mars_a_dot;
    mars.l_dot = mars_l_dot;
    mars.e_dot = mars_e_dot ;
    mars.i_dot = mars_i_dot;
    mars.an_dot = mars_pa_dot;
    mars.pa_dot = mars_an_dot;
end

