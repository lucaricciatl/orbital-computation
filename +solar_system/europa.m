function europa = europa()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    solar_system.solar_system_parameters();
    europa = solar_system.body();
    europa.name="europa";
    europa.mass = europa_mass;
    europa.diameter = europa_diameter;
    europa.a = europa_a;
    europa.i = europa_i;
    europa.e= europa_e;
    europa.an = europa_an;
    europa.pa = europa_pa;
    europa.l = europa_l;
end

