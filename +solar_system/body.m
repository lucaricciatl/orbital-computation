
%p0 = x0,y0,z0
%v0 = vx0,vy0,vz0
%---
%---ORBITAL PARAMETERS---
%here the animation : https://www.youtube.com/watch?v=QZrYaKwZwhI
%a -> major semiaxis (Km)
%e -> eccentricity 
%i -> inclination [deg]
%an -> ascending node longitude [deg] (big omega)
%pa -> Argument of periapsis [deg] (little omega)
%ma -> mean anomaly (M)
%attractor -> body class , main attractor of the aster

classdef body < handle
    %celestial body class with 
    properties
        name 
        color
        t0 = 0; %il tempo iniziale è il primo gennaio del 2000 
        t = 0;
        %posizioni relative all'attrattore
        x = 0;
        y = 0;
        z = 0;
        %posizione assoluta J2000 rispetto all'attrattore principale
        x0 = 0;
        y0 = 0;
        z0 = 0;
        %posizione assoluta rispetto all'attrattore principale
        xg = 0;
        yg = 0;
        zg = 0;
        %velocità relative all'atrattore
        vx = 0;
        vy = 0;
        vz = 0;
        a %req
        e %req
        i %req
        an %req
        pa % req
        l %req
        %---
        l_dot = 0;
        e_dot = 0;
        i_dot = 0;
        a_dot = 0;
        pa_dot= 0;
        an_dot = 0;
        %---
        mass %req
        diameter %req
        attractor %---body---%
        orbit
        influence_radius
        periapsis = [0,0,0]
        apoapsis = [0,0,0]
        ascending_node = [0,0,0]
        descending_node = [0,0,0]
        T
        G = 6.67*10^-11;
        h %Specific relative angular momentum
        E %eccentric anomaly
        coe 
       
    end
    
    %----methods---%  
    methods     
        function compute_orbital_elements(self)
        	self.compute_influence_radius();
            self.compute_apoapsis();
            self.compute_periapsis();
            self.compute_period();
            self.compute_ascending_node();
            self.compute_descending_node();
            self.compute_angolar_moment();
            self.compute_eccenric_anomaly();
            self.setpos();
            self.compute_initial_position();
            self.extract_coe()
        end
        
        function extract_coe(self)
            a = self.a;
            e = self.e;
            i = self.i;
            an = self.an;
            pa = self.pa;
            l = self.l;
            self.coe = [a ,e , i ,an ,pa ,l];
        end
        
        function compute_angolar_moment(self)
            if isempty(self.attractor)== true
                hm = 0;
            else
                hm = power((1-power(self.e,2)*self.G*self.attractor.mass),1/2);
            end
            self.h = hm;
        end
        function setpos(self)
            if isempty(self.attractor)== true
                xc = self.x;
                yc = self.y;
                zc = self.z;
            else
                xc = self.x + self.attractor.x;
                yc = self.y + self.attractor.y;
                zc = self.z + self.attractor.z;
            end
            self.xg = xc;
            self.yg = yc;
            self.zg = zc;
        end
        function [xa,ya,za] = compute_apoapsis(self)
            apoapsis_readius = self.a*(1-self.e^2)/(1+self.e*cos(self.pa-pi)); 
            Tr = T_matrix.Tzxz(-self.an,-self.i,-(self.pa+pi));
            Ta = T_matrix.Ttr(apoapsis_readius,0,0) ;
            T_orb = Tr*Ta;
            self.apoapsis = T_orb(1:3,4);
            xa = self.apoapsis(1);
            ya = self.apoapsis(2);
            za = self.apoapsis(3);
        end
        function [xp,yp,zp] = compute_periapsis(self)
            periapsis_readius = self.a*(1-self.e^2)/(1+self.e*cos(self.pa));
            Tr = T_matrix.Tzxz(-self.an,-self.i,-self.pa);
            Ta = T_matrix.Ttr(periapsis_readius,0,0) ;
            T_orb = Tr*Ta;
            self.periapsis = T_orb(1:3,4);
            xp = self.periapsis(1);
            yp = self.periapsis(2);
            zp = self.periapsis(3);
        end
        function [xa,ya,za] = compute_ascending_node(self)
            Tr = T_matrix.Tzxz(-self.an,-self.i,0);
            r = self.a*(1-self.e^2)/(1+self.e*cos(0));
            Ta = T_matrix.Ttr(r,0,0) ;
            T_orb = Tr*Ta;
            self.ascending_node = T_orb(1:3,4);
            xa = self.ascending_node(1);
            ya = self.ascending_node(2);
            za = self.ascending_node(3);
        end
        function [xa,ya,za] = compute_descending_node(self)
            Tr = T_matrix.Tzxz(-self.an,-self.i,-pi);
            r = self.a*(1-self.e^2)/(1+self.e*cos(pi));
            Ta = T_matrix.Ttr(r,0,0) ;
            T_orb = Tr*Ta;
            self.descending_node = T_orb(1:3,4);
            xa = self.descending_node(1);
            ya = self.descending_node(2);
            za = self.descending_node(3);
        end
        function r = compute_influence_radius(self)
            if isempty(self.attractor)== true
                r = 0;
            else 
                r = (self.a+self.a*self.e) * power(self.mass/self.attractor.mass,2/5);
            end
            self.influence_radius = r;
        end
        function compute_period(self)
            if isempty(self.attractor)== true
            else
            period = 2*pi*power(power(self.a/((self.attractor.mass+self.mass)*self.G),3),1/2);
            self.T = period;
            end
        end  
        function [xs,ys,zs] = compute_influence_sphere(self)
            [xs,ys,zs] = sphere(100);
            self.compute_influence_radius();
            xs = self.influence_radius * xs + self.x;
            ys = self.influence_radius * ys + self.y;
            zs = self.influence_radius * zs + self.z;
        end
        function E = compute_eccenric_anomaly(self)
            e = self.e;
            M = self.l-self.pa;     %mean anomaly
            E = M + e*sin(M);
            delta_E = 1;
            tol = 1e-06; %errore max consentito

            while(abs(delta_E)>tol)
                delta_M = M-(E-e*sin(E));
                delta_E = delta_M/(1-e*cos(E));
                E = E + delta_E;
            end
            self.E = E;
        end
        
        function pos = compute_initial_position(self)
            om_bar = self.pa;
            om_big = self.an;
            L = self.l;
            e = self.e;
            I = self.i;
            a = self.a;

            om = om_bar - om_big;   %periheion
            %M = L - om_bar;                         %mean anomaly
            E = self.E;

            position_first = [ a*(cos(E)-e);
                               a*sqrt(1-e^2)*sin(E);
                               0 ];
            Trot = T_matrix.Tzxz(-om_big,-I,-om);
            pos = Trot(1:3,1:3) * position_first;
            self.x0 = pos(1);
            self.y0 = pos(2);
            self.z0 = pos(3);
        end
       
        function orbit_pos = orbitf(self,s)
           %calculate orbit parameter in s 
           %t = t*360/(2*pi);
           Tr = T_matrix.Tzxz(-self.an,-self.i,-s);
           %mu = G*self.attractor.mass;
           %theta ->angle with periapsis
           %l = self.m*r^2*theta_dot
           r = norm(self.a*(1-self.e^2)/(1+self.e*cos(s)));
           %r = self.a*(1-self.e*cos(self.E))
           Ta = T_matrix.Ttr(r,0,0) ;
           T_orb = Tr*Ta;
           orbit_pos = T_orb(1:3,4);
           %---
           %x=self.a*(cos(self.E)-self.e);
           %y= self.a.*sqrt(1-self.e^2)*sin(self.E);
           %z=0;
           %T_orb = Tr*Ta;
           %p = [x y z 0]
           %orbit_pos = Tr.*p
           %orbit_pos = orbit_pos(1:3)
        end

        function [x,y,z] = orbitpath(self)
                compute_eccenric_anomaly(self)
                s = linspace(0,2*pi,10000);
                ts = size(s);
                track = zeros(ts(2),3);
                for k=1:ts(2)
                    track(k,:) = self.orbitf(s(k));
                end
                x = track(:,1);
                y = track(:,2);
                z = track(:,3);
        end
        
        function [x,y,z] = bodysphere(self)
            [X,Y,Z] = sphere(100);
            x = self.diameter/2 * X;
            y = self.diameter/2 * Y;
            z = self.diameter/2 * Z;
        end
        
        function pos = stepf(self,T)
        pa = self.pa %+ self.pa * T;
        an = self.an %+ self.an * T;
        l = self.l + self.l_dot * T;
        e = self.e %+ self.e_dot * T;
        i = self.i %+ self.i_dot * T;
        a = self.a %+ self.a_dot * T;
        
        om = -(pa - an);
        
        e = self.e;
        M = l-pa;     %mean anomaly
        E = M + e*sin(M);
        delta_E = 1;
        tol = 1e-06; %errore max consentito

        while(abs(delta_E)>tol)
            delta_M = M-(E-e*sin(E));
            delta_E = delta_M/(1-e*cos(E));
            E = E + delta_E;
        end

        position = [ a*(cos(E)-e);
                               a*sqrt(1-e^2)*sin(E);
                               0 ];
                           
        r = a*(1-e^2)/(1+e*cos(om));
         
        Ta = T_matrix.Ttr(r,0,0);
        Trot = T_matrix.Tzxz(-an,-self.i,-om);
        T_orb = Trot * Ta;
        pos = Trot(1:3,1:3) * position;

        self.x = pos(1);
        self.y = pos(2);
        self.z = pos(3);
      
        end
        
        function ta = compute_time_anomaly(self,T)
        pa = self.pa %+ self.pa * T;
        an = self.an %+ self.an * T;
        l = self.l + self.l_dot * T;
        e = self.e %+ self.e_dot * T;
        i = self.i %+ self.i_dot * T;
        a = self.a %+ self.a_dot * T;
        
        om = -(pa - an);
        
        e = self.e;
        M = l-pa;     %mean anomaly
        E = M + e*sin(M);
        delta_E = 1;
        tol = 1e-06; %errore max consentito

        while(abs(delta_E)>tol)
            delta_M = M-(E-e*sin(E));
            delta_E = delta_M/(1-e*cos(E));
            E = E + delta_E;
        end

        position = [ a*(cos(E)-e);
                               a*sqrt(1-e^2)*sin(E);
                               0 ];
                           
        r = a*(1-e^2)/(1+e*cos(om));
         
        Ta = T_matrix.Ttr(r,0,0);
        Trot = T_matrix.Tzxz(-an,-self.i,-om);
        T_orb = Trot * Ta;
        pos = Trot(1:3,1:3) * position;

        ta = 2*atan( sqrt((1+e)/(1-e)) * tan(E/2));
        end
           %%%%%importante --> metodo per calcolare l'orbita ne ltempo ,
           %%%%%iterativo
%         function En = newton_iteraction(self)
%             E = self.E; %
%             e = self.E; %
%             Mm = self.Mm %
%             M0 = self.M0 %
%             %%%iterazioni????
%             while abs(En-Ea) < 0.001
%               Ea = En;
%               En = Ea -(Ea-e*sin(e)-(Mm-M0))/(1-e*cos(Ea))); 
%               end
%             self.E = En; %anomalia eccentrica
%       end

%         function [xn,yn,zn] = timestep(self,dt)
%             
%             end
        end
end
