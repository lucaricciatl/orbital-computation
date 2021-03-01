function body_sphere(obj_id,obj_pos)
% BODY_SPHERE(obj_id,rr) plots a sphere at coordinates OBJ_POS
%   and applies the image of OBJ_ID body to its surface.
%   Available: all solar planets, Pluto, Vesta, Ceres, Sun.
%
%   obj_id   - numeric identifier of the target body (1-12)
%
%   obj_pos  - (x,y,z) coordinates of target body, wrt the Sun
%

    %% Constants
    body = ["venus.jpg" 
            "earth.jpg"  
            "jupiter.jpg"
            "Io.jpg"
            "Io.jpg"];
        
    radii = [6051.8 %venere
             6371   %terra
             69911  %jupiter
             1822  %Sun
             1822];  %Io 
         
	R = radii(obj_id); %[km]
    
    %% Sphere creation
    [xx,yy,zz] = sphere(100);
    sp_hand = surface(obj_pos(1)+R*xx,obj_pos(2)+R*yy,obj_pos(3)+R*zz);
    
    %% Surface change
    img = imread(body(obj_id));

	%img = imrotate(img,180);
% 	tform = affine2d([1 0 0; 0 -1 0; 0 0 1]);
	tform = affine3d(...
			[[1 0 0; 0 -1 0; 0 0 1], [0 0 0]'; 0 0 0 1] );
		
	img = imwarp(img,tform);
	
    set(sp_hand,'facecolor','texture',...
        'cdata',im2double(img),...
        'edgecolor','none');
    
end