function P = symstep(vessel,ss,dt)
%symstep : function to perform simulation step of solar system
%   Takes in input vessel,solar system and delta time of simulation and 
%   give as a result the logging position of all bodies and the vessel.
fn = fieldnames(ss);
num_struct = numel(fn);
%qui variabile di logging per le posizioni di tutti i corpi
P = zeros(4,num_struct+1);
for k=1:num_struct
    body = ss.(fn{k});
    %muovi corpo
    [x,y,z] = body.step_2body(dt); %posizioni al passo successivo del corpo

    %salva le posizioni
    P(1,k) = x ;
    P(2,k) = y ;
    P(3,k) = z ;
end
%muovi astronave
[xv,yv,zv] = vessel.step(); %posizione navicella al passo successivo
%vessel_check--->   controlla le azioni di influenza gravitazionale
%                   velocità,posizione, ed eventuali manovre da effettuare
P(1,num_struct+1) = xv;
P(2,num_struct+1) = yv;
P(3,num_struct+1) = zv;
end

