%  Grav is vector with (1xn), and n is the number of 
%    Degrees of freedom of(DOF) the robotic Manipulators
function Grav=Grav(q,grav)  
% defualt grav=[0,0,-9.81]
    n = length(q);
    qdot = zeros(1,n);
    qddot = zeros(1,n);
    Grav = rnedyn(q,qdot,qddot,grav);
end