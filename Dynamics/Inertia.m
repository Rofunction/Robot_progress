%  M is vector with (nxn), and n is the number of 
%    Degrees of freedom of(DOF) the robotic Manipulators
function M=Inertia(q) % (n*n)
grav=[0,0,0];
n=length(q);
M=zeros(n,n);
qdot=zeros(1,n);
for i=1:n
    qddot=zeros(1,n);
    qddot(i)=1;  %% qddot最后变成一个单位阵
    M(i,:) = rnedyn (q,qdot,qddot,grav);
end
end

