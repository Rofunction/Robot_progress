%{
        Note: the function is used to computes 10 parameters of each
    joints, which contains 6X10 parameters.
%}
function full_params=full_paramts()
ur5_mode;
Ip(:,1)=[Ix1,0,0,Iy1,0,Iz1];
Ip(:,2)=[Ix2,0,0,Iy2,0,Iz2];
Ip(:,3)=[Ix3,0,0,Iy3,0,Iz3];
Ip(:,4)=[Ix4,0,0,Iy4,0,Iz4];
Ip(:,5)=[Ix5,0,0,Iy5,0,Iz5];
Ip(:,6)=[Ix6,0,0,Iy6,0,Iz6];
    for i=1:DOF
        m(i)=ur5.links(i).m;
        mr(:,i)=ur5.links(i).m .* ur5.links(i).r;
        full_params(:,i)=[Ip(:,i);mr(:,i);m(i)];
    end
end