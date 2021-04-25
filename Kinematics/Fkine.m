%{
 >>File name: Robot.fkinematics/Robot.fkine.m
 >>Author: 
 >>Creat Time : 23/11/2020
 Base on 机器人工具箱的tr2rpy.m文件  The FK.m is optimized
 Note:RPY-XYZ, Roll and TAW is in [-pi,pi], i.e. Alpha,Gamma; Pitch (Beta) is in [-pi/2,pi/2]
 Rolat.m is stored in the folder Kinematics.
 XYZ->Rx_y,Ry_p,Rz_r; ZYX->Rz_y,Ry_p,Rx_r
 XYZ是欧拉角的内旋（绕自身坐标系，右乘；旋转顺序是Rx->yaw,Ry->pitch,Rz->roll），故，返回的角的顺序应该是[y,p,r]
 the XYZ is used in the Robot manipulator （RTB工具箱内的所有旋转都是内旋）
 Reference:机器人视觉英文版（p251）
%}
function PoseRef=Fkine(q,opt)
if length(q)~=6
    error('Input error, input is not 6 aixes')  
end
    R=cell(1,6);T=cell(1,6);
    dh.a=[0,-0.425,-0.39225,0,0,0];
    dh.d=[0.089459,0,0,0.10915,0.09465,0.0823];
    dh.alpha=[pi/2,0,0,pi/2,-pi/2,0];
    dh.theta=q;
for i=1:6
    R{i}=Rolat(i,dh);
    if i>1
        T{i}=T{i-1}*R{i};
    else
        T{i}=R{i};
    end
end
    Fkm=T{6};                                   %% Fkm=fkinemartix  
if strcmp(opt,'xyz')
if abs(abs(Fkm(1,3)) -1)<eps(single(1/2))
        R_r=0;                                  %%  there is a singularity case where Pitch=+-pi/2 
                                                %  in which case ROLL is arbitraily set to zero
    if  Fkm(1,3)>0
        R_y=atan2(Fkm(3,2),Fkm(2,2));   % R+Y0
    else
        R_y=-atan2(Fkm(2,1),Fkm(3,1));  % R-Y
    end
        R_p=asin(Fkm(1,3));                                 %% Beta=pi/2
else
        R_r=atan2(-Fkm(1,2),Fkm(1,1));                      %% R_Z 
        R_y=atan2(-Fkm(2,3),Fkm(3,3));                      %% R_X
        R_p=atan2(Fkm(1,3),sqrt(Fkm(3,3)^2+Fkm(2,3)^2));    %% R_Y
end 
end

if strcmp(opt,'zyx')
if abs(abs(Fkm(3,1)) -1)<eps(single(1/2))
        R_r=0;                                  %%  there is a singularity case where Pitch=+-pi/2 
                                                %  in which case ROLL is arbitraily set to zero
    if  Fkm(3,1)<0
        R_y=-atan2(Fkm(1,2),Fkm(1,3));     % R-Y
    else
        R_y=atan2(-Fkm(1,2),-Fkm(1,3));    % R+Y
    end
        R_p=-asin(Fkm(3,1));                               %% Beta=pi/2
else
        R_r=atan2(Fkm(3,2),Fkm(3,3));                      %% R _x
        R_y=atan2(Fkm(2,1),Fkm(1,1));                      %% Y _z 
        R_p=atan2(-Fkm(3,1),sqrt(Fkm(3,3)^2+Fkm(3,2)^2));  %% P _y 
end 
end

 PoseRef=[Fkm(1,4), Fkm(2,4), Fkm(3,4), R_r, R_p, R_y];
%  fkinemartix=Fkm;
end