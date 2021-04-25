%% Kinematics parameters
DOF=6;
    sdh.theta = [pi,0,0,0,pi/2,0];  %% 处于qr状态，即就绪状态；机械臂伸直切垂直
    sdh.a=[0,-0.425,-0.39225,0,0,0];
    sdh.d=[0.089459,0,0,0.10915,0.09465,0.0823];
    sdh.alpha=[pi/2,0,0,pi/2,-pi/2,0]; 
    for i=1:DOF
       ur5link(i)=Link( [sdh.theta(i),sdh.d(i),sdh.a(i),sdh.alpha(i)],'standard' );
    end
    ur5=SerialLink( ur5link,'name','ur5');  
%% Dynanacis parameters
% 连杆参数 i=1:6
% I(:,:,i)表示连杆i的惯量矩阵
m=[3.7, 8.393, 2.33, 1.219, 1.219, 0.1897]; 
% I(:,:,i)表示连杆i的惯量矩阵
[Ix1,Iy1,Iz1]=deal(0.00375, 0.00765, 0.00765);
[Ix2,Iy2,Iz2]=deal(0.0085, 0.208, 0.208);
[Ix3,Iy3,Iz3]=deal(2.46e-4, 0.00719, 0.00719);
[Ix4,Iy4,Iz4]=deal(9.09e-4, 0.00119, 0.00119);
[Ix5,Iy5,Iz5]=deal(9.09e-4, 0.00119, 0.00119);
[Ix6,Iy6,Iz6]=deal(1.22e-4, 8.21e-5, 8.21e-5);
I(:,:,1)=diag([Ix1,Iy1,Iz1]);
I(:,:,2)=diag([Ix2,Iy2,Iz2]);
I(:,:,3)=diag([Ix3,Iy3,Iz3]);
I(:,:,4)=diag([Ix4,Iy4,Iz4]);
I(:,:,5)=diag([Ix5,Iy5,Iz5]);
I(:,:,6)=diag([Ix6,Iy6,Iz6]);
% sdh.r(i)表示连杆质心位置
sdh.r(:,1)=[0,-0.02561,0.00193]';         
sdh.r(:,2)=[0.2125,0,0.1134]';
sdh.r(:,3)=[0.15,0,0.0265]';
sdh.r(:,4)=[0, -0.0018, 0.01634]';
sdh.r(:,5)=[0, -0.0018, 0.01634]';
sdh.r(:,6)=[0, 0, -0.001159]';

for i=1:DOF
   ur5.links(i).qlim=[-2*pi,2*pi];
   ur5.links(i).m  =m(i); 
   ur5.links(i).r  =sdh.r(:,i); 
   ur5.links(i).I  =I(:,:,i);
   ur5.links(i).Jm =0;
%    ur5.links(i).Bm =0;      %% Bm 表示电机的粘性摩擦系数 
   ur5.links(i).Tc =0;      %% Tc 表示库伦摩擦力，非线性；因此 (0+,0-)(这里是不考虑摩擦力)
   ur5.links(i).G  =0;      %% G  表示齿轮箱的转数比
end