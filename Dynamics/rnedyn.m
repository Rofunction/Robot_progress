 % ************************************************************************
% >>File name: Robot.fkinematics/Robot.fkine.m
% >>Author: 
% >>Creat Time : 10/1/2021
% >>Note base on the dynamacis paramaters of Puma560
% grav=[0;0;-9.81];(Defualt)
% *************************************************************************
function tau=rnedyn(q,dq,ddq,grav) 
%% SDH参数
sdh.a=[0,-0.425,-0.39225,0,0,0];
sdh.d=[0.089459,0,0,0.10915,0.09465,0.0823];
sdh.alpha=[pi/2, 0, 0, pi/2, -pi/2, 0];
sdh.theta=q; 
%% 连杆参数 i=1:6
m=[3.7, 8.393, 2.33, 1.219, 1.219, 0.1897];    %% I(:,:,i)表示连杆i的惯量矩阵，
[Ix1,Iy1,Iz1]=deal(0.00375, 0.00765, 0.00765);
[Ix2,Iy2,Iz2]=deal(0.0085, 0.208, 0.208);
[Ix3,Iy3,Iz3]=deal(2.46e-4, 0.00719, 0.00719);
[Ix4,Iy4,Iz4]=deal(9.09e-4, 0.00119, 0.00119);
[Ix5,Iy5,Iz5]=deal(9.09e-4, 0.00119, 0.00119);
[Ix6,Iy6,Iz6]=deal(1.22e-4, 8.21e-5, 8.21e-5);
I(:,:,1)=diag([Ix1,Iy1,Iz1]);                  %% I(:,:,i)表示连杆i的惯量矩阵，
I(:,:,2)=diag([Ix2,Iy2,Iz2]);
I(:,:,3)=diag([Ix3,Iy3,Iz3]);
I(:,:,4)=diag([Ix4,Iy4,Iz4]);
I(:,:,5)=diag([Ix5,Iy5,Iz5]);
I(:,:,6)=diag([Ix6,Iy6,Iz6]);
sdh.r(:,1)=[0,-0.02561,0.00193]';              %% sdh.r(i)表示连杆质心位置，
sdh.r(:,2)=[0.2125,0,0.1134]';
sdh.r(:,3)=[0.15,0,0.0265]';
sdh.r(:,4)=[0, -0.0018, 0.01634]';
sdh.r(:,5)=[0, -0.0018, 0.01634]';
sdh.r(:,6)=[0, 0, -0.001159]'; 
z0=[0;0;1];
%% Forward recursion
N=length(q);
for i=1:N
    T=DH(i,sdh);
    R(:,:,i)=T(1:3,1:3)';      %% R(i)=R_{i,i-1}^T
    r(:,i)=R(:,:,i)*T(1:3,4);  %% r(i)=r_{i-1,i}
%   r(:,i)=[sdh.a(i);sdh.d(i)*sin(sdh.alpha(i));sdh.d(i)*cos(sdh.alpha(i))];  %% r(i)=r_{i-1,i}
    if i>1 
        w(:,i)=R(:,:,i)*(w(:,i-1)+dq(i)*z0);
        wdot(:,i)=R(:,:,i)*(wdot(:,i-1)+ddq(i)*z0+cross(w(:,i-1),dq(i)*z0));
        vdot(:,i)=R(:,:,i)*vdot(:,i-1)+cross(wdot(:,i),r(:,i))+cross(w(:,i),cross(w(:,i),r(:,i)));
    else
        w(:,i)=R(:,:,i)*dq(i)*z0; 
        wdot(:,i)=R(:,:,i)*ddq(i)*z0;
        vdot(:,i)=-R(:,:,i)*grav'+cross(wdot(:,i),r(:,i))+cross(w(:,i),cross(w(:,i),r(:,i)));
    end
    vcdot(:,i)=vdot(:,i)+cross(wdot(:,i),sdh.r(:,i))+cross(w(:,i),cross(w(:,i),sdh.r(:,i)));
end
%% Backward recursion
for j=N:-1:1
    if j<N
        f(:,j)=R(:,:,j+1).'*f(:,j+1)+m(j)*vcdot(:,j);
        n(:,j)=-cross(f(:,j), r(:,j)+sdh.r(:,j) )+R(:,:,j+1).'*n(:,j+1)+cross( R(:,:,j+1).'*f(:,j+1),sdh.r(:,j))+...
                I(:,:,j)*wdot(:,j)+cross(w(:,j),I(:,:,j)*w(:,j) ); 
    else
        f(:,j)=m(j)*vcdot(:,j);
        n(:,j)=-cross(f(:,j), (r(:,j)+sdh.r(:,j)) )+I(:,:,j)*wdot(:,j)+...
                    cross(w(:,j), I(:,:,j)*w(:,j) );
                
    end 
    tau(j)=n(:,j).'*R(:,:,j)*z0;
end
end
%%
function T=DH(i,sdh)
cth = cos(sdh.theta(i));
sth = sin(sdh.theta(i));
ca  = cos(sdh.alpha(i));
sa  = sin(sdh.alpha(i));
d = sdh.d(i);
a = sdh.a(i);
   T = [cth, -sth*ca,  sth*sa, a*cth;
        sth,  cth*ca, -cth*sa, a*sth;
         0 ,   sa   ,    ca  ,   d  ; 
         0 ,   0    ,    0   ,   1 ]; 
end