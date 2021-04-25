%{
  Note:  OnlyGrav is Only gravity is a kind of tau equal to 0, 
    the robot arm is only affected by its own gravity, 
    and the result is that the robot arm is in an unstable state. 
    The shoulder and elbow joints will hang down and swing back and forth; 
    the waist joint will rotate due to the Coriolis coupling torque. 
    改变ode45fcn.m可以实现在不同的情况下的结果，如在重力和摩檫力作用下
%}
%%
    clear all;
    clc;
    if isempty(gcp('nocreate'))
        parpool;
    end  
    tic
    q0=zeros(1,6);
    dq0=zeros(1,6);
    endTime=5;
%     disp('runing at OnlyGrav/friction');
% ode45 is used in the second method
    [t,qlist,qdlist]=ode45fcn(endTime,q0,dq0);
    toc
    this.t       =t;
    this.qlist   =qlist;
    this.qdlist  =qdlist;
    this.endTime =endTime;
    PlotFcn(this,'pose');
    save('F:\Robot progress\filedata\OnlyGrav.mat','this');  
%% ode45fcn
    function [t,q,dq]=ode45fcn(tf,q0,qd0)
    n=length(q0);                                    %%  initialize q0\dq0 is column with (1x6)
    tspan=[0 tf];
    q0=q0(:); qd0=qd0(:);
    q0=[q0; qd0]; 
    opts=odeset('RelTol',1e-4,'AbsTol',1e-5);
    [t,y]=ode45(@(t,y) fdy2(t,y,tf), tspan, q0, opts);
    q  =y(:,1:n);                                    %% q =[q1,q2,...,q6] is column
    dq =y(:,n+1:2*n);                                %% dq=[dq1,dq2,...,dq6] is column
    end
function xd=fdy2(t,x,tf)                                %% ~ =t is Independent variable must be exist
    n=6;
    q=x(1:n);                                        %% q ,dq are column
    dq=x(n+1:2*n);
    M_q=(Inertia(q));   
    Cor=Coriolis(q,dq);     
    G_q=(Grav(q,[0,0,-9.81])).';                     %% grav=[0,0,-9.81];
    
    disp(['runing at = ', num2str(t/tf*100),'%']);
    tau=zeros(6,1);                                  %% OnlyGrav/Friction is tau=zeros 
    ddq=M_q\( tau -Cor*dq -G_q- fric(0.5,0.1,dq) );  %% fric is friction and it saved in the fric.m
    xd=[x(n+1:2*n); ddq]; 
end
%%
% first method 
%     tau=zeros(1,6);
%     steps=20000;
%     dt=endTime/steps;
%     [qlist,qdlist] = deal(zeros(length(0:dt:endTime),6));
%     p560_model;
%     for i=1:steps
%         ddq = V(q0,dq0,tau);
%         dq = dq0 + ddq * dt;  %% 脱胎于积分中值定理
%         q = q0 + dq * dt;
%         qlist(i+1,:)  = q;
%         qdlist(i+1,:) =dq;
%         disp(['run at =', num2str(i/steps*100),'%']);
%     end
    
    