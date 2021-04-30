 %{
    e=qd-q
%}
clear all;  
    clc;
    if isempty(gcp('nocreate'))
        parpool;
    end  
    tic
    q0=zeros(1,6);
    dq0=zeros(1,6);
    tau0=zeros(1,6);
    deta=0;
    e0=[-2.8414,-0.9626,1.5597,-2.1679,1.5708,-1.8710];        
    endTime=10;
    
    disp('runing at Control simulation');
%   Simulation results [t, qlist, qdlist]
    [t,qlist,qdlist,tau,e]=ode45fcn(endTime,q0,dq0,tau0,e0,deta);
    toc
    this.t       =t;       
    this.qlist   =qlist;           
    this.qdlist  =qdlist;
    this.endTime =endTime;
    this.tau     =tau;
    this.e       =e ;  
    PlotFcn(this,'pose');
%     save ('C:\Users\Administrator\Desktop\Robot_progress\Result11.mat','this');

    save('F:\Robot progress\filedata\Result.mat','this'); 
    
    