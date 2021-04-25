 %{ 
 Note : The friction tauque is not considered // fodyn is a abbreviation of forward dynamics 
 q/dq are vectors with (1xn)
 由于输出结果y=[q,dq]是行向量，因此，状态变量中 x 是列向量。
%} 
function [t,q,dq,tau,e]=ode45fcn(tf,q0,qd0,tau0,e0)
    n=length(q0);       %%  initialize q0\dq0 is column with (1x6)
    tspan=[0 tf]; 
    q0=q0(:); qd0=qd0(:); tau0=tau0(:);e0=e0(:);
    q0=[q0; qd0;tau0;e0]; 
    
    ts=linspace(0,tf,10000);
    [qd,dqd,ddqd]=DesirTrj_gd(ts);
    
    opts=odeset('MaxStep',1e-3,'RelTol',1e-4,'AbsTol',1e-5);
%     opts=odeset('RelTol',1e-4,'AbsTol',1e-5);
    [t,y]=ode45(@(t,x) fdy2(t,x,tf,ts,qd,dqd,ddqd), tspan, q0, opts);
    q  =y(:,1:n);                                  %% q =[q1,q2,...,q6] is column
    dq =y(:,n+1:2*n);                              %% dq=[dq1,dq2,...,dq6] is column
    tau=y(:,2*n+1:3*n);
    e  =y(:,3*n+1:4*n);
end
function xd=fdy2(t,x,tf,ts,qd,dqd,ddqd)            %% ~ =t is Independent variable must be exist
    n=6;  
    q=x(1:n);                                      %% q ,dq are column
    dq=x(n+1:2*n);
    tau=x(2*n+1:3*n);
    e=x(3*n+1:4*n);
    
    M_q=(Inertia(q));   
    Cor=Coriolis(q,dq);     
    G_q=(Grav(q,[0,0,-9.81])).';                   %% grav=[0,0,-9.81];
    disp(['runing at = ', num2str(t/tf*100),'%']);
    
    q_d   =interp1(ts,qd,  t,'spline');            %% q_d  =[qd1,qd2,...,qd6] is column
    dq_d  =interp1(ts,dqd, t,'spline');            %% dq_d =[dq_d1,...,dq_d6] is column
    ddq_d =interp1(ts,ddqd,t,'spline');            %% ddq_d=[ddq_d1,...,ddq_d6] is column
    
    e=q_d.'-q; de=dq_d.'-dq;
    k1=300; k2=100; k3=50; epi=50;
    s= de + k1*e; 
    tau=M_q*(ddq_d.' + k1*de + k3*sign(s) + k2*s ) + Cor * dq+ G_q + epi ;   %% M\fric < epi
 %   tau=k1*de+k2*e;
    
    ddq=M_q\( tau -Cor*dq -G_q -fric(0.5,0.1,dq) );
    xd=[x(n+1:2*n); ddq; tau-x(2*n+1:3*n); e-x(3*n+1:4*n)]; 
end