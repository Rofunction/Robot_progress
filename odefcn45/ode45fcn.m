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
    n  =6;  
    q  =x(1:n);                                      %% q ,dq are column
    dq =x(n+1:2*n);
    tau=x(2*n+1:3*n);
    e =x(3*n+1:4*n);
    
    M_q=(Inertia(q));   
    Cor=Coriolis(q,dq);     
    G_q=(Grav(q,[0,0,-9.81])).';                   %% grav=[0,0,-9.81];
    q_d   =interp1(ts,qd,  t,'spline');            %% q_d  =[qd1,qd2,...,qd6] is column
    dq_d  =interp1(ts,dqd, t,'spline');            %% dq_d =[dq_d1,...,dq_d6] is column
    ddq_d =interp1(ts,ddqd,t,'spline');            %% ddq_d=[ddq_d1,...,ddq_d6] is column
    
    fai1=400; fai2=200; alpha1=0.55; alpha2=1.75; mu=0.05; L1=((alpha2-alpha1)/(alpha2-1))*mu^(alpha1-1);
    L2=((alpha1-1)/(alpha2-1))*mu^(alpha1-alpha2); gama1=100; gama2=10; epi=30;
%     k1=10; k2=500; k3=50; epi=50;  alpha=800;    % k1=1500 k2=1500 ,k3=5,
    e=q_d.'-q; de=dq_d.'-dq; s_bat=de + fai1*e + fai2*(abs(e).^alpha1).*sign(e);                     % e2=-de-alpha*e;
   
    % theta
    for i=1:n
    if ( s_bat(i)==0 ) || ( (s_bat(i)~=0)&&(abs(e(i))>=mu) )
        theta(i,:)=(abs(e(i))^alpha1)*sign(e(i)); 
    elseif ((s_bat(i)~=0) && (abs(e(i))<mu))
        theta(i,:)=L1*e(i) + L2 * ( abs(e(i)) ^alpha2) *sign(e(i) );
    end
    end
    % d_theta
    for i=1:n
    if (s_bat(i)==0) || ( (s_bat(i)~=0)&&(abs(e(i))>=mu) )
        d_theta(i,:)=alpha1*(abs(e(i))^(alpha1-1))*de(i);
    elseif ((s_bat(i)~=0) && (abs(e(i))<mu))
        d_theta(i,:)=L1*de(i) + L2*alpha2*( abs(e(i))^(alpha2-1) ) *de(i);
    end
    end
    s= de + fai1 * e + fai2 * theta;
%     s =k1 * e + e2; （backstepping）
%     s= de1 + k1*e1; 
%     tau=M_q*(ddq_d.' - (k1-alpha)*de - k2*s - k3*sign(s)  ) + Cor * dq+ G_q + epi ;   %% M\fric < epi
%     tau=k1*de+k2*e;
    tau=M_q*(ddq_d.' - fai1*de -fai2*d_theta - gama1*s -gama2*sign(s)) + Cor * dq+ G_q + epi ;
    ddq=M_q\( tau - Cor*dq - G_q - fric(0.5,0.1,dq) - noise(t,1,0.1,10) );
    xd=[x(n+1:2*n); ddq; tau-x(2*n+1:3*n); e-x(3*n+1:4*n)];
    
    disp(['runing at = ', num2str(t/tf*100),'%']);
end