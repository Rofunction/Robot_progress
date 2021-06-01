%{
   ----------------------------------------------------------------------------------------------
    Symbol describtion:
        w - fundamental frequency
        A - sine coefficient
        B - cosine coefficient
        C - coefficient of fifth order polynomail 
        N - number of harmonics
       q0 - initial offset (nX1)
       qf - the computation result of Fourier series trajectory
       qp - the computation result of fifth order polynomail trajectory
    Note:
        A,B,C are optimized parameters by patter search meothod or genetic
     algorithm (GA), which to make the trajectory sustainably excited, also, to make the velocity
     and acceleration zero at the initial and final.
    ----------------------------------------------------------------------------------------------
%}
function [q,dq,ddq]=Fourier_series_trj(t,wf,A,B,q0,N,C)
    qf=q0; % q0 = zeros(6,1);
    qp = C(:,1); dqp = zeros(size(qp));  ddqp = zeros(size(qp));
    dqf=zeros(size(q0));  ddqf=zeros(size(q0));
    for k=1:N
        qf   =  qf  + A(:,k)/(wf*k).*sin(wf*k*t) - B(:,k)/(wf*k)*cos(wf*k*t);
        dqf  =  dqf + A(:,k)*cos(wf*k*t) + B(:,k)*sin(wf*k*t);
        ddqf = ddqf - A(:,k)*wf*k*sin(wf*k*t) + B(:,k)*wf*k*cos(wf*k*t);
    end
    for i=1:5
        qp   =  qp  + C(:,i+1)*t.^i;
        dqp  = dqp  + C(:,i+1)*i*t.^(i-1);
        ddqp = ddqp + C(:,i+1)*i*(i-1)*t.^abs(i-2);
    end
    q   = qf + qp;
    dq  = dqf + dqp;
    ddq = ddqf + ddqp;
end
