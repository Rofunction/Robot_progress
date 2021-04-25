%{
    Desired trajectory ([px_d, py_d, pz_d]) is the trajectoties of
  end-effector. The trajectory is Rose pattern. this method is only solve
  qd in the time domain. dqd ,ddqd are suitable for this method.
%}
function [buff,qd]=DesiredTrj_q(buff,t)
         q_d=[]; 
         trj=DesiredTrj(t);
         q_d=ur5ik(trj);
         qd=buff;                        %% 这个buff在t=0时，无意义  
         
        if t==0
            qd=q_d(1,:);                 %% 期望轨迹的起始点   
        elseif t>0
            dim_n=size(q_d,1);
            D_e=zeros(dim_n,6);          %% 每次重新定义误差矩阵维度
            D_a = zeros(dim_n,1);
            for i=1:dim_n
                D_e(i,:)=q_d(i,:)-qd;
                D_a(i,:)=sum ( abs(D_e(i,:)) );
            end
            [~,sort_i]=sort(D_a);
            qd=q_d(sort_i(1),:);
        end
        buff=qd;
end

function    Trj=DesiredTrj(t)            %% [px_d, py_d, pz_d, r, p, y]=DesiredTrj_q(t)
            x0=0.3 ; y0=0.3; z0=0.3; eta=0.3;
            px_d= eta * cos(4*pi*t) * cos(2*pi*t) + x0;
            py_d= eta * cos(4*pi*t) * sin(2*pi*t) + y0;
            pz_d= eta * cos(4*pi*t) * sin(2*pi*t) + z0;
            [r,p,gamma]=deal(0);
            Trj=[px_d, py_d, pz_d, r, p, gamma];  
            
end 