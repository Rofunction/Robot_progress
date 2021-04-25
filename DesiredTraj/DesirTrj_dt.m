%{
    Note: dqd=\delta (qd)/ \delta t;
          ddqd=\delta (dqd)/ \delta t;
    the precision of this approach is lower than the Gradient derivation.
    this desired trajectories is solved in advance.
%}
function [qd,dqd,ddqd]=DesirTrj_dt(ts)
%         ts=linspace(0,10,10000);
        qd=zeros(10000,6); q_d={};
        
        for i=1:length(ts)
            t=ts(i);
         trj=DesiredTrj(t);
         q_d1={ur5ik(trj)};
         q_d(i,1)=q_d1;                      
        end
        qd(1,:)=q_d{1,1}(1,:);                                %% 期望轨迹的起始点 
        dqd(1,:)=zeros(1,6); ddqd(1,:)=zeros(1,6);       
        for j=2:length(ts)
            dim_n=size(q_d{j,:},1);
            D_e=zeros(dim_n,6);                               %% 每次重新定义误差矩阵维度
            D_a = zeros(dim_n,1);
            for i=1:dim_n
                D_e(i,:)=q_d{j,1}(i,:)-qd(j-1,:);
                D_a(i,:)=sum ( abs(D_e(i,:)) );
            end
            [~,sort_i]=sort(D_a);
            qd(j,:)=q_d{j,1}(sort_i(1),:); 
            dqd(j,:)=(qd(j,:)-qd(j-1,:))/(ts(j)-ts(j-1));
            ddqd(j,:)=(dqd(j,:)-dqd(j-1,:))/(ts(j)-ts(j-1));
            disp(['running at =',num2str(j/10000*100),'%']);
        end
end
function    Trj=DesiredTrj(t)            %% [px_d, py_d, pz_d, r, p, y]=DesiredTrj_q(t)
            x0=0.3 ; y0=0.3; z0=0.3; eta=0.3;
            px_d= eta * cos(4*pi*t) * cos(2*pi*t) + x0;
            py_d= eta * cos(4*pi*t) * sin(2*pi*t) + y0;
            pz_d= eta * cos(4*pi*t) * sin(2*pi*t) + z0;
            [r,p,gamma]=deal(0);
            Trj=[px_d, py_d, pz_d, r, p, gamma];  
            
end 