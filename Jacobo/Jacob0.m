%{
    'df' represent Geometric Jacobian; 'rpy' represent Analytical Jacobian
%}
function Ja0=Jacob0(q,opt) 
        dh.a=[0,-0.425,-0.39225,0,0,0];
        dh.d=[0.089459,0,0,0.10915,0.09465,0.0823];
        dh.alpha=[pi/2,0,0,pi/2,-pi/2,0];
        dh.theta=q;
        R=cell(1,6);T=cell(1,6); P=zeros(1,6); R_w=zeros(3,3);
      for i=1:6
        R{i}=Rolat(i,dh);
          if i>1
            T{i}=T{i-1}*R{i};
          else
            T{i}=R{i};
          end
           z(:,i)=T{i}(1:3,3);
           p_e(:,i)=T{i}(1:3,4);
      end 
       z0=[0;0;1]; p0=[0;0;1];
       B(:,1)=cross(z0,p_e(:,6)-p0); 
       Ja0(:,1)=[B(:,1);z0];
      for i=2:6
          B(:,i)=cross(z(:,i-1),p_e(:,6)-p_e(:,i-1)); 
          Ja0(:,i)=[B(:,i);z(:,i-1)];
      end
      
      if strcmp(opt,'df')
          Ja0=Ja0;
      end 
      if strcmp(opt,'rpy')
          P=Fkine(q,'xyz');
          [p,y]=deal(P(:,5),P(:,6));
          R_w=[sin(p), 0, 1; -cos(p)*sin(y), cos(y), 0;
               cos(p)*cos(y), sin(y), 0];
           Ja0 = blkdiag( eye(3,3),inv(R_w) ) * Ja0;
          if rcond(R_w)< eps
               error('Jacobo is singularity');
          end     
      end
end