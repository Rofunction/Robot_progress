%{
Note: J_dot(q,dq) is Analytical Jacobian martix;
      jac_dot=J_dot(q,dq) is the product (6x1) of the derivative of the 
    Jacobian (in the world frame) and the joint rates. i.e. jac_dot=J_dot*dq
      R{i} represents the transformation martix of the adjacent position
    and pose. 
      Q{i} represents rolation matrix of transformation martix R{i}.
      a{i} represents the translation vector of the transformation martix
   R{i}, namely R{i}(1:3,4).
       Rolat.m is stored in the folder Kinematics.
Reference: Fundamentals of Robotic Mechanical Systems Theory, Methods, and
Algorithms by Jorge Angeles. p229 and RTB jaco_dot.m
%}
function jac_dot=J_dot(q,dq,opt)
        n=length(q);
        dh.a=[0,-0.425,-0.39225,0,0,0];
        dh.d=[0.089459,0,0,0.10915,0.09465,0.0823];
        dh.alpha=[pi/2,0,0,pi/2,-pi/2,0];
        dh.theta=q;
        R=cell(1,6);
        for i=1:6
            R{i}=Rolat(i,dh);
            Q{i}=R{i}(1:3,1:3);
            a{i}=R{i}(1:3,4);
        end
        P{1}=Q{1};
        e{1}=[0,0,1]';
        for i=2:n
            P{i}=P{i-1}*Q{i};
            e{i}=P{i}(:,3);
        end
% step 1
        w{1}=dq(1)*e{1};
        for i=1:n-1
            w{i+1}=dq(i+1)*[0,0,1]'+Q{i}'*w{i};
        end
% step 2
        ed{1}=[0,0,0]';
        for i=2:n
        ed{i}=cross(w{i},e{i});
        end
% step 3 
        rd{n}=cross(w{n},a{n});
        for i=(n-1):-1:1
            rd{i}=cross(w{i},a{i})+Q{i}*rd{i+1};
        end
        
        r{n}=a{n};
        for i=(n-1):-1:1
            r{i}=a{i}+Q{i}*r{i+1};
        end
% step 4
       ud{1}=cross(e{1},rd{1});
       for i=2:n
           ud{i}=cross(ed{i},r{i}) + cross(e{i},rd{i});
       end
   if strcmp(opt,'Jdot')
% Solve J_dot(q)
       jac_dot=zeros(6,6);
       for i=1:6
       jac_dot(:,i)=[ud{i};ed{i}];
       end 
   end
   if strcmp(opt,'dJdq')
% Solve J_dot*dq
       v{n} = dq(n)*[ud{n}; ed{n}];
    for i=(n-1):-1:1
        Ui = blkdiag(Q{i}, Q{i});
        v{i} = dq(i)*[ud{i}; ed{i}] + Ui*v{i+1};
    end
    jac_dot = v{1};
   end
end