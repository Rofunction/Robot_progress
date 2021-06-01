%{
  --------------------------------------------------------------------------------
    The precedure is used for compute full regressor paramaters of Robot by
  Newton-Euler method;
    the Newton-Euler is refer to Christopher G. Atkeson.etc's article "Estimation 
  of Inertial Parameters of Manipulator Loads and Links" and specific
  derivation refers to "Newton-Euler identified parameters.pdf".
    Note: although the results of the method is different from the result of
    Euler-Lagrangian method, the effect is the same.
  --------------------------------------------------------------------------------
%}
function Y=full_regressor(q,dq,ddq)
% --------------------------------------------------------------------------------
% Default variable size  
% J_f(i) = Jf'(i); J_n(i) = Jn'(i); Jz(i)=z0'*R(i)*Jn(i)
% --------------------------------------------------------------------------------
    w=zeros(3,6); w_dot=zeros(3,6); v_dot=zeros(3,6); r=zeros(3,6);
    R=cell(1,6); T=cell(1,6); z0=[0,0,1]'; 
    J_f=cell(1,6); J_n=cell(1,6); Jf=cell(1,6); Jn=cell(1,6); Jz=cell(1,6);
    
% --------------------------------------------------------------------------------
% kinematics parameters of robot (Stadard D-H)
% --------------------------------------------------------------------------------
    dh.a=[0,-0.425,-0.39225,0,0,0];
    dh.d=[0.089459,0,0,0.10915,0.09465,0.0823];
    dh.alpha=[pi/2,0,0,pi/2,-pi/2,0];
    dh.theta=q;
% --------------------------------------------------------------------------------
% Solve J'(i)f and J'(i)n
% --------------------------------------------------------------------------------
    for i=1:6
        R{i}=Rolat(i,dh);
        r(:,i)=R{i}(1:3,1:3)' *R{i}(1:3,4);  % r(:,i) represents r_i(superscript)_(i-1,i)=p_i-p_i-1
        if i>1
            T{i}=R{i}*T{i-1}; 
            w(:,i) = R{i}(1:3,1:3)' *( w(:,i-1) + dq(i)*z0 );
            w_dot(:,i) = R{i}(1:3,1:3)' *( w_dot(:,i-1) + ddq(i)*z0 + dq(i)*cross(w(:,i),z0) );
            v_dot(:,i) = v_dot(:,i-1) + cross(w_dot(:,i),r(:,i)) + cross( w(:,i),cross( w(:,i),r(:,i) ) );
        else
            T{i}=R{i};
            w(:,i) = R{i}(1:3,1:3)' *dq(i)*z0;
            w_dot(:,i) = R{i}(1:3,1:3)' *( ddq(i)*z0 + dq(i)*cross(w(:,i),z0) );
            v_dot(:,i) = cross(w_dot(:,i),r(:,i)) + cross( w(:,i),cross( w(:,i),r(:,i) ) );
        end
        wX=func1( w(:,i) ); w_dotX=func1(w_dot(:,i)); v_dotX=func1(v_dot(:,i)); 
        w_dot_d = func2(w_dot(:,i)); w_d=func2( w(:,i) );
        J_f{i}  = [v_dot(:,i), w_dotX + wX*wX , zeros(3,6)];
        J_n{i}  = [zeros(3,1), -v_dotX, w_dot_d + wX*w_d];
    end
    
% --------------------------------------------------------------------------------
% Solve Jf(i) and Jn(i)
% --------------------------------------------------------------------------------
    for j=6:-1:1
        rX=func1(r(:,i));
        if j<6
            Jf{j}=[ J_f{j}, R{j+1}(1:3,1:3)*Jf{j+1} ];
            Jn{j}=[ J_n{j}, R{j+1}(1:3,1:3)*Jn{j+1} ] + rX *Jf{j};
        else
             Jf{6}=J_f{6};
             Jn{6}=J_n{6} + rX * Jf{6};
        end
    end
% --------------------------------------------------------------------------------
% Get full regressor 
% --------------------------------------------------------------------------------
    for i=1:6
        Jz{i}=z0' * R{i}(1:3,1:3) * Jn{i};
    end
    Y=[ Jz{1}; zeros(1,10),Jz{2}; zeros(1,20),Jz{3}; zeros(1,30),Jz{4};...
            zeros(1,40),Jz{5};zeros(1,50),Jz{6} ];
    
end

% --------------------------------------------------------------------------------
% Define a specific matrix (3X3)
% --------------------------------------------------------------------------------
function hX=func1(x)
[x1,x2,x3]=deal(x(1),x(2),x(3));
hX=[0,-x3,x2; x3,0,-x1; -x2,x1,0];
end

% --------------------------------------------------------------------------------
% Define a specific matrix (6X6)
% --------------------------------------------------------------------------------
function h_d=func2(x)
[x1,x2,x3]=deal(x(1),x(2),x(3));
h_d=[x1,x2,x3,zeros(1,3); 0,x1,0,x2,x3,0; 0,0,x1,0,x2,x3];
end