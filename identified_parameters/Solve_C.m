%{
    Note: Coefficient C is a function about A,B,q_init, i.e. C=F(A,B,q_init);
    the specific solution of C is obtained based on the constraints conditions of the
   robotic manipulators. i.e.:
    1、q(0)   =q_init , q(T)   =q_init;
    2、q'(0)  =0      , q'(T)  =0;
    3、q''(0) =0      , q''(T) =0;
    T= 2Pi/wf;
    wf - fundamental frequecncy
    the specific trajector refers to the follow articles:
    1、Closed-loop Dynamic Parameter Identification of Robot Manipulators Using Modified Fourier Series,
 International Journal of Advanced Robotic Systems; DOI: 10.5772/45818
    2、Practical Aspects of Model-Based Collision Detection, doi: 10.3389/frobt.2020.571574
%}
function C=Solve_C(T,A,B,N,wf,q_init)
    % q_init = traj.q0;
    Coe_M =zeros(36,36); 
    q0   =  sum(B./ (1:N)*wf ,2);
    dq0  = -sum(A,2);
    ddq0 = -sum(B.* (1:N)*wf ,2);
    O_6  = zeros(6); I_6= eye(6);
    Coe_M(1:6,  :)=[ I_6, O_6, O_6, O_6, O_6, O_6 ];
    Coe_M(7:12, :)=[ O_6, I_6, O_6, O_6, O_6, O_6 ];
    Coe_M(13:18,:)=[ O_6, O_6, 2*I_6, O_6, O_6, O_6 ];
    Coe_M(19:24,:)=[ I_6, T*I_6, T^2*I_6, T^3*I_6, T^4*I_6, T^5*I_6 ];
    Coe_M(25:30,:)=[ O_6, I_6, T*2*I_6, 3*T^2*I_6, 4*T^3*I_6, 5*T^4*I_6 ];
    Coe_M(31:36,:)=[ O_6, O_6, 2*I_6, 6*T*I_6, 12*T^2*I_6, 20*T^3*I_6 ];

    SolveC = Coe_M \ [ q_init+q0; dq0; ddq0; q_init+q0; dq0; ddq0;];
    C = reshape (SolveC,[6,6]);
end