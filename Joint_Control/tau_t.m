%{
    Note:
%}
function tau=tau_t(t,q,dq,M_q,Cor,G_q)
    q_d   = sin(t)*ones(6,1) ;
    dq_d  = cos(t)*ones(6,1) ;
    ddq_d =-sin(t)*ones(6,1) ;
    e=q-q_d; de=dq-dq_d;
% Control parameters
    k=100; epi=5;
    s= de + k*e;
    tau=M_q*( ddq_d- k*de-epi*sign(s) ) + Cor * dq+ G_q;
end 
 