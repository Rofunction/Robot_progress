%{
  Note : The friction tauque is not considered // fodyn is a abbreviation of forward dynamics 
         q/dq are vectors with (1xn)
%}
function ddq=Vq(q,dq,tau)
    dq=dq(:);
    grav=[0,0,-9.81];
    M=Inertia(q);
    G=Grav(q,grav).';     %(1xn)
    Cor=Coriolis(q,dq);
    
    ddq=M\(tau.'- Cor*dq- G);
    ddq=ddq.';
end
    
    
    