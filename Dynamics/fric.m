%{
   Bm is the coefficient of viscous friction and assumed that 
  the friction of the end effector is 
%}
function friction=fric(Bm,fc,dq)     %% Simple friction model
           dq=dq(:);     
           friction=Bm*dq+fc.*sign(dq);
           friction(6)=0;
end