%{
    Velocity ([dx,dy,dz]) is the Velocities of end-effector
%}
function V = Velocity(t)
     eta=0.3;
       dx=-2*pi*eta*(2*sin(4*pi*t) *cos(2*pi*t)+cos(4*pi*t) *sin(2*pi*t));
       dy=-2*pi*eta*(2*sin(4*pi*t) *cos(2*pi*t)-cos(4*pi*t) *cos(2*pi*t));
       dz=-2*pi*eta*(2*sin(4*pi*t) *cos(2*pi*t)-cos(4*pi*t) *cos(2*pi*t));
       [dr,dp,d_gamma]=deal(0);
       V=[dx,dy,dz,dr,dp,d_gamma];
end 