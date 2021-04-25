%{
    Note: R=Rolation(i,dh) represents the transformation martix of the adjacent positon
    and pose.
%}
function R=Rolat(i,dh)  
cth=cos(dh.theta(i));
sth=sin(dh.theta(i));
ca=cos(dh.alpha(i));
sa=sin(dh.alpha(i));
d=dh.d(i); a=dh.a(i);
     R  =[cth, -sth*ca, sth*sa, a*cth;
          sth,  cth*ca,-cth*sa, a*sth;
           0 ,   sa   ,   ca  ,   d  ;
           0 ,   0    ,   0   ,   1 ];
end