function [mom]=f_john_mom(xx,type,ka3t,ka4t)
%
% Computes the Johnson distribution function (ul is the upper limit of the integral)
%
a=xx(1); b=xx(2); c=xx(3); d=xx(4);

if type==1
   ll=c; ul=50;
elseif type==2
   ll=-50; ul=50;
elseif type==3
   ll=c; ul=c+d;
elseif type==4
   ll=-10; ul=10;
end

mu =quad(@(x) x   .*f_john_dens(x,a,b,c,d,type),ll,ul);
si2=quad(@(x) x.^2.*f_john_dens(x,a,b,c,d,type),ll,ul);
ka3=quad(@(x) x.^3.*f_john_dens(x,a,b,c,d,type),ll,ul);
ka4=quad(@(x) x.^4.*f_john_dens(x,a,b,c,d,type),ll,ul);

mom=[mu-0; sqrt(si2); ka3; ka4;];