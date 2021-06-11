function [fv]=f_john_dens(x,a,b,c,d,type)
%
%  Computes the Johnson density
%

u=(x-c)./d;
gp=fgp(u,type); 
g=fg(u,type);

fv=b./(d.*sqrt(2*pi)).*gp.*exp(-0.5.*((a+b.*g).^2));

function gp=fgp(u,type)
if type==1
   gp=1./u;
elseif type==2
   gp=1./sqrt(u.^2+1);
elseif type==3
   gp=-1./(u.*(u-1));
elseif type==4
   gp=ones(size(u));
end

function g=fg(u,type)
if type==1
   g=log(u);
elseif type==2
   g=log(u+sqrt(u.^2+1));
elseif type==3
   g=log(u./(1-u));
elseif type==4
   g=u;
end