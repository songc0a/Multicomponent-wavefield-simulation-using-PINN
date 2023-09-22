function [u,v] = green_elastic(w,vp0,vs0,rho,x,z,xs,zs);

    F = 10;
    r = @(zz,xx)(zz.^2+xx.^2).^0.5;
   
    G1=-1i*pi/2*(1/vp0^2*besselh(0,2,(w*r(z-zs,x-xs)/(vp0)))+ ...
        1./(w*r(z-zs,x-xs)*vs0).*besselh(1,2,(w*r(z-zs,x-xs)/(vs0)))- ...
        1./(w*r(z-zs,x-xs)*vp0).*besselh(1,2,(w*r(z-zs,x-xs)/(vp0))));
    
    G2=1i*pi/2*(1/vs0^2*besselh(0,2,(w*r(z-zs,x-xs)/(vs0)))- ...
        1./(w*r(z-zs,x-xs)*vs0).*besselh(1,2,(w*r(z-zs,x-xs)/(vs0)))+ ...
        1./(w*r(z-zs,x-xs)*vp0).*besselh(1,2,(w*r(z-zs,x-xs)/(vp0))));
    
    u=F/(2*pi*rho)*((x-xs).*(z-zs)./r(z-zs,x-xs).^2).*(G1+G2);
    v=F/(2*pi*rho)*(1./r(z-zs,x-xs).^2).*((z-zs).^2.*G1-(x-xs).^2.*G2);
    
end