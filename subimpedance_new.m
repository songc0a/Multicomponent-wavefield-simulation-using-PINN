% evaluate the impedance matrix
function [Asparse]=subimpedance_new(vp,vs,npmlx,npmlz,w,fm,h)
[nz,nx]=size(vp);
npml=npmlx;
[Vp,Vs]=AddBC(nz,nx,npml,vp,vs);
% miu=Rho.*Vs.^2;
% lamda=Rho.*Vp.^2-2*Rho.*Vs.^2;
Nz=nz+2*npmlz;
Nx=nx+2*npmlx;


dx=h;
dz=h;
Fw=w;
Fw=Fw*2*pi;
%--------------设置PML吸收边界------------
%衰减系数设置
a0=1.79;
% a0=0.1;
ddx=zeros(1,Nx);
ddx(1:npml)=2*pi*a0*fm*((npml:-1:1)/npml).^2;
ddx(Nx-npml+1:Nx)=2*pi*a0*fm*((1:1:npml)/npml).^2;
ex=1./(1+1i*ddx/Fw);

ddz=zeros(Nz,1);
ddz(1:npml)=2*pi*a0*fm*((npml:-1:1)/npml).^2;
ddz(Nz-npml+1:Nz)=2*pi*a0*fm*((1:1:npml)/npml).^2;
ez=1./(1+1i*ddz/Fw);

ex=ones(Nz,1)*ex;
ez=ez*ones(1,Nx);

a1=0.4798;
a2=0.1517;
a3=-0.0249;
b1=0.5580;
b2=0.2178;
N=0;
for iz=1:Nz%方程一
    for ix=1:Nx
        row=(iz-1)*Nx+ix;
        column=row;
        N=N+1;
        Row(N,1)=row;
        Column(N,1)=column;
        A(N,1)=Fw^2*a1+Vp(iz,ix)*Vp(iz,ix)*(-2)*b1/(ex(iz,ix)^2*dx^2)+Vs(iz,ix)*Vs(iz,ix)*(-2)*b1/(ez(iz,ix)^2*dz^2);
        
        if ix+1<=Nx
            N=N+1;
            column=(iz-1)*Nx+ix+1;
            Row(N,1)=row;
            Column(N,1)=column;
            A(N,1)=Fw^2*a2+Vp(iz,ix)*Vp(iz,ix)*b1/(ex(iz,ix)^2*dx^2)+Vs(iz,ix)*Vs(iz,ix)*(-2)*b2/(ez(iz,ix)^2*dz^2);
        end
        
        if iz-1>0
            N=N+1;
            column=(iz-1-1)*Nx+ix;
            Row(N,1)=row;
            Column(N,1)=column;
            A(N,1)=Fw^2*a2+Vp(iz,ix)*Vp(iz,ix)*(-2)*b2/(ex(iz,ix)^2*dx^2)+Vs(iz,ix)*Vs(iz,ix)*b1/(ez(iz,ix)^2*dz^2);
        end
        
        if ix-1>0
            N=N+1;
            column=(iz-1)*Nx+ix-1;
            Row(N,1)=row;
            Column(N,1)=column;
            A(N,1)=Fw^2*a2+Vp(iz,ix)*Vp(iz,ix)*b1/(ex(iz,ix)^2*dx^2)+Vs(iz,ix)*Vs(iz,ix)*(-2)*b2/(ez(iz,ix)^2*dz^2);
        end
        
        if iz+1<=Nz
            N=N+1;
            column=(iz-1+1)*Nx+ix;
            Row(N,1)=row;
            Column(N,1)=column;
            A(N,1)=Fw^2*a2+Vp(iz,ix)*Vp(iz,ix)*(-2)*b2/(ex(iz,ix)^2*dx^2)+Vs(iz,ix)*Vs(iz,ix)*b1/(ez(iz,ix)^2*dz^2);
        end
        
        if ix+1<=Nx && iz+1<=Nz
            N=N+1;
            column=(iz-1+1)*Nx+ix+1;
            Row(N,1)=row;
            Column(N,1)=column;
            A(N,1)=Fw^2*a3+Vp(iz,ix)*Vp(iz,ix)*b2/(ex(iz,ix)^2*dx^2)+Vs(iz,ix)*Vs(iz,ix)*b2/(ez(iz,ix)^2*dz^2);
            
            N=N+1;
            column=(iz-1+1)*Nx+ix+1+Nz*Nx;
            Row(N,1)=row;
            Column(N,1)=column;
            A(N,1)=(Vp(iz,ix)*Vp(iz,ix)-Vs(iz,ix)*Vs(iz,ix))/(ex(iz,ix)*ez(iz,ix)*4*dx*dz);
        end
        
        if ix+1<=Nx && iz-1>0
            N=N+1;
            column=(iz-1-1)*Nx+ix+1;
            Row(N,1)=row;
            Column(N,1)=column;
            A(N,1)=Fw^2*a3+Vp(iz,ix)*Vp(iz,ix)*b2/(ex(iz,ix)^2*dx^2)+Vs(iz,ix)*Vs(iz,ix)*b2/(ez(iz,ix)^2*dz^2);
            
            N=N+1;
            column=(iz-1-1)*Nx+ix+1+Nz*Nx;
            Row(N,1)=row;
            Column(N,1)=column;
            A(N,1)=-(Vp(iz,ix)*Vp(iz,ix)-Vs(iz,ix)*Vs(iz,ix))/(ex(iz,ix)*ez(iz,ix)*4*dx*dz);
        end
        
        if ix-1>0 && iz-1>0
            N=N+1;
            column=(iz-1-1)*Nx+ix-1;
            Row(N,1)=row;
            Column(N,1)=column;
            A(N,1)=Fw^2*a3+Vp(iz,ix)*Vp(iz,ix)*b2/(ex(iz,ix)^2*dx^2)+Vs(iz,ix)*Vs(iz,ix)*b2/(ez(iz,ix)^2*dz^2);
            
            N=N+1;
            column=(iz-1-1)*Nx+ix-1+Nz*Nx;
            Row(N,1)=row;
            Column(N,1)=column;
            A(N,1)=(Vp(iz,ix)*Vp(iz,ix)-Vs(iz,ix)*Vs(iz,ix))/(ex(iz,ix)*ez(iz,ix)*4*dx*dz);
        end
        
        if ix-1>0 && iz+1<=Nz
            N=N+1;
            column=(iz-1+1)*Nx+ix-1;
            Row(N,1)=row;
            Column(N,1)=column;
            A(N,1)=Fw^2*a3+Vp(iz,ix)*Vp(iz,ix)*b2/(ex(iz,ix)^2*dx^2)+Vs(iz,ix)*Vs(iz,ix)*b2/(ez(iz,ix)^2*dz^2);
            
            N=N+1;
            column=(iz-1+1)*Nx+ix-1+Nz*Nx;
            Row(N,1)=row;
            Column(N,1)=column;
            A(N,1)=-(Vp(iz,ix)*Vp(iz,ix)-Vs(iz,ix)*Vs(iz,ix))/(ex(iz,ix)*ez(iz,ix)*4*dx*dz);
        end
    end
end

for iz=1:Nz%方程二
    for ix=1:Nx
        row=(iz-1)*Nx+ix+Nz*Nx;
        column=row;
        N=N+1;
        Row(N,1)=row;
        Column(N,1)=column;
        A(N,1)=Fw^2*a1+Vp(iz,ix)*Vp(iz,ix)*(-2)*b1/(ez(iz,ix)^2*dz^2)+Vs(iz,ix)*Vs(iz,ix)*(-2)*b1/(ex(iz,ix)^2*dx^2);
        
        if ix+1<=Nx
            N=N+1;
            column=(iz-1)*Nx+ix+1+Nz*Nx;
            Row(N,1)=row;
            Column(N,1)=column;
            A(N,1)=Fw^2*a2+Vp(iz,ix)*Vp(iz,ix)*(-2)*b2/(ez(iz,ix)^2*dz^2)+Vs(iz,ix)*Vs(iz,ix)*b1/(ex(iz,ix)^2*dx^2);
        end
        
        if iz-1>0
            N=N+1;
            column=(iz-1-1)*Nx+ix+Nz*Nx;
            Row(N,1)=row;
            Column(N,1)=column;
            A(N,1)=Fw^2*a2+Vp(iz,ix)*Vp(iz,ix)*b1/(ez(iz,ix)^2*dz^2)+Vs(iz,ix)*Vs(iz,ix)*(-2)*b2/(ex(iz,ix)^2*dx^2);
        end
        
        if iz+1<=Nz
            N=N+1;
            column=(iz-1+1)*Nx+ix+Nz*Nx;
            Row(N,1)=row;
            Column(N,1)=column;
            A(N,1)=Fw^2*a2+Vp(iz,ix)*Vp(iz,ix)*b1/(ez(iz,ix)^2*dz^2)+Vs(iz,ix)*Vs(iz,ix)*(-2)*b2/(ex(iz,ix)^2*dx^2);
        end
        
        if ix-1>0
            N=N+1;
            column=(iz-1)*Nx+ix-1+Nz*Nx;
            Row(N,1)=row;
            Column(N,1)=column;
            A(N,1)=Fw^2*a2+Vp(iz,ix)*Vp(iz,ix)*(-2)*b2/(ez(iz,ix)^2*dz^2)+Vs(iz,ix)*Vs(iz,ix)*b1/(ex(iz,ix)^2*dx^2);
        end
        
        if ix+1<=Nx && iz+1<=Nz
            N=N+1;
            column=(iz-1+1)*Nx+ix+1+Nz*Nx;
            Row(N,1)=row;
            Column(N,1)=column;
            A(N,1)=Fw^2*a3+Vp(iz,ix)*Vp(iz,ix)*b2/(ez(iz,ix)^2*dz^2)+Vs(iz,ix)*Vs(iz,ix)*b2/(ex(iz,ix)^2*dx^2);
            
            N=N+1;
            column=(iz-1+1)*Nx+ix+1;
            Row(N,1)=row;
            Column(N,1)=column;
            A(N,1)=(Vp(iz,ix)*Vp(iz,ix)-Vs(iz,ix)*Vs(iz,ix))/(4*dx*dz*ex(iz,ix)*ez(iz,ix));
        end
        
        if ix+1<=Nx && iz-1>0
            N=N+1;
            column=(iz-1-1)*Nx+ix+1+Nz*Nx;
            Row(N,1)=row;
            Column(N,1)=column;
            A(N,1)=Fw^2*a3+Vp(iz,ix)*Vp(iz,ix)*b2/(ez(iz,ix)^2*dz^2)+Vs(iz,ix)*Vs(iz,ix)*b2/(ex(iz,ix)^2*dx^2);
            
            N=N+1;
            column=(iz-1-1)*Nx+ix+1;
            Row(N,1)=row;
            Column(N,1)=column;
            A(N,1)=-(Vp(iz,ix)*Vp(iz,ix)-Vs(iz,ix)*Vs(iz,ix))/(4*dx*dz*ex(iz,ix)*ez(iz,ix));
        end
        
        if ix-1>0 && iz+1<=Nz
            N=N+1;
            column=(iz-1+1)*Nx+ix-1+Nz*Nx;
            Row(N,1)=row;
            Column(N,1)=column;
            A(N,1)=Fw^2*a3+Vp(iz,ix)*Vp(iz,ix)*b2/(ez(iz,ix)^2*dz^2)+Vs(iz,ix)*Vs(iz,ix)*b2/(ex(iz,ix)^2*dx^2);
            
            N=N+1;
            column=(iz-1+1)*Nx+ix-1;
            Row(N,1)=row;
            Column(N,1)=column;
            A(N,1)=-(Vp(iz,ix)*Vp(iz,ix)-Vs(iz,ix)*Vs(iz,ix))/(4*dx*dz*ex(iz,ix)*ez(iz,ix));
        end
        
        if ix-1>0 && iz-1>0
            N=N+1;
            column=(iz-1-1)*Nx+ix-1+Nz*Nx;
            Row(N,1)=row;
            Column(N,1)=column;
            A(N,1)=Fw^2*a3+Vp(iz,ix)*Vp(iz,ix)*b2/(ez(iz,ix)^2*dz^2)+Vs(iz,ix)*Vs(iz,ix)*b2/(ex(iz,ix)^2*dx^2);
            
            N=N+1;
            column=(iz-1-1)*Nx+ix-1;
            Row(N,1)=row;
            Column(N,1)=column;
            A(N,1)=(Vp(iz,ix)*Vp(iz,ix)-Vs(iz,ix)*Vs(iz,ix))/(4*dx*dz*ex(iz,ix)*ez(iz,ix));
        end
    end
end

Asparse=sparse(Row,Column,A,2*Nz*Nx,2*Nz*Nx);
clear Row Column A ex ez

% [L,U]=lu(Asparse);