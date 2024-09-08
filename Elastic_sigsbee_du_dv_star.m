clear all
close all

load('bluewhitered.mat')

load('vp.mat')
load('vs.mat')
n = size(vp);

dx = 0.025; dz = 0.025;
h  = [dz dx];
nz = n(1); nx = n(2);
z  = [0:n(1)-1]'*h(1);
x  = [0:n(2)-1]*h(2);

figure;
pcolor(x,z,vp);
shading interp
axis ij
xlabel('Distance (km)','FontSize',12)
ylabel('Depth (km)','FontSize',12);
colormap(bluewhitered)
colorbar
set(gca,'FontSize',14)
caxis([1.5 3.0])

figure;
pcolor(x,z,vs);
shading interp
axis ij
xlabel('Distance (km)','FontSize',12)
ylabel('Depth (km)','FontSize',12);
colormap(bluewhitered)
colorbar
set(gca,'FontSize',14)
caxis([0.8 2])

vp0 = 1.83; vs0 = vp0/sqrt(3); 
rho = 1;
npmlz = 50; npmlx = 50; 
Nz = nz + 2*npmlz;
Nx = nx + 2*npmlx;
NN = (Nx)*(Nz);

src_x = 11:10:91;
src_z = 2; 
ns = length(src_x);

Ps1 = elastic_getP_new(n,npmlz,npmlx,src_z,src_x);Ps1 = Ps1'*16500;
Ps2 = zeros(size(Ps1));
Pss = [Ps2; Ps1];
Ps = -sparse(Pss);

fm = 8;
f = 5;
Ak  = subimpedance_new(vp,vs,npmlx,npmlz,f,fm,dx);
Uv  = Ak\(Ps);

vp0 = vp0*ones(n);
vs0 = vs0*ones(n);

A0  = subimpedance_new(vp0,vs0,npmlx,npmlz,f,fm,dx);
Uv0  = A0\(Ps);

for is = 1:ns
    
    u_star = Uv(1:NN,is);
    v_star = Uv(1+NN:end,is);
    
    U = reshape(u_star,[Nx,Nz]);
    V = reshape(v_star,[Nx,Nz]);
    
    U = U.'; V = V.';
    
    u = U(npmlz+1:npmlz+nz,npmlx+1:npmlx+nx);
    v = V(npmlz+1:npmlz+nz,npmlx+1:npmlx+nx);
    
    u0_star = Uv0(1:NN,is);
    v0_star = Uv0(1+NN:end,is);
    
    U0 = reshape(u0_star,[Nx,Nz]);
    V0 = reshape(v0_star,[Nx,Nz]);
    
    U0 = U0.'; V0 = V0.';
    
    u0 = U0(npmlz+1:npmlz+nz,npmlx+1:npmlx+nx);
    v0 = V0(npmlz+1:npmlz+nz,npmlx+1:npmlx+nx);
    
    du = u-u0; dv = v - v0; 
    
    du_star( ((is-1)*nz*nx+1) : (is*nz*nx) ,1) = full(du(:));
    dv_star( ((is-1)*nz*nx+1) : (is*nz*nx) ,1) = full(dv(:));
    
end

figure;
pcolor(x,z,real(du));
shading interp
axis ij
xlabel('Distance (km)','FontSize',12)
ylabel('Depth (km)','FontSize',12);
colormap(bluewhitered)
colorbar
caxis([-0.5 0.5])
set(gca,'FontSize',14)

figure;
pcolor(x,z,real(dv));
shading interp
axis ij
xlabel('Distance (km)','FontSize',12)
ylabel('Depth (km)','FontSize',12);
colormap(bluewhitered)
colorbar
caxis([-0.5 0.5])
set(gca,'FontSize',14)

save du_star_mar_5Hz_sx.mat du_star
save dv_star_mar_5Hz_sx.mat dv_star
