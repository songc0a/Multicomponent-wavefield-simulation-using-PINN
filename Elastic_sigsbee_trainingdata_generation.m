clear all
close all

load('bluewhitered.mat')

str2 =['vpt'];
filename2=['' str2 '.rsf@'];
fid2=fopen(filename2,'rb');
vel=fread(fid2,[143,465],'float');
fclose(fid2);

vp=vel(1:1:101,1:1:101)/1000 ;
vs = vp/sqrt(3);

n  = size(vp); N = n(1)*n(2);

vp = imgaussfilt(vp,1);
vs = imgaussfilt(vs,1);

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

N_train = 50000; %% number of the random points

z_train = 2.45*rand(N_train,1)+0.025;
x_train = 2.45*rand(N_train,1)+0.025;
sx_train = 2.475*rand(N_train,1)+0.025;

vp0 = 1.83; vs0 = vp0/sqrt(3); 
rho = 1;

src_z = 2; 
zs = (src_z-1)*dz; %% depth of the source

z  = [0:n(1)-1]'*h(1);
x  = [0:n(2)-1]*h(2);
[X,Y] = meshgrid(x,z);
x1 = [0:2501-1]*0.001;
z1 = [0:2501-1]'*0.001;
[Xq,Yq] = meshgrid(x1,z1);
vp_in = interp2(x,z,vp,Xq,Yq);
vs_in = interp2(x,z,vs,Xq,Yq);

xx_in = round(x_train/0.001)+1;
zz_in = round(z_train/0.001)+1;

alpha_train = zeros(N_train,1);
beta_train = zeros(N_train,1);

dxxu0_train = zeros(N_train,1);
dzzu0_train = zeros(N_train,1);
dxzu0_train = zeros(N_train,1);

dzzv0_train = zeros(N_train,1);
dxxv0_train = zeros(N_train,1);
dxzv0_train = zeros(N_train,1);

F = 1;
f = 4;

for is = 1:N_train
  
    w = 1*2*pi*f;
    
    x0 = x_train(is); z0 = z_train(is); xs0 = sx_train(is); 
    x01 = x_train(is)+dx; z01 = z_train(is)+dz; 
    x_1 = x_train(is)-dx; z_1 = z_train(is)-dz; 
    
    % position x, z   
    [u0,v0] = green_elastic(w,vp0,vs0,rho,x0,z0,xs0,zs);
    
    % position x-1, z
    [u0_x_1,v0_x_1] = green_elastic(w,vp0,vs0,rho,x_1,z0,xs0,zs);    
    
    % position x+1,z
    [u0_x01,v0_x01] = green_elastic(w,vp0,vs0,rho,x01,z0,xs0,zs);      
    
    % position x, z-1
    [u0_z_1,v0_z_1] = green_elastic(w,vp0,vs0,rho,x0,z_1,xs0,zs);        
    
    % position x, z+1
    [u0_z01,v0_z01] = green_elastic(w,vp0,vs0,rho,x0,z01,xs0,zs);         

    % position x-1, z-1
    [u0_x_1_z_1,v0_x_1_z_1] = green_elastic(w,vp0,vs0,rho,x_1,z_1,xs0,zs);  
    
    % position x+1, z-1
    [u0_x01_z_1,v0_x01_z_1] = green_elastic(w,vp0,vs0,rho,x01,z_1,xs0,zs);  
 
    % position x-1, z+1
    [u0_x_1_z01,v0_x_1_z01] = green_elastic(w,vp0,vs0,rho,x_1,z01,xs0,zs);  
     
    % position x+1, z+1
    [u0_x01_z01,v0_x01_z01] = green_elastic(w,vp0,vs0,rho,x01,z01,xs0,zs);  
    
    dxxu0_train(is,1) = (u0_x_1 - 2*u0 + u0_x01)/dx/dx;
    dzzu0_train(is,1) = (u0_z_1 - 2*u0 + u0_z01)/dz/dz;
    dxzu0_train(is,1) = (u0_x_1_z_1 - u0_x01_z_1 - u0_x_1_z01 + u0_x01_z01)/(4*dx*dz);
    
    dzzv0_train(is,1) = (v0_z_1 - 2*v0 + v0_z01)/dz/dz;
    dxxv0_train(is,1) = (v0_x_1 - 2*v0 + v0_x01)/dx/dx;
    dxzv0_train(is,1) = (v0_x_1_z_1 - v0_x01_z_1 - v0_x_1_z01 + v0_x01_z01)/(4*dx*dz);
    
    alpha_train(is,1) = vp_in(zz_in(is),xx_in(is))^2;
    beta_train(is,1) = vs_in(zz_in(is),xx_in(is))^2;
    
end

alpha0_train = vp0^2 * ones(N_train,1);
beta0_train = vs0^2 * ones(N_train,1);

z  = [0:n(1)-1]'*h(1);
x  = [0:n(2)-1]*h(2);
sx = (10:10:90)*dx; ns = length(sx);

[zz,xx] = ndgrid(z,x);
sx = repmat(sx,nx*nz,1);

x1 = xx(:); x_star = (repmat(x1,ns,1)); 
z1 = zz(:); z_star = (repmat(z1,ns,1));
sx_star = sx(:);

save Elastic_sigsbee_4Hz_sx_F10.mat x_star z_star sx_star ...
    x_train sx_train z_train  ...
    alpha_train alpha0_train beta_train beta0_train ...
    dzzu0_train dxxu0_train dxzu0_train dzzv0_train dxxv0_train dxzv0_train 

