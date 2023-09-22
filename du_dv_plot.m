clear all
close all

load('bluewhitered.mat')
load('misfit_fre_elastic_4hz_sx.mat');
load('misfit1_fre_elastic_4hz_sx.mat');

mis = [misfit misfit1(3:end)];
nn = length(mis);

figure;
semilogy(1:nn, mis,'LineWidth',1.5);
xlabel('Epoch','FontSize',12)
ylabel('Loss','FontSize',12);
set(gca,'FontSize',14)

nz = 101; nx = 101; n = [nz,nx];
dx = 25; dz = dx; h  = [dz dx];
z  = [0:n(1)-1]'*h(1)/1000;
x  = [0:n(2)-1]*h(2)/1000;

amp = 0.5;

load('du_pred_atan_4hz_fre_elastic_sx.mat')
load('dv_pred_atan_4hz_fre_elastic_sx.mat')

is = 3;

du_pred = imag(du_pred( ((is-1)*nz*nx+1) : (is*nz*nx) ,1));
dv_pred = imag(dv_pred( ((is-1)*nz*nx+1) : (is*nz*nx) ,1));

figure;
pcolor(x,z,reshape(du_pred,n));
shading interp
axis ij
colorbar; colormap(bluewhitered)
% xlim([0 2]);ylim([0 2])
caxis([-amp amp]);
xlabel('Distance (km)','FontSize',12)
ylabel('Depth (km)','FontSize',12);
set(gca,'FontSize',14)

figure;
pcolor(x,z,reshape(dv_pred,n))
shading interp
axis ij
colorbar; colormap(bluewhitered)
caxis([-amp amp]);
xlabel('Distance (km)','FontSize',12)
ylabel('Depth (km)','FontSize',12);
set(gca,'FontSize',14)

load('du_star_sigsbee_4Hz_sx1.mat')
load('dv_star_sigsbee_4Hz_sx1.mat')

du_star = imag(du_star( ((is-1)*nz*nx+1) : (is*nz*nx) ,1));
dv_star = imag(dv_star( ((is-1)*nz*nx+1) : (is*nz*nx) ,1));

figure;
pcolor(x,z,reshape(du_star,n));
shading interp
axis ij
colorbar; colormap(bluewhitered)
% xlim([0 2]);ylim([0 2])
caxis([-amp amp]);
xlabel('Distance (km)','FontSize',12)
ylabel('Depth (km)','FontSize',12);
set(gca,'FontSize',14)

figure;
pcolor(x,z,reshape(dv_star,n))
shading interp
axis ij
colorbar; colormap(bluewhitered)
% xlim([0 2]);ylim([0 2])
caxis([-amp amp]);
xlabel('Distance (km)','FontSize',12)
ylabel('Depth (km)','FontSize',12);
% title('True scattered wavefield');
set(gca,'FontSize',14)

% amp = 0.01;
figure;
pcolor(x,z,reshape(du_star-du_pred,n));
shading interp
axis ij
colorbar; colormap(bluewhitered)
% xlim([0 2]);ylim([0 2])
xlabel('Distance (km)','FontSize',12)
ylabel('Depth (km)','FontSize',12);
% title('Wavefield difference');
caxis([-amp amp]);
set(gca,'FontSize',14)

figure;
pcolor(x,z,reshape(dv_star-dv_pred,n))
shading interp
axis ij
colorbar; colormap(bluewhitered)
% xlim([0 2]);ylim([0 2])
xlabel('Distance (km)','FontSize',12)
ylabel('Depth (km)','FontSize',12);
% title('Wavefield difference');
caxis([-amp amp]);
set(gca,'FontSize',14)