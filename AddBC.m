function [Vp,Vs]=AddBC(nz,nx,bn,vp,vs)
%给已经建好的地质模型镶边
%vp镶边：
tmpv=zeros(nz+2*bn,nx+2*bn);
tmpv(bn+1:nz+bn,bn+1:nx+bn)=vp;
tmpv(1:bn,bn+1:nx+bn)=ones(bn,1)*vp(1,:);
tmpv(nz+bn+1:nz+2*bn,bn+1:nx+bn)=ones(bn,1)*vp(nz,:);
tmpv(:,1:bn)=tmpv(:,bn+1)*ones(1,bn);
tmpv(:,nx+bn+1:nx+2*bn)=tmpv(:,nx+bn)*ones(1,bn);
Vp=tmpv;
%vs镶边：
tmpv=zeros(nz+2*bn,nx+2*bn);
tmpv(bn+1:nz+bn,bn+1:nx+bn)=vs;
tmpv(1:bn,bn+1:nx+bn)=ones(bn,1)*vs(1,:);
tmpv(nz+bn+1:nz+2*bn,bn+1:nx+bn)=ones(bn,1)*vs(nz,:);
tmpv(:,1:bn)=tmpv(:,bn+1)*ones(1,bn);
tmpv(:,nx+bn+1:nx+2*bn)=tmpv(:,nx+bn)*ones(1,bn);
Vs=tmpv;
