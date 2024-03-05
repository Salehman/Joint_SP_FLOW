function [Pwell, Pwell2]=FLOW_3D_1(NX,NY,NZ,L,W,H,dt,NT,Q,Well_L,K_x,K_y,K_z,phi,Cc,Mu,ob)
%for parallel

TN=NX*NY*NZ; %Total number of grids

Kxx=permute(reshape(K_x,[NY,NX,NZ]),[2 1 3]);
Kyy=permute(reshape(K_y,[NY,NX,NZ]),[2 1 3]);
Kzz=permute(reshape(K_z,[NY,NX,NZ]),[2 1 3]);
    Kyy=Kxx;
    Kzz=Kxx;
% Zdat=0;   %The depth of datum

dx=L/NX;
dz=H/NZ;
dy=W/NY;
% dt=1;
% Kxx=Kxx*milli*darcy;
% Kyy=Kyy*milli*darcy;
% Kzz=Kzz*milli*darcy;


CF=0.001127;
 CF=1;



DeltaX=ones(NX,NY,NZ)*dx;
DeltaY=ones(NX,NY,NZ)*dy;
DeltaZ=ones(NX,NY,NZ)*dz;


% x=dx:dx:L;z=0:dz:H-dz;y=dy:dy:W; 

% 
% [X,Y,Z] = meshgrid(x,y,z);
% 
%  X=permute(X,[2 1 3]);
%  Y=permute(Y,[2 1 3]);
% Z=permute(Z,[2 1 3]);
% Zdiff=Z-Zdat;


% 
%  Kxx=ones(NX,NY,NZ)*10;
%  Kyy=ones(NX,NY,NZ)*10;
%  Kzz=ones(NX,NY,NZ)*50;
% Kxx=SPEKx(1:NX,1:NY,1:NZ);
% Kyy=SPEKy(1:NX,1:NY,1:NZ);
% Kzz=SPEKz(1:NX,1:NY,1:NZ);
%K(:,:,2:end)=00.0000001;
KX=zeros(NX+1,NY,NZ);
KY=zeros(NX,NY+1,NZ);
KZ=zeros(NX,NY,NZ+1);



% phi=ones(NX,NY,NZ)*0.5;
%  phi=SPEPhi(1:NX,1:NY,1:NZ);
Vp=(dx*dy*dz).*phi./1;
% Vp=(dx*dy*dz).*phi./1;
VPV=reshape(permute(Vp,[2 1 3]), [TN,1,1]);
 Ct=spdiags((Cc*ones(1,TN).*(VPV'/dt))',0,TN,TN);
%%
% % La=K./Mu;
% 
% Zall=zeros(NX+2,NY+2,NZ+2);
% 
% Zall(2:end-1,2:end-1,2:end-1)=Z;
% Zall(:,1,:)=Zall(:,2,:);Zall(:,end,:)=Zall(:,end-1,:);Zall(1,:,:)=Zall(2,:,:);
% Zall(end,:,:)=Zall(end-1,:,:);Zall(:,:,1)=Zall(:,:,2);Zall(:,:,end)=Zall(:,:,end-1);
% 
% ZX=zeros(NX+1,NY,NZ);
% ZY=zeros(NX,NY+1,NZ);
% ZZ=zeros(NX,NY,NZ+1);
% 
% 
% 
% 
% ZX(:,:,:)=Zall(2:end,2:end-1,2:end-1)-Zall(1:end-1,2:end-1,2:end-1)+Zdat; 
% 
% ZY(:,:,:)=Zall(2:end-1,2:end,2:end-1)-Zall(2:end-1,1:end-1,2:end-1)+Zdat;                                                                                                                                                                                      ;
% 
% ZZ(:,:,:)=Zall(2:end-1,2:end-1,2:end)-Zall(2:end-1,2:end-1,1:end-1)+Zdat;                                                                                                                                                                                      ;
% 
% 
% 
% 
% ZXm(:,:,:)=Zall(2:end,2:end-1,2:end-1)-Zall(1:end-1,2:end-1,2:end-1)+Zdat; 
% 
% ZY(:,:,:)=Zall(2:end-1,2:end,2:end-1)-Zall(2:end-1,1:end-1,2:end-1)+Zdat;                                                                                                                                                                                      ;
% 
% ZZ(:,:,:)=Zall(2:end-1,2:end-1,2:end)-Zall(2:end-1,2:end-1,1:end-1)+Zdat;  



KX(2:NX,:,:)=((DeltaX(2:end,:,:)+DeltaX(1:end-1,:,:)).*(Kxx(1:end-1,:,:).*Kxx(2:end,:,:)))...
    ./(DeltaX(2:end,:,:).*Kxx(1:end-1,:,:)+DeltaX(1:end-1,:,:).*Kxx(2:end,:,:));
    
KY(:,2:NY,:)=((DeltaY(:,2:end,:)+DeltaY(:,1:end-1,:)).*(Kyy(:,1:end-1,:).*Kyy(:,2:end,:)))...
    ./(DeltaY(:,2:end,:).*Kyy(:,1:end-1,:)+DeltaY(:,1:end-1,:).*Kyy(:,2:end,:));


KZ(:,:,2:NZ)=((DeltaZ(:,:,2:end)+DeltaZ(:,:,1:end-1)).*(Kzz(:,:,1:end-1).*Kzz(:,:,2:end)))...
    ./(DeltaZ(:,:,2:end).*Kzz(:,:,1:end-1)+DeltaZ(:,:,1:end-1).*Kzz(:,:,2:end));

TX=zeros(NX+1,NY,NZ);
TY=zeros(NX,NY+1,NZ);
TZ=zeros(NX,NY,NZ+1);




TX(2:NX,:,:)=(CF./Mu).*KX(2:NX,:,:).*DeltaY(2:end,:,:).*DeltaZ(2:end,:,:)...
    ./((DeltaX(2:end,:,:)+DeltaX(1:end-1,:,:))./2);

TY(:,2:NY,:)=(CF./Mu).*KY(:,2:NY,:).*DeltaX(:,2:end,:).*DeltaZ(:,2:end,:)...
    ./((DeltaY(:,2:end,:)+DeltaY(:,1:end-1,:))./2);


TZ(:,:,2:NZ)=(CF./Mu).*KZ(:,:,2:NZ).*DeltaX(:,:,2:end).*DeltaY(:,:,2:end)...
    ./((DeltaZ(:,:,2:end)+DeltaZ(:,:,1:end-1))./2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%Boundaries non-zeros%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  TBoun=5;
 TXB=TX;TYB=TY;
%  TZB=TZ;
%   TXB(1,:,:)=TXB(2,:,:)*2;TXB(end,:,:)=TXB(end-1,:,:)*2;TYB(:,1,:)=TYB(:,2,:)*2;
%   TYB(:,end,:)=TYB(:,end-1,:)*2;%;TZB(:,:,1)=5;TZB(:,:,end)=5;
%   
   TXB(1,:,:)=(CF./Mu).*Kxx(1,:,:).*DeltaY(1,:,:).*DeltaZ(1,:,:)...
    ./((DeltaX(1,:,:)./2));
   TXB(end,:,:)=(CF./Mu).*Kxx(end,:,:).*DeltaY(end,:,:).*DeltaZ(end,:,:)...
    ./((DeltaX(end,:,:)./2));

   TYB(:,1,:)=(CF./Mu).*Kyy(:,1,:).*DeltaX(:,1,:).*DeltaZ(:,1,:)...
    ./((DeltaY(:,1,:))./2);

   TYB(:,end,:)=(CF./Mu).*Kyy(:,end,:).*DeltaX(:,end,:).*DeltaZ(:,end,:)...
    ./((DeltaY(:,end,:))./2);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    TXB(1,:,:)=ones(1,NY)*000;TXB(end,:,:)=ones(1,NY)*000;TYB(:,1,:)=ones(NX,1)*000;
%    TYB(:,end,:)=ones(NX,1)*000;%;TZB(:,:,1)=5;TZB(:,:,end)=5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%   Gravity terms    %%%%%%%%%%%
% De=55/144;
% De=1;
% 
% GTX=De.*TX.*ZX;
% GTY=De.*TY.*ZY;
% GTZ=De.*TZ.*ZZ;
% 
% GTXs=GTX(2:end,:,:)+GTX(1:end-1,:,:);
% GTYs=GTY(:,2:end,:)+GTY(:,1:end-1,:);
% GTZs=GTZ(:,:,2:end)+GTZ(:,:,1:end-1);
% 
%     Gxv=reshape(GTXs(:,:,1)',NX*NY,1);
%     Gyv=reshape(GTYs(:,:,1)',NX*NY,1);
%     Gzv=reshape(GTZs(:,:,1)',NX*NY,1);
% GV=Gxv+Gyv+Gzv;
% for w=2:NZ
%     Gxv=reshape(GTXs(:,:,w)',NX*NY,1);
%     Gyv=reshape(GTYs(:,:,w)',NX*NY,1);
%     Gzv=reshape(GTZs(:,:,w)',NX*NY,1);
%     GV=[GV ;Gxv+Gyv+Gzv];
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TXe1=reshape(TX,NX*NY,1);
TXe=reshape(permute(TX(1:NX,:,:),[2 1 3]),[NX*NY,1,NZ]);
TXw=reshape(permute(TX(2:end,:,:),[2 1 3]),[NX*NY,1,NZ]);
TYn=reshape(permute(TY(:,1:NY,:),[2 1 3]),[NX*NY,1,NZ]);
TYs=reshape(permute(TY(:,2:end,:),[2 1 3]),[NX*NY,1,NZ]);


TZu=reshape(permute(TZ(:,:,1:NZ),[2 1 3]),[TN,1,1]);
TZd=reshape(permute(TZ(:,:,2:end),[2 1 3]),[TN,1,1]);






%%%%%%%%%%%%%%%%  The Big Matrix %%%%%%%%%%%%%
T=zeros(TN,TN);
% Tm=zeros(9,9,3);
h=1;
for g=1:NZ
    
    Tm=(spdiags([TXw(:,:,g) TYs(:,:,g) -(TXe(:,:,g)+TXw(:,:,g)+TYn(:,:,g)+TYs(:,:,g)+TZu(h:h+NX*NY-1)+TZd(h:h+NX*NY-1))...
    TYn(:,:,g) TXe(:,:,g)], [-NY -1 0 1 NY], NX*NY,NX*NY));

 T(h:h+NX*NY-1,h:h+NX*NY-1)=Tm;
 h=h+NX*NY;
end
TZu=[TZu ;zeros(NX*NY,1)];
TZd=[zeros(NX*NY,1);TZd];
[B, ~]  =spdiags(T);
T=spdiags([-TZu(NX*NY+1:end) -B(:,1) -B(:,2) -B(:,3) -B(:,4) -B(:,5) -TZd(1:end-(NX*NY))], [-NX*NY -NY -1 0 1 NY NX*NY], TN,TN);
% TZdm=diag(TZd,NX*NY);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  TT=full(T);

P=ones(NX*NY*NZ,1)*0;
% P(50)=5000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%Boundries%%%%%%%%%%%%%%%%
b=zeros(NX,NY,NZ);

PB=2;

PNorth=PB;
PSouth=PB;
PEast=PB;
PWest=PB;


PBx=zeros(NX+1,NY,NZ);
PBx(1,:,:)=PWest;
PBx(end,:,:)=PEast;
PBy=zeros(NX,NY+1,NZ);
PBy(:,1,:)=PSouth;
PBy(:,end,:)=PNorth;
%%%%%%%%%%%%%%%%%%Corner Boudaries%%%



b(1,1,:)=TXB(1,1,:).*PBx(1,1,:)+TYB(1,1,:).*PBy(1,1,:);
b(1,end,:)=TXB(1,end,:).*PBx(1,end,:)+TYB(1,end,:).*PBy(1,end,:);
b(end,end,:)=TXB(end,end,:).*PBx(end,end,:)+TYB(end,end,:).*PBy(end,end,:);
b(end,1,:)=TXB(end,1,:).*PBx(end,1,:)+TYB(end,1,:).*PBy(end,1,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%Inside Boudaries%%%%

b(1,2:end-1,:)=TXB(1,2:end-1,:).*PBx(1,2:end-1,:);
b(end,2:end-1,:)=TXB(end,2:end-1,:).*PBx(end,2:end-1,:);

b(2:end-1,1,:)=TYB(2:end-1,1,:).*PBy(2:end-1,1,:);
b(2:end-1,end,:)=TYB(2:end-1,end,:).*PBy(2:end-1,end,:);

bt=(reshape(permute(b(:,:,:),[2 1 3]),[TN,1,1]));
%%%%%%%%%%%%
%%%%%%%%%%%%
BOUT=zeros(NX,NY,NZ);

BOUT(1,1,:)=TXB(1,1,:)+TYB(1,1,:);
BOUT(1,end,:)=TXB(1,end,:)+TYB(1,end,:);
BOUT(end,end,:)=TXB(end,end,:)+TYB(end,end,:);
BOUT(end,1,:)=TXB(end,1,:)+TYB(end,1,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%Boudaries TRANS%%%%

BOUT(1,2:end-1,:)=TXB(1,2:end-1,:);
BOUT(end,2:end-1,:)=TXB(end,2:end-1,:);

BOUT(2:end-1,1,:)=TYB(2:end-1,1,:);
BOUT(2:end-1,end,:)=TYB(2:end-1,end,:);
  bm=(reshape(permute(BOUT(:,:,:),[2 1 3]),[TN,1,1]));
  bmt=spdiags(bm,0,TN,TN);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  
  
  
  
 qs1=Q;

% wello=sub2ind(size(Kxx),ceil(NX/2),ceil(NY/2));
 wello=Well_L;
bt(wello,1)=bt(wello,1)-qs1;
% bt(1,1)=0;



% A=(T+Ct+bmt);
 A=(T+Ct+bmt);
% figure 
% set(gca, 'Zdir', 'reverse')
Pwell(:,1)=P(ob);
%%
for i=1:NT
             
    Bf=Ct*P+bt;
    P=A\Bf;
    
    
%     dt=dt+dt*0.3;
%      Ct=spdiags((Cc*ones(1,TN).*(VPV'/dt))',0,TN,TN);
%      A=(T+Ct+bmt);
    
%     Pall(:,i)=P;
%     j=1;
%     for r=1:NZ
%         P3d(:,:,r)=(reshape(P(j:NX*NY+j-1),NX,NY));
%         j=j+NX*NY;
%     end 
% %    
%     Pf=reshape(Pall(:,i),NX,NY);
%      imagesc(P3d(:,:,2))
%      hold on
%     
%    [fo,vo] = isosurface(X,Y,Z,P3d,5); 
%    [fe,ve,ce] = isocaps(X,Y,Z,flip(P3d),5); 
%    
% %    figure
% % p1 = patch('Faces', fo, 'Vertices', vo);       % draw the outside of the volume
% % p1.FaceColor = 'red';
% % p1.EdgeColor = 'none';
% % hold on
% 
% p2 = patch('Faces', fe, 'Vertices', ve, ...    % draw the end caps of the volume
%    'FaceVertexCData', ce);
% p2.FaceColor = 'interp';
% p2.EdgeColor = 'none';
% view(-40,24)
% caxis([1000 3000])
% box on
% hold on

    Pwell(:,i)=P(ob);
end
Pwell2=reshape(Pwell',[length(ob)*(NT),1]);


%%
% 
% figure(2)
% 
%  hslice1=slice(X,Y,Z,(Kxx),[],[],0);
%  hold on
%   hslice2=slice(X,Y,Z,(Kxx),700,[],[]);
%   set(hslice1,'EdgeColor','none','FaceColor','interp');
%   set(hslice2,'EdgeColor','none','FaceColor','interp');
%   fv = isosurface(X,Y,Z,P3d,55); hIsoSurf1=patch(fv);
%   set(hIsoSurf1,'CData',55,'FaceColor','flat','EdgeColor','none');
% %   hstreamslice=streamslice(X,Y,Z,-u,-v,-w,[],[],0);
% %   hstreamslice2=streamslice(X,Y,Z,-u,-v,-w,700,[],[]);
% %   set(hstreamslice,'Color','w');
% %     set(hstreamslice2,'Color','w');
%   camlight;
%   camlight(-90,0);
%   lighting gouraud
% %%
% 
% figure(3)
% [u,v,w]=gradient(P3d);
%  hslice1=slice(X,Y,Z,(P3d),[],[],0);
%  hold on
%   hslice2=slice(X,Y,Z,(P3d),700,[],[]);
%   set(hslice1,'EdgeColor','none','FaceColor','interp');
%   set(hslice2,'EdgeColor','none','FaceColor','interp');
%   fv = isosurface(X,Y,Z,P3d,5500); hIsoSurf1=patch(fv);
%   set(hIsoSurf1,'CData',5500,'FaceColor','flat','EdgeColor','none');
%   hstreamslice=streamslice(X,Y,Z,-u,-v,-w,[],[],0);
%   hstreamslice2=streamslice(X,Y,Z,-u,-v,-w,700,[],[]);
%   set(hstreamslice,'Color','w');
%     set(hstreamslice2,'Color','w');
%   camlight;
%   camlight(-90,0);
%   lighting gouraud
%   
% %%
% figure(4)
%     [fo,vo] = isosurface(X,Y,Z,Kxx,5); 
%    [fe,ve,ce] = isocaps(X,Y,Z,Kxx,5); 
%    
% %    figure
% % p1 = patch('Faces', fo, 'Vertices', vo);       % draw the outside of the volume
% % p1.FaceColor = 'red';
% % p1.EdgeColor = 'none';
% % hold on
% 
% p2 = patch('Faces', fe, 'Vertices', ve, ...    % draw the end caps of the volume
%    'FaceVertexCData', ce);
% p2.FaceColor = 'interp';
% p2.EdgeColor = 'none';
% view(-40,24)
%  caxis([0 400])
% box on
% 
% %  plot(1:dt:NT,Pwell,'LineWidth',1.5)
% %  ylabel('Pressure (Psi)', 'FontWeight','bold','Color','k');
% % xlabel('Time (day)','FontWeight','bold','Color','k')
% % title('Pressure at Production Well');
