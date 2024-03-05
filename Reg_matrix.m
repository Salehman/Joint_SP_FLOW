function RM=Reg_matrix(NX,NY,NZ)


TN=NX*NY*NZ;
TX=zeros(NX+1,NY,NZ);
TY=zeros(NX,NY+1,NZ);
TZ=zeros(NX,NY,NZ+1);




TX(2:NX,:,:)=1;

TY(:,2:NY,:)=1;


TZ(:,:,2:NZ)=1;

   

   TXe=reshape(permute(TX(1:NX,:,:),[2 1 3]),[NX*NY,1,NZ]);
TXw=reshape(permute(TX(2:end,:,:),[2 1 3]),[NX*NY,1,NZ]);
TYn=reshape(permute(TY(:,1:NY,:),[2 1 3]),[NX*NY,1,NZ]);
TYs=reshape(permute(TY(:,2:end,:),[2 1 3]),[NX*NY,1,NZ]);


TZu=reshape(permute(TZ(:,:,1:NZ),[2 1 3]),[TN,1,1]);
TZd=reshape(permute(TZ(:,:,2:end),[2 1 3]),[TN,1,1]);

RM=zeros(TN,TN);
% Tm=zeros(9,9,3);
h=1;
for g=1:NZ
    
    Tm=(spdiags([TXw(:,:,g) TYs(:,:,g) -(TXe(:,:,g)+TXw(:,:,g)+TYn(:,:,g)+TYs(:,:,g)+TZu(h:h+NX*NY-1)+TZd(h:h+NX*NY-1))...
    TYn(:,:,g) TXe(:,:,g)], [-NY -1 0 1 NY], NX*NY,NX*NY));

 RM(h:h+NX*NY-1,h:h+NX*NY-1)=Tm;
 h=h+NX*NY;
end
TZu=[TZu ;zeros(NX*NY,1)];
TZd=[zeros(NX*NY,1);TZd];
[B, ~]  =spdiags(RM);
RM=spdiags([-TZu(NX*NY+1:end) -B(:,1) -B(:,2) -B(:,3) -B(:,4) -B(:,5) -TZd(1:end-(NX*NY))], [-NX*NY -NY -1 0 1 NY NX*NY], TN,TN);
