function [rmse, F_image]=FLOW_inversion2(TRUE_DATA,TN,MaxI,OB,...
    NX,NY,NZ,L,W,H,fluid,p_init,dt,NT,wc,...
    phi,CF,Q,Frac_l,Frac_a,Frac_poro,Frac_k)

NDS=(length(TRUE_DATA));

Wd=ones(length(TRUE_DATA),1);

KV_I=10;
K_Ix=KV_I*ones(TN,1);
K_Ix(1)=1;
%  K_Ix=K_x+K_x.*rand(TN,1);
% K_Ix=
% K_Ix(in)=K_x(in);
K_Iy=KV_I*ones(TN,1);
K_Iz=KV_I*ones(TN,1);

K_MIN=0.001;
K_MAX=1000;
%%%%%%%%%%%%%Inversion%%%%%%%%%%%%
stp=1e-20;
PAR=5; % Percentage increase in sensitivity
PAR=PAR/100;


iteration=1;
maxiteration=MaxI;
alphaalpha=0.1;
dfit=1;

% set(figure(4), 'Position', [500 600 400 300]);
% clf
% plotCellData(G, rock.perm(:, 1)/(milli*darcy))
% view(50, 50), axis tight
% colorbar
% title('K_x [mD]')
% caxis([-2 2])


RM=Reg_matrix(NX,NY,NZ);
DATA1=zeros(1,length(OB)*(NT-1));
DATA2=zeros(1,length(OB)*(NT-1));
fig1=figure(1)
fig2=figure(2)
fig3=figure(3)
set(fig1, 'Position', [0 50 400 300]);
set(fig2, 'Position', [500 50 400 300]);
set(fig3,'Position', [1000 50 400 300]);
while iteration<maxiteration
    tic

    jj=iteration

    [states1, G_1]=FR_MRST_FRAC(NX,NY,NZ,L,W,H,fluid,p_init,dt,NT,wc,...
 K_Ix,K_Ix,K_Ix,phi,CF,Q,Frac_l,Frac_a,Frac_poro,Frac_k,bc);
 
     
DATA1=Vector_data2(OB,states1,NT); 
    EE=TRUE_DATA-DATA1;
    E1=reshape(EE',NDS,1); % NUMBERING FROM TOP LEFT TO TOP RIGHT, THEN NEXT ROW ETC
    EE2=Wd.*E1; % For real data
    dd=E1;
  
    misfit1=(EE2'*EE2);
    misfit(jj)=sqrt(misfit1/NDS);

NorMi(jj)=misfit(end)/(mean(DATA1));
if jj>1
    
    NorMi2(jj-1)=abs((misfit(end)-misfit(end-1))/(misfit(end-1)));
end
    
    
set(0,'CurrentFigure',fig2);
    axis square
%     plot(jj,log(misfit(jj)),'k*','linewidth',2)
    plot(jj,log(NorMi(jj)),'k*','linewidth',2)
     xlabel('Iteration');
ylabel('NRMSE');
    xlim([0 jj+1])
    drawnow

    hold on;
    
%     set(fig3);
    axis square
    set(0,'CurrentFigure',fig3);
%     plot(TRUE_DATA(end-NT:end),'b*','linewidth',2)
    plot(TRUE_DATA,'b*','linewidth',2,'DisplayName','Data')
    
    hold on
%     plot(DATA1(end-NT:end),'r','linewidth',2)
    plot(DATA1,'r','linewidth',2,'DisplayName','Model')
    legend('Location','southeast')
    hold off
%     xlim([0 jj+1])
    

    
          if jj>1
    prmse=abs((misfit(jj-1)-misfit(jj))/1);
 
    if prmse < stp
        break
    end
          end

if jj>1
if misfit(jj-1)>misfit(jj)
    
alphaalpha=alphaalpha/3;
else 
    alphaalpha=alphaalpha*4;
end
end
    J=zeros(NDS,TN);
%     KK=0;
    
    K_Ixj=repmat(K_Ix,1,TN);
    K_Ixj2=K_Ixj;
for j2=1:TN
        
        K_Ixj2(j2,j2)=K_Ixj(j2,j2)*(1+PAR);  % 10% perturbatn of each layer resistivity
   
%         K_Ixj2=K_Ixj;        % PX_K2 is initialized to PX_K so that only d next model parameter is perturbed.
  
end

KK=TN;
tic
        [PA]=FR_FRAC_In_PACK(NX,NY,NZ,L,W,H,G_1,fluid,p_init,dt,NT,wc,...
   K_Ixj2,K_Ix,K_Ix,phi,CF,Q,KK,OB);

%         [PALL, PA(:,KK)]=FLOW_3D_1(NX,NY,NZ,L,W,H,dt,NT,Q,Well_L,K_Ixj2(:,KK),K_y,K_z,phi,Cc,Mu,Ob);
%         PA(:,KK)=Vector_data2(OB,states1,NT)';
toc

%  PA=PA./barsa;
    
    parfor LL=1:TN
%          J(:,LL)=((PA(:,LL)-DATA1')./(K_Ix(LL)*PAR));
         J(:,LL)=((PA(:,LL)-DATA1')./(K_Ix(LL)*PAR)).*(K_Ix(LL)./DATA1');
   
    end 
    
    



logD=log(TRUE_DATA./DATA1);

        
  
       j=1;
       

 dm=(J'*J+(alphaalpha^2)*(RM'*RM))\(J'*logD');



        K_Ix=K_Ix.*exp(dm);
%            K_Ix=K_Ix+dm;
         for hh=1:TN
    
                if(K_Ix(hh)<K_MIN)
                K_Ix(hh)=K_MIN;
                end
                if(K_Ix(hh)>K_MAX)
                K_Ix(hh)=K_MAX;
                end 
              
         
        end

             iteration=iteration+1;
    toc         
end
F_image=K_Ix;
rmse=NorMi;
end