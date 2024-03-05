%%fracture parameter inversion 13-Feb-2021 (most updated version)
%% Water injection into a 3D fractured porous media
% Two-phase example with vertical producer/injector pair simulating
% water injection in a 3-dimensional fractured porous media using the EDFM
% method
clear Q wc src
close all

%% Load necessary modules
mrstModule add hfm;             % hybrid fracture module
mrstModule add mrst-gui;        % plotting routines
mrstModule add ad-blackoil;     % AD blackoil solver
mrstModule add ad-core;         % AD core module
mrstModule add ad-props;        % AD properties

%% Set up a matrix grid
% We first set up a Cartesian matrix grid with dimensions 500m x 200m x
% 100m with grid blocks that are 10m x 10m x 10m. Matrix permeability is
% 100mD and porosity is 30%.
NX=21;NY=21; NZ=3;
 Half_i=11;
 Half_z=2;
TN=NX*NZ*NY;
celldim = [NX, NY, NZ];
L=500;W=500;H=100;
% dx=L/NX;
% dz=H/NZ;
% dy=W/NY;
physdim = [L, W, H];
G = cartGrid(celldim, physdim); 
G = computeGeometry(G);
Kxx=ones(NX,NY,NZ)*30;
K=reshape(permute(Kxx,[2 1 3]), [TN,1,1]);
 phi=ones(TN,1)*0.3;
pR  = 200*barsa;
% CF=1e-5/barsa;




%% Set up fracture planes
% Fracture planes are set up by defining their vertices. Additionally,
% the aperture, porosity and permeability of the fractures are provided.
% Fracture planes 1 and 3 will be vertical while fracture 2 is slanted.
 a=30; b=470;
%  d=0;f=40;
 e=0;g=100;
 frac_numb=1;
%  fz=[0 100];
% fz=zeros(frac_numb,2);
fy=zeros(frac_numb,2);
fx=zeros(frac_numb,2);
% fz(1,:)=[randdim(d,f) randdim(e,g)];
fz(1,:)=[e g];
fy(1,:)=[randdim(a,b) randdim(a,b)];
fx(1,:)=[randdim(a,b) randdim(a,b)];
for i=1:frac_numb-1
fz(1+i,:)=[e g]; 
fy(1+i,:)=[randdim(a,b) randdim(a,b)];
fx(1+i,:)=[randdim(a,b) randdim(a,b)];
end

fx(1,:)=[0 500];
fx(2,:)=[100 400];
fy(2,:)=[100 400];
fy(1,:)=[400 100];
% 
Frac_l1=zeros(4,3,frac_numb);
Frac_a1= zeros(frac_numb,1);
Frac_poro1=zeros(frac_numb,1);
Frac_k1=zeros(frac_numb,1);
% 

fr_K=[10,200,300,1000000,1000,50000,100,5000,50,400000];
apr=[0.0001,0.0001,0.0001,0.2,0.0005,0.006,0.0001,0.001,0.0001,0.02];
for i_f=1:frac_numb
% Fracture plane 1
Frac_l1(:,:,i_f) = [fx(i_f,1) fy(i_f,1) fz(i_f,1);
                         fx(i_f,2) fy(i_f,2) fz(i_f,1);
                        fx(i_f,2) fy(i_f,2) fz(i_f,2);
                        fx(i_f,1) fy(i_f,1) fz(i_f,2)];
%  Frac_a1(i_f)= 1/100;       %original
%  Frac_poro1(i_f)=0.5;        %original
%  Frac_k1(i_f)=10000*darcy; % original
%%%%%%%%%%%%%%%%%%  k=ap^2/12
Frac_a1(i_f)= 1/50;       
Frac_poro1(i_f)=0.5;        
Frac_k1(i_f)=10000*darcy; 
% Frac_k1(i_f)=0.1*milli*darcy;

end
%    load('frac_locations2')
%%%%%%% for the one fractures
 xloca=250
 
 Frac_l1 = [xloca 0 0;
                         xloca 500 0;
                         xloca 500 100;
                            xloca 0 100;]

Frac_a1= 0.0011;       %original
 Frac_poro1=0.5;        %original
 Frac_k1=100000*darcy; % original
%  Frac_k1=0.1*milli*darcy;
% for i_f=1:frac_numb
% % Fracture plane 1
% Frac_l1(:,3,i_f) = [ fz(1);
%                           fz(1);
%                          fz(2);
%                        fz(2)];
% end

% Fracture plane 1
% Frac_l(:,:,2) = [250 500 250;
%                         850 500 250;
%                         850 500 850;
%                         250 500 850];
% Frac_a(2)= 1;
% Frac_poro(2)=1;https://acsess.onlinelibrary.wiley.com/doi/full/10.1002/vzj2.20157
% Frac_k(2)=3.3e-7;
% 
% Frac_l(:,:,3) = [250 700 250;
%                         850 700 250;
%                         850 700 850;
%                         250 700 850];
% Frac_a(3)= 1;
% Frac_poro(3)=1;
% Frac_k(3)=3.3e-7;

% Frac_l(:,:,4) = [250 400 250;
%                         850 400 250;
%                         850 400 850;
%                         250 400 850];
% Frac_a(4)= 1;
% Frac_poro(4)=1;
% Frac_k(4)=3.3e-7;

% Frac_l(:,:,5) = [250 800 250;
%                         850 800 250;
%                         850 800 850;
%                         250 800 850];
% Frac_a(5)= 1;
% Frac_poro(5)=1;
% Frac_k(5)=3.3e-7;




%% Define fluid properties
% Define a three-phase fluid model without capillarity.

% fluid = initSimpleADIFluid('mu' , [   1,  5, 0.2] .* centi*poise     , ...
%                            'rho', [1000, 700, 250] .* kilogram/meter^3, ...
%                            'n'  , [   2,   2, 2]);
                       
% fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1000*kilogram/meter^3);
CR=1e-4/barsa;
CF=1e-4/barsa;
fluid = initSimpleADIFluid('phases','W',           ... % Fluid phase: water
                           'mu',  1*centi*poise,   ... % Viscosity
                           'rho', 1000,...   % Surface density [kg/m^3]
                            'pRef',    0*barsa, ... 
                           'c',   CF,      ... % Fluid compressibility
                           'cR',  CR       ... % Rock compressibility
                           );


%%



%% Add BC
% Set boundary condition P = 10 bars on the left face and P = 1 bar on the
% right face of the domain.

%% Add wells
% An injector/producer pair is added. 1PV of water is injected over the
% course of 5 years.

 total_time = 5*day;

wc(1) = sub2ind(G.cartDims, 1,1,Half_z);
Q(1)=100/day;
Q(2)=300/day;
Q(3)=-400/day;

wc(2) = sub2ind(G.cartDims, NX,NY,Half_z);

 wc(3) = sub2ind(G.cartDims, 11,11,Half_z);

% o=[wc(1) wc(2) wc(3)];

o=[sub2ind(G.cartDims, 1,11,Half_z);...
    sub2ind(G.cartDims, 9,11,Half_z);...
    sub2ind(G.cartDims, 10,11,Half_z);...
    sub2ind(G.cartDims, 12,11,Half_z);...
    sub2ind(G.cartDims, 13,11,Half_z);...
    sub2ind(G.cartDims, 21,11,Half_z);]



%% Set up initial state and schedule
% We set up a initial state with the reference pressure and a mixture of
% water and oil initially. We also set up a simple-time step strategy that
% ramps up gradually towards 30 day time-steps.



dt = rampupTimesteps(total_time, 0.5*day, 8);
p_init=pR;
NT=length(dt);

% schedule = simpleSchedule(dt, 'W', W);

%% Simulate problem
tic

bc = [];
% bc  = pside(bc, G, 'LEFT', 200*barsa);
% bc  = pside(bc, G, 'RIGHT', 200*barsa);
% bc  = pside(bc, G, 'North', 200*barsa);
% bc  = pside(bc, G, 'South', 200*barsa);
[states_t, G_t]=FR_MRST_FRAC(NX,NY,NZ,L,W,H,fluid,p_init,dt,NT,wc,...
K,K,K,phi,CF,Q,Frac_l1,Frac_a1,Frac_poro1,Frac_k1,bc);
toc
OB=o;
figure(4)
TRUE_DATA=Vector_data2(OB,states_t,NT);
plot(TRUE_DATA);
figure(5)
%    [i j k] = ind2sub(G_t.cartDims, 1:G_t.cells.num);
%    I = false(NX,1); I([1 end])=true;
%    J = false(NY,1); J(6)=true;
%    K = false(NZ,1); K([1 end]) = true;
plotCellData(G_t, states_t{end}.pressure./barsa,...
    'EdgeColor','none')
    figure(6)
plotEDFMgrid(G_t);
% hold on;
% plotWell(G,W);
%%
clear G fracplanes
G = cartGrid(celldim, physdim); 
G = computeGeometry(G);
G.rock=makeRock(G,5e-4,1);

for i_f=1:frac_numb;
% Fracture plane 1

Frac_a(i_f)= 1/100;
Frac_poro(i_f)=1;
Frac_k(i_f)=5e-2;
end


TF=length(Frac_a);
for f=1:TF

fracplanes(f).points = Frac_l1(:,:,f);
fracplanes(f).aperture = Frac_a(f);
fracplanes(f).poro=Frac_poro(f);
fracplanes(f).perm=Frac_k(f);
end 
 re=ones(NX,NY,NZ)*5e-2;
% re(:,:,1:2)=2e-2;
res=reshape(permute(re,[2 1 3]), [TN,1,1]);
tol=1e-5;
[G,fracplanes]=EDFMgrid(G,fracplanes,...
    'Tolerance',tol,'fracturelist',1:TF);

tol=1e-5;
G=fracturematrixNNC3D(G,tol);

tol=1e-5;
[G,fracplanes]=fracturefractureNNCs3D(G,fracplanes,tol,'Verbose',true);

% Plot fracture gridblocks that have Fracture-Fracture NNCs
fraccells = G.nnc.cells(strcmp(G.nnc.type,'fracfrac'),:);

bc = [];
% bc  = pside(bc, G, 'LEFT', 200*barsa);
% bc  = pside(bc, G, 'RIGHT', 200*barsa);
% bc  = pside(bc, G, 'North', 200*barsa);
% bc  = pside(bc, G, 'South', 200*barsa);

src = addSource([], wc(1), Q(1));
for i=2:length(wc)
    
src = addSource(src, wc(i), Q(i));
end



fluid = initSingleFluid('mu', 1, 'rho', 1);
T = computeTrans(G, G.rock);
cf = G.cells.faces(:,1);
nf = G.faces.num;
T  = 1 ./ accumarray(cf, 1./T, [nf, 1]);
T = [T;G.nnc.T];
state  = initResSol (G, 1*barsa);
[Ab, ~] = getSystemIncompTPFA(state, G, T, fluid,'use_trans', true,  'bc', bc,'src',src);
%%%%
L_factor=8;
clear G fracplanes
G = cartGrid(celldim, physdim); 
G = computeGeometry(G);
G.rock=makeRock(G,5e-4*L_factor/barsa,1);




for i_f=1:frac_numb;
% Fracture plane 1
Frac_a(i_f)= 1/100;
Frac_poro(i_f)=1;
Frac_k(i_f)=5e-2*L_factor/barsa;
end


TF=length(Frac_a);
for f=1:TF

fracplanes(f).points = Frac_l1(:,:,f);
fracplanes(f).aperture = Frac_a(f);
fracplanes(f).poro=Frac_poro(f);
fracplanes(f).perm=Frac_k(f);
end 


tol=1e-5;
[G,fracplanes]=EDFMgrid(G,fracplanes,...
    'Tolerance',tol,'fracturelist',1:TF);

tol=1e-5;
G=fracturematrixNNC3D(G,tol);

tol=1e-5;
[G,fracplanes]=fracturefractureNNCs3D(G,fracplanes,tol,'Verbose',true);

% Plot fracture gridblocks that have Fracture-Fracture NNCs
fraccells = G.nnc.cells(strcmp(G.nnc.type,'fracfrac'),:);


fluid = initSingleFluid('mu', 1, 'rho', 1);
T = computeTrans(G, G.rock);
cf = G.cells.faces(:,1);
nf = G.faces.num;
T  = 1 ./ accumarray(cf, 1./T, [nf, 1]);
T = [T;G.nnc.T];
state  = initResSol (G, 1*barsa);
[Ac, rhs] = getSystemIncompTPFA(state, G, T, fluid,'use_trans', true,  'bc', bc,'src',src);
% Ac=Ac(1:end-1,1:end-1);
% state_sp={};
% for N=1:length(dt)
% pf2=Ab\-Ac*states_t{end}.pressure;
%  state_sp{N}=0-pf2;
% end
OB_sp=linspace(1,NX*NY,NX*NY);
state_sp=cell(NT,1);
state_sp{1}=zeros(TN,1);
for tt=1:NT
    
pf2=Ab\-Ac*states_t{tt}.pressure;
state_sp{tt}=pf2(1)-pf2;
end
TRUE_DATA_sp=Vector_data_sp2(OB_sp,state_sp,NT);
%  figure(6)
%  plotCellData(G, state_sp,'EdgeAlpha',.1)
 %%
 

% TN=NX*NZ*NY;
% celldim = [NX, NY, 1];
% physdim = [L, W, H];
% Gs = cartGrid(celldim, physdim); 
% Gs = computeGeometry(Gs);
% RM=Reg_matrix(NX,NY,NZ);
% A=Ab;
% source=A(1:NX*NY,1:NX*NY)*state_sp(1:NX*NY);
% figure(7)
% plotCellData(Gs,source)

%% inversion 
%%%% choose cases 1 and 2 for only fluid flow, 3 for the joint%%%%
n = 3;

switch n
    case 1
        bcFaces1 = find(all([G.faces.centroids(:,2) == W],2));
        bcFaces2 = find(all([G.faces.centroids(:,2) == 0],2));
        bcFaces3 = find(all([G.faces.centroids(:,1) == L],2));
        bcFaces4 = find(all([G.faces.centroids(:,1) == 0],2));
        bc = [];

        % bc = addBC(bc, bcFaces1, 'pressure', 200*barsa);
        % bc = addBC(bc, bcFaces2, 'pressure', 200*barsa);
        % bc = addBC(bc, bcFaces3, 'pressure', 200*barsa);
        % bc = addBC(bc, bcFaces4, 'pressure', 200*barsa);
        fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1000*kilogram/meter^3);
        MaxI=40;
        [output]=FLOW_inversion(TRUE_DATA,TN,MaxI,o,...
            NX,NY,NZ,L,W,H,fluid,p_init,dt,NT,wc,...
            phi,CR,CF,Q,bc);
        celldim = [NX, NY, NZ];
        physdim = [L, W, H];
        Gi = cartGrid(celldim, physdim); 
        Gi = computeGeometry(Gi);
        figure(8)
        plotCellData(Gi,output{1})
    case 2
        
        fluid = initSimpleADIFluid('phases','W',           ... % Fluid phase: water
                           'mu',  1*centi*poise,   ... % Viscosity
                           'rho', 1000,...   % Surface density [kg/m^3]
                            'pRef',    0*barsa, ... 
                           'c',   1e-4/barsa,      ... % Fluid compressibility
                           'cR',  1e-4/barsa       ... % Rock compressibility
                           );
           MaxI=10;
        [rmse, F_image]=FLOW_inversion2(TRUE_DATA,TN,MaxI,o,...
            NX,NY,NZ,L,W,H,fluid,p_init,dt,NT,wc,...
            phi,CF,Q,Frac_l1,Frac_a1,Frac_poro1,Frac_k1);
        celldim = [NX, NY, NZ];
        physdim = [L, W, H];
        Gi = cartGrid(celldim, physdim); 
        Gi = computeGeometry(Gi);
        figure(8)
        plotCellData(Gi,F_image)
        
            case 3
        bcFaces1 = find(all((G.faces.centroids(:,2) == W),2));
        bcFaces2 = find(all((G.faces.centroids(:,2) == 0),2));
        bcFaces3 = find(all((G.faces.centroids(:,1) == L),2));
        bcFaces4 = find(all((G.faces.centroids(:,1) == 0),2));
        bc = [];

%         bc = addBC(bc, bcFaces1, 'pressure', 200*barsa);
%         bc = addBC(bc, bcFaces2, 'pressure', 200*barsa);
%         bc = addBC(bc, bcFaces3, 'pressure', 200*barsa);
%         bc = addBC(bc, bcFaces4, 'pressure', 200*barsa);
        fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1000*kilogram/meter^3);
        MaxI=80;
        [output]=FLOW_SP_inversion(TRUE_DATA,TRUE_DATA_sp,TN,MaxI,o,...
            OB_sp,NX,NY,NZ,L,W,H,fluid,p_init,dt,NT,wc,...
            phi,CR,CF,Q,bc,res,L_factor);
        celldim = [NX, NY, NZ];
        physdim = [L, W, H];
        Gi = cartGrid(celldim, physdim); 
        Gi = computeGeometry(Gi);
        figure(8)
        plotCellData(Gi,output{1})
                       
                       
                       
    
end
%%  Ploting
% 100mD and porosity is 30%.

% NX=20;NY=20; NZ=2;
TN=NX*NZ*NY;
celldim = [NX, NY, NZ];
L=500;W=500;H=100;
physdim = [L, W, H];
Gp = cartGrid(celldim, physdim); 
Gp = computeGeometry(Gp);
TF=length(Frac_a);
for f=1:TF

fracplanes(f).points = Frac_l1(:,:,f);
fracplanes(f).aperture = Frac_a1(f);
fracplanes(f).poro=Frac_poro1(f);
fracplanes(f).perm=Frac_k1(f);
end 

Wells = addWell([], G, G.rock,  wc(1),     ...
        'Type', 'rate', 'Val', Q(1), ...
        'Radius', 0.1, 'Name', 'I_1');
for i=2:length(wc)-1
    
    Wells = addWell(Wells, G, G.rock,  wc(i),     ...
        'Type', 'rate', 'Val', Q(i), ...
        'Radius', 0.1, 'Name', 'I_2');
end
Wells2 = addWell([], G, G.rock,  sub2ind(G.cartDims, Half_i,Half_i,Half_z),     ...
        'Type', 'rate', 'Val', Q(1), ...
        'Radius', 0.1, 'Name', 'P');


figure(9);
 plotfracongrid(Gp,fracplanes)
% plotfracongrid_notext(Gp,fracplanes); % visualize to check before pre-process
view(30,45)
hold on
plotWell(G,Wells,'color','k','radius',0.5);
plotWell(G,Wells2,'color','k','radius',0.5);
xlabel('x(m)', 'fontweight', 'bold')
ylabel('y(m)', 'fontweight', 'bold')
zlabel('z(m)', 'fontweight', 'bold')
ZMin=0;
dx=L/NX;
dz=H/NZ;
dy=W/NY;
 x=dx:dx:L;z=0:dz:H-dz;y=dy:dy:W;
for i=1:length(x)
    for j=1:length(y)
        plot3(x(i)-12.5,y(j)-12.5,ZMin,'kv','LineWidth',1,'MarkerSize',5)
    end
end
figure(10)
plotfracongrid_notext(Gp,fracplanes)
hold on
h=plotCellData(Gp,output{1});
view(30,45)
set(h, 'EdgeColor','none', 'FaceColor','flat')
caxis([10 200])
xlabel('x(m)', 'fontweight', 'bold')
ylabel('y(m)', 'fontweight', 'bold')
zlabel('z(m)', 'fontweight', 'bold')
colormap cool
a = colorbar;
a.Label.String = 'k(mD)';
figure(11)
   [i_1 ,j_1, k_1] = ind2sub(Gp.cartDims, 1:Gp.cells.num);
   I_1 = false(NX,1); I_1([ ])=true;
   J_1 = false(NY,1); J_1([3 9])=true;
   K_1 = false(NZ,1); K_1([8]) = true;
   plotCellData(Gp, output{1}, I_1(i_1) | J_1(j_1) | K_1(k_1),'EdgeAlpha',.1);
view(30,45)
%% Ploting2
to=length(dt)
figure(6)
dx=L/NX;
dz=H/NZ;
dy=W/NY;
x=dx:dx:L;z=0:dz:H-dz;y=dy:dy:W;
subplot(1,2,1)
sp_1=reshape(permute(TRUE_DATA_sp(5:to:end),[2 1 3]), [NX,NY,1]);
imagesc(x,y,sp_1')
%  caxis([-5 0])
axis square
set(gca,'YDir','normal')
xlabel('x(m)')
ylabel('y(m)')
c=colorbar
c.Label.String = 'Self-Potential (mV)';
subplot(1,2,2)
sp_1=reshape(permute(TRUE_DATA_sp(to:to:end),[2 1 3]), [NX,NY,1]);
imagesc(x,y,sp_1')
colormap jet
%  caxis([-5 0])
axis square
xlabel('x(m)')
ylabel('y(m)')
c=colorbar;
c.Label.String = 'Self-Potential (mV)';
set(gca,'YDir','normal')

%%


% RM=Reg_matrix(NX,NY,NZ);
epsilon=output{6}
% epsilon=0.002

TN=21*21*3;
z=output{3}(:,1:TN);

si=size(z)

%  M = z'*inv((z*z'+epsilon*(eye(si(1),si(1)))))*z;

 M=output{4};
 

 d1=M;
 
 d2=eye(TN,TN);
  sp=norm(d1-d2, "fro")/(63)
lg=TN;

i=NX*NY;
W=zeros(TN,TN);
d=size(diag(W,i),1);
e=0.5*ones(d,1);
W(i*size(W,1)+1:size(W,1)+1:end)=e;
W(i+1:size(W,1)+1:1+size(W,1)*min(size(W,1)-i,size(W,2)))=e;
i2=2*i
d=size(diag(W,i2),1);
e=0.25*ones(d,1);
W(i2*size(W,1)+1:size(W,1)+1:end)=e;
W(i2+1:size(W,1)+1:1+size(W,1)*min(size(W,1)-i2,size(W,2)))=e;

SpM=0;
I=eye(lg);
for i=1:lg
    for j=1:lg
        SpM=SpM+W(i,j)*(M(i,j)-I(i,j))^2;
    end
end
sp=SpM;
figure(13)
% imagesc(M)
spy(M>0.005);
colormap('jet')
caxis([0.0001 0.1])
xlabel('m')
ylabel('m')
% str=num2str(SpM);
    str=num2str((sp));
title({
    ['Model Resolution Matrix' ] 
    ['Spread = ' str  ] 
%     ['IRMSE = ' str2 ]
    });

%  cov_m=M*M';
%  cov_m=cov_m./max(max(cov_m));
% str=trace(cov_m)
%  imagesc(cov_m)

%%
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontName','helvetica','DefaultAxesTickLength',[.02 .02],'DefaultAxesLineWidth',2) 

set(figure(1), 'Position', [0 50 500 600]);
PS = PLOT_STANDARDS();
clear alltime
colors=[PS.DYellow4; PS.DOrange4; PS.DBlue4; PS.LGreen4; PS.Yellow5]
alltime(1)=dt(1)/day;
for i=2:length(dt)
    alltime(i)=alltime(i-1)+dt(i)/day;
end
n=1;
subplot(2,1,2)
for i=1:3
    
semilogx(alltime,TRUE_DATA(n:n+NT-1),'.','LineWidth',3, 'MarkerSize', 20,'DisplayName','DATA','Color',colors(i,:))
hold on

n=n+NT
end
grid on
xlabel('Time (days)');
ylabel('Pressure(bar)');
%%




 function ra=randdim(a,b)
ra=a+rand(1)*(b-a);
end 