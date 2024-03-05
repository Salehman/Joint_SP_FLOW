%% Example: Depletion of a closed or open reservoir compartment
% In this tutorial, we will show how to set up a simulator from scratch in
% the automatic differentiation, object-oriented (AD-OO) framework without
% the use of input files. As an example we consider a 2D rectangular
% reservoir compartment with homogeneous properties, drained by a single
% producer at the midpoint of the top edge. The compartment is either
% closed (i.e., sealed by no-flow boundary conditions along all edges), or
% open with constant pressure support from an underlying, infinite aquifer,
% which we model as a constant-pressure boundary condition.
close all
clear all
mrstModule add ad-props  ad-core ad-blackoil add optimization
mrstModule add incomp
tic
set(0,'DefaultAxesFontSize',12,'DefaultAxesFontName','helvetica','DefaultAxesTickLength',[.02 .02],'DefaultAxesLineWidth',2) 

%% Grid, petrophysics, and fluid objects
% To create a complete model object using the AD-OO framework, we first
% need to define three standard MRST structures representing the grid and
% the rock and fluid properties
close all
clear all
% The grid and rock model
L=1000
Lz=10
[nx,ny,nz] = deal(51,51,1);

TN_b=nx*ny;
xx= ceil(nx/2),yy= ceil(ny/2);

Gb    = computeGeometry(cartGrid([nx ny nz],[L L Lz]));
r1 = sqrt(sum(bsxfun(@minus,Gb.cells.centroids,[L/2 L/2 5]).^2,2));
% r2 = sum(bsxfun(@minus,G.cells.centroids,[0.5 0]).^2,2);
clf, plotCellData(Gb, double((r1<L/2) ) );

G = extractSubgrid(Gb,(r1<L/2 ));
wl=sub2ind([nx,ny],ceil(nx/2),ceil(ny/2))

[w v]=find(G.cells.indexMap==wl)





% rock = makeRock(G, 500*milli*darcy, 1);

TN=G.cells.num;



K=ones(TN,1)*10;
phi=ones(TN,1)*0.1;
rock = makeRock(G, [K, K, K]*milli*darcy,0.3);
% 
% [IX,IY,IZ]=ndgrid(1:G.cartDims(1),1:G.cartDims(2),1:G.cartDims(3));
% slice_ind=5;
% cellind=sub2ind(G.cartDims,IX(:,:,slice_ind),IY(:,:,slice_ind),IZ(:,:,slice_ind));
% cellind=cart2active(G,cellind(:));

op_r=setupOperatorsTPFA(G,rock);
% Fluid properties

gravity off
mu=1*centi*poise ;
fluid     = initSingleFluid('mu' ,    mu     , ...
                            'rho', 1014*kilogram/meter^3);
Kcon=(mu)/rock.perm(1,1);

resSol = initResSol(G, 0.0);
display(resSol);
%   bc = fluxside([], G, 'LEFT',  1*meter^3/day());
%  bc.value(2:end)=0;
% bc = fluxside([], G, 'LEFT', -10000);
% bc = fluxside(bc, G, 'RIGHT',  -01000);
% bc = fluxside(bc, G, 'NORTH',  -010000);
% bc = fluxside(bc, G, 'SOUTH',  -01000);
bcc=100*barsa;
bc = pside   ([], G, 'LEFT', bcc);
bc = pside   (bc, G, 'RIGHT', bcc);
bc = pside   (bc, G, 'NORTH', bcc);
bc = pside   (bc, G, 'SOUTH', bcc);
% bc = pside   (bc, G, 'top', 0);
% bc = pside   (bc, G, 'bottom', 0);
src=addSource([],w,1*meter^3/day());
% src.cell=[1:1:nx]';
% src.rate=G.cells.centroids(:,1);
display(bc);
rock.cr = 1e-4/barsa

T = computeTrans(G, rock, 'Verbose', true);

resSol = incompTPFA(resSol, G, T, fluid, ...
                   'bc', bc,'src',src ,'MatrixOutput', true);
display(resSol)
clf

W = verticalWell([], G, rock,ceil(nx/2), ceil(ny/2),1, ...
            'Type', 'bhp' , 'Val', 0*barsa(), ...
            'Radius', 0.1, 'InnerProduct', 'ip_tpf', ...
            'Comp_i', [1, 0]);
 W.name='Injecting Well'

   [i ,j, k] = ind2sub(G.cartDims, 1:G.cells.num);
   I = false(nx,1); I([1 end])=true;
   J = false(ny,1); J([11])=true;
   K2 = false(nz,1); K2([1 end]) = true;
fig1=figure(1)
fig1.Position=[0 0 500 500]
% plotCellData(G, convertTo(resSol.pressure(1:G.cells.num), barsa()),...
%     I(i) | J(j) | K2(k), ...
%              'EdgeColor', 'k');
title('Cell Pressure [bar]')
xlabel('x'), ylabel('y'), zlabel('Depth');
view(3); shading faceted; camproj perspective; axis tight;
colorbar
figure(2)
 plot(resSol.pressure)
 hold on
%% Drive mechansims and schedule
% The second thing we need to specify is the mechanisms that will drive
% flow in the reservoir, i.e., the wells and boundary conditions. These may
% change over time and MRST therefore uses the concept of a schedule that
% describes how the drive mechansims change over time. In our case, we use
% the same setup for the whole simulation. The schedule also enables us to
% specify the time steps we wish the simulator to use, provided that these
% give convergent and stable computations. (If not, the simulator may cut
% the time step).

% Well: at the midpoint of the south edge






rock = makeRock(G, 1e-3, 1);
mrstModule add ad-core
nc=G.cells.num
op1=setupOperatorsTPFA(G,rock)


rock2 = makeRock(G, 1e-6, 1);
op2=setupOperatorsTPFA(G,rock2)

pv = initVariablesADI(zeros(nc,1));
p2 = initVariablesADI(resSol.pressure);



q=zeros(G.cells.num,1);
q(w)=src.rate
Source_V=rock2.perm(1,1)*-Kcon*q(w)*10e-5;
sigm=rock.perm(1,1)

rad=sqrt(sum([(G.cells.centroids(:,:)-G.cells.centroids(w,:)).^2],2));
rs=rad(w)
rs=0
RR=max(rad)
r=linspace(-L/2,L/2,nx)
dv=G.cells.volumes(1);
def=abs(sqrt(r.^2+rs^2));
def=abs(rad);
% def=abs(rad)
def(def<0.0001)=1;

pp=-q(w)*Kcon./(2*pi*Lz)*(log(def)-log(RR))+bcc
% pp=q(w)*Kcon./(4*pi*def)-(1/(2*pi))*log(RR);

 eq    = op1.Div(op1.T.*op1.Grad(pv))+op2.Div(op2.T.*op2.Grad(p2)); % equation

% eq=op1.Div(op1.Grad(p))+(rock2.perm(1,1)*-Kcon*q);
% eq(1) = eq(1) + p(1); % make solution unique
eq(1) = eq(1) + pv(1); 
Jac=eq.jac{1};
ep=0.01;
% pv     = -(Jac+eye(G.cells.num)*ep)\eq.val; 
 pv=-Jac\eq.val;
% solve equation
% clf, plotCellData(G, p);
rad=sqrt(sum([(G.cells.centroids(:,:)-G.cells.centroids(1,:)).^2],2));
dr=(rad(2)-rad(1));





plot(pp,'o')


  bc = pside   ([], G, 'LEFT', 0);
bc = pside   (bc, G, 'RIGHT', 0);
bc = pside   (bc, G, 'NORTH', 0);
bc = pside   (bc, G, 'SOUTH', 0);

T = computeTrans(G, rock, 'Verbose', true);
volSol = initResSol(G, 0.0);
src.rate=rock2.perm(1,1)*-Kcon*q(w)*10e-5;
fluid     = initSingleFluid('mu' ,    1     , ...
                            'rho', 1);
volSol = incompTPFA(volSol, G, T, fluid, ...
                   'bc', bc,'src',src ,'MatrixOutput', true);
dv=G.cells.volumes(1);
% for ix=1:length(rad)
%     g=0;
%     g2=0;
%     
%     for is=1:length(rad)
% %         if rad(ix)>=rad(is)
% %             h=1;
% % %             green=rad(ix)-1;
% %             green=(rad(ix)-L)*(rad(is)/L);
% %         else
% %             h=0;
% % %             green=rad(is)-1;
% %             green=(-rad(ix)+(rad(ix)*rad(is)/L));
% %         end
%          green=(1/(4*pi))*1/abs(rad(ix)-rad(is));
%          green(isnan(green)|isinf(green))=0;
%          g=g+green*q(is)*dv;
%          g2=g2+green*Source_V(is)*dr;
% %          g=g+green*pp(is)*1e-3/(50*milli*darcy)*dr
% %          g=g+(((rad(is)-rad(ix))*h+rad(is)*(rad(ix)-1))*q(is)*1e-3/(50*milli*darcy));
%     end
%     
%     pressure_analytical(ix)=-Kcon*g;
%     voltage_analytical(ix)=-(1/sigm)*g2;
%     pressure_analytical(ix)=-g+rad(ix)*pp(1);
% end
%%
p_analytical=-Kcon*((rad.^3)/6)+Kcon*10000*rad/6;
% SPALL=p(1)-p;
Source_V=rock2.perm(1,1)*q(w)*Kcon*10e-5;
vv=Source_V./(2*sigm*pi*Lz)*(log(def)-log(RR));



% figure(2)
% 
% plotCellData(G, convertTo(pressure_analytical(1:G.cells.num)', barsa()), ...
%              'EdgeColor', 'k');
% title('Cell Pressure [bar]')
% xlabel('x'), ylabel('y'), zlabel('Depth');
% view(3); shading faceted; camproj perspective; axis tight;
% colorbar
% fig2=figure(2)
% fig2.Position=[0 1000 500 500]
% % plot(pressure_analytical)
% % plot(p_analytical)
% plot(pressure_analytical)
% hold on
% % plot(pressure_analytical,'x')
% plot(resSol.pressure, 'o')


% fig3=figure(3)
% fig3.Position=[1000 1000 1000 500]

% subplot(1,2,1)
%    [i ,j, k] = ind2sub(G.cartDims, 1:G.cells.num);
%    I = false(nx,1); I([1 end])=true;
%    J = false(ny,1); J([11])=true;
%    K = false(nz,1); K([1 end]) = true;
%    plotCellData(G, volSol.pressure,'EdgeAlpha',.1);
%    title(' Potential [mV]')
%    view(-40,70), camproj perspective, axis tight, box on
%     set(gca,'XTick',[],'YTick',[],'ZTick',[]);
%     colorbar
% %    subplot(1,2,2)
% %    [i ,j, k] = ind2sub(G.cartDims, 1:G.cells.num);
% %    I = false(nx,1); I([1 end])=true;
% %    J = false(ny,1); J([11])=true;
% %    K = false(nz,1); K([1 end]) = true;
% % %    plotCellData(G, p_analytical, I(i) | J(j) | K(k),'EdgeAlpha',.1);
% % [nx,ny,nz] = deal(100,1,1);
% % G2    = computeGeometry(cartGrid([nx ny nz],[L 1 1]));
% %       plotCellData(G, voltage_analytical','EdgeAlpha',.1);
% %    title(' Potential [mV]')
% %    view(-40,70), camproj perspective, axis tight, box on
% %    set(gca,'XTick',[],'YTick',[],'ZTick',[]);
% %     colorbar
toc
fig4=figure(4)
fig4.Position=[1500 0 500 500]


plot(vv,'o')
hold on
% plot(pressure_analytical,'x')
 plot(volSol.pressure)
 plot(pv,'.')
%%
set(0,'DefaultAxesFontSize',12,'DefaultAxesFontName','helvetica','DefaultAxesTickLength',[.02 .02],'DefaultAxesLineWidth',2) 
PS = PLOT_STANDARDS();
fig5=figure(5)
fig5.Position=[0 500 500 500]
[i j k] = ind2sub(G.cartDims, 1:G.cells.num);
clf;
% plotGrid(G, 'FaceAlpha', 0.0, 'EdgeAlpha', 0.1)
% hold on
c=[0.4941    0.1843    0.5569]
face=find(G.faces.normals(:,3)>0);
[IX,IY,IZ]=ndgrid(1:G.cartDims(1),1:G.cartDims(2),1:G.cartDims(3));
slice_ind=1;
cellind=sub2ind(G.cartDims,IX(:,:,slice_ind),IY(:,:,slice_ind),IZ(:,:,slice_ind));
cellind=cart2active(G,cellind(:));

% plotFaces(G, face(cellind),'EdgeAlpha', 0, 'FaceColor', PS.Blue1,'FaceAlpha', 0.5);
plotGrid(G,'EdgeAlpha', 0, 'FaceColor', PS.Blue1,'FaceAlpha', 0.5)
% Plot the wells
plotWell(G, W,'color','k','radius',    1);
xlabel('x'), ylabel('y'), zlabel('Depth');
% shading faceted; camproj perspective; axis tight;
view(-14,18)
% plotSlice(G,resSol.pressure,5,2)
zlim([0 10])
pbaspect([1 1 0.1])
title('Reservoir Model')
fig5=figure(6)
fig5.Position=[0 500 800 800]
[X,Y] = meshgrid(Gb.nodes.coords(1:nx,1),Gb.nodes.coords(1:nx+1:nx*ny,2));

P_temp=zeros(TN_b,1)+bcc;
P_temp(G.cells.indexMap)=resSol.pressure;
pp_temp=zeros(TN_b,1)+bcc;
pp_temp(G.cells.indexMap)=pp;
Pf_s=convertTo(reshape(P_temp,[nx,ny]),barsa);
Pf_a=convertTo(reshape(pp_temp,[nx,ny]),barsa);
subplot(2,2,1)
contour(X,Y,Pf_a,'LineWidth',3,'ShowText','on','color',PS.Yellow4)

hold on
contour(X,Y,Pf_s,'k--','LineWidth',2,'color',PS.Blue3)

legend('Analytical','Simulation')

% plotSlice(G,resSol.pressure,5,3)
title( 'Fluid-flow')
grid on
axis tight
axis square
xlabel('x'), ylabel('y')
colormap(winter)
subplot(2,2,2)

V_temp=zeros(TN_b,1);
V_temp(G.cells.indexMap)=volSol.pressure;
vv_temp=zeros(TN_b,1);
vv_temp(G.cells.indexMap)=vv;

V_s=reshape(V_temp,[nx,ny]);

V_a=reshape(vv_temp,[nx,ny]);
V_a(V_a==0)=V_a(5,12);
V_s(V_s==0)=V_s(5,12);
% contour(X,Y,V_s,'k--','LineWidth',2)
contour(X,Y,V_a*1000,'LineWidth',3,'ShowText','on','color',PS.Blue5)
hold on

contour(X,Y,V_s*1000,'k--','LineWidth',2,'color',PS.Orange1)
legend('Analytical','Simulation')
grid on
% clim([-20 -.50])
% caxis([-20 1])
axis tight
axis square
title( 'Self-potential')
colormap(winter)
xlabel('x'), ylabel('y')


[corsin ~]=find(G.cells.centroids(:,2)==500);

subplot(2,2,3)

plot(X(1,:),convertTo(pp(corsin),barsa), 'color',PS.Yellow4,'lineWidth',2)
hold on
plot(X(1,:),convertTo(resSol.pressure(corsin),barsa),'o','color',PS.Blue3,'MarkerSize',6,...
    'lineWidth',2)
xlabel('x')
ylabel('Pressure(bar)')
legend('Analytical','Simulation')
grid on
axis tight
axis square
subplot(2,2,4)

plot(X(1,:),vv(corsin)*1000,'color',PS.Blue5,'lineWidth',2)
hold on
plot(X(1,:),volSol.pressure(corsin)*1000,'o','color',PS.Orange1,'MarkerSize',6,...
    'lineWidth',2)
xlabel('x')
ylabel('Potential(mV)')
legend('Analytical','Simulation','location','southeast')

grid on
axis tight
axis square
