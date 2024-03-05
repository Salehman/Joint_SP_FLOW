function [PA,state_sp]=FR_F_SP_MATRIX2(NX,NY,NZ,L,W,H,fluid,p_init,dt,NT,wc,...
K_x,K_y,K_z,phi,CR,CF,Q,BC,res,L_factor)


[nx,ny,nz] = deal(NX,NY,NZ);
G = computeGeometry(cartGrid([nx ny nz],[L W H]));

TN=nx*ny*nz;




rock = makeRock(G, [K_x, K_x, K_x]*milli*darcy, phi);
rock.cr=CR;

rock1 = makeRock(G, res, 1);
T1 = computeTrans(G, rock1);
cf = G.cells.faces(:,1);
nf = G.faces.num;
T1  = 1 ./ accumarray(cf, 1./T1, [nf, 1]);
rock2 = makeRock(G, res.*L_factor./barsa, 1);

T2 = computeTrans(G, rock2);
cf = G.cells.faces(:,1);
nf = G.faces.num;
T2  = 1 ./ accumarray(cf, 1./T2, [nf, 1]);

fluid.cf=CF;
T = computeTrans(G, rock);
pR  = 200*barsa;
Wells = addWell([],G,rock,wc(1),'type','rate','val',Q(1));
src = addSource([], wc(1), Q(1));
for i=2:length(wc)
    
Wells = addWell(Wells,G,rock,wc(i),'type','rate','val',Q(i));
src = addSource(src, wc(i), Q(i));
end
fluid2 = initSingleFluid('mu', 1, 'rho', 1);
state2  = initResSol (G, 1*barsa);
% src = addSource([], wc, Q);

[Ab, ~] = getSystemIncompTPFA(state2, G, T1, fluid2,'use_trans', true,  'bc', BC,'src',src);
[Ac, ~] = getSystemIncompTPFA(state2, G, T2, fluid2,'use_trans', true,  'bc', BC,'src',src);
initState.pressure=p_init;

state = initResSol(G, initState.pressure, 1);
    

pv_ref=poreVolume(G,rock);
pRef=0;
PA=cell(NT,1);
PA{1}=state;
state_sp=cell(NT,1);
state_sp{1}=zeros(TN,1);
tic
 for N=2:NT 
     pv = pv_ref .* (1 + (PA{N-1}.pressure-pRef)*CR);
   PA{N} = lincompTPFA(dt(N-1), PA{N-1}, G, T,pv, fluid,...
       rock, 'MatrixOutput', true,'wells', Wells,'bc',BC);
%    
%      pf2=Ab\(-Ac*PA{N}.pressure);
%   state_sp{N}=pf2(1)-pf2;
%    
 end
% 
 parfor N=2:NT
         pf2=Ab\(-Ac*PA{N}.pressure);
 state_sp{N}=pf2(1)-pf2;
 
 end
end 
