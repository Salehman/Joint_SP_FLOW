function [PA]=FR_MRST_MATRIX2(NX,NY,NZ,L,W,H,fluid,p_init,dt,NT,wc,...
K_x,K_y,K_z,phi,CR,CF,Q,BC)


[nx,ny,nz] = deal(NX,NY,NZ);
G = computeGeometry(cartGrid([nx ny nz],[L W H]));

TN=nx*ny*nz;



rock = makeRock(G, [K_x, K_x, K_x]*milli*darcy, phi);
rock.cr=CR;

fluid.cf=CF;
T = computeTrans(G, rock);
pR  = 200*barsa;
Wells = addWell([],G,rock,wc(1),'type','rate','val',Q(1));
for i=2:length(wc)
    
Wells = addWell(Wells,G,rock,wc(i),'type','rate','val',Q(i));
end
% src = addSource([], wc, Q);
initState.pressure=p_init;


state = initResSol(G, initState.pressure, 1);
    

pv_ref=poreVolume(G,rock);
pRef=0;
PA=cell(NT,1);
PA{1}=state;
 for N=2:NT 
     pv = pv_ref .* (1 + (PA{N-1}.pressure-pRef)*CR);
   PA{N} = lincompTPFA(dt(N-1), PA{N-1}, G, T,pv, fluid,...
       rock, 'MatrixOutput', true,'wells', Wells,'bc',BC);
   
 end
