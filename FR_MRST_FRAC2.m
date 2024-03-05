function [PA,G]=FR_MRST_FRAC2(NX,NY,NZ,L,W,H,fluid,p_init,dt,NT,wc,...
K_x,K_y,K_z,phi,CR,CF,Q,BC,Frac_l,Frac_a,Frac_poro,Frac_k)



[nx,ny,nz] = deal(NX,NY,NZ);
G = computeGeometry(cartGrid([nx ny nz],[L W H]));

TN=nx*ny*nz;



rock = makeRock(G, K_x*milli*darcy, phi);
rock.cr=CR;

TF=length(Frac_a);
G.faces.tags=zeros(G.faces.num,1);
apt = zeros(G.faces.num,1);
for f=1:TF


G.faces.tags(Frac_l(:,f)) = 1;
aperture = Frac_a(f);
apt(Frac_l(:,f)) = aperture;
    
end 
G = addhybrid(G,G.faces.tags > 0,apt);

hybridInd = find(G.cells.hybrid);

nCells = G.cells.num;


% The fracture permeability is computed from the parallel plate assumption,
% which states that the permeability is aperture squared divided by 12.
rock.perm(hybridInd,:) = Frac_k(1);
rock.poro(hybridInd) = Frac_poro(1);


T = computeTrans_DFM(G,rock,'hybrid',true);
[G,T2] = computeHybridTrans(G,T);



gravity reset off

Wells = addWell([],G,rock,wc(1),'type','rate','val',Q(1));
for i=2:length(wc)
    
Wells = addWell(Wells,G,rock,wc(2),'type','rate','val',Q(2));
end

state = initState(G,Wells,p_init);

pv_ref= poreVolume(G, rock);
pRef=0;
PA{1}=state;
NT=length(dt);

fluid.cf=CF;

PA{1}=state;
% plotCellData(G, states{end}.pressure./barsa,'EdgeColor','none')
for i=2:length(dt)+1
 pv = pv_ref .* (1 + (PA{i-1}.pressure-pRef)*CR);
PA{i} = lincompTPFA_mod_Frac(dt(i-1), PA{i-1}, G, T, pv, fluid, rock,...
    'MatrixOutput', true, 'wells', Wells, 'c2cTrans',T2,'bc', BC);
end
end
