function [PA,G]=FR_MRST_FRAC(NX,NY,NZ,L,W,H,fluid,p_init,dt,NT,wc,...
K_x,K_y,K_z,phi,CF,Q,Frac_l,Frac_a,Frac_poro,Frac_k,bc)



[nx,ny,nz] = deal(NX,NY,NZ);
G = computeGeometry(cartGrid([nx ny nz],[L W H]));

TN=nx*ny*nz;



G.rock = makeRock(G, K_x*milli*darcy, phi);
G.rock.cr=CF;

TF=length(Frac_a);
for f=1:TF

fracplanes(f).points = Frac_l(:,:,f);
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

TPFAoperators = setupEDFMOperatorsTPFA(G, G.rock, tol);


initState.pressure=p_init;


state = initResSol(G, initState.pressure, 1);
    


gravity reset off

model = WaterModel(G,G.rock, fluid);
 model = model.validateModel();
model.operators = TPFAoperators;
Wells = addWell([], G, G.rock,  wc(1),     ...
        'Type', 'rate', 'Val', Q(1), ...
        'Radius', 0.1, 'Name', 'P1');
for i=2:length(wc)
    
    Wells = addWell(Wells, G, G.rock,  wc(i),     ...
        'Type', 'rate', 'Val', Q(i), ...
        'Radius', 0.1, 'Name', 'P1');
end
% plotEDFMgrid(G);
% hold on;
% plotWell(G,Wells);



schedule = simpleSchedule(dt, 'bc', bc, 'W', Wells);
PA={};
PA{1}=state;
src = addSource([], wc(1), Q(1));
PA{1}=state;

tic
[ws, states, report] = simulateScheduleAD(state, model, schedule)
toc
figure(3)
% plotCellData(G, states{end}.pressure./barsa,'EdgeColor','none')
for i=2:length(states)+1
PA{i}=states{i-1};
end
end
