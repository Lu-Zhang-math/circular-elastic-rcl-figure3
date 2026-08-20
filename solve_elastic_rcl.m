function result=solve_elastic_rcl(omega,L,params)

rho0   = params.rho;
lambda = params.lambda;
mu     = params.mu;
eta0   = rho0*omega^2;

r0 = params.r0;
a0 = params.a;
b0 = params.b;

M   = params.M;
nPhysicalLayers=params.nPhysicalLayers;
N1  = params.N1;
N2  = params.N2;
N3  = params.N3;

kp = omega*sqrt(rho0/(lambda+2*mu));
ks = omega*sqrt(rho0/mu);

tau0=L/(b0-a0);

%% Mesh and connectivity
nLayers=nPhysicalLayers+1;
Npxvec=(N3+1)*ones(1,nLayers*M);
Npyvec=[(N1+1)*ones(1,nPhysicalLayers*M),(N2+1)*ones(1,M)];
[Vert,Ele,xml,yml]=circular_mesh( ...
    r0,a0,b0,nPhysicalLayers,M,Npxvec,Npyvec);

Nel = size(Ele,1);
Nv  = numel(Vert.x);

ElasticElements=1:(nPhysicalLayers*M);
PotentialElements=nPhysicalLayers*M+(1:M);

InterfaceElasticElements=(nPhysicalLayers-1)*M+(1:M);
InterfacePotentialElements=nPhysicalLayers*M+(1:M);

[rx,sx,ry,sy,xr,xs,yr,ys,J,Js]=geometric_factors_2d(xml,yml);

GnuEdge=element_connectivity(Ele);
IfInv = edge_orientations(Ele);

[GnuEdgeSplit,EleSplit,NvSplit] = ...
    split_interface( ...
    GnuEdge,Ele,Nv,PotentialElements);

[Ngrid,Node_Edge,G] = assemble_global_nodes( ...
    GnuEdgeSplit,EleSplit,NvSplit, ...
    Npxvec,Npyvec,IfInv);

GNodeIn = Node_Edge(1);
Ydim = Ngrid-GNodeIn;
Nub = 2*Ydim;

%% Element matrices
A11 = cell(1,Nel);
A12 = cell(1,Nel);
A21 = cell(1,Nel);
A22 = cell(1,Nel);

for iel = 1:Nel
    if mod(iel-1,M)~=0
        continue;
    end

    Npx = Npxvec(iel);
    Npy = Npyvec(iel);

    [zx,wx]=gll_nodes_weights(Npx-1);
    [zy,wy]=gll_nodes_weights(Npy-1);
    Dr=differentiation_matrix(zy);
    Ds=differentiation_matrix(zx);

    W  = kron(diag(wy),diag(wx));
    Is = eye(Npx);
    Ir = eye(Npy);

    if ismember(iel,ElasticElements)
        [A11{iel},A12{iel},A21{iel},A22{iel}] = ...
            elastic_matrix( ...
            lambda,mu,rho0,omega, ...
            rx{iel},sx{iel},ry{iel},sy{iel}, ...
            J{iel},Dr,Ds,W,Is,Ir);
    else
        A11{iel}=rcl_matrix( ...
            iel,kp,eta0,a0,b0,tau0, ...
            xml,yml,xr,xs,yr,ys,J,Js, ...
            Npxvec,Npyvec);

        A22{iel}=rcl_matrix( ...
            iel,ks,eta0,a0,b0,tau0, ...
            xml,yml,xr,xs,yr,ys,J,Js, ...
            Npxvec,Npyvec);

    end
end

%% Static condensation
FBD = complex(zeros(Nub,1));

recoveryInfo=cell(1,Nel);
condensationCache=cell(1,nLayers);
schurRows=cell(1,Nel);
schurColumns=cell(1,Nel);
schurValues=cell(1,Nel);

for iel = 1:Nel
    Npx = Npxvec(iel);
    Npy = Npyvec(iel);

    [~,~,linearIndIN,~,~,linearIndBD] = ...
        local_node_indices(Npx,Npy);

    Xdim = numel(linearIndIN);
    Bdim = numel(linearIndBD);
    GlobalSub=G{iel};
    GlobalIndBD=GlobalSub(linearIndBD)-GNodeIn;
    globalBoundaryDofs=reshape( ...
        [2*GlobalIndBD(:).'-1;2*GlobalIndBD(:).'],[],1);

    radialLayer=ceil(iel/M);
    angularElement=mod(iel-1,M)+1;

    if angularElement==1
        alpha=0;
        if ismember(iel,PotentialElements)
            KpII=A11{iel}(linearIndIN,linearIndIN);
            KpIB=A11{iel}(linearIndIN,linearIndBD);
            KpBI=A11{iel}(linearIndBD,linearIndIN);
            KpBB=A11{iel}(linearIndBD,linearIndBD);

            KsII=A22{iel}(linearIndIN,linearIndIN);
            KsIB=A22{iel}(linearIndIN,linearIndBD);
            KsBI=A22{iel}(linearIndBD,linearIndIN);
            KsBB=A22{iel}(linearIndBD,linearIndBD);

            recoveryP=KpII\KpIB;
            recoveryS=KsII\KsIB;
            localSchur=zeros(2*Bdim,2*Bdim);
            localSchur(1:2:end,1:2:end)=KpBB-KpBI*recoveryP;
            localSchur(2:2:end,2:2:end)=KsBB-KsBI*recoveryS;
            condensationCache{radialLayer}=struct( ...
                'recoveryP',recoveryP,'recoveryS',recoveryS, ...
                'schur',localSchur);
        else
            HII=zeros(2*Xdim,2*Xdim);
            HIB=zeros(2*Xdim,2*Bdim);
            HBI=zeros(2*Bdim,2*Xdim);
            HBBLocal=zeros(2*Bdim,2*Bdim);

            HII(1:2:end,1:2:end)=A11{iel}(linearIndIN,linearIndIN);
            HII(1:2:end,2:2:end)=A12{iel}(linearIndIN,linearIndIN);
            HII(2:2:end,1:2:end)=A21{iel}(linearIndIN,linearIndIN);
            HII(2:2:end,2:2:end)=A22{iel}(linearIndIN,linearIndIN);

            HIB(1:2:end,1:2:end)=A11{iel}(linearIndIN,linearIndBD);
            HIB(1:2:end,2:2:end)=A12{iel}(linearIndIN,linearIndBD);
            HIB(2:2:end,1:2:end)=A21{iel}(linearIndIN,linearIndBD);
            HIB(2:2:end,2:2:end)=A22{iel}(linearIndIN,linearIndBD);

            HBI(1:2:end,1:2:end)=A11{iel}(linearIndBD,linearIndIN);
            HBI(1:2:end,2:2:end)=A12{iel}(linearIndBD,linearIndIN);
            HBI(2:2:end,1:2:end)=A21{iel}(linearIndBD,linearIndIN);
            HBI(2:2:end,2:2:end)=A22{iel}(linearIndBD,linearIndIN);

            HBBLocal(1:2:end,1:2:end)=A11{iel}(linearIndBD,linearIndBD);
            HBBLocal(1:2:end,2:2:end)=A12{iel}(linearIndBD,linearIndBD);
            HBBLocal(2:2:end,1:2:end)=A21{iel}(linearIndBD,linearIndBD);
            HBBLocal(2:2:end,2:2:end)=A22{iel}(linearIndBD,linearIndBD);

            localRecovery=HII\HIB;
            localSchur=HBBLocal-HBI*localRecovery;
            condensationCache{radialLayer}=struct( ...
                'recovery',localRecovery,'schur',localSchur, ...
                'angle',atan2(mean(yml{iel}(:)),mean(xml{iel}(:))));
        end
    else
        cached=condensationCache{radialLayer};
        if ismember(iel,PotentialElements)
            alpha=0;
            localSchur=cached.schur;
        else
            elementAngle=atan2(mean(yml{iel}(:)),mean(xml{iel}(:)));
            alpha=angle(exp(1i*(elementAngle-cached.angle)));
            rotation=[cos(alpha),-sin(alpha);sin(alpha),cos(alpha)];
            transformBoundary=kron(speye(Bdim),rotation);
            localSchur=transformBoundary*cached.schur* ...
                transformBoundary';
        end
    end

    recoveryInfo{iel}=struct('dofs',globalBoundaryDofs, ...
        'layer',radialLayer,'alpha',alpha);

    [rowIndex,columnIndex]=ndgrid(globalBoundaryDofs);
    schurRows{iel}=rowIndex(:);
    schurColumns{iel}=columnIndex(:);
    schurValues{iel}=localSchur(:);
end

HBB=sparse(vertcat(schurRows{:}),vertcat(schurColumns{:}), ...
    vertcat(schurValues{:}),Nub,Nub);
clear A11 A12 A21 A22 schurRows schurColumns schurValues
clear KpII KpIB KpBI KpBB KsII KsIB KsBI KsBB recoveryP recoveryS
clear HII HIB HBI HBBLocal localRecovery localSchur

%% Interface terms
for kk=1:numel(InterfaceElasticElements)
    ielE=InterfaceElasticElements(kk);
    ielP=InterfacePotentialElements(kk);

    NpxE=Npxvec(ielE);
    NpyP=Npyvec(ielP);

    jE=1;
    jP=NpyP;

    xEedge=xml{ielE}(:,jE);
    yEedge=yml{ielE}(:,jE);

    edgeE=(1:NpxE).';
    edgeP=edgeE+(jP-1)*NpxE;

    [Ux,Uy]=boundary_trace(ielE,edgeE,G,GNodeIn,Nub);
    [Vp,Vs]=boundary_trace(ielP,edgeP,G,GNodeIn,Nub);

    thetaEdge=atan2(yEedge,xEedge);
    nx=cos(thetaEdge);
    ny=sin(thetaEdge);
    tx=-sin(thetaEdge);
    ty=cos(thetaEdge);

    Nx=spdiags(nx,0,NpxE,NpxE);
    Ny=spdiags(ny,0,NpxE,NpxE);
    T1mat=spdiags(tx,0,NpxE,NpxE);
    T2mat=spdiags(ty,0,NpxE,NpxE);

    Un=Nx*Ux+Ny*Uy;
    Ut=T1mat*Ux+T2mat*Uy;

    [zeta,edgeWeights]=gll_nodes_weights(NpxE-1);
    Dedge=differentiation_matrix(zeta);
    edgeJacobian=Js{ielE}(:,jE);
    Dt=spdiags(1./edgeJacobian,0,NpxE,NpxE)*Dedge;
    Wedge=spdiags( ...
        edgeWeights.*edgeJacobian,0,NpxE,NpxE);

    Kappa=(1/a0)*speye(NpxE);
    Qn=Dt*Ut+Kappa*Un;
    Qt=Dt*Un-Kappa*Ut;

    Wp=exp(1i*kp*a0);
    Ws=exp(1i*ks*a0);
    Phi=Wp*Vp;
    Psi=Ws*Vs;
    DtPhi=Wp*(Dt*Vp);
    DtPsi=Ws*(Dt*Vs);

    HBB=HBB ...
        +eta0*(Un'*Wedge*Phi-Ut'*Wedge*Psi) ...
        +2*mu*(Un'*Wedge*Qn-Ut'*Wedge*Qt) ...
        +eta0*(Phi'*Wedge*(Un-DtPsi) ...
        +Psi'*Wedge*(DtPhi-Ut));
end

for iel=1:M
    Npx=Npxvec(iel);
    Npy=Npyvec(iel);
    jInner=Npy;

    for i=1:Npx
        id=G{iel}(i,jInner)-GNodeIn;
        x=xml{iel}(i,jInner);
        y=yml{iel}(i,jInner);

        [uxBoundary,uyBoundary]=exact_field( ...
            x,y,kp,ks,params.source);

        HBB(2*id-1,:)=0;
        HBB(2*id,:)=0;
        HBB(2*id-1,2*id-1)=1;
        HBB(2*id,2*id)=1;

        FBD(2*id-1)=uxBoundary;
        FBD(2*id)=uyBoundary;
    end
end

Araw=HBB;
braw=FBD;

rowScale=full(max(abs(Araw),[],2));
rowScale(rowScale==0)=1;

Ascaled=spdiags(1./rowScale,0,Nub,Nub)*Araw;
bscaled=braw./rowScale;

ub=Ascaled\bscaled;

scaledResidual=norm(Ascaled*ub-bscaled)/max(norm(bscaled),realmin);
denominator=norm(Araw,1)*norm(ub,1)+norm(braw,1);
backwardError=norm(Araw*ub-braw,1)/max(denominator,realmin);

u=cell(1,Nel);

for iel=1:Nel
    Npx=Npxvec(iel);
    Npy=Npyvec(iel);

    u{iel}=zeros(2,Npx,Npy);

    [rowsubIN,colsubIN,linearIndIN, ...
     rowsubBD,colsubBD,linearIndBD]= ...
        local_node_indices(Npx,Npy);

    info=recoveryInfo{iel};
    cache=condensationCache{info.layer};
    boundaryValues=ub(info.dofs);
    if ismember(iel,PotentialElements)
        uinner=zeros(2*numel(linearIndIN),1);
        uinner(1:2:end)=-cache.recoveryP*boundaryValues(1:2:end);
        uinner(2:2:end)=-cache.recoveryS*boundaryValues(2:2:end);
    else
        rotation=[cos(info.alpha),-sin(info.alpha); ...
            sin(info.alpha),cos(info.alpha)];
        transformInterior=kron(speye(numel(linearIndIN)),rotation);
        transformBoundary=kron(speye(numel(linearIndBD)),rotation);
        uinner=-transformInterior*(cache.recovery* ...
            (transformBoundary'*boundaryValues));
    end

    for q=1:numel(linearIndIN)
        u{iel}(1,rowsubIN(q),colsubIN(q))=uinner(2*q-1);
        u{iel}(2,rowsubIN(q),colsubIN(q))=uinner(2*q);
    end

    for q=1:numel(linearIndBD)
        id=G{iel}(rowsubBD(q),colsubBD(q))-GNodeIn;
        u{iel}(1,rowsubBD(q),colsubBD(q))=ub(2*id-1);
        u{iel}(2,rowsubBD(q),colsubBD(q))=ub(2*id);
    end
end

relativeH1=relative_h1_error( ...
    u,ElasticElements,Npxvec,Npyvec,xml,yml,rx,sx,ry,sy,J, ...
    kp,ks,params.source);

result=struct();
result.omega=omega;
result.L=L;
result.uRelH1=relativeH1;
result.scaledResidual=scaledResidual;
result.backwardError=backwardError;
end

function A=rcl_matrix( ...
        iel,k,eta0,a0,b0,tau0, ...
        xml,yml,xr,xs,yr,ys,J,Js, ...
        Npxvec,Npyvec)

Npx=Npxvec(iel);
Npy=Npyvec(iel);

[zx,wx]=gll_nodes_weights(Npx-1);
[zy,wy]=gll_nodes_weights(Npy-1);

Dr=differentiation_matrix(zy);
Ds=differentiation_matrix(zx);

Is=eye(Npx);
Ir=eye(Npy);

DrFull=kron(Dr,Is);
DsFull=kron(Ir,Ds);

W=kron(diag(wy),diag(wx));

rr=hypot(xml{iel},yml{iel});
theta=atan2(yml{iel},xml{iel});

tau=exponential_map(rr,a0,tau0);

Arr=1./(rr*tau0);
Att=rr*tau0;

gamma=rr./tau;
gamma_r=gamma.*(1./rr-tau0);

ct=cos(theta);
st=sin(theta);
rEta=ct.*xr{iel}+st.*yr{iel};
thetaXi=(-st.*xs{iel}+ct.*ys{iel})./rr;
elementJ=J{iel};

Mradial=diag(reshape( ...
    gamma.*Arr.*elementJ./rEta.^2,[],1))*W;
Mtangential=diag(reshape( ...
    gamma.*Att.*elementJ./(rr.^2.*thetaXi.^2),[],1))*W;
Mgradient=diag(reshape( ...
    Arr.*gamma_r.*elementJ./rEta,[],1))*W;
Mphase=diag(reshape(elementJ./rEta,[],1))*W;
Mzero=diag(reshape((1./rr-tau0).*elementJ,[],1))*W;

K=DrFull'*Mradial*DrFull+DsFull'*Mtangential*DsFull;
P=Mgradient*DrFull;
C=1i*k*(DrFull'*Mphase-Mphase*DrFull);
Z=1i*k*Mzero;

A=eta0*(K+P+C+Z);

jOuter=1;
gammaB=gamma(:,jOuter);
mR=tau(:,jOuter)/b0;
beta=ones(Npx,1);
boundaryWeight=gammaB.*mR.*beta;

MRobin=diag(wx(:).*Js{iel}(:,jOuter).*boundaryWeight(:));
outer=1:Npx;
A(outer,outer)=A(outer,outer)-eta0*1i*k*MRobin;
end

function tau=exponential_map(r,a0,tau0)
tau=a0*exp(tau0*(r-a0));
end

function [GnuEdgeNew,EleNew,NvNew]= ...
    split_interface(GnuEdge,Ele,Nv,PotentialElements)

GnuEdgeNew=GnuEdge;
EleNew=Ele;

edgeLocal=2;
vertLocal=[3,4];

oldEdgeIDs=GnuEdge(PotentialElements,edgeLocal);
uniqueOldEdges=unique(oldEdgeIDs(:));

NedOld=max(GnuEdge(:));
newEdgeIDs=NedOld+(1:numel(uniqueOldEdges));

for k=1:numel(PotentialElements)
    iel=PotentialElements(k);
    pos=find(uniqueOldEdges==GnuEdge(iel,edgeLocal),1);
    GnuEdgeNew(iel,edgeLocal)=newEdgeIDs(pos);
end

oldVertIDs=Ele(PotentialElements,vertLocal);
uniqueOldVerts=unique(oldVertIDs(:));
newVertIDs=Nv+(1:numel(uniqueOldVerts));

for k=1:numel(PotentialElements)
    iel=PotentialElements(k);

    for lv=vertLocal
        pos=find(uniqueOldVerts==Ele(iel,lv),1);
        EleNew(iel,lv)=newVertIDs(pos);
    end
end

NvNew=Nv+numel(uniqueOldVerts);
end

function [T1,T2]=boundary_trace( ...
        iel,localIndices,G,GNodeIn,Nub)

globalNodes=G{iel}(localIndices)-GNodeIn;
rows=(1:numel(localIndices)).';
values=ones(size(rows));
T1=sparse(rows,2*globalNodes-1,values,numel(rows),Nub);
T2=sparse(rows,2*globalNodes,values,numel(rows),Nub);
end

function [Dx,Dy]=physical_derivatives( ...
        iel,Npxvec,Npyvec,rx,sx,ry,sy)

Npx=Npxvec(iel);
Npy=Npyvec(iel);

[zx,~]=gll_nodes_weights(Npx-1);
[zy,~]=gll_nodes_weights(Npy-1);

Dr=differentiation_matrix(zy);
Ds=differentiation_matrix(zx);

Is=eye(Npx);
Ir=eye(Npy);

DrFull=kron(Dr,Is);
DsFull=kron(Ir,Ds);

Dx=diag(rx{iel}(:))*DrFull+diag(sx{iel}(:))*DsFull;
Dy=diag(ry{iel}(:))*DrFull+diag(sy{iel}(:))*DsFull;
end

function [ux,uy,ux_x,ux_y,uy_x,uy_y]= ...
    exact_field(x,y,kp,ks,zsrc)

dx=x-zsrc(1);
dy=y-zsrc(2);
r=hypot(dx,dy);

nx=dx./r;
ny=dy./r;

phi=besselh(0,1,kp*r);
psi=besselh(0,1,ks*r);

phi_r=-kp*besselh(1,1,kp*r);
psi_r=-ks*besselh(1,1,ks*r);

phi_rr=kp*besselh(1,1,kp*r)./r-kp^2*phi;
psi_rr=ks*besselh(1,1,ks*r)./r-ks^2*psi;

phi_x=phi_r.*nx;
phi_y=phi_r.*ny;
psi_x=psi_r.*nx;
psi_y=psi_r.*ny;

phi_xx=phi_rr.*nx.^2+(phi_r./r).*(1-nx.^2);
phi_xy=(phi_rr-phi_r./r).*nx.*ny;
phi_yy=phi_rr.*ny.^2+(phi_r./r).*(1-ny.^2);

psi_xx=psi_rr.*nx.^2+(psi_r./r).*(1-nx.^2);
psi_xy=(psi_rr-psi_r./r).*nx.*ny;
psi_yy=psi_rr.*ny.^2+(psi_r./r).*(1-ny.^2);

ux=phi_x+psi_y;
uy=phi_y-psi_x;

ux_x=phi_xx+psi_xy;
ux_y=phi_xy+psi_yy;
uy_x=phi_xy-psi_xx;
uy_y=phi_yy-psi_xy;
end

function [A11,A12,A21,A22]=elastic_matrix( ...
        lambda,mu,rho0,omega,rx,sx,ry,sy,J,Dr,Ds,W,Is,Ir)

AG11=((lambda+2*mu)*(rx.*rx)+mu*(ry.*ry)).*J;
AG12=((lambda+2*mu)*(rx.*sx)+mu*(ry.*sy)).*J;
AG21=AG12;
AG22=((lambda+2*mu)*(sx.*sx)+mu*(sy.*sy)).*J;

BG11=(mu*(ry.*rx)+lambda*(ry.*rx)).*J;
BG12=(lambda*(rx.*sy)+mu*(ry.*sx)).*J;
BG21=(lambda*(sx.*ry)+mu*(sy.*rx)).*J;
BG22=(mu*(sy.*sx)+lambda*(sy.*sx)).*J;

CG11=((lambda+2*mu)*(ry.*ry)+mu*(rx.*rx)).*J;
CG12=((lambda+2*mu)*(ry.*sy)+mu*(rx.*sx)).*J;
CG21=CG12;
CG22=((lambda+2*mu)*(sy.*sy)+mu*(sx.*sx)).*J;

DG11=(mu*(ry.*rx)+lambda*(ry.*rx)).*J;
DG12=(lambda*(ry.*sx)+mu*(rx.*sy)).*J;
DG21=(lambda*(sy.*rx)+mu*(sx.*ry)).*J;
DG22=(mu*(sy.*sx)+lambda*(sy.*sx)).*J;

massMatrix=kron(Ir',Is')* ...
    diag(reshape(rho0.*J,[],1))*W*kron(Ir,Is);

A11= ...
    kron(Dr',Is')*diag(AG11(:))*W*kron(Dr,Is)+ ...
    kron(Dr',Is')*diag(AG12(:))*W*kron(Ir,Ds)+ ...
    kron(Ir',Ds')*diag(AG21(:))*W*kron(Dr,Is)+ ...
    kron(Ir',Ds')*diag(AG22(:))*W*kron(Ir,Ds)- ...
    omega^2*massMatrix;

A12= ...
    kron(Dr',Is')*diag(BG11(:))*W*kron(Dr,Is)+ ...
    kron(Dr',Is')*diag(BG12(:))*W*kron(Ir,Ds)+ ...
    kron(Ir',Ds')*diag(BG21(:))*W*kron(Dr,Is)+ ...
    kron(Ir',Ds')*diag(BG22(:))*W*kron(Ir,Ds);

A22= ...
    kron(Dr',Is')*diag(CG11(:))*W*kron(Dr,Is)+ ...
    kron(Dr',Is')*diag(CG12(:))*W*kron(Ir,Ds)+ ...
    kron(Ir',Ds')*diag(CG21(:))*W*kron(Dr,Is)+ ...
    kron(Ir',Ds')*diag(CG22(:))*W*kron(Ir,Ds)- ...
    omega^2*massMatrix;

A21= ...
    kron(Dr',Is')*diag(DG11(:))*W*kron(Dr,Is)+ ...
    kron(Dr',Is')*diag(DG12(:))*W*kron(Ir,Ds)+ ...
    kron(Ir',Ds')*diag(DG21(:))*W*kron(Dr,Is)+ ...
    kron(Ir',Ds')*diag(DG22(:))*W*kron(Ir,Ds);

end

function relativeH1=relative_h1_error( ...
        field,ElasticElements,Npxvec,Npyvec,xml,yml, ...
        rx,sx,ry,sy,J,kp,ks,zsrc)

uErrorL2=0;
uExactL2=0;
uErrorGrad=0;
uExactGrad=0;
for iel=ElasticElements
    Npx=Npxvec(iel);
    Npy=Npyvec(iel);

    [~,wx]=gll_nodes_weights(Npx-1);
    [~,wy]=gll_nodes_weights(Npy-1);
    weight=kron(wy(:),wx(:)).*J{iel}(:);

    uxNum=squeeze(field{iel}(1,:,:));
    uyNum=squeeze(field{iel}(2,:,:));

    [uxExact,uyExact,uxXExact,uxYExact,uyXExact,uyYExact]= ...
        exact_field( ...
        xml{iel},yml{iel},kp,ks,zsrc);

    [Dx,Dy]=physical_derivatives( ...
        iel,Npxvec,Npyvec,rx,sx,ry,sy);

    uxXNum=reshape(Dx*uxNum(:),Npx,Npy);
    uxYNum=reshape(Dy*uxNum(:),Npx,Npy);
    uyXNum=reshape(Dx*uyNum(:),Npx,Npy);
    uyYNum=reshape(Dy*uyNum(:),Npx,Npy);

    ex=uxNum-uxExact;
    ey=uyNum-uyExact;

    exX=uxXNum-uxXExact;
    exY=uxYNum-uxYExact;
    eyX=uyXNum-uyXExact;
    eyY=uyYNum-uyYExact;

    uErrorL2=uErrorL2+sum(weight.*(abs(ex(:)).^2+abs(ey(:)).^2));
    uExactL2=uExactL2+sum(weight.*( ...
        abs(uxExact(:)).^2+abs(uyExact(:)).^2));

    uErrorGrad=uErrorGrad+sum(weight.*( ...
        abs(exX(:)).^2+abs(exY(:)).^2+ ...
        abs(eyX(:)).^2+abs(eyY(:)).^2));

    uExactGrad=uExactGrad+sum(weight.*( ...
        abs(uxXExact(:)).^2+abs(uxYExact(:)).^2+ ...
        abs(uyXExact(:)).^2+abs(uyYExact(:)).^2));

end

relativeH1=sqrt((uErrorL2+uErrorGrad)/ ...
    max(uExactL2+uExactGrad,realmin));
end

function [Vert,Ele,xml,yml]=circular_mesh( ...
    r0,a,b,nPhysicalLayers,M,Npxvec,Npyvec)

radii=[linspace(r0,a,nPhysicalLayers+1),b];
theta=(0:M-1)*(2*pi/M);
nRadii=numel(radii);
Vert.x=zeros(M*nRadii,1);
Vert.y=zeros(M*nRadii,1);
for j=1:nRadii
    ids=(j-1)*M+(1:M);
    Vert.x(ids)=radii(j)*cos(theta);
    Vert.y(ids)=radii(j)*sin(theta);
end

Ele=zeros((nRadii-1)*M,4);
for j=1:nRadii-1
    for m=1:M
        mp=mod(m,M)+1;
        e=(j-1)*M+m;
        inner=(j-1)*M+[m,mp];
        outer=j*M+[m,mp];
        Ele(e,:)=[outer,fliplr(inner)];
    end
end

Nel=size(Ele,1);
xml=cell(1,Nel);
yml=cell(1,Nel);
for e=1:Nel
    j=ceil(e/M);
    m=mod(e-1,M)+1;
    [xi,~]=gll_nodes_weights(Npxvec(e)-1);
    [eta,~]=gll_nodes_weights(Npyvec(e)-1);
    angle0=(m-1)*2*pi/M;
    angle1=m*2*pi/M;
    angles=0.5*((1-xi)*angle0+(1+xi)*angle1);
    radius=0.5*((1-eta)*radii(j+1)+(1+eta)*radii(j));
    [Theta,R]=ndgrid(angles,radius);
    xml{e}=R.*cos(Theta);
    yml{e}=R.*sin(Theta);
end
end

function [z,w]=gll_nodes_weights(p)
n=p+1;
if p==1
    z=[-1;1];
    w=[1;1];
    return;
end

z=cos(pi*(0:p)'/p);
P=zeros(n,n);
for iteration=1:100
    zOld=z;
    P(:,1)=1;
    P(:,2)=z;
    for k=2:p
        P(:,k+1)=((2*k-1)*z.*P(:,k)-(k-1)*P(:,k-1))/k;
    end
    z=zOld-(z.*P(:,p+1)-P(:,p))./((p+1)*P(:,p+1));
    if max(abs(z-zOld))<5e-15
        break;
    end
end
P(:,1)=1;
P(:,2)=z;
for k=2:p
    P(:,k+1)=((2*k-1)*z.*P(:,k)-(k-1)*P(:,k-1))/k;
end
w=2./(p*(p+1)*P(:,p+1).^2);
[z,order]=sort(z,'ascend');
w=w(order);
end

function D=differentiation_matrix(x)
x=x(:);
n=numel(x);
lambda=ones(n,1);
for j=1:n
    lambda(j)=1/prod(x(j)-x([1:j-1,j+1:n]));
end
D=zeros(n,n);
for i=1:n
    for j=1:n
        if i~=j
            D(i,j)=lambda(j)/(lambda(i)*(x(i)-x(j)));
        end
    end
    D(i,i)=-sum(D(i,:));
end
end

function GnuEdge=element_connectivity(Ele)
Nel=size(Ele,1);
GnuEdge=zeros(Nel,4);
localPairs=[1,2;4,3;1,4;2,3];
edgeMap=containers.Map('KeyType','char','ValueType','int32');
nextEdge=0;
for e=1:Nel
    for le=1:4
        vertices=Ele(e,localPairs(le,:));
        key=sprintf('%d_%d',min(vertices),max(vertices));
        if isKey(edgeMap,key)
            edgeId=double(edgeMap(key));
        else
            nextEdge=nextEdge+1;
            edgeId=nextEdge;
            edgeMap(key)=int32(edgeId);
        end
        GnuEdge(e,le)=edgeId;
    end
end
end

function IfInv=edge_orientations(Ele)
Nel=size(Ele,1);
IfInv=false(Nel,4);
localPairs=[1,2;4,3;1,4;2,3];
orientationMap=containers.Map('KeyType','char','ValueType','any');
for e=1:Nel
    for le=1:4
        pair=Ele(e,localPairs(le,:));
        key=sprintf('%d_%d',min(pair),max(pair));
        if ~isKey(orientationMap,key)
            orientationMap(key)=pair;
        else
            IfInv(e,le)=isequal(pair,fliplr(orientationMap(key)));
        end
    end
end
end

function [Ngrid,Node_Edge,G]=assemble_global_nodes( ...
    GnuEdge,Ele,Nv,Npxvec,Npyvec,IfInv)

Nel=size(GnuEdge,1);
Ned=max(GnuEdge(:));
Ninside=0;
for e=1:Nel
    Ninside=Ninside+(Npxvec(e)-2)*(Npyvec(e)-2);
end

edgeNodeCount=zeros(Ned,1);
for e=1:Nel
    for le=1:4
        edge=GnuEdge(e,le);
        if le<=2
            count=Npxvec(e)-2;
        else
            count=Npyvec(e)-2;
        end
        if edgeNodeCount(edge)==0
            edgeNodeCount(edge)=count;
        end
    end
end

Node_Edge=zeros(Ned+1,1);
Node_Edge(1)=Ninside;
for edge=1:Ned
    Node_Edge(edge+1)=Node_Edge(edge)+edgeNodeCount(edge);
end

Ngrid=Node_Edge(end)+Nv;
G=cell(1,Nel);
nextInterior=0;
for e=1:Nel
    Npx=Npxvec(e);
    Npy=Npyvec(e);
    G{e}=zeros(Npx,Npy);
    count=(Npx-2)*(Npy-2);
    if count>0
        G{e}(2:Npx-1,2:Npy-1)=reshape( ...
            nextInterior+(1:count),Npx-2,Npy-2);
    end
    nextInterior=nextInterior+count;
end

for e=1:Nel
    Npx=Npxvec(e);
    Npy=Npyvec(e);
    for le=1:4
        edge=GnuEdge(e,le);
        ids=Node_Edge(edge)+1:Node_Edge(edge+1);
        if IfInv(e,le)
            ids=fliplr(ids);
        end
        switch le
            case 1
                G{e}(2:Npx-1,1)=reshape(ids,Npx-2,1);
            case 2
                G{e}(2:Npx-1,Npy)=reshape(ids,Npx-2,1);
            case 3
                G{e}(1,2:Npy-1)=ids;
            case 4
                G{e}(Npx,2:Npy-1)=ids;
        end
    end
end

vertexOffset=Node_Edge(end);
for e=1:Nel
    Npx=Npxvec(e);
    Npy=Npyvec(e);
    for lv=1:4
        id=vertexOffset+Ele(e,lv);
        switch lv
            case 1
                G{e}(1,1)=id;
            case 2
                G{e}(Npx,1)=id;
            case 3
                G{e}(Npx,Npy)=id;
            case 4
                G{e}(1,Npy)=id;
        end
    end
end
end

function [rx,sx,ry,sy,xr,xs,yr,ys,J,Js]=geometric_factors_2d(xml,yml)
Nel=size(xml,2);
rx=cell(1,Nel);
sx=cell(1,Nel);
ry=cell(1,Nel);
sy=cell(1,Nel);
xr=cell(1,Nel);
xs=cell(1,Nel);
yr=cell(1,Nel);
ys=cell(1,Nel);
J=cell(1,Nel);
Js=cell(1,Nel);
for iel=1:Nel
    [Npx,Npy]=size(xml{iel});
    [zx,~]=gll_nodes_weights(Npx-1);
    [zy,~]=gll_nodes_weights(Npy-1);
    Dr=differentiation_matrix(zy);
    Ds=differentiation_matrix(zx);
    xs{iel}=Ds*xml{iel};
    xr{iel}=xml{iel}*Dr';
    ys{iel}=Ds*yml{iel};
    yr{iel}=yml{iel}*Dr';
    detJ=xr{iel}.*ys{iel}-xs{iel}.*yr{iel};
    J{iel}=abs(detJ);
    rx{iel}=ys{iel}./detJ;
    sx{iel}=-yr{iel}./detJ;
    ry{iel}=-xs{iel}./detJ;
    sy{iel}=xr{iel}./detJ;
    Js{iel}=hypot(xs{iel},ys{iel});
end
end

function [rowsubIN,colsubIN,linearIndIN, ...
    rowsubBD,colsubBD,linearIndBD]=local_node_indices(Npx,Npy)

middleRows=repmat([1,Npx],1,Npy-2);
middleColumns=repelem(2:Npy-1,2);
rowsubBD=[1:Npx,middleRows,1:Npx];
colsubBD=[ones(1,Npx),middleColumns,Npy*ones(1,Npx)];
linearIndBD=sub2ind([Npx,Npy],rowsubBD,colsubBD);

[rowsubIN,colsubIN]=ndgrid(2:Npx-1,2:Npy-1);
rowsubIN=rowsubIN(:).';
colsubIN=colsubIN(:).';
linearIndIN=sub2ind([Npx,Npy],rowsubIN,colsubIN);
end
