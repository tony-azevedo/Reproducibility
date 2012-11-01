function Results = rod_diffusion_PDE_WT(nsim, NumSteps, InitialShutoffRate, DecayFact, GainFact, DeterministicModel, FigIndex)

% Revision history:
% FMR 11/11
%   moved some parameters to command line so we could explore them more
%   systematically
%
%   modified several parameters to make more consistent with our other
%   simulations.  Likely these reflect changes I made in rhodopsin model:
% kRK same as before 10.5 * kRK_mult; %s^(-1)
% nu_RG=100 * nu_RG_mult;  was 330 s^(-1)
% dtfactor = 1; was (1/60)*(1/0.07); adjustment to give correct time course
% of SPR variance, again probably reflects changes in Rh model (as
% described in that function)

%D  - cGMP diffusion constant
%D_Ca - Ca2+ diffusion coefficient
%D_PDE - PDE* diffusion coefficient
%R  - disk radius
%br  - annlular radius (between disk rim and OS membrane)
%dl - height is inter-diskal compartment (between disks)
%     and thickness of disk as well
%L - full length of the rod os
%nr - number of divisions of disk radius
      %(nr - 1/2)*dr=R
%     each discretized element is an annulus except for central circle
%nd - number of inter-diskal compartments to be computed (extended by
%     symmetry)
%alpha_dark - dark rate of synthesis of cGMP
         %(dark conc. of cGMP per s)
%%[kRK,      beta_dark, GC_gain_factor, beta_sub, tau_PDE, gammaCa, k2,
%a, b,kf, kGAPf, kGAPr, eta, k_G_PDE];
%x=[525.4191    2.8362    7.7000    2.3988    0.1691   22.2336    5.8895 ...
%    0.6431   13.4928   35.8120    0.0053    0.0180    0.7735   27.1804];

%-------------------------
% experimental target data
%-------------------------
load WTSingles36

%nsim=50;%number of simulated SPRs
tau_PDE=0.2;%time constant for PDE* decay
%tau_PDE=0.05; %for RGS9 OX simulation
nu_PDE=1/tau_PDE;

%-------------------------
% parameters controlling calcium feedback to cyclase
%-------------------------

GC_gain_factor=7.7;%factor by which GC activity is increaed at minimum Ca2+
%The following 4 parameters are for Ca buffering and pumpimg
a=0.6431;
b=13.4928;
gamma_Ca=22.2336;
k2=5.8895;


Ca_dark=0.250; % (mu M) dark Ca2+ concentration, [Ca2+] (nominal)
Ca_min=0.020; % (mu M) minimum value to which [Ca2+] can fall
%m_Ca=4; %Hill coefficient for Ca dependence of gunylate cyclase rate
m_Ca=2.2; %Hill coefficient for Ca dependence of gunylate cyclase rate
K_GC=Ca_dark*( (1-GC_gain_factor*(Ca_min/Ca_dark)^m_Ca)/...
               (GC_gain_factor-1) )^(1/m_Ca); % (mu M) parameter
     %determining dependence of guanlyate cylclase activity on [Ca2+]
% FMR: check above - should it be 1/m_Ca in both places (from inverting Hill relation)

y_star=K_GC/Ca_dark; %cyclase activation parameter, dimensionless
y_min=Ca_min/Ca_dark;% dimensionless
f_Ca=0.15;%fraction of light-sensitive current carried by Ca2+
z_Ca=2;% valence of Ca


%-------------------------
% parameters controlling diffusion
%-------------------------

%D=150; %diffusion constant for cGMP, sq. micron/s
D=100; 
D_Ca=15;%%diffusion constant for Ca, sq. micron/s
%D_PDE=5;
D_PDE=2.5;%diffusion constant for PDE, sq. micron/s

R=0.7; %disk radius (microsn) as given for mouse rodes in Bisegna et al.
dl=0.015; % thickness of disk and spacing between disks (microns)
L=20; % lenght of rod OS, microns
br=dl;%radial distance between disk edge and cell membrane (microns)

nr=50;%number of discretized values of radial disk coordinate r
nd=200;%number of disks to be modeled; only the active one is in detatil
dr=R/(nr-1/2); % distance between each simulated disk compartment
r=( (1:nr)-1/2 )*dr; %radii of concentric circles
rplot=r-dr/2; % values of r at the center of each discretized disk element
xd=(0:(nd-1))'*2*dl;


alpha_dark=2.8362; %same as beta_dark - dark PDE activity in sec-1
Pdark=alpha_dark;%rate of hydrolysis of normalized cGMP in dark


%-------------------------
% simulated concentrations of diffusing components
%-------------------------

cr=ones(nr,1); %cGMP conc. as function of r in active disk compartment
yr=cr;         %Ca2+ ...
ur=cr;         %bound Ca2+ buffer (spatially fixed) ...
P=Pdark*cr;    %PDE*. All PDE* emanates from a centrlal small spot on disc (cGMP hydrolysis rate)
alpha_vec_r=alpha_dark*ones(nr,1); % PDE enzymatic activity (not hydrolysis rate b/c only # active enzymes)
cd=ones(nd,1); % cGMP
yd=cd;       % Ca2+
ud=cd;       % bound Ca2+
c=[cr;cd];
y=[yr;yd];
u=[ur;ud];

dt=8*(dr^2/D);%time step
dt=4*dt; % 4 
T=1.2;% simulation time

nt=round(T/dt);     % number of time points

t=(0:(nt-1))*dt;

Mr=zeros(nr); %diffusion matrix operator for radial cGMP
Md=zeros(nd); %... for longitudingal cGMP
Myr=Mr;       %... for radial Ca
Myd=Md;       %... for longitudingal Ca
MPr=zeros(nr);%... for  radial PDE*

n=nr+nd;        % total number of compartments

M=zeros(n); %overall diffusion operator for cGMP, active disc
            %   and the other lumped interdiscal compartments
My=M;       %.... for free Ca2+
C=ones(n,nt);
Y=C;
U=C;

%-------------------------
% parameters controlling stochastic Rh activity
%-------------------------
RhCondition.EpochPts = nt;
RhCondition.InitialShutoffRate = InitialShutoffRate*dt;
RhCondition.NumSteps = NumSteps;
RhCondition.RhDecayFact = DecayFact;
RhCondition.GainFact = GainFact;
RhCondition.Verbose = 0;

Mr(1,1:2)=(4)*[-1 1];
for j=2:(nr-1)
    Mr(j,(j-1):(j+1))=[r(j-1) -(r(j)+r(j-1)) r(j)]/(r(j-1)+dr/2);
end
Mr(nr,(nr-1):nr)=[r(nr-1) -( r(nr-1)+r(nr)*2*dr/(dr+br) )]/(r(nr-1)+dr/2);
Myr=Mr;
Mr=Mr*D/dr^2;
Myr=Myr*D_Ca/dr^2;

for jd=2:(nd-1)
    Md(jd,(jd-1):(jd+1))=[1 -2 1];
end
Md(nd,(nd-1):nd)=[1 -1];
Myd=Md;
Md = ( ( (br*D/R)/dl^2 )/(1+4*br/R) )*Md;
Myd = ( ( (br*D_Ca/R)/dl^2 )/(1+4*br/R) )*Myd;

Md(1,1:2)=(D/dl^2)*[-( dl^2/(br*(dr+br)) + 1/2 ) 1/2];
Myd(1,1:2)=(D_Ca/dl^2)*[-( dl^2/(br*(dr+br)) + 1/2 ) 1/2];

alpha_vec_d=alpha_dark*[0;ones(nd-1,1)];
P_vec_d=Pdark*[0;ones(nd-1,1)];



M(1:nr,1:nr)=Mr;
M((nr+1):n,(nr+1):n)=Md;
M(nr,nr+1)=((D/dr^2)*r(nr)*(2*dr/(dr+br))/(r(nr-1)+dr/2));
M(nr+1,nr)=D/(br*(dr+br));
M=sparse(M);

My(1:nr,1:nr)=Myr;
My((nr+1):n,(nr+1):n)=Myd;
My(nr,nr+1)=((D_Ca/dr^2)*r(nr)*(2*dr/(dr+br))/(r(nr-1)+dr/2));
My(nr+1,nr)=D_Ca/(br*(dr+br));
My=sparse(My);

In=speye(n);

Pfixed=Pdark*cr;

PfixD=spdiags([Pfixed;P_vec_d/(1+4*br/R)],0,n,n);

b_vec=[alpha_vec_r;alpha_vec_d/(1+4*br/R)];

C(:,1)=c;

MPr=zeros(nr);


MPr(1,1:2)=(4)*[-1 1];
for j=2:(nr-1)
    MPr(j,(j-1):(j+1))=[r(j-1) -(r(j)+r(j-1)) r(j)]/(r(j-1)+dr/2);
end
%For no-flux boundary condition at r=R:
MPr(nr,(nr-1):nr)=[r(nr-1) -r(nr-1)]/(r(nr-1)+dr/2);

MPr=MPr*D_PDE/dr^2;

PR=zeros(nr,nt);
nuPR_vec=zeros(nr,1);
nuPR_vec(1)=Pdark*R^2/(r(1)^2);
APr=( MPr - nu_PDE*eye(nr) );
BPr=sparse( eye(nr)-dt*APr );
Pdiag=spdiags([PR(:,1);zeros(n-nr,1)],0,n,n);

Ydiag=spdiags((a*k2/(b-1))*ones(n,1),0,n,n);
c_vec_y=[zeros(nr,1);gamma_Ca*[( R/(4*br) +1 );ones(nd-1,1)]];
YfixD=spdiags( c_vec_y,0,n,n );

ind=sub2ind(size(M),1:nr,1:nr);
jsub=ind2sub([n n],ind);

indY=sub2ind(size(M),1:n,1:n);
jsubY=ind2sub([n n],indY);

indU=indY;
jsubU=jsubY;

Udiag=spdiags(ones(n,1),0,n,n);



%%%%HAMER MODEL
%kRK=525.4191;
%np=6;
%omega=0.6;
%kRK=kRK*exp((np-6)*0.6);
%kRK_n=kRK*exp(-omega*(0:np));
%mu_T_n=1./kRK_n;
%%%%%%%%

%nu_Arr=60; %s^(-1)  % NOTE: this is no longer doing anything with changes FMR made to Rh model
%nu_Arr=nu_Arr*1e-6; %for simulating Arr-/-
%kRK=10.5 * kRK_mult; %s^(-1)
%kRK=15 * kRK_mult; %s^(-1)
%kRK = kRK * 0.5;    % for simulating RK+/-
%kRK=kRK*1e-8;  %for simulating RK-/-
%nu_RG=330 * nu_RG_mult; %330 s^(-1) - can turn this down to reduce compression
%nu_RG=100 * nu_RG_mult; %330 s^(-1) - can turn this down to reduce compression
%np=4;

%dtfactor = 1;
%dtfactor=(1/60)*(1/0.07);%Dan's time scale adjustment to give correct time
                         %course of SPR variance
            
% nu_Arr=nu_Arr*dtfactor;
% kRK=kRK*dtfactor;
% nu_RG=nu_RG*dtfactor;
% 
% Results.np = np;
% Results.kRK = kRK;
% Results.nu_RG = nu_RG;

Ihat_mat=zeros(nsim,nt);%matrix of normalized photocurrents

tb=cputime;

% simulate each response
tic

for ksim=1:nsim
    if (mod(ksim, 5) == 0)
        fprintf(1, '.');
        pause(0.1);
    end
    
    %T_random=exprnd(mu_T_n); %for random Hamer R* model
    %nuRP=2*nu_R_PDE_random(t,T_random,kRK,np); %for random Hamer R* model
    %nuRP=2*nu_R_PDE(t,kRK,np); %for deterministic Hamer R*

    %nuRPfactor=1.4*8;
    %nuRP=nuRPfactor*nu_R_PDE_random_caruso(t,nu_RG,kRK,nu_Arr,np);%random caruso mod
    %nuRPfactor=1.4*14;%slow random caruso mod
    %nuRP=nuRPfactor*nu_R_PDE_random_caruso(t,nu_RG,kRK,nu_Arr,np);%slow random caruso mod
    if (DeterministicModel)
        nuRP = zeros(1, RhCondition.EpochPts);
        for epoch = 1:100
            nuRP = nuRP + RhTrajectoryGenerator(RhCondition);
        end
    else
        nuRP = RhTrajectoryGenerator(RhCondition);
    end
    %nuRP=1.4*10*nu_R_PDE_caruso(t,nu_RG,kRK,nu_Arr,np);%determ caruso mod
    %nuRP=1.4*14*nu_R_PDE_caruso(t,nu_RG,kRK,nu_Arr,np);%slow determ caruso mod
    
    g_vec=ones(n,1); %w/o feedback
    
    % step across time
    for j=2:nt
        PR(:,j)=BPr\(PR(:,j-1) + dt*nuRP(j)*nuPR_vec);
        Ydiag(jsubY)=(a*k2/(b-1))*(b-U(:,j-1));
        Ay=My - YfixD - Ydiag;
        By=In-dt*Ay;
        c_vec_y_var=c_vec_y;
        c_vec_y_var((nr+1):n)=c_vec_y_var((nr+1):n).* ...
            ( (1-y_min)*C((nr+1):n,j-1).^3 + y_min );
        Y(:,j)=By\( Y(:,j-1) + dt*( a*k2*U(:,j-1) + c_vec_y_var) );
        Udiag(jsubU)=k2*Y(:,j)/(b-1);
        Au=-Udiag-k2*In;
        Bu=In-dt*Au;
        U(:,j)=Bu\( U(:,j-1) + dt*k2*b*Y(:,j)/(b-1) );
        Pdiag(jsub)=PR(:,j);
        A = M - PfixD - Pdiag;
        B=In-dt*A;
        g_vec=( 1+(1/y_star)^m_Ca )./( 1+(Y(:,j)/y_star).^m_Ca ); % w feedback
        C(:,j)=B\( C(:,j-1) + dt*b_vec.*g_vec );
    end

    CR=C(1:nr,:); %cGMP conc as function of r (rows) and t (columns)
    cr0=CR(1,:); %%cGMP conc at r=0 in active compartment
    CD=C((nr+1):n,:);%cGMP conc as function of compartment postition x (rows) and t (columns)
    cd0=CD(1,:); %cGMP conc at cell membrane in active compartment
    ID=CD.^3;
    I=( 2*sum(ID,1) + L/(2*dl) - 2*nd )/(L/(2*dl));
    Ihat=1-I;
    [Imax,Jmax]=max(Ihat);

    Ihat_mat(ksim,:)=Ihat;

end % for simulations, indexed by ksim
toc
tc=cputime-tb

if (DeterministicModel)
    figure(FigIndex); clf;
    plot(Ihat_mat(1, :));
    [MaxVal, MaxLoc] = max(Ihat_mat(1, :));
    fprintf(1, 'time to peak = %d\n', MaxLoc * dt);
else
    figure(FigIndex); clf;
    subplot(1, 2, 1); hold on
    for resp = 1:min(50, ksim)
        plot(Ihat_mat(resp, :));
    end

    v=var(Ihat_mat,0,1);  %time-dependent variance of the photocurrent
    m=mean(Ihat_mat,1); %time-dependent mean photocurrent
    subplot(1, 2, 2);
    plot(t,m.^2,'b-',t,v,'r-')
    [z,jj]=max(m);
    if(ksim>1)
        cv_t_peak=sqrt(v(jj))/m(jj);
    end
    A=(sum(Ihat_mat,2)-0.5*Ihat_mat(:,1)-0.5*Ihat_mat(:,end))*dt; %SPR area
    mA=mean(A); %mean SPR area
    vA=var(A); %variance of SPR are
    cvA=sqrt(vA)/mA; %coefficient of variation of SPR area

    fprintf(1, 'steps = %d T gain fact = %d CV = %d\n', NumSteps, GainFact, cvA);

    %Sub sample viariable timecourses to save on storage
    nn=floor(length(t)/floor(length(t)/1000));
    mm=floor(length(t)/nn);
    kk=(0:mm:nn*mm)+1;
    Ihat_mat=Ihat_mat(:,kk);
    t=t(kk);
    CR=CR(:,kk);
    CD=CD(:,kk);
    PR=PR(:,kk);
    m=m(kk);
    v=v(kk);

    % ttt=round(clock);
    % filename=['I_np5_rnd_normalPDE_Rslow_',date,'_',num2str(ttt(4)),...
    %     '_',num2str(ttt(5)),'_',num2str(ttt(6)),'.mat'];
    % 
    % save filename Ihat_mat...
    %     A t dt v m mA cvA D D_Ca D_PDE CR CD PR np nu_Arr kRK nu_RG m_Ca...
    %     nuRPfactor tau_PDE

    % likelihood parameters
    [MaxAmp, TPeakMean] = max(m);
    [MaxVar, TPeakVar] = max(v);
    Indices = find(v > 0.5*MaxVar);
    FirstTime = Indices(min(find(Indices < TPeakVar)));
    SecondTime = Indices(max(find(Indices > TPeakVar)));
    VarWidth = SecondTime - FirstTime;

    if (1)
        ParamsVec(1) = cvA - Target.MeanCVArea; 
        ParamsVec(2) = TPeakVar/TPeakMean - Target.MeanTPeakRatio;
        ParamsVec(3) = VarWidth/TPeakMean - Target.MeanVarWidth;
        ParamsProjs = ParamsVec * Target.ParamsEigVec;
        Likelihood = normpdf(abs(ParamsProjs(1)/Target.ParamsEigSEM(1)));
        Likelihood= Likelihood * normpdf(abs(ParamsProjs(2)/Target.ParamsEigSEM(2)));;
        Likelihood = Likelihood * normpdf(abs(ParamsProjs(3)/Target.ParamsEigSEM(3)));
        Results.Likelihood = Likelihood;

        fprintf(1, 'CV = %d TPeakVar = %d Width Var = %d Amp = %d\n', cvA, TPeakVar / TPeakMean, VarWidth / TPeakMean, MaxAmp);

        fprintf(1, 'Likelihood = %d\n', Likelihood);

        Likelihood = normpdf((abs(cvA - Target.MeanCVArea) / Target.SEMCVArea));
        fprintf(1, '\tCV Likelihood = %d\n', Likelihood);
        Likelihood = normpdf((abs(TPeakVar/TPeakMean - Target.MeanTPeakRatio) / Target.SEMTPeakRatio));
        fprintf(1, '\tTPeakVar Likelihood = %d\n', Likelihood);
        Likelihood = normpdf((abs(VarWidth/TPeakMean - Target.MeanVarWidth) / Target.SEMVarWidth));
        fprintf(1, '\tVarWidth Likelihood = %d\n', Likelihood);
    end

    Results.CVArea = cvA;
    Results.TPeakVar = TPeakVar / TPeakMean;
    Results.VarWidth = VarWidth / TPeakMean;
    Results.PeakAmp = MaxAmp;
end


