function dp_plus=nu_R_PDE_random_caruso(t,nu_RG,kRK,nu_Arr,np)

% Revision history:
%   11/11 FMR
% modified to make more similar to our other multistep model. Forced
% arrestin to always be final step after all phosphorylation completed,
% make activity and shutoff rate closely linked to product of activity and
% mean lifetime of a step are constant.
% Specific changes were redefining kG_n, introducing ShutoffRate
% and altering nu_Arr_n

%T_random - row vector of random dwell times in R_n state (n=0,..,np)


omega=0.5;
kRK_n=kRK*(np-(0:np));
kG_n = nu_RG * (1 - (0:np)/(np+1));
ShutoffRate_n = kRK * np * (1 - (0:np)/(np+1)); % FMR - introduced to create equal steps

%kG_n=nu_RG*exp(-omega*(0:np));%rate of G* production for each R_n state
%Arrestin binding competes with phosphorylation
%nu_Arr_n=nu_Arr*[0, 0, 0, ones(1,np-2)];%rate of arrestin binding, eqaul to
                                       %zero for 3 or fewer phosporylations
nu_Arr_n=nu_Arr*[zeros(1, np) 1];  % only bind arrestin after all phosphorylation completed
%mu_T_n=1./(nu_Arr_n + kRK_n);%mean lifetime of each phosphorylated state
mu_T_n=1./(ShutoffRate_n);%mean lifetime of each phosphorylated state
T_random=exprnd(mu_T_n);%random lifetime of each phosphorylated state
T=[0,cumsum(T_random)];%event times: R* poduction, phosphorylations, and
                       %arrestin binding
p_bernoulli=kRK_n./(kRK_n+nu_Arr_n);%probability of phosphorylation instead
                                    %of arrestin binding
N=binornd(1,p_bernoulli);%Bernoulli random tirals to decide between
                         %phosphorylation (N=1) or arrestin binding (N=0)
n_quench=sum(cumsum(N)==(1:(np+1))) + 1;%total number of R* states that
                                       %contribute to G* (PDE*) production

dp_plus=zeros(size(t));
for j=1:(n_quench)
    %dp_plus=dp_plus+kG_n(j)*( (t>T(j)*ones(size(t)))&(t<=T(j+1)*ones(size(t))) );
    dp_plus=dp_plus+kG_n(j)*( (t>T(j))&(t<=T(j+1)) ); %rate of PDE* production
end
