%
% get mean single photon response right
% time to peak should be 200-210 ms.
%rod_diffusion_PDE_WT(3, 7, 0.8, 0.2, 1, 1) % steep dependence on RK rate


%%
% rod_diffusion_PDE_WT(100, 3, 1, 0.5, 1)
% rod_diffusion_PDE_WT(100, 3, 1, 1, 1)
% rod_diffusion_PDE_WT(100, 3, 1, 1.5, 1)
% rod_diffusion_PDE_WT(100, 3, 1, 2, 1)
% 
% %%
%rod_diffusion_PDE_WT(200, 3, 2, 2, 1)
%rod_diffusion_PDE_WT(200, 3, 2, 2, 1)
%rod_diffusion_PDE_WT(200, 3, 2, 2.5, 1)
% 
% %%

%%
%rod_diffusion_PDE_WT(30, 3, 0.5, 1.75, 2)
%rod_diffusion_PDE_WT(30, 3, 1e-8, 1.75, 3)

%%
%rod_diffusion_PDE_WT(30, 5, 0.5, 0.5, 2)
%rod_diffusion_PDE_WT(30, 5, 1e-8, 0.5, 3)

%%

% TRateMult = [0.125 0.25 0.5 1 2]; 
% % 0.5 and 1 quite unlikely (5e-3 and 1e-4, compared to 0.4 or so for small TRate values), others ok for 5 or 7 steps, 0.5 corresponds to AmpRatio of 1.35, 1 to AmpRatio of 1.2
% % for 3 steps 1 and 2 unlikely, 0.5 ok
% load AveSingles;
% RKPeakAmp = max(AveSingles.RKAveSingle);
% 
% AmpRatio = AveSingles.CSMPeakAmp / RKPeakAmp;
% AmpRatioSEM = AveSingles.CSMPeakAmpSEM / RKPeakAmp;
% 
% for cond = 1:length(TRateMult)
%     RKResults = rod_diffusion_PDE_WT(60, 5, 0.5, TRateMult(cond), 2);
%     CSMResults = rod_diffusion_PDE_WT(60, 5, 1e-8, TRateMult(cond), 3);
%     TestAmpRatio = CSMResults.PeakAmp / RKResults.PeakAmp;
%     Likelihood(cond) = normpdf(abs(TestAmpRatio - AmpRatio) / AmpRatioSEM);
%     fprintf(1, 'TRate %d TargetAmpRatio %d AmpRatio %d Likelihood %d\n', TRateMult(cond), AmpRatio, TestAmpRatio, Likelihood(cond));
% end
% 
% figure(1);
% plot(TRateMult, Likelihood, 'o');

%% 
% these look reasonable for CV and variance
%rod_diffusion_PDE_WT(300, 3, 2, 2, 1)
%Results = rod_diffusion_PDE_WT(200, 5, 1, 0.5, 1) % steep dependence on RK rate
%rod_diffusion_PDE_WT(300, 7, 0.75, 0.15, 1)
% time to peak of variance increases as saturation increases (as expected)

% variance width also appears to increase with increasing saturation- for 7
% steps this is noticeable (and hurts fits) when activation rate above 0.3.

% time to peak of variance can be brought down by increasing RK rate - for
% 7 steps something around 0.7 seems best

%Results_02 = rod_diffusion_PDE_WT(300, 7, 0.75, 0.15, 1) % steep dependence on RK rate
%Results_03 = rod_diffusion_PDE_WT(300, 7, 0.75, 0.3, 1) % steep dependence on RK rate
%Results_03 = rod_diffusion_PDE_WT(300, 7, 0.75, 0.45, 1) % steep dependence on RK rate

%%
%Results_03 = rod_diffusion_PDE_WT(300, 7, 0.4, 0.2, 1) % steep dependence on RK rate
%Results_03 = rod_diffusion_PDE_WT(300, 7, 0.6, 0.2, 1) % steep dependence on RK rate
%Results_03 = rod_diffusion_PDE_WT(300, 7, 0.8, 0.2, 1) % steep dependence on RK rate


%Results_03 = rod_diffusion_PDE_WT(3, 7, 0.8, 0.2, 1) % steep dependence on RK rate


%Results_1 = rod_diffusion_PDE_WT(300, 7, 30, 1, 1, 0, 1) % steep dependence on RK rate ***
%Results_2 = rod_diffusion_PDE_WT(300, 7, 50, 1, 1, 0, 1) % steep dependence on RK rate
%Results_3 = rod_diffusion_PDE_WT(300, 7, 70, 1, 1, 0, 1) % steep dependence on RK rate

%Results_4 = rod_diffusion_PDE_WT(300, 7, 50, 0, 1, 0, 1) % steep dependence on RK rate
%Results_5 = rod_diffusion_PDE_WT(300, 7, 50, 0.2, 1, 0, 1) % steep dependence on RK rate
%Results_6 = rod_diffusion_PDE_WT(300, 7, 50, 0.4, 1, 0, 1) % steep dependence on RK rate
%Results_7 = rod_diffusion_PDE_WT(300, 7, 50, 0.8, 1, 0, 1) % steep dependence on RK rate
%Results_8 = rod_diffusion_PDE_WT(300, 7, 50, 1, 1, 0, 1) % steep dependence on RK rate ***

%Results_1 = rod_diffusion_PDE_WT(300, 7, 20, 1, 1, 0, 1) % 
%Results_2 = rod_diffusion_PDE_WT(300, 7, 30, 1, 3, 0, 1) % 
%Results_3 = rod_diffusion_PDE_WT(300, 7, 30, 1, 10, 0, 1) %  **
