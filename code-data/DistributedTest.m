%%
% scan across several parameters and assess likelihood of each given
% properties of experimental responses

NumSteps =  [8 10];
RhTCon = [200 250];
RhDecayFact = [0.6 1];
load GCAPSingles36
SimParams.NumResponses = 10000;       % with 5000 resps simulation std ~1%, about 10% std in likelihood
SimParams.Resp = Target.AveSingle;
SimParams.FreqCutoff = 10;                                
NumLoops = 1;
Verbose = 1;

tic
clear SimReturns;
for loop = 1:NumLoops
    
    for step = 1:length(NumSteps)

        for tcon = 1:length(RhTCon)

            for decay = 1:length(RhDecayFact)

                SimReturns(step, tcon, decay) = DoOneRhCascadeModel(NumSteps(step), RhDecayFact(decay), RhTCon(tcon), SimParams)
                                
            end

        end

    end

end
toc

%%