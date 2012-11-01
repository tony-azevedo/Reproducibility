function features =  DoOneCascadeModelTA(sim)
% see also DoOneRhCascadeModel  DoOneRhCascadeModel(NumSteps, RhDecayFact,RhTCon, SimParams)

% R*(t) decay time constant - units ms
sim.InitialShutoffRate = sim.SamplingInterval * sim.RhSteps / sim.RhShutoffFact;
if isfield(sim,'TransRate')
    sim.BaseTransRate = sim.TransRate * sim.SamplingInterval;          % cannot exceed 1 - decrease sampling interval to get higher rate
    sim.TransDecayTimeConst = sim.TransTConst / sim.SamplingInterval;
end

% run RhCascadeModelTA or RhPlusTransducinCascadeModelTA
% Both functions generate stochastic elements, average, produce filter,
% then filter responses
sim = sim.SimFunc(sim, [], sim.NumResponses);

AverageResponse = mean(sim.EpochData);
VarianceResponse = var(sim.EpochData);

CV = std(sum(sim.EpochData,2)) / mean(sum(sim.EpochData,2));

[MaxVal, SimMeanTPeak] = max(AverageResponse);
[MaxVal, SimVarTPeak] = max(VarianceResponse);
TPeakRatio = SimVarTPeak / SimMeanTPeak;

Indices = find(VarianceResponse > 0.5*MaxVal);
FirstTime = min(Indices);
SecondTime = max(Indices);
VarWidth = SecondTime - FirstTime;

t = sim.target;

Likelihood = normpdf((CV - t.MeanCVArea) / t.SEMCVArea);
Likelihood = Likelihood * normpdf((TPeakRatio - t.MeanTPeakRatio) / t.SEMTPeakRatio);
Likelihood = Likelihood * normpdf((VarWidth - t.MeanVarWidth) / t.SEMVarWidth);

features.CV = CV;
features.TPeakRatio = TPeakRatio;
features.VarWidth = VarWidth;
features.Likelihood = Likelihood;
features.SimulateCondition = stripSimulateCondition(sim);

% use this to return examples for storing
if isfield(sim,'storeExamples') && sim.storeExamples
    sim.AverageResponse = AverageResponse;
    sim.VarianceResponse = VarianceResponse;
    sim.Filter = Filt;
    features.SimulateCondition = sim;
end

function sim = stripSimulateCondition(sim)
sim = rmfield(sim,'EpochData');
field = fieldnames(sim);
temp = logical(size(field));
% for rmnames = 1:length(field)
%     temp(rmnames) = ~isempty(strfind(field{rmnames},'Examples'));
% end
for rmnames = 1:length(field)
    temp(rmnames) = length(sim.(field{rmnames}))>1 || ~isnumeric(sim.(field{rmnames}));
end
sim = rmfield(sim,field(temp));

 