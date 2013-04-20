function RhTimeCourse = RhTrajectoryGenerator(Condition)

% RhTrajectoryGenerator
%
% Generate stochastic trajectory for Rhodopsin activity.  Activity is
% assumed to decline through a specified number of shutoff steps,
% each modeled as a first-order process.  The activity level of each step
% is determined by a recursion relation
%       Rh(n) = Rh(0) * DecayFact^n
% where DecayFact is a user-defined variable.
% 
% Output is vector of Rh activity
% Input is structure with following fields:
%   Condition.EpochPts: total points
%   Condition.InitialShutoffRate: sets rate of initial shutoff step -
%           steps have rate of InitialShutoffRate/Rh(n)
%   Condition.NumSteps: number of shutoff steps
%   Condition.RhDecayFact: decay factor for above recursion relation
%   Condition.Verbose: generate error messages as appropriate
%
%   Created 11/11 FMR

% initial settings for rhodopsin activity
RhTimeCourse(1:Condition.EpochPts) = 0;
CurrentStep = 0;
CurrentCatalyticActivity = 1;
ShutoffRate = Condition.InitialShutoffRate;
	
% for each time point decide whether shutoff reaction has occurred.
% if it has, update catalytic activity and shutoff rate.
for cnt=1:Condition.EpochPts
    RhTimeCourse(cnt) = CurrentCatalyticActivity;
    % generate random number between 0 and 1 (uniform dist) and compare to shutoff rate
    if (rand(1) < ShutoffRate)
        CurrentStep = CurrentStep + 1;
        CurrentCatalyticActivity = 1 - CurrentStep * Condition.RhDecayFact/Condition.RhSteps;
        ShutoffRate = Condition.InitialShutoffRate * CurrentCatalyticActivity;
    end
    % are we done yet?
    if (CurrentStep == Condition.RhSteps)
        break;
    end
end	
RhTimeCourse = RhTimeCourse * Condition.GainFact;

if (Condition.Verbose)
    if (cnt == Condition.EpochPts)
        fprintf(1, 'Warning - Rh trajectory too long\n');
    end
end