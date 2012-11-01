%%
% scan across several parameters and assess likelihood of each given
% properties of experimental responses

% find the rieke-server job manager
jm = findResource('scheduler', 'configuration', 'rieke-server');

% create the an empty job on the job manager
% FileDependencies are files or folders that each of the computers on the 
% cluster will need to execute your function. Ideally you would keep these 
% small and to a minimum to reduce network chatter, but during the
% execution of your job they are only sent to each worker once even if that
% worker is used multiple times.
job = createJob(jm, 'FileDependencies', {'DoOneRhCascadeModel.m','RhCascadeModel.m'});

NumSteps =  [1 2 3 4 5 6 7 8 9 10 11 12 13]; %[8 10];
RhTCon = [200 250 240 210]; %[200 250];
RhDecayFact = [0.6 1 1.2 0.7]; %[0.6 1];
load WTSingles36
SimParams.NumResponses = 5000;       % with 5000 resps simulation std ~1%, about 10% std in likelihood
SimParams.Resp = Target.AveSingle;
SimParams.FreqCutoff = 10;                                
NumLoops = 1;
Verbose = 1;

tic
clear SimReturns;
i = 0;
for loop = 1:NumLoops
    
    for step = 1:length(NumSteps)

        for tcon = 1:length(RhTCon)

            for decay = 1:length(RhDecayFact)
                
                % create a task inside our job
                % each task created will be sent to one worker (core) on
                % the cluster. Ideally you would create a number of tasks
                % >= the number of workers available on the cluster (26).
                createTask(job, ...                  % the job to add the task to
                           @DoOneRhCascadeModel, ... % function to be run
                           1, ...                    % number of output arguments returned
                           {NumSteps(step), RhDecayFact(decay), RhTCon(tcon), SimParams}); % input arguments
                           
                %SimReturns(step, tcon, decay) = DoOneRhCascadeModel(NumSteps(step), RhDecayFact(decay), RhTCon(tcon), SimParams);
                
                i = i + 1;
            end

        end

    end

end

% submit the job to the cluster
submit(job);

% wait for the job to be complete
% normally you wouldn't do this, but we want to time to job with tic/toc
% you can also get the executation time of a job without waiting for it
waitForState(job);

toc
disp(num2str(i));

% retrieve the results from the job
% this will return a cell array (or matrix) of all the output arguments
SimReturns = getAllOutputArguments(job);

% destroy the job
destroy(job);



%%