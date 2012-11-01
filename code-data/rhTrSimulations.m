classdef rhTrSimulations 
    properties
        parameters
        features
    end
    properties (Hidden = true)
        filename
        featureindices
        pnames
        fnames
        verbose
        saved
    end

    
    methods

        function self = rhTrSimulations(simStruct,simfunc)
            switch func2str(simfunc)
                case 'RhCascadeModelTA'
                    % find the file, point to file if not there
                    self.filename = 'RhSimulationDB';
                case 'RhPlusTransducinCascadeModelTA'
                    % find the file, point to file if not there
                    self.filename = 'RhTrSimulationDB';
            end
            self.verbose = true;
            self = makeDataBase(self,simStruct);
            self.verbose = false;
        end
    
        
        function [self,status] = writeSimulations(self)
            if exist(sprintf('%s_temp.mat',self.filename),'file');
                db = load(sprintf('%s_temp',self.filename));
                db = db.db;
            else 
                db = load(sprintf('%s',self.filename));
                db = db.db;
            end
            for li_sim = 1:numel(self.features.Likelihood)
                il_sim = linearIndex(self,li_sim);
                il_db = il_sim;
                for p = 1:length(self.pnames)
                    il_db(p) = find(db.parameters.(self.pnames{p})==self.parameters.(self.pnames{p})(il_sim(p)));
                end
                li_db = rhTrSimulations.dblinearIndex(db.features.Likelihood,il_db);
                for f = 1:length(self.fnames)
                    db.features.(self.fnames{f})(li_db) = self.features.(self.fnames{f})(li_sim);
                end
            end
            status = sum(isnan(db.features.Likelihood(:)));
            self.saved = true;
            if self.verbose
                fprintf('Saving status (0 is good): %d\n',status);
            end
            if ~status
                save(self.filename,'db');
                if exist(sprintf('%s_temp.mat',self.filename),'file')
                    delete(sprintf('%s_temp.mat',self.filename));
                end
            end
        end
        
        function l = simulationExists(self,sim)
            indices = simIndex(self,sim);
            li = linearIndex(self,indices);
            l = ~isnan(self.features.CV(li)) && sim.NumResponses <= self.features.NumResponses(li);
            if self.verbose && l
                fprintf('Simulation exists, Likelihood: %g\n', self.features.Likelihood(li));
            elseif self.verbose && sim.NumResponses > self.features.NumResponses(li)
                fprintf('Rerun simulation with # Respon = %d vs. %d\n',sim.NumResponses, self.features.NumResponses(li));
            end
        end
        
        function self = putSimulation(self,sim,features)
            indices = simIndex(self,sim);
            li = linearIndex(self,indices);
            features.NumResponses = sim.NumResponses;
            
            %keep the values with largest number of responses
            if isnan(self.features.NumResponses(li))...
                    || features.NumResponses > self.features.NumResponses(li)
                for f = 1:length(self.fnames);
                    self.features.(self.fnames{f})(li) = features.(self.fnames{f});
                end
                if (self.verbose)
                    fprintf('\tLikelihood = %d CV = %d TPeakRatio = %d VarWidth = %d\n', features.Likelihood, features.CV, features.TPeakRatio, features.VarWidth);
                end
            end
        end

        function features = getFeatures(self,sim)
            indices = simIndex(self,sim);
            li = linearIndex(self,indices);
            for f = 1:length(self.fnames);                
                features.(self.fnames{f}) = self.features.(self.fnames{f})(li);
            end
            if self.verbose
                dbstr = '';
                for p = 1:length(self.pnames)
                    dbstr = sprintf('%s%s: %g\n',dbstr,self.pnames{p},self.parameters.(self.pnames{p})(indices(p)));
                end
                fprintf(dbstr)
            end
        end
        
        function il_p = mostLikely(self)
            [val,li] = max(self.features.Likelihood(:));
            il_p = linearIndex(self,li);
            if (self.verbose)
                s = sprintf('Most Likely parameter: Likelihood = %d\n',self.features.Likelihood(li));
                for p = 1:length(self.pnames)
                    % , features.CV, features.TPeakRatio, features.VarWidth);
                    s = sprintf('%s\t%d %s: %g\n',s,il_p(p),self.pnames{p},self.parameters.(self.pnames{p})(il_p(p)));
                end
                fprintf(s);
            end
            
        end
        
        function features = featuresByIndices(self,indices)
            li = linearIndex(self,indices);
            for f = 1:length(self.fnames);
                features.(self.fnames{f}) = self.features.(self.fnames{f})(li);
            end
        end

        function sim = simByIndices(self,indices)
            for p = 1:length(self.pnames);
                sim.(self.pnames{p}) = self.parameters.(self.pnames{p})(indices(p));
            end
        end

        function lims = getFeatureLims(self,feature)
            lims = [min(self.features.(feature)(:)), max(self.features.(feature)(:))];
        end
        
        function self = verboseOn(self)
            self.verbose = true;
        end
        function self = verboseOff(self)
            self.verbose = false;
        end
        function l = isVerbose(self)
            l = self.verbose;
            if self.verbose
                fprintf('Database is printing\n');
            end
        end
        
    end
    
    methods (Hidden = true)
                
        function indices = simIndex(self,sim)
            indices = zeros(1,length(self.pnames));
            for p = 1:length(self.pnames)
                indices(p) = find(self.parameters.(self.pnames{p})==sim.(self.pnames{p}));
            end
            if sum(indices==0)
                error('Simulation is not in the simulation object database')
            end
        end
              
        function li = linearIndex(self,il)
            if length(il)>1
                li = il(1);
                d = size(self.features.CV);
                %   (l-1)(d3)(d2)(d1)+(k-1)(d2)(d1)+(j-1)(d1)+i
                for i = 2:length(il)
                    li = li + (il(i)-1)*prod(d(1:i-1));
                end
            else
                d = size(self.features.CV);
                d = d(1:ndims(self.features.CV));
                li = zeros(size(d));
                li(1) = mod(il,d(1));
                if li(1)==0
                    li(1)=d(1);
                end
                il = il-li(1);
                il = il/d(1);
                for i = 2:length(li)
                    li(i) = rem(il,d(i))+1;
                    il = floor((il-(li(2)-1))/d(i));
                end
                if length(li) < length(self.pnames)
                   li = [li ones(1, length(self.pnames)-length(li))];
                end
            end
        end
        
        function self = makeDataBase(self,sim)
            % load database of parameters and features
            if exist(sprintf('%s.mat',self.filename),'file');
                db = load(self.filename);
                db = db.db;
            else
                db.features.NumResponses = [];
                db.features.TPeakRatio = [];
                db.features.VarWidth = [];
                db.features.CV = [];
                db.features.Likelihood = [];
                
                fn = sort(fieldnames(sim));
                fn = fn(~strcmp(fn,'FreqCutoff'));
                fn = fn(~strcmp(fn,'NumResponses'));
                for n = 1:length(fn)
                    db.parameters.(fn{n}) = [];
                end
            end
            
            % expand the database
            % expand the parameters
            self.pnames = fieldnames(db.parameters);
            self.fnames = fieldnames(db.features);
            dbpdims = zeros(size(self.pnames));
            if size(dbpdims,1)>size(dbpdims,2)
                dbpdims = dbpdims';
            end
            newp = dbpdims;

            dbstr = 'Old database: ';
            dbstr2 = 'New database: ';
            for p = 1:length(self.pnames)
                dbpdims(p) = length(db.parameters.(self.pnames{p}));
                newp(p) = length(setdiff(sim.(self.pnames{p}),db.parameters.(self.pnames{p})));
                
                dbstr = sprintf('%s%dx',dbstr,dbpdims(p));
                dbstr2 = sprintf('%s%dx',dbstr2,newp(p));
            end
            if (self.verbose)
                fprintf('%s\n%s\n',dbstr(1:end-1),dbstr2(1:end-1));
            end
            
            % expand the database
            if sum(newp)                
                expsize = dbpdims+newp;
                % expsize = expsize(1:find(expsize-1,1,'last'));
                for p = 1:length(self.pnames)
                    inds{p} = (1:dbpdims(p));
                end
                % inds = inds(1:length(expsize));
                for f = 1:length(self.fnames)
                    expansion = nan(expsize);
                    if sum(dbpdims)
                        expansion(inds{1:ndims(expansion)}) = db.features.(self.fnames{f});
                    end
                    db.features.(self.fnames{f}) = expansion;
                end
                
                % make string 'order,:,:,:,...'
                orderstrcell{1} = 'order,';
                for i = 2:length(self.pnames)%ndims(expansion)
                    orderstrcell{i} = ':,';
                end
                % (2,3,4,...,1) to last singleton dimension
                pervec = circshift(1:length(orderstrcell),[0 1]);
                
                % rearange feature matrices
                for p = 1:length(self.pnames)
                    % add param values
                    db.parameters.(self.pnames{p}) = ...
                        [db.parameters.(self.pnames{p}),...
                        setdiff(sim.(self.pnames{p}),db.parameters.(self.pnames{p}))];
                    [db.parameters.(self.pnames{p}),order] = sort(db.parameters.(self.pnames{p}));
                    
                    orderstr = strcat(orderstrcell{:});
                    orderstr = orderstr(1:end-1);
                    if sum(order~=(1:length(db.parameters.(self.pnames{p}))))
                        reorderstr = sprintf('db.features.(self.fnames{n}) = db.features.(self.fnames{n})(%s)',orderstr);
                        for f = 1:length(self.fnames)
                            eval(reorderstr);
                        end
                    end
                    orderstrcell = orderstrcell(pervec);
                end
                save(sprintf('%s_temp',self.filename),'db');
                self.saved = false;
            end
            
            dbstr = sprintf('Running simulations for n values of:\n');
            % save those that are of importance now
            for p = 1:length(self.pnames);
                self.parameters.(self.pnames{p}) = sim.(self.pnames{p});
                self.featureindices{p} = ismember(db.parameters.(self.pnames{p}),self.parameters.(self.pnames{p}));
                
                % debug statement
                dbstr = [dbstr, sprintf(' %d - %s\n',sum(self.featureindices{p}),self.pnames{p})];
            end
            if (self.verbose)
                fprintf('%s\n',dbstr);
            end
            
            % self.featureindices = self.featureindices(1:length(size(db.features.(self.fnames{1}))));
            for f = 1:length(self.fnames)
                self.features.(self.fnames{f}) = db.features.(self.fnames{f})(self.featureindices{:});
            end
            
        end

    end
    
    methods (Static)
        
        function li = dblinearIndex(feature,il)
            if length(il)>1
                li = il(1);
                d = size(feature);
                %   (l-1)(d3)(d2)(d1)+(k-1)(d2)(d1)+(j-1)(d1)+i
                for i = 2:length(il)
                    li = li + (il(i)-1)*prod(d(1:i-1));
                end
            else
                d = size(feature);
                li = zeros(size(d));
                li(1) = mod(il,d(1));
                if li(1)==0
                    li(1)=d(1);
                end
                il = il-li(1);
                il = il/d(1);
                for i = 2:length(li)
                    li(i) = rem(il,d(i))+1;
                    il = (il-(li(2)-1))/d(2);
                end
            end
        end

    end
end

