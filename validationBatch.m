classdef validationBatch
    %validationBatch    Validates given batch of data.
    %   Primarily developed for the purpose of finding outlier candidates
    %   via ellipsoidal pealing in data. This object does also invalid data
    %   search and frozen value detection.


    %% Authors:
    % Markus Neuvonen (MNe)
    % University of Oulu
    % email: markus.neuvonen@oulu.fi
    %
    %% Revision history:
    % 15th of September 2022, v1, MNe: Initial version.
    %
    %% BEGIN CODE


    properties
        numericBatch %Only numeric values of given data
        metaBatch %Table of all other info given data
    end %Visible properties

    properties (Hidden = true)
        validMethods = {'remove','zoh'} %List of valid fault methods
        faultMethod = 'remove' %How to treat invalid data (NaN, text etc.)
                               %'remove': Delete invalid data rows
                               %'zoh': Keep latest acceptable value
        faultIndices
        faultSamples
        nSameEntries = 10 %How many values in row means "stuck"?
        relTolEntries = 1e-6 %How close values are considered "same"?
        stuckIndices
    end %Hidden properties


    methods (Hidden = true)

        function obj = validationBatch(givenData)
            %Constructor method. Check & initialize system.
            assert(nargin==1,'Provide exactly one input.')
            assert(istable(givenData),'Batch data must be a table.')

            %CVX installation is required for outlier detection. Ensure an
            %installation is found:
            obj.initializeCvx();

            %Find faulty instances of numeric data:
            obj.numericBatch = givenData(:,vartype('numeric'));
            assert(~isempty(obj.numericBatch),...
                'Provided table contains no numeric data!')
            assert(width(obj.numericBatch)==width(table2array(obj.numericBatch)),...
                'Variables with multiple columns in data table NOT supported.')
            obj.faultIndices = isnan(obj.numericBatch.Variables) |...
                isinf(obj.numericBatch.Variables) ;
            obj.faultSamples = any(obj.faultIndices,2);
            assert(~all(obj.faultSamples),'All samples contain some faulty data!')

            %Find the metadata of given table:
            metaVariables = givenData.Properties.VariableNames;
            metaVariables(ismember(metaVariables,...
                obj.numericBatch.Properties.VariableNames)) = [];
            %Add also initialization of status variables:
            obj.metaBatch = [givenData(:,metaVariables) ...
                table(zeros(height(givenData),1),...
                'VariableNames', {'nOutlier'})];
            
            %Do data checks with initial parameters:
            obj = obj.setFaultyEntriesParams(obj.faultMethod);
            obj.stuckIndices = zeros(size(obj.numericBatch)); %Initialize...
            obj = obj.setStuckValueParams(obj.nSameEntries,obj.relTolEntries);
        end

    end %Hidden methods

    methods (Hidden = true, Static = true)

        function initializeCvx()
            %Make sure cvx is installed. Automatize running of cvx setup if
            %needed (if no admin-rights to computer, for example).
            % //20.9.22 MNe: Switched to method from subfunction to enable
            % more agile unittest of class.
            try
                %Try a "random" low-level cvx command
                cvx_clear
            catch
                cvxFolderLocation = 'C:/Program Files/MATLAB/cvx';
                prompt = ['Please enter CVX installation folder. \n Default folder is ',...
                    cvxFolderLocation, '. \n Press enter to use default folder. \n'];
                x = input(prompt,"s");
                if ~isempty(x)
                    cvxFolderLocation = x;
                end
                [oldFolderLocation] = attemptFolderChange(cvxFolderLocation);
                try
                    cvx_setup; %Set up the cvx then!
                catch ME
                    clc
                    cd(oldFolderLocation);
                    rethrow(ME);
                end
                clc
                cd(oldFolderLocation);
            end
        end

    end %Hidden static methods

    methods

        function [obj] = setFaultyEntriesParams(obj,desiredMethod)
            %Set the method determines what to do with corrupt data.
            if any(strcmp(desiredMethod,obj.validMethods))
                if ~strcmp(obj.faultMethod,desiredMethod)
                    obj.faultMethod = desiredMethod;
                    %Change of method causes recalculation of all ellisoids
                    obj.metaBatch.nOutlier = zeros(size(obj.metaBatch.nOutlier));
                end
            else
                disp('Invalid fault method given. Valid methods are:')
                disp(obj.validMethods)
            end
        end

        function [validData] = getFaultFreeData(obj)
            %Do the changes to numeric data according to method:
            switch obj.faultMethod
                case 'remove'
                    validData = obj.numericBatch(~obj.faultSamples,:);
                case 'zoh'
                    validData = obj.numericBatch; %Initialize
                    for iRow = reshape(find(obj.faultSamples),1,[]) %For-loop
                                          %needs row vector to work expectedly!
                        if iRow > 1 %Normal ZOH-procedure
                            validData(iRow,obj.faultIndices(iRow,:)) = ...
                                validData(iRow-1,obj.faultIndices(iRow,:));
                        else %Special if starting faulty
                            for jColumn = find(obj.faultIndices(1,:))
                                %Overwrite with first valid value
                                validData(1,jColumn) = ...
                                    obj.numericBatch(find(...
                                    ~obj.faultIndices(:,jColumn),1),jColumn);
                            end
                        end
                    end
            end
        end

        function [obj] = setStuckValueParams(obj,n,tol)
            %Update stuck value parameters (OF ORIGINAL numeric data)
            assert(nargin==3,...
                'Provide amount of numbers and relative tolerance as inputs.')
            assert(isnumeric(n) && isnumeric(tol),'Provide numeric inputs!')
            assert(mod(n,1)==0,'Amount must be integer!')

            obj.nSameEntries = n;
            obj.relTolEntries = tol;

            %Calculate which entries fail to change more than tolerance:
            sameValues = abs(diff(obj.numericBatch.Variables))./...
                min(abs(obj.numericBatch{1:end-1,:}),...
                    abs(obj.numericBatch{2:end,:}))    < tol;
            if any(sameValues(:))
                %re-initialize if changes found:
                obj.stuckIndices = zeros(size(obj.numericBatch));
                %Now find the locations of 0->1->0 changes at each variable:
                for iVariable = 1:width(obj.numericBatch)
                    f = find(diff([0 ; sameValues(:,iVariable) ; 0 ] == 1));
                    p = f(1:2:end-1);  % Start indices
                    y = f(2:2:end)-p;  % Consecutive ones counts
                    for jViolation = find(y>n) %Mark all violations:
                        obj.stuckIndices(p(jViolation):p(jViolation)+y(jViolation)-1)...
                            = ones(y(jViolation),1);
                    end
                end
            end
        end

        function T = getStuckValues(obj)
            %Return binary table (1 = stuck, 0 = OK) of numeric data.
            T = array2table(obj.stuckIndices,...
                'VariableNames',obj.numericBatch.Properties.VariableNames);
        end

        function [v,obj] = getOutlierCandidates(obj,n,absTol,plotOptions)
            %Apply (recursive) ellipsoidal pealing to detect possible
            %outliers
            arguments
                obj validationBatch
                n (1,1) double {mustBePositive,mustBeInteger}
                absTol (1,1) double {mustBeNonnegative} = 1e-5
                plotOptions.doPlot (1,1) logical = false
            end

            %Check if results are already saved in input object:
            histMax = max(obj.metaBatch.nOutlier);
            if n > histMax %Have to do new routine

                %Create vector for index mapping:
                indVec = reshape(1:length(obj.metaBatch.nOutlier),[],1);
                %And another internal helper-vector
                tempVec = obj.metaBatch.nOutlier;
                %Get the data:
                numericData = table2array(obj.getFaultFreeData());
                %Check if data has been removed:
                if length(indVec) ~= size(numericData,1)
                    indVec(obj.faultSamples) = [];
                    tempVec(obj.faultSamples) = [];
                    assert(length(indVec) == size(numericData,1),...
                        'Problem with removed numeric data')
                end
                %Do the required amount of ellipsoidal pealing:
                for iLayer = histMax+1:1:n
                    numericData(tempVec > 0,:) = [];
                    indVec(tempVec > 0) = [];
                    tempVec(tempVec > 0) = [];
                    assert(~isempty(numericData),'All data pealed already!')
                    
                    %Do the actual ellipsoidal pealing:
                    x = numericData'; %Detect COLUMN outliers, so transpose.
                    [nVars,nSamples] = size(x); %#ok<ASGLU> 
                    % Create and solve the model
                    cvx_begin quiet
                    variable A(nVars,nVars) symmetric
                    variable b(nVars)
                    maximize( det_rootn( A ) )
                    subject to
                    norms( A * x + b * ones( 1, nSamples ), 2 ) <= 1; %#ok<VUNUS>
                    cvx_end
                    %Find the candidates of outlying:
                    outlierIndices = find(norms( A*x+b*ones(1,nSamples),2 )...
                        > 1-absTol);
                    assert(~isempty(outlierIndices),['Outlier detection failed.' ...
                        'Try increasing tolerance or ' ...
                        'check your data for anomalities...']);
                    %Update internal helpers accordingly:
                    tempVec(outlierIndices) = 1;
                    %Map back to original information vector:
                    obj.metaBatch.nOutlier(indVec(outlierIndices)) = iLayer;
                end
            end
            %Write the results:
            v = (obj.metaBatch.nOutlier >  0) &...
                (obj.metaBatch.nOutlier <= n);
            %Plot the results:
            if plotOptions.doPlot
                switch nVars
                    case 2
                        figure();
                        noangles = 200;
                        angles   = linspace( 0, 2 * pi, noangles );
                        ellipse  = A \ [ cos(angles) - b(1) ; sin(angles) - b(2) ];
                        plot( x(1,:), x(2,:), 'ro', ellipse(1,:), ellipse(2,:), 'b-' );
                        plot( x(1,v), x(2,v), 'r*', x(1,~v), x(2,~v), 'bo', ...
                            ellipse(1,:), ellipse(2,:), 'k-' );
                        axis equal
                        legend('Outlier candidates','Not candidates', ...
                            'Found ellipsoid')
                        xlabel(obj.numericBatch.Properties.VariableNames{1})
                        ylabel(obj.numericBatch.Properties.VariableNames{2})
                    case 3
                        [X,Y,Z] = sphere; %Unit sphere coordinates
                        ellipse  = A \ [ (X(:))' - b(1) ;...
                                         (Y(:))' - b(2) ;...
                                         (Z(:))' - b(3) ]; %Map into ellipsoid
                        X = reshape(ellipse(1,:),[],21);%Reformulate
                        Y = reshape(ellipse(2,:),[],21);% into
                        Z = reshape(ellipse(3,:),[],21);% matrices.
                        figure();plot3(x(1,v),  x(2,v),  x(3,v), 'r*');hold on;
                                 plot3(x(1,~v), x(2,~v), x(3,~v), 'bo')
                           surf(X,Y,Z,'FaceColor','none');axis equal;hold off
                        legend('Outlier candidates','Not candidates', ...
                            'Found ellipsoid')
                        xlabel(obj.numericBatch.Properties.VariableNames{1})
                        ylabel(obj.numericBatch.Properties.VariableNames{2})
                        zlabel(obj.numericBatch.Properties.VariableNames{3})
                    otherwise
                        disp('Plotting possible only in 2D and 3D.')
                end
            end


        end

    end %non-static public methods


    methods (Static)

    end %public static methods


    methods (Access = private)

    end %private methods
end

%% Subfunctions

function [oldFolderLocation] = attemptFolderChange(newFolderLocation)
%attemptFolderChange   Change folder with custom error message.
%   Try to change to given folder and save the original folder as output.
try
    oldFolderLocation = cd(newFolderLocation);
catch ME
    if (strcmp(ME.identifier,'MATLAB:cd:NonExistentFolder'))
        msg = ['The folder location provided is incorrect. ', ...
            num2str(newFolderLocation), ' is not found.'];
    elseif (strcmp(ME.identifier,'MATLAB:string'))
        msg = ['The folder location provided is of incorrect type! ', ...
            num2str(newFolderLocation), ' is not valid string!'];
    end
    causeException = MException('MATLAB:attemptFolderChange:location',msg);
    ME = addCause(ME,causeException);
    rethrow(ME);
end
end
