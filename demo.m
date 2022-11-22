%Demo - Run this code to have a demonstration how to use Outlier Detector
%-tool. Two simple tables are created, and available functionality of the
%tool displayed.

%NOTE: Installation of cvx convex optimization tool is required to run!

clc;
%Create example table 1:
times = reshape(datetime(date())+(seconds(0):seconds(5):seconds(20)),[],1);
gdata = randn(size(times,1),size(times,2));
bdata = gdata;bdata(1)=NaN;bdata(end)=inf;bdata(floor(length(bdata)/2))="lol";
T = table(times,gdata,bdata) %display the created table.
myObject = validationBatch(T);

%% Faulty data detection.
% Display the data from which faulty data is taken away of:
myObject.getFaultFreeData()
%Change the behaviour, use zero-order-hold for faulty data:
myObject = myObject.setFaultyEntriesParams('zoh');
%Display results with new settings:
myObject.getFaultFreeData()

%% Stuck value detection:
%Modify the table; update the object to use this modified table:
T(2:3,3) = T(4,3);
T(2,2) = T(3,2)
myObject = validationBatch(T);
% Set detection of stuck values when more than 2 repetitions and with zero
% tolerance (i.e. exactly same value required):
myObject = myObject.setStuckValueParams(2,0);
% Show the indices of stuck values (binary):
myObject.getStuckValues()

% Create example table 2 and update the object to use it:
x = randn(100,3);
T = table(x(:,1),x(:,2),x(:,3));
myObject = validationBatch(T);
%% Find the first layer of outlier candidates and also visualize the results:
myObject.getOutlierCandidates(1,'doPlot',true)
%When making successive pealing and not wanting to re-calculate previous steps,
%the object must be saved between calls:
[~,myObject] = myObject.getOutlierCandidates(1,'doPlot',false);
myObject.getOutlierCandidates(2,'doPlot',true);
%Some deeper layer (shown if points still remain only):
myObject.getOutlierCandidates(10,'doPlot',true)
