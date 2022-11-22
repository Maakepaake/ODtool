classdef testValidationBatch < matlab.unittest.TestCase

    properties
        oldfolder
    end

    properties (TestParameter)
        tTable = struct(...
            "small",getMyTestData(1), ...
            "medium",getMyTestData(2))
        Answers = struct(...
            "small",getMyTestAnswers(1), ...
            "medium",getMyTestAnswers(2))
    end

    methods (TestClassSetup)
        function initializeTests(testCase)
            testCase.oldfolder = cd('..');
            validationBatch.initializeCvx();
        end
    end
    methods (TestClassTeardown)
        function goBackToStart(testCase)
            cd(testCase.oldfolder);
        end
    end

    methods (Test,ParameterCombination="sequential")

        function testInitialization(testCase,tTable)
            myObject = validationBatch(tTable);

            testCase.fatalAssertClass(myObject,"validationBatch")
            testCase.assertClass(myObject.numericBatch,"table")
        end

        function testSetFaultyEntries(testCase,tTable)
            myObject = validationBatch(tTable);
            try 
                myObject.setFaultyEntriesParams('LoLz');
            catch
                testCase.verifyTrue('false','Should not fail due to invalid input')
            end
            testCase.verifyTrue(strcmp(myObject.faultMethod,'remove'));

            myObject = myObject.setFaultyEntriesParams('zoh');
            testCase.verifyTrue(strcmp(myObject.faultMethod,'zoh'));
        end

        function testGetFaultFreeData(testCase,tTable)
            myObject = validationBatch(tTable);
            myResult = myObject.getFaultFreeData();
            testCase.verifyEqual(myResult,myObject.numericBatch(~myObject.faultSamples,:));

            myObject = myObject.setFaultyEntriesParams('zoh');
            myResult = myObject.getFaultFreeData();
            testCase.verifySize(myResult,size(myObject.numericBatch));
            disp(myResult);
            prompt = 'Is the table above valid zoh? y/n \n';
            x = input(prompt,"s");
            testCase.verifyTrue(strcmp(x,'y'));
        end

        function testSetStuckValueParams(testCase,tTable,Answers)
            %Also getStuckValues gets tested here.
            myObject = validationBatch(tTable);
            myResult = myObject.getStuckValues();
            testCase.verifyEqual(myResult.Variables,Answers.stuckData_1);
            testCase.verifyEqual(myResult.Properties.VariableNames,...
                myObject.numericBatch.Properties.VariableNames)
            myObject = myObject.setStuckValueParams(1,0);
            myResult = myObject.getStuckValues();
            testCase.verifyEqual(myResult.Variables,Answers.stuckData_2);
            testCase.verifyEqual(myResult.Properties.VariableNames,...
                myObject.numericBatch.Properties.VariableNames)
            myObject = myObject.setFaultyEntriesParams('zoh');
            myObject = myObject.setStuckValueParams(1,0);
            myResult = myObject.getStuckValues();
            testCase.verifyEqual(myResult.Variables,Answers.stuckData_3);
            testCase.verifyEqual(myResult.Properties.VariableNames,...
                myObject.numericBatch.Properties.VariableNames)
        end

        function testGetOutlierCandidates(testCase,tTable,Answers)
            myObject = validationBatch(tTable);
            myObject = myObject.setFaultyEntriesParams('zoh');
            myResult = myObject.getOutlierCandidates(1,'doPlot',true);
            testCase.verifyEqual(myResult,Answers.ellipsoid);
        end

    end %Test methods


    methods (Static)

    end %Static methods

    methods

    end %non-static methods


end %Classdef

%Subfunctions:
function T = getMyTestData(n)
switch n
    case 1
        %Create some table to begin with:
        times = reshape(datetime(date())+(seconds(0):seconds(5):seconds(20)),[],1);
        gdata = getRepeatableRandomNumber(size(times,1),size(times,2));
        bdata = gdata;bdata(1)=NaN;bdata(end)=inf;bdata(floor(length(bdata)/2))="lol";
        T = table(times,gdata,bdata);
        %Modify to get some repetitions going:
        T(2,2) = T(1,2);
        T(3,3) = T(4,3);
    case 2
        %Some 3D data for plotting with stuck values:
        x = getRepeatableRandomNumber(10,3);
        T = table(x(:,1),x(:,2),x(:,3));
end
end
function T = getMyTestAnswers(n)
switch n
    case 1
        T.stuckData_1 = zeros(5,2);
        T.stuckData_2 = zeros(5,2);
        T.stuckData_2(1:2,1) = [1;1];T.stuckData_2(3:4,2) = [1;1];
        T.stuckData_3 = T.stuckData_2;

        T.ellipsoid   = [0 ; 0 ; 0 ; 1 ; 0];
    case 2
        T.stuckData_1 = zeros(10,3);
        T.stuckData_2 = zeros(10,3);
        T.stuckData_3 = zeros(10,3);

        T.ellipsoid   = [0 ; ones(4,1) ; 0 ; 1 ; 1 ; 0 ; 1];
end
end
function randNums = getRepeatableRandomNumber(n,m)
s = rng;
rng(1);
randNums = randn(n,m);
rng(s)
end