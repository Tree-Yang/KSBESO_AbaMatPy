%Main Script for soft-kill bidirection evolution method
%*******************************************************
%By J S Yang
%Date: 2018-12-04
%*******************************************************

clear; close all; clc;

%Parameters
EvoluationRate = 0.02;      %Evolution rate
VolFracCon = 0.4;           %Volume constraints
IterNumMax = 1000;          %Maximum iteration number
ConvergeCondition = 1e-5;   %Converge Condition
IterNum = 0;                %Initial iternation number
Rmin = 4;                   %Filter Radium

%+++++++++++++++++++++++++++++++++++++++++++++++++++++
%Finish in python file InpGenerater.py
% DV = ones(6400,1);          %Initial design variables
% save('./DesignVariables/DV_Iter0.dat','DV','-ascii')
%+++++++++++++++++++++++++++++++++++++++++++++++++++++

%Create folder for Design variables
if exist('DesignVariables','dir')==0
    mkdir('DesignVariables');
end 

%Create file to transfer iteration number to python script
save('IterNum.dat', 'IterNum', '-ascii');
ObjChange10 = 1;

while ObjChange10 > ConvergeCondition
    
    %Change of objective values in recent 10 loops
    if IterNum >= 10
        Obj = load('ModelTotExtWork.dat');
        ObjChange10 = abs(sum(Obj(end-9 : end-5)) - sum(Obj(end-4 : end))) / abs(sum(Obj(end-9 : end-5)));
    end

    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % This is for high version Matlab with interface to python
    % %Call InpGenerater0.py to generate new *.inp file
    % import py.InpGenerater0.*
    % InpGenerater0(IterNum)
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % This is for low version Matlab without interface to python
    %Call InpGenerater.py to generate new *.inp file
    dos('python InpGenerater.py', '-echo');
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    %Run Abaqus in command window from *.inp file
    JobName = ['Job-Iter', num2str(IterNum)];
    JobRun = ['abaqus job=',JobName,' int'];
    % JobRun = 'abaqus job=Job-Iter1';
    % Attention: NO SPACE before and after '=' to avoid Abaqus Error
    % 'int' displays the Abaqus process: Run pre.exe ...
    dos(JobRun, '-echo');

    %Get element elastic strain energy from *.odb file
    ScriptName = 'OdbReader.py';
    ScriptRun = ['abaqus cae noGUI=',ScriptName];
    dos(ScriptRun, '-echo');

    if IterNum == 0
        %Get all nodes and all elements
        ElementF  = load('ElementsFull.dat');
        NodesF = load('NodesFull.dat');
        CoordNum = length((NodesF(1,:))) - 1;   %Number of coordinates per node
        NodePerEle = length(ElementF(1,:)) - 1;   %Number of nodes per element
        EleCentersF = zeros(length(ElementF(:,1)), CoordNum + 1);

        if CoordNum == 2                        %2d model
            NodesXF = NodesF(:, 2);
            EleNodesXF = NodesXF(ElementF(:, 2:end));
            NodesYF = NodesF(:, 3);
            EleNodesYF = NodesYF(ElementF(:, 2:end));
            EleCentersF(:,1) = ElementF(:,1);
            EleCentersF(:,2) = sum(EleNodesXF,2) / NodePerEle;
            EleCentersF(:,3) = sum(EleNodesYF,2) / NodePerEle;
            clear NodesXF NodesYF EleNodesXF EleNodesYF;
        elseif CoordNum == 3                    %3d model
            NodesXF = NodesF(:, 2);
            EleNodesXF = NodesXF(ElementF(:, 2:end));
            NodesYF = NodesF(:, 3);
            EleNodesYF = NodesYF(ElementF(:, 2:end));
            NodesZF = NodesF(:, 4);
            EleNodesZF = NodesZF(ElementF(:, 2:end));
            EleCentersF(:,1) = ElementF(:, 1);
            EleCentersF(:,2) = sum(EleNodesXF,2) / NodePerEle;
            EleCentersF(:,3) = sum(EleNodesYF,2) / NodePerEle;
            EleCentersF(:,4) = sum(EleNodesZF,2) / NodePerEle;
            clear NodesXF NodesYF NodesZF EleNodesXF EleNodesYF EleNodesZ
        end

        %Get distance between centers of elements
        %This function is so fast!!
        %distance between element jj and all the other elements in row jj
        EleCenDispF = pdist2(EleCentersF(:,2:end),EleCentersF(:,2:end));
        %elements in the range of filter and corresponding element ids
        FilterIndex = (EleCenDispF < Rmin) * 1.0;
        DispInFilter = (Rmin - EleCenDispF) .* FilterIndex;
        TotalDisp = sum(DispInFilter,2);
        TotalDispM = repmat(TotalDisp,1,size(DispInFilter,2));
        %filter matrix
        FilterMat = DispInFilter ./ TotalDispM;
    end

    %Target volume fraction in current loop
    TargetVolFrac = (1-EvoluationRate) ^ (IterNum+1);
    TargetVolFrac = max(TargetVolFrac, VolFracCon);

    %Get Element elastic energy density
    EleElasEnerDen = load('EleElasEnerDen.dat');

    %Element sensitivity
    if IterNum == 0
        DV = ones(size(EleElasEnerDen));
        save('./DesignVariables/DV_Iter0.dat','DV','-ascii');
    end
    EleSensitivity = EleElasEnerDen ./ DV;

    %Initialize filtered element sensitivity
    if IterNum == 0
        EleSenFiltered = zeros(size(EleSensitivity));
    end

    %Filterate the element sensitivity
    for jj = 1:1:size(DispInFilter,1)
        EleSenFiltered(jj) = FilterMat(jj,:) * (EleSensitivity .* (FilterIndex(jj,:))');
    end

    %Average the filtered element sensitivity bewteen the last and current loop for algorithm staibility.
    if IterNum >= 1
        EleSenFiltered = (EleSenFiltered + EleSenFilteredO) / 2.0;
    end

    %Bi-section methods for critical element sensitivity
    EleSenH = max(EleSenFiltered);
    EleSenL = min(EleSenFiltered);
    while abs((EleSenH - EleSenL) / EleSenH) > 10e-5
        EleSenCr = (EleSenH + EleSenL) / 2;
        DV = 0.999 * (EleSenFiltered > EleSenCr) + 0.001;
        if sum(DV) > TargetVolFrac * length(DV)
            EleSenL = EleSenCr;
        else
            EleSenH = EleSenCr;
        end
    end

    % Save design variables to files
    DVName = ['./DesignVariables/DV_Iter',num2str(IterNum + 1),'.dat'];
    save(DVName,'DV','-ascii');

    %Save the last filterd element sensitivity
    EleSenFilteredO = EleSenFiltered;

    %Update the iteration counter
    IterNum = IterNum + 1;
    save('IterNum.dat', 'IterNum', '-ascii');
    
    % Max Iteration number
    if IterNum > IterNumMax
        break;
    end
end

%get topology
dos('abaqus cae noGUI=PostProcessor.py','-echo')

%History curve plot
[DVHis,ObjValHis] = HistoryPlot;