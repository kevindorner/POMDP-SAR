close all; clear all;
file2='ExtendedCaseStudy.policy';
load('ExtendedCaseStudy.mat');
struct=parseXML(file2);
set(0,'DefaultAxesColor','none')
%Parse policy
i=2;j=1;
zerovector=zeros(1,Nr_victims+NrRobots*SizeRobotT);
alphas=zeros(str2double(struct.Children(2).Attributes(2).Value),str2double(struct.Children(2).Attributes(3).Value));
while j<str2double(struct.Children(2).Attributes(2).Value)+1
    alphas(j,:)=str2num(struct.Children(2).Children(i).Children.Data);
    actionvector(j)=str2num(struct.Children(2).Children(i).Attributes(1).Value);
    i=i+2;
    j=j+1;
end
%Initial belief
belief=zeros(str2double(struct.Children(2).Attributes(3).Value),1);
startingstate=0;
for power=1:Nr_victims
    startingstate=startingstate+2^(power-1);
end
belief(startingstate+1,1)=1;

for i=1:63 %Number of timesteps
    beliefcumsum(i,:)=cumsum(belief);
    randcheck=rand>beliefcumsum(i,:);
    statesim(i)=sum(randcheck);
    statesimvec(i,:)=dec2bin(statesim(i)+1,zerovector);
    locsim(i)=bin2dec(statesimvec(i,1:SizeRobotT))-1;
    locdecode(i,2)=locate(mapsize,locsim(i),1); locdecode(i,1)=locate(mapsize,locsim(i),2);
    [action,reward(i)]=bestaction(belief,alphas,actionvector);
    actionbackup(i)=action;
    %Belief updates
    switch action
        case 0
            newbelief=MovementTransitionP0*belief;
        case 1
            newbelief=MovementTransitionP1*belief;
        case 2
            newbelief=MovementTransitionP2*belief;
        case 3
            newbelief=MovementTransitionP3*belief;
        case 4
            newbelief=MovementTransitionP4*belief;
        otherwise
            newbelief=belief;
    end
    belief=newbelief;
end

[B,I] = sort(sum(locdecode==[1 1],2),'desc');
XAxisMax=length(reward);

for x=1:Nr_victims-1
    sp(x)=subplot(2,Nr_victims-1,x); 
    pbaspect(sp(x),[1 1 1])
    set(gca,'YDir','reverse'); xlim([0.5 mapsize+0.5]); ylim([0.5 mapsize+0.5]);
    set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]); hold on;
    %Plot waypoints, starting state, victims     
    for i=1:mapsize; for j=1:mapsize; scatter(i,j,1,'.','MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1); end; end
    scatter(1,1,35,'x','MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1)
    for i=1:Nr_victims
        plot(locate(mapsize,VictimsLocation(i)-1,2),locate(mapsize,VictimsLocation(i)-1,1),'-s','MarkerSize',3,'MarkerEdgeColor','b','MarkerFaceColor','b')
    end
    %Plot inaccessible squares
    rectangle('Position',[0.5,0.5,mapsize,mapsize],'EdgeColor','k')
    matrix= [0 1 1 5; 1 5 2 1; 3 0 1 2; 4 1 1 1; 5 6 4 1; 6 0 1 4; 7 0 1 2; 7 5 1 2; 9 2 1 1; 10 1 1 4; 11 4 1 1; 10 6 2 1; 0 7 2 1; 1 9 3 1; 3 7 3 1; 3 9 1 2; 5 8 1 1; 5 10 1 2; 7 8 4 1; 7 9 1 3; 8 10 1 2; 10 9 1 1; 10 11 1 1]; 
    for i=1:length(matrix)
        rectangle('Position',[matrix(i,1)+0.5,matrix(i,2)+0.5,matrix(i,3),matrix(i,4)],'FaceColor','k','EdgeColor','k')
    end
    sp(x+Nr_victims-1)=subplot(2,Nr_victims-1,x+Nr_victims-1);
    pbaspect(sp(x+Nr_victims-1),[1 0.3 1])
    set(gca,'XTickLabel',[]);
    hold on;
    %Plot arrowheads, labels
    xticks([XAxisMax]);xticklabels({'t'});
    axis([0 XAxisMax 0 round(max(reward),2)*1.2])
    text(XAxisMax-1.9,-0.1,'$\rightarrow$','Interpreter','latex','FontSize',8)
    text(-0.93,round(max(reward),2)*1.2,'$\uparrow$','Interpreter','latex','FontSize',8)

end
x=1; victimsaved=[];
for i=1:length(locdecode)-1
    if x<Nr_victims 
        subplot(2,Nr_victims-1,x); 
        y=1;
        while y<x
            plot(victimsaved(y,1),victimsaved(y,2),'-s','MarkerSize',4,'MarkerEdgeColor','w','MarkerFaceColor','w');
            scatter(victimsaved(y,1),victimsaved(y,2),1,'.','MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1);
            y=y+1;
        end
        if locdecode(i+1,1)==locdecode(i+1,2)
            if locdecode(i+1,1)==1
                %Plot rover, exp.reward before rescue, current timestep
                plot(locdecode(i,1),locdecode(i,2),'-p','MarkerSize',5,'MarkerEdgeColor',[0 0.6 0.3],'MarkerFaceColor',[0 0.6 0.3])
                victimsaved=[victimsaved; locdecode(i,1) locdecode(i,2)];
                subplot(2,Nr_victims-1,x+Nr_victims-1);
                line([i i],[0 round(max(reward),2)*1.2],'Color','r','LineStyle','-')
                text(i+1,round(max(reward),2)-2,'r=100','Interpreter','latex','FontSize',8)
                x=x+1;
            else
                %Plot rescue path
                subplot(2,Nr_victims-1,x);
                plot([locdecode(i,1) locdecode(i+1,1)],[locdecode(i,2) locdecode(i+1,2)],'Color',[0.4660 0.6740 0.1880],'LineWidth',1)
            end
        else
            %Plot rescue path
            subplot(2,Nr_victims-1,x);
            plot([locdecode(i,1) locdecode(i+1,1)],[locdecode(i,2) locdecode(i+1,2)],'Color',[0.4660 0.6740 0.1880],'LineWidth',1)
        end
    end
end

%Plot expected future reward
for x=Nr_victims:2*Nr_victims-2
    subplot(2,Nr_victims-1,x);
    axis([0 XAxisMax 0 round(max(reward),2)*1.2]);
    for i=1:length(locdecode)-1
        plot([i-1 i],[reward(i) reward(i)],'b','LineWidth',1)
        plot([i i],[reward(i) reward(i+1)],'b','LineWidth',1)
    end
    text(XAxisMax-1.9,-0.1,'$\rightarrow$','Interpreter','latex','FontSize',8)
    text(-1,round(max(reward),2)*1.2,'$\uparrow$','Interpreter','latex','FontSize',8)
    text(0,round(max(reward),2)*1.2+5,'$V_{\pi}(b_t)$','Interpreter','latex','FontSize',8)
%     text(XAxisMax,0,'$t$','Interpreter','latex','FontSize',8)
    text(0,round(max(reward),2),num2str(round(max(reward),2)),'Interpreter','latex','FontSize',8)
end

%% Functions
%Select best action to take
function [action, reward]=bestaction(belief,alphavectors,actionvector)
    Ereward=alphavectors*(belief);
    [m, i]=max(Ereward);
    action=actionvector(i);
    reward=m;
end

%XML Parser
function theStruct = parseXML(filename)
    % PARSEXML Convert XML file to a MATLAB structure.
    try
       tree = xmlread(filename);
    catch
       error('Failed to read XML file %s.',filename);
    end
    % Recurse over child nodes. This could run into problems 
    % with very deeply nested trees.
    try
       theStruct = parseChildNodes(tree);
    catch
       error('Unable to parse XML file %s.',filename);
    end
    % ----- Local function PARSECHILDNODES -----
    function children = parseChildNodes(theNode)
    % Recurse over node children.
    children = [];
    if theNode.hasChildNodes
       childNodes = theNode.getChildNodes;
       numChildNodes = childNodes.getLength;
       allocCell = cell(1, numChildNodes);
       children = struct(             ...
          'Name', allocCell, 'Attributes', allocCell,    ...
          'Data', allocCell, 'Children', allocCell);
        for count = 1:numChildNodes
            theChild = childNodes.item(count-1);
            children(count) = makeStructFromNode(theChild);
        end
    end
    end
    % ----- Local function MAKESTRUCTFROMNODE -----
    function nodeStruct = makeStructFromNode(theNode)
    % Create structure of node info.
    nodeStruct = struct(                        ...
       'Name', char(theNode.getNodeName),       ...
       'Attributes', parseAttributes(theNode),  ...
       'Data', '',                              ...
       'Children', parseChildNodes(theNode));
    if any(strcmp(methods(theNode), 'getData'))
       nodeStruct.Data = char(theNode.getData); 
    else
       nodeStruct.Data = '';
    end
    end
    % ----- Local function PARSEATTRIBUTES -----
    function attributes = parseAttributes(theNode)
    % Create attributes structure.
    attributes = [];
    if theNode.hasAttributes
       theAttributes = theNode.getAttributes;
       numAttributes = theAttributes.getLength;
       allocCell = cell(1, numAttributes);
       attributes = struct('Name', allocCell, 'Value', ...
                           allocCell);
       for count = 1:numAttributes
          attrib = theAttributes.item(count-1);
          attributes(count).Name = char(attrib.getName);
          attributes(count).Value = char(attrib.getValue);
       end
    end
    end
end

%Decimal to binary converter
function binVec=dec2bin(decNumber,zerovector)
    binVec=zerovector;
    decNumber=decNumber-1;
    if decNumber<=0
        return
    end
    for i=length(zerovector):-1:1
        decNumber2 = floor(decNumber/2);
        binVec(i) = decNumber-2*decNumber2;
        decNumber=decNumber2;
    end   
end

%Binary to decimal convertor
function decNumber=bin2dec(binVec)
    decNumber=0; j=1;
    for i=length(binVec):-1:1
        decNumber=decNumber+binVec(j).*2^(i-1);
        j=j+1;
    end
    decNumber=decNumber+1;
end

%Translate location state into coordinate
function location=locate(mapsize,locationencoding,coord)
    if mod(locationencoding,mapsize)==0
        xcoord=floor(locationencoding/mapsize)+1;
        ycoord=1;
    else
        xcoord=floor(locationencoding/mapsize)+1;
        ycoord=mod(locationencoding,mapsize)+1;
    end
    if coord==1
        location=xcoord;
    else
        location= ycoord;
    end
end