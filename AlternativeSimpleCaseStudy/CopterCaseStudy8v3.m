clear all; close all;
% profile on
%% Parameters
P_survival=[1 1 1]; %Survival probabilities
Discount=0.9; Reward=100;
VictimsLocation=[5 34 49]; %Location of each victim (starting from index 1)
NrRobots=1; Nr_victims=3; mapsize=8;
startingstate=0;
for power=1:Nr_victims
    startingstate=startingstate+2^(power-1);%Initiate starting state to index 0 and all victims alive
end
%Array of inacessible waypoints
exclude=[];
exclude=[4 7:8 9 12:13 15:16 17 23 25 31 33 41:43 48 54:56];

%% Construct Transition Matrices
%Determine number of states
SizeRobotT=ceil(log(mapsize^2)/log(2));
zerovector=zeros(1,Nr_victims+NrRobots*SizeRobotT);
NrStatesV=[];
mapsizeV=dec2bin(mapsize^2,zerovector);
mapsizeV(:,1:Nr_victims)=[];
for i=1:NrRobots
    NrStatesV=cat(2,NrStatesV,mapsizeV);
end
for i=1:Nr_victims
    NrStatesV=cat(2,NrStatesV,[1]);
end
NrStates=bin2dec(NrStatesV);

%Allocate sparse matrices
MovementT1=spalloc(mapsize^2,mapsize^2,1);
MovementT2=spalloc(mapsize^2,mapsize^2,1);
MovementT3=spalloc(mapsize^2,mapsize^2,1);
MovementT4=spalloc(mapsize^2,mapsize^2,1);
MovementTransitionP0=spalloc(NrStates,NrStates,1);
MovementTransitionP1=spalloc(NrStates,NrStates,1);
MovementTransitionP2=spalloc(NrStates,NrStates,1);
MovementTransitionP3=spalloc(NrStates,NrStates,1);
MovementTransitionP4=spalloc(NrStates,NrStates,1);

%Build finding transition matrices
x=1; action=0;
for initialState=1:NrStates %loop over all intial states
    initialStateV=dec2bin(initialState,zerovector);R1initS=bin2dec(initialStateV(1:SizeRobotT));                   
    for endState=1:NrStates %loop over all end states
        endStateV=dec2bin(endState,zerovector);
        SearchTransitionPx(x)=initialState; SearchTransitionPy(x)=endState; %register current states+1           
        %robot 1 movement
        R1endS=bin2dec(endStateV(1:SizeRobotT));
        SearchTransitionPv(x)=RobotsT(R1initS,action,R1endS,mapsize); %register transition probability
        %victim transitions
        for i=1:Nr_victims
            if action==0
                    if (R1initS== VictimsLocation(i)||...
                        (locate(R1initS-1,1,mapsize)==locate(VictimsLocation(i)-1,1,mapsize)-1&&locate(R1initS-1,2,mapsize)==locate(VictimsLocation(i)-1,2,mapsize))||...
                        (locate(R1initS-1,1,mapsize)==locate(VictimsLocation(i)-1,1,mapsize)+1&&locate(R1initS-1,2,mapsize)==locate(VictimsLocation(i)-1,2,mapsize))||...
                        (locate(R1initS-1,2,mapsize)==locate(VictimsLocation(i)-1,2,mapsize)-1&&locate(R1initS-1,1,mapsize)==locate(VictimsLocation(i)-1,1,mapsize))||...
                        (locate(R1initS-1,2,mapsize)==locate(VictimsLocation(i)-1,2,mapsize)+1&&locate(R1initS-1,1,mapsize)==locate(VictimsLocation(i)-1,1,mapsize)))
                        rescuedvictim=1;
                    else
                        rescuedvictim=0;
                    end
            else
                rescuedvictim=0;
            end
            SearchTransitionPv(x)=SearchTransitionPv(x)*VictimsT(i,initialStateV(NrRobots*SizeRobotT+i),action,endStateV(NrRobots*SizeRobotT+i),rescuedvictim);
        end
        x=x+1;
    end
end     
MovementTransitionP0=sparse(SearchTransitionPy,SearchTransitionPx,SearchTransitionPv);

%Create transition matrices
A=[1 1 0; 0 0 1; 0 0 0]; B=[0 0 0; 1 0 0; 0 1 1];
Av=zeros(mapsize-1,1);Av(1,1)=1;
A=cat(2,Av,eye(mapsize-1));A=cat(1,A,zeros(1,length(A)));
Bv=zeros(mapsize-1,1);Bv(end,1)=1;
B=cat(2,eye(mapsize-1),Bv);B=cat(1,zeros(1,length(B)),B);

MovementT1=kron(sparse(A),eye(mapsize));
MovementT3=kron(sparse(B),eye(mapsize));
MovementT2=kron(eye(mapsize),sparse(B));
MovementT4=kron(eye(mapsize),sparse(A));
VictimTr=[1]; VictimT=[1];
VictimTr0=sparse([1 1; 0 0]);
VictimT0=sparse([1 1-P_survival(1); 0 P_survival(1)]);
for i=1:Nr_victims
    VictimTr=kron(VictimTr,VictimTr0);
    VictimT=kron(VictimT,VictimT0);
end

%Build movement transition matrices
MovementTransitionP1=kron(MovementT1,VictimT);
MovementTransitionP2=kron(MovementT2,VictimT);
MovementTransitionP3=kron(MovementT3,VictimT);
MovementTransitionP4=kron(MovementT4,VictimT);    

%% Transition probabilities
% Transfer matrix [T: <action> : <start-state> : <end-state> %f]
nrStates=NrStates;
file = fopen('CopterCaseStudy.POMDP','w');
fprintf(file,'discount: %f \n',Discount);
fprintf(file,'values: reward \n');
fprintf(file,'actions: rescue N E S W \n');
fprintf(file,'observations: not_found found rescued\n');
fprintf(file,'states: %d \n\n',nrStates);
fprintf(file,'start: \n');
for i=1:NrStates
     if i==startingstate+1 
        fprintf(file,'1.0 ');
    else
        fprintf(file,'0.0 ');
    end
end
fprintf(file,'\n\n');        

[r0, c0, v0] = find(MovementTransitionP0);
[r1, c1, v1] = find(MovementTransitionP1);
[r2, c2, v2] = find(MovementTransitionP2);
[r3, c3, v3] = find(MovementTransitionP3);
[r4, c4, v4] = find(MovementTransitionP4);

for i=1:length(v0)
    if v0(i)==1
        fprintf(file,'T: rescue : %d : %d   1.0 \n',c0(i)-1,r0(i)-1);    
    end
end
for i=1:length(v1)
    if v1(i)==1
        fprintf(file,'T: N : %d : %d   1.0 \n',c1(i)-1,r1(i)-1);    
    else
        fprintf(file,'T: N : %d : %d   %s \n',c1(i)-1,r1(i)-1,num2str(round(v1(i),4)));
    end
end
for i=1:length(v2)
    if v2(i)==1
        fprintf(file,'T: E : %d : %d   1.0 \n',c2(i)-1,r2(i)-1);    
    else
        fprintf(file,'T: E : %d : %d   %s \n',c2(i)-1,r2(i)-1,num2str(round(v2(i),4)));
    end
end
for i=1:length(v3)
    if v3(i)==1
        fprintf(file,'T: S : %d : %d   1.0 \n',c3(i)-1,r3(i)-1);    
    else
        fprintf(file,'T: S : %d : %d   %s \n',c3(i)-1,r3(i)-1,num2str(round(v3(i),4)));
    end
end
for i=1:length(v4)
    if v4(i)==1
        fprintf(file,'T: W : %d : %d   1.0 \n',c4(i)-1,r4(i)-1);    
    else
        fprintf(file,'T: W : %d : %d   %s \n',c4(i)-1,r4(i)-1,num2str(round(v4(i),4)));
    end
end

%% Observations and Rewards
% O : <action> : <end-state> : <observation> %f
% R: <action> : <start-state> : <end-state> : <observation> %f
fprintf(file,'O: rescue : * : rescued   1.0 \n');
ObservationP=zeros(5,NrStates);
for a=2:5 %Loop over all movement actions (a=1:r; a=2:N; a=3:E; a=4:S; a=5:W)
    action=a-1; %Correct offset
    for initialState=1:NrStates %Loop over all intial states
            initialStateV=dec2bin(initialState,zerovector);
        for endState=1:NrStates %Loop over all end states
            endStateV=dec2bin(endState,zerovector);AliveVictimFound=0;R1endS=bin2dec(endStateV(1:SizeRobotT));
            for victim=1:Nr_victims  
                if bin2dec(endStateV(SizeRobotT+victim:SizeRobotT+victim))-1==1
                    %Check victim location in 1 waypoint radius
                    if (R1endS== VictimsLocation(victim)||...
                        (locate(R1endS-1,1,mapsize)==locate(VictimsLocation(victim)-1,1,mapsize)-1&&locate(R1endS-1,2,mapsize)==locate(VictimsLocation(victim)-1,2,mapsize))||...
                        (locate(R1endS-1,1,mapsize)==locate(VictimsLocation(victim)-1,1,mapsize)+1&&locate(R1endS-1,2,mapsize)==locate(VictimsLocation(victim)-1,2,mapsize))||...
                        (locate(R1endS-1,2,mapsize)==locate(VictimsLocation(victim)-1,2,mapsize)-1&&locate(R1endS-1,1,mapsize)==locate(VictimsLocation(victim)-1,1,mapsize))||...
                        (locate(R1endS-1,2,mapsize)==locate(VictimsLocation(victim)-1,2,mapsize)+1&&locate(R1endS-1,1,mapsize)==locate(VictimsLocation(victim)-1,1,mapsize)))
                        AliveVictimFound=1; %Victim found
                    end
                end
            end
            if AliveVictimFound 
                    %Print found observation and allocate reward
                    ObservationP(action,endState)=1;                                          
                    fprintf(file,'O: %d : %d : found   1.0 \n',action,endState-1);
                    fprintf(file,'R: rescue : %d : * : *   %d \n',endState-1,Reward);
            end
        end
        if a==5
            R1initS=bin2dec(initialStateV(1:SizeRobotT));
            if ismember(R1initS,exclude)
                fprintf(file,'R: * : * : %d : * -10000 \n',initialState-1); %Punish impossible move
            end
        end
    end
end

%Print not found observations
for action=1:4
    for endState=1:NrStates
        endStateV=dec2bin(endState,zerovector);
        if ObservationP(action,endState)==0
            fprintf(file,'O: %d : %d : not_found   1.0 \n',action,endState-1);
        end
    end
end

fclose(file);
save('CopterCaseStudy.mat');
% profile viewer

%% Functions
%decimal to binary convertor
function binVec=dec2bin(decNumber,zerovector)
    binVec=zerovector;len=length(zerovector);
    decNumber=decNumber-1;
    if decNumber<=0
        return
    end
    for i=len:-1:1
        decNumber2 = floor(decNumber/2);
        binVec(i) = decNumber-2*decNumber2;
        decNumber=decNumber2;
    end   
end

%binary to decimal convertor
function decNumber=bin2dec(binVec)
    decNumber=0; len=length(binVec);
    for i=len:-1:1
        decNumber=decNumber+binVec(len-i+1).*2^(i-1);
    end
    decNumber=decNumber+1;
end

%Translate location state into coordinate
function location=locate(locationencoding,coord,mapsize)
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

%Survival transfer function
function TransferP=VictimsT(victimIndex,initialState,action,endState,recuedvictim,P_survival)
TransferP=0;
    if action==0 %Rescue
        if endState==initialState %No health state change
            TransferP=1;
        else
            TransferP=0;
        end
        if recuedvictim==1
            if (initialState==1 && endState==0) || (initialState==0 && endState==0)
                TransferP=1;
            elseif (initialState==1 && endState==1) || (initialState==0 && endState==1)
                TransferP=0;
            end
        end
    elseif initialState==0 && endState==0
        TransferP=1;
    elseif initialState==1 && endState==0
        TransferP=1-P_survival(victimIndex);
    elseif initialState==1 && endState==1
        TransferP= P_survival(victimIndex);
    end
end

%Movement transfer function 
function TransferP=RobotsT(initialState,action,endState,mapsize)
TransferP=0;
    switch action
        case 0
            if endState==initialState
                TransferP=1;
            end
        case 1
            if locate(initialState-1,1,mapsize)==1 
                if endState==initialState %north boundary
                    TransferP=1;
                end
            elseif endState==initialState-mapsize
                TransferP=1;
            end
        case 2
            if locate(initialState-1,2,mapsize)== mapsize 
                if endState==initialState %east boundary
                    TransferP=1;
                end
            elseif endState==initialState+1
                TransferP=1;
            end
        case 3
            if locate(initialState-1,1,mapsize)== mapsize
                if endState==initialState %south boundary
                    TransferP=1;
                end
            elseif endState==initialState+mapsize
                TransferP=1;
            end
        case action==4
            if locate(initialState-1,2,mapsize)==1 
                if endState==initialState %west boundary
                    TransferP=1;
                end
            elseif endState==initialState-1
                TransferP=1;
            end  
        otherwise
            if initialState==endState %Victim health decay witout movement not possible
                TransferP=0; 
            elseif(initialState>mapsize^2) || (endState>mapsize^2)
                TransferP=0;
            end
    end
end