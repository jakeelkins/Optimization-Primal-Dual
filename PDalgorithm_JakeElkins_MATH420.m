clc
clear
format short
%%
%make sure filenames are correct

%A = importdata('a.txt', ' ');
%b = importdata('b.txt', ' ');
%c = importdata('c.txt', ' ');
%lambda = importdata('y.txt', ' ');

%b = b'; %I transpose here to fit the way I wrote the algorithm
%c=c';
%lambda = lambda';
%%
%This is where I was inputting my variables for testing. 

%A = [1,1,2;2,1,3];
%b = [3;5];
%c = [2;1;4];
%lambda = [0;0];

%A = [2,1,2;3,3,1];
%b = [4; 3]; 
%c = [4; 1; 1];
%lambda = [0;0];

A=[1, -2,  1,  3,  2, -3,  4,  1,  0,  1,  2,  4,  4, -1,  0; 0,  2,  4, -1,  0,  1,  2, -1,  1,  2,  0,  2,  1, -2,  4; 2,  1,  0,  1,  2, -1, -2,  2,  1,  0, -1, -2,  1,  0,  1; -2, -4,  2,  1, -2,  1,  0,  1,  2, -1,  0,  0,  2,  1, -1; 1,  0,  2,  4,  1,  0,  1, -1,  1,  2,  0,  1, -1,  2, -2];
b=[37;  30;  8;  -4;  8];
c= [1;  -10;   15;  20;  0;   0;  12;  3;   10;   8;  4;   12;  3;  10;  -10];              
lambda=[1; 0; -1; 2; 4];
%%
%------------------algorithm start---------------------

P = []; %initializing P
Pc = []; %initialing P-complement
[lenA,widthA] = size(A);

%making tableau for ARP
tableau = [A eye(lenA) b];

%ARP cost coefficients
cARP = [zeros([1,widthA]) ones([1,lenA]) 0];

%ARP relative cost coefficients
rARP = [-sum(A) zeros([1, lenA]) -sum(b)];

%keep track of basic variable index
bvi = widthA+1:(widthA+lenA);
    
%c-lambda*a used in algorithm, 
cla = c' - (lambda.'*A);

counter = 0; %iteration counter

while true %forever loop will break at stop conditions
    counter = counter+1;
    fprintf('iteration #%d:\n',counter);
    P = [];
    Pc=[];
    
    %this loop checks c-lambda*a to update P and P-complement.
    %if cla is 0, the column index goes in P. if not, it goes in Pc.
    for i=1:(length(cla))
        if cla(i) <= 0
            P(1, length(P)+1) = i;
        else
            Pc(1,length(Pc)+1) = i;
        end
    end
    
    %this loop checks to see if columns need updating (or if P isn't
    %empty), if so, does Gaussian elimination on the tableau to solve the
    %associated restricted primal problem. if not, we have optimal solution
    %for ARP, and we move on.
    if ~isempty(P)
        for i=1:length(P)
            %SOLVING THE TABLEAU (ARP)
            %pivot columns are only indices of P
            pivcol = tableau(:, P(i));
            %getting min rowratio, index
            rowratios = tableau(:, end)./pivcol; %calc row ratios
            ratioindex = find(rowratios>(5*eps)); %find all positive row ratios
            [minratio, minindex] = min(rowratios(ratioindex)); %find minimum ratio -> get pivot row index
            rowindex = ratioindex(minindex);
            %now do gaussian elimination with just this column
            tableau(rowindex,:) = tableau(rowindex,:)./tableau(rowindex,P(i)); %normalizing
            %this loop does the elementary row operations. leaves the row
            %alone if it's the pivot row or if the element is already zero
            for j=1:lenA
                if (j ~= rowindex) && ((abs(tableau(j, P(i))))>(5*eps))
                    tableau(j, :) = (tableau(rowindex, :).*(-tableau(j, P(i))))+tableau(j, :);
                end
            end
            %then update rARP (relative cost coefficients of associated restricted
            %primal problem
            if abs(rARP(P(i))) > 5*eps
                rARP = (tableau(rowindex, :)).*(-rARP(P(i)))+rARP;
            end
        end
    end
    
    
    %update basic variable index by comparing each column to column of an identity matrix
    tableaucompare = tableau(:,1:end-1);
    [lenT, widthT] = size(tableaucompare);
    eyetest = eye(lenT);
    bvi=[];
    
    for i = 1:lenT
        for j = 1:widthT
            if abs((tableaucompare(:,j))-(eyetest(:,i))) < (5*eps)
                bvi(1, length(bvi)+1) = j;
                break
            end
        end
    end
    
    bvi = sort(bvi);
    
    %filling in solution of (x,y) by reading off the tableau 
    xy = zeros([1,(widthA+lenA)]); %initialize where to update
    bvivals = (tableau(:, bvi)*tableau(:, end))'; %multiply our permutation matrix by b to get basic values
    bvcount = 1;
    %fill in basic variable values in correct places in (x,y)
    for i = 1:length(xy)
        if ismember(i,bvi)
            xy(i) = bvivals(bvcount);
            bvcount = bvcount+1;
        end
    end
    
    %separating x from y
    x = xy(1:widthA); y = xy(widthA+1:end);
    
    %from here, (x,y) optimal to current ARP
    disp('current x: '); disp(num2str(x));

    if all(y < 5*eps) %first stop condition, optimality condition
        %BREAK THE LOOP, we have opt soln, else continue
        disp('-----------------------------------------------------------');disp(' ');
        fprintf('OPT SOLUTION FOUND'); disp(' '); disp(' ');
        disp('optimal x: '); disp(num2str(x));
        z = c'*x'; %cost
        disp('optimal z: '); disp(z);
        break
    end
    
    if all((rARP(1:widthA))>=0) %second stop condition, infeasibility
        %BREAK THE LOOP, print out (D)=unbounded and (P)=infeasible.
        fprintf('THIS PROBLEM IS INFEASIBLE. The dual problem is unbounded. Check inputs and try again.')
        break
    end
    
    %calculating epsilon
    cla_Pc = cla(Pc); %cla contained in P-complement
    rARP_Pc = rARP(Pc); %ARP relative costs contained in P-complement

    negRindex = (find(rARP_Pc<0)); %get all rARP that are negative

    eps_ratios = (cla_Pc(negRindex))./(-rARP_Pc(negRindex)); %calculate ratios

    epsilon = min(eps_ratios); %get minimum ratio for epsilon to use for updating lambda, new cla

    cla = cla + (epsilon*(rARP(1:widthA))); %updating new cla
    
    disp('current epsilon: '); disp(epsilon);
    
    %repeat, where next loop will use this loop's epsilon and tableau.
end   