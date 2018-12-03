clc
clear
%what do we do about a b(i) < 0? or degeneracy?
A = [1, 1, 2; 2, 1, 3];
b = [3; 5]; 
c = [2; 1; 4];

P = [];
Pc = [];
[lenA,widthA] = size(A);

%making tableau for ARP
tableau = [A eye(lenA) b];

%ARP cost coefficients
cARP = [zeros([1,widthA]) ones([1,lenA]) 0];

%ARP relative cost coefficients
rARP = [-sum(A) zeros([1, lenA]) -sum(b)];

%keep track of basic variable index
bvi = widthA+1:(widthA+lenA);


if all(b) >= 0
    lambda = [0; 0]; %for initial feasibility to (dLP)
end
    

cla = c' - (lambda.'*A);

counter = 0;
while true
    counter = counter+1
    P = [];
    
    for i=1:(length(cla))
        if cla(i) <= 0
            P(1, length(P)+1) = i;
        else
            Pc(1,length(Pc)+1) = i;
        end
    end
    
    
    if ~isempty(P)
        for i=1:length(P)
            %SOLVING THE TABLEAU (ARP)
            %pivot columns only indices of P
            pivcol = tableau(:, P(i));
            %getting min rowratio, index
            rowratios = tableau(:, end)./pivcol;
            [minratio, rowindex] = min(rowratios);
            %update bvi
            bvi(1, i) = P(i)
            %now do gaussian elimination with just this column
            tableau(rowindex,:) = tableau(rowindex,:)./tableau(rowindex,P(i)); %normalize
            for j=1:lenA
                if j ~= rowindex
                    tableau(j, :) = (tableau(rowindex, :).*(-tableau(j, P(i))))+tableau(j, :);
                end
            end
            
            %then update rARP
            rARP = (tableau(rowindex, :)).*(-rARP(P(i)))+rARP;
        end
    end
    %tableau
    %rARP
    %filling in solution
    xy = zeros([1,(widthA+lenA)]);
    bvivals = [tableau(:, bvi)*tableau(:, end)]'
    bvcount = 1;
    for i = 1:length(xy)
        if ismember(i,bvi)
            bvcount
            xy(i) = bvivals(bvcount);
            bvcount = bvcount+1;
        end
    end
    xy
    x = [xy(1:widthA)]; y = [xy(widthA+1:end)];
    %from here, (x,y) optimal to ARP

    if all(y < 5*eps) %how big does this need to be?
        %BREAK THE LOOP, we have opt soln, else continue
        fprintf('OPT SOLUTION FOUND:');
        disp(x)
        break
    end
    
    if all((rARP(1:widthA))>=0) %how big does this need to be?
        %BREAK THE LOOP, print out (D)=unbounded and (P)=infeasible.
        fprintf('THIS PROBLEM IS INFEASIBLE')
        break
    end

    %add in update lambda later?

    ratios = cla./(-rARP(1:widthA)); %calc ratios

    epsilon = min(ratios); %get epsilon to use for updating lambda, new cla

    cla = cla + (epsilon*(rARP(1:widthA)));
    

end   