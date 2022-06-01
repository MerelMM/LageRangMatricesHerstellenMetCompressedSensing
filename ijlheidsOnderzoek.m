%bekijk het effect van verschillende ijlheidsgraden van de metingsmatrix

drempel = 10e-5;
rijen = 7;
kolommen = 6;
rang = 2;
start = (rijen+kolommen-rang)*rang;
aantalMetingen = 2*start;
aantalKeer = 20;
y = zeros(aantalKeer*(aantalMetingen-start+1),1);
k=0;
m=0;
lambda=10;

spaars2 = zeros(aantalKeer*(aantalMetingen-start+1),1);
spaars3 = zeros(aantalKeer*(aantalMetingen-start+1),1);
spaars4 = zeros(aantalKeer*(aantalMetingen-start+1),1);
spaars5 = zeros(aantalKeer*(aantalMetingen-start+1),1);
aantalMetingen
for i= start:aantalMetingen
    i
    for j=1:aantalKeer
        k=k+1;
        [spaars2(k),spaars3(k),spaars4(k),spaars5(k)] = losOp(rijen,kolommen,rang,i,lambda);        
    end
end

%%gemiddelde
g=0;
k=0;
sr2 = zeros((aantalMetingen-start+1),1);
sr3 = zeros((aantalMetingen-start+1),1);
sr4 = zeros((aantalMetingen-start+1),1);
sr5 = zeros((aantalMetingen-start+1),1);

for j=start:aantalMetingen
    g=g+1;
    s2=0;
    s3=0;
    s4=0;
    s5=0;

    for i = 1:aantalKeer
        k=k+1;
        if (spaars2(k)<=drempel)
            s2 = s2 + 1;
        end
        if (spaars3(k)<=drempel)
            s3 = s3 + 1;
        end
        if (spaars4(k)<=drempel)
            s4 = s4 + 1;
        end
        if (spaars5(k)<=drempel)
            s5 = s5 + 1;
        end
    end
    sr2(g) = s2/aantalKeer;
    sr3(g) = s3/aantalKeer;
    sr4(g) = s4/aantalKeer;
    sr5(g) = s5/aantalKeer;
end
x=start/start:1/start:aantalMetingen/start;
fig1 = figure(1)
plot(x,sr2,'--ks')
title('~1/2 ijlheid',fontsize=16)
xlabel('oversamplingsfactor \phi [aantal metingen / (m+n-rang)*rang]',fontsize=16)
ylabel('succesvol hersteld [%]',fontsize=16)
grid on

fig2 = figure(2)
plot(x,sr3,'--ks')
title('~1/10 ijlheid',fontsize=16)
xlabel('oversamplingsfactor \phi [aantal metingen / (m+n-rang)*rang]',fontsize=16)
ylabel('succesvol hersteld [%]',fontsize=16)
grid on

fig3 = figure(3)
plot(x,sr4,'--ks')
title('~1/20 ijlheid',fontsize=16)
xlabel('oversamplingsfactor \phi [aantal metingen / (m+n-rang)*rang]',fontsize=16)
ylabel('succesvol hersteld [%]',fontsize=16)
grid on

fig4 = figure(4)
plot(x,sr5,'--ks')
title('~1/50 ijlheid',fontsize=16)
xlabel('oversamplingsfactor \phi [aantal metingen / (m+n-rang)*rang]',fontsize=16)
ylabel('succesvol hersteld [%]',fontsize=16)
grid on

function [x2,x3,x4,x5]= losOp(rijen,kolommen,rang,aantalMetingen,lambda)

%maak originele matrix en metingen
A = randn(rijen,rang);
B = randn(rang, kolommen);
matrix = A*B; %maak matrix van rang rang
vector = matrix(:);
[C2,meting2] = MaakMetingenMetSpaarseMatrix2(2,aantalMetingen, rijen*kolommen, vector);
[C3,meting3] = MaakMetingenMetSpaarseMatrix2(10,aantalMetingen, rijen*kolommen, vector);
[C4,meting4] = MaakMetingenMetSpaarseMatrix2(20,aantalMetingen, rijen*kolommen, vector);
[C5,meting5] = MaakMetingenMetSpaarseMatrix2(50,aantalMetingen, rijen*kolommen, vector);

%herstelling met meegegeven penalty
options = optimoptions(@fminunc,'Algorithm','trust-region','SpecifyObjectiveGradient',true, 'MaxFunctionEvaluations',3500, 'MaxIterations', 1e6,...
    'OptimalityTolerance', 1e-7, 'StepTolerance', 1e-7, 'FunctionTolerance', 1e-7 );

A = randn(rijen,rang);
B = randn(rang,kolommen);
z0 = [A(:);B(:)];

%herstel 
objfun = @(z)recoverWithGradientEnPenalty(C2,z, kolommen, rijen, rang, meting2, lambda)
[z] = fminunc(objfun,z0,options);
g1 = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
x2 = norm(matrix-g1) /norm(matrix);

objfun = @(z)recoverWithGradientEnPenalty(C3,z, kolommen, rijen, rang, meting3, lambda)
[z] = fminunc(objfun,z0,options);
g1 = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
x3 = norm(matrix-g1) /norm(matrix);

objfun = @(z)recoverWithGradientEnPenalty(C4,z, kolommen, rijen, rang, meting4, lambda)
[z] = fminunc(objfun,z0,options);
g1 = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
x4 = norm(matrix-g1) /norm(matrix);


objfun = @(z)recoverWithGradientEnPenalty(C5,z, kolommen, rijen, rang, meting5, lambda)
[z] = fminunc(objfun,z0,options);
g1 = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
x5 = norm(matrix-g1) /norm(matrix);
        
end
function[f,gradientt] = recoverWithGradientEnPenalty(C,z, kolommen, rijen, rang, meting, lambda)
f = 0.5*norm(C * reshape((reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen)),rijen*kolommen,1)-meting)^2 + ...
    lambda * (norm(reshape(z(1:rijen*rang,:), rijen, rang), "fro")^2 - norm(reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen),"fro")^2)^2;
A = reshape(z(1:rijen*rang,:), rijen, rang);
B  = reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen);
%met A% = eiejT
pos = 1;
gradientt  = zeros(rijen*rang+kolommen*rang,1);
%eerst voor A mxr
for j=1:rang
    for i=1:rijen
        ei = zeros(rijen,1);
        ej = zeros(rang,1);
        ei(i) = 1;
        ej(j)=1;
        E = ei*ej';
          gradientt(pos) = dot(C*reshape((A*B), rijen*kolommen, 1)-meting, C*reshape(E*B, rijen*kolommen, 1)) + ...
            4*lambda*(norm(A,"fro")^2-norm(B,"fro")^2)*trace(A'*E);
        pos = pos+1;
    end
end
%dan voor B rxn
for j=1:kolommen
    for i=1:rang
        ei = zeros(rang,1);
        ej = zeros(kolommen,1);
        ei(i) = 1;
        ej(j)=1;
        E = ei*ej';
        gradientt(pos) = dot(C*reshape((A*B), rijen*kolommen, 1)-meting, C*reshape(A*E, rijen*kolommen, 1)) - ...
           4*lambda*(norm(A,"fro")^2-norm(B,"fro")^2)*trace(B'*E);        
       pos = pos+1;
    end
end
end



