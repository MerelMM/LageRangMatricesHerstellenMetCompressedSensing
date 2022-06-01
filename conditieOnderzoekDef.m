%meer rijen dan kolommen anders matrix'
drempel = 10e-2;
rijen = 7;
rang = 3;
kolommen =6;
lambda =10;



eigenwaardes1 =[0.99 0.8 0.7999 0 0 0]
eigenwaardes2 =[0.99 0.8 0.7 0 0 0]
eigenwaardes3 =[0.99 0.8 0.2 0 0 0]
conditie1 = 1/(1-(0.7999/0.8))
conditie2 = 1/(1-(0.7/0.8))
conditie3 = 1/(1-(0.2/0.8))
%conditie2
%conditie3

U= orth(randn(rijen,kolommen));
S1 = diag(eigenwaardes1);
S2 = diag(eigenwaardes2);
S3 = diag(eigenwaardes3);
V = orth(randn(kolommen));
matrix1 = U*S1*V;
matrix2 = U*S2*V;
matrix3 = U*S3*V;


start = (rijen+kolommen-rang)*rang;
aantalMetingen = 2*start;
aantalKeer = 6;
Res1 = zeros(aantalKeer*(aantalMetingen-start+1),1);
Res2 = zeros(aantalKeer*(aantalMetingen-start+1),1);
Res3 = zeros(aantalKeer*(aantalMetingen-start+1),1);

k=0;
fouteRang = rang-1;
for i = start: aantalMetingen
    for j = 1:aantalKeer
        k=k+1;
        Res1(k) = losOp(matrix1, rijen, kolommen, rang-1,i, lambda);
        Res2(k) = losOp(matrix2, rijen, kolommen, rang-1,i, lambda);
        Res3(k) = losOp(matrix3, rijen, kolommen, rang-1,i, lambda);
    end
end

g=0;
k=0;
Rec1 = zeros((aantalMetingen-start+1),1);
Rec2 = zeros((aantalMetingen-start+1),1);
Rec3 = zeros((aantalMetingen-start+1),1);
for j=start:aantalMetingen
    g=g+1;
    fout1=0;
    fout2=0;
    fout3=0;
    for i = 1:aantalKeer
        k=k+1;
        fout1 = fout1 + Res1(k);
        fout2 = fout2 + Res2(k);
        fout3 = fout3 + Res3(k);
    end
    Rec1(g) = fout1/aantalKeer;
    Rec2(g) = fout2/aantalKeer;
    Rec3(g) = fout3/aantalKeer;
end

x=start/start:1/start:aantalMetingen/start;
fig1 = figure(1)
plot(x,Rec1,'--ks')
title('verschil 0.001, conditiegetal =  8.0000e+03',fontsize=16)
xlabel('oversamplingsfactor \phi [aantal metingen / (m+n-rang)*rang]',fontsize=16)
ylabel('gemiddelde fout',fontsize=16)
grid on

fig2 = figure(2)
plot(x,Rec2,'--ks')
title('verschil 1, conditiegetal= 8.0000',fontsize=16)
xlabel('oversamplingsfactor \phi [aantal metingen / (m+n-rang)*rang]',fontsize=16)
ylabel('gemiddelde fout',fontsize=16)
grid on

fig3 = figure(3)
plot(x,Rec3,'--ks')
title('verschil 100, conditiegetal = 1.3333',fontsize=16)
xlabel('oversamplingsfactor \phi [aantal metingen / (m+n-rang)*rang]',fontsize=16)
ylabel('gemiddelde fout',fontsize=16)
grid on

conditie1
conditie2
conditie3

function [x] = losOp(matrix, rijen, kolommen, rang, aantalMetingen, lambda)
A = randn(rijen, rang);
B = randn(rang, kolommen);
z0 = [A(:) ; B(:)];
vector = matrix(:);
[C,meting] = maakMetingen(aantalMetingen,rijen*kolommen, vector);

options = optimoptions(@fminunc,'Algorithm','trust-region','SpecifyObjectiveGradient',true, 'MaxFunctionEvaluations',999999, 'MaxIterations', 1e6,...
'OptimalityTolerance', 1e-6, 'StepTolerance', 1e-6, 'FunctionTolerance', 1e-6 );
objfun = @(z)recoverWithGradientEnPenalty(C,z, kolommen, rijen, rang, meting, lambda);
[z] = fminunc(objfun,z0,options);
g = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
x = norm(matrix-g) /norm(matrix);

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
