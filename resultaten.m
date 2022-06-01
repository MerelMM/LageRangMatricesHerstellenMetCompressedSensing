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
itgew = zeros(aantalKeer*(aantalMetingen-start+1),1);
itgrad = zeros(aantalKeer*(aantalMetingen-start+1),1);
itpen = zeros(aantalKeer*(aantalMetingen-start+1),1);
itGradEnPen = zeros(aantalKeer*(aantalMetingen-start+1),1);
gewoon = zeros(aantalKeer*(aantalMetingen-start+1),1);
gradient = zeros(aantalKeer*(aantalMetingen-start+1),1);
penalty = zeros(aantalKeer*(aantalMetingen-start+1),1);
gradEnPen = zeros(aantalKeer*(aantalMetingen-start+1),1);
for i= start:aantalMetingen
    for j=1:aantalKeer
        k=k+1;
        [gewoon(k),gradient(k),penalty(k),gradEnPen(k),itgew(k),itgrad(k),itpen(k),itGradEnPen(k)]=allRecs(rijen,kolommen,rang,i,lambda);
    end
end

%%gemiddelde
g=0;
k=0;
succesRecoveryGew = zeros((aantalMetingen-start+1),1);
succesRecoveryGrad = zeros((aantalMetingen-start+1),1);
succesRecoveryPen = zeros((aantalMetingen-start+1),1);
gemItGew = zeros((aantalMetingen-start+1),1);
gemItGrad = zeros((aantalMetingen-start+1),1);
gemItPen = zeros((aantalMetingen-start+1),1);
for j=start:aantalMetingen
    g=g+1;
    sgew=0;
    sgrad=0;
    spen=0;
    igew=0;
    igrad=0;
    ipen=0;
    for i = 1:aantalKeer
        k=k+1;
        igew = igew + itgew(k);
        igrad= igrad+ itgrad(k);
        ipen=ipen+itpen(k);
        if (gewoon(k)<=drempel)
            sgew = sgew + 1;
        end
        if (gradient(k)<=drempel)
            sgrad = sgrad + 1;
        end
        if (penalty(k)<=drempel)
            spen = spen + 1;
        end
    end
    succesRecoveryGew(g) = sgew/aantalKeer;
    succesRecoveryGrad(g) = sgrad/aantalKeer;
    succesRecoveryPen(g) = spen/aantalKeer;
    gemItGew(g) = igew/aantalKeer;
    gemItGrad(g) = igrad/aantalKeer;
    gemItPen(g) = ipen/aantalKeer;
end
x=start/start:1/start:aantalMetingen/start;
fig1 = figure(1)
plot(x,succesRecoveryGew,'--ks')
title('0,5‖𝐶∗𝑣𝑒𝑐(𝑥) −𝑦‖2^2',fontsize=16)
xlabel('oversamplingsfactor \phi [aantal metingen / (m+n-rang)*rang]',fontsize=16)
ylabel('succesvol hersteld [%]',fontsize=16)
grid on

fig2 = figure(2)
plot(x,gemItGew,'--ks')
title('gemiddeld aantal iteraties',fontsize=16)
xlabel('oversamplingsfactor \phi [aantal metingen / (m+n-rang)*rang]',fontsize=16)
ylabel('succesvol hersteld [%]',fontsize=16)
grid on

fig3 = figure(3)
plot(x,succesRecoveryGrad,'--ks')
title('met meegegeven gradiënt',fontsize=16)
xlabel('oversamplingsfactor \phi [aantal metingen / (m+n-rang)*rang]',fontsize=16)
ylabel('succesvol hersteld [%]',fontsize=16)
grid on

fig4 = figure(4)
plot(x,gemItGrad,'--ks')
title('gemiddeld aantal iteraties',fontsize=16)
xlabel('oversamplingsfactor \phi [aantal metingen / (m+n-rang)*rang]',fontsize=16)
ylabel('succesvol hersteld [%]',fontsize=16)
grid on

fig5 = figure(5)
plot(x,succesRecoveryPen,'--ks')
title('met penalty',fontsize=16)
xlabel('oversamplingsfactor \phi [aantal metingen / (m+n-rang)*rang]',fontsize=16)
ylabel('succesvol hersteld [%]',fontsize=16)
grid on

fig6 = figure(6)
plot(x,gemItPen,'--ks')
title('gemiddeld aantal iteraties',fontsize=16)
xlabel('oversamplingsfactor \phi [aantal metingen / (m+n-rang)*rang]',fontsize=16)
ylabel('succesvol hersteld [%]',fontsize=16)
grid on
