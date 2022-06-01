
drempel = 10e-5;
aantalKeer = 10;
x = zeros(200,1);
k=0;
m=0;
z= zeros(200,1);

rijen= 6;
kolommen= 5;
rang= 3;
aantalMetingen= 20;
%
for i=1:aantalMetingen
    for j=1:aantalKeer
        k=k+1;
        [x(k)] = gradientEnPenalty(i);
    end
end

%%gemiddelde
g=0;
k=0;
succesRecovery = zeros(20,1);

for j=1:aantalMetingen
    g=g+1;
    succes=0;
    for i = 1:aantalKeer
        k=k+1;
        if (x(k)<=drempel)
            succes = succes + 1;
        end

    end
    succesRecovery(g) = succes/aantalKeer;

end
x= 1:1:20;
fig1= figure(1)
plot(x,succesRecovery, '--ks')
title('random meting')
xlabel('aaltal metingen')
ylabel ('aantal iteraties')
grid on
