%% BESSLER/HARTZ

%%Autocorrelation

close all, clear, clc
%% Q1
% fonction d'autocorrelation entre t1 et t2: Rxx(t1,t2) = sigma^2*delta(t2-t1) + mu^2
% densite de proba : (1/sigma*sqrt(2*pi))*exp(-(x-mu)^2/(2*sigma^2))
% DSP : sigma*sigma
%% Q2
X1 = randn(1,1000);
X2 = randn(1,1000)*2;
X3 = randn(1,1000)*3;
% commentaire : La variance influe sur la plage de valeurs du BBG mais le
% BBG est toujours centré en 0 (moyenne nulle)

%% Q3
axe_r = -length(X1)+1:length(X1)-1;
r11 = xcorr(X1, 'unbiased');
figure, plot(axe_r, r11)
r31 = xcorr(X3, 'unbiased');
figure, plot(axe_r,r31)
% commentaire: Les fluctuations augmentent plus l'on s'éloigne du centre
% car le nombre d'echantillons considérés pour l'estimateur non-biaisé
% diminue. On observe un pic en 0 d'amplitude valant la variance du BBG.

%% Q4
r12 = xcorr(X1, 'biased');
figure, plot(axe_r, r12)
r32 = xcorr(X3, 'biased');
figure, plot(axe_r,r32)
%commentaire : Idem question 3 excepté le fait que les fluctuations
%diminuent plus on s'éloigne du centre. Cela est dû au fait que le facteur
%de division dans l'estimateur biaisé reste constant égale au nombre total
%d'échantillons (N = 1000) tandis que le nombre d'échantillons considérés dans
%la somme diminue.

%% Q5
T11 = toeplitz(r11(1000:1009));
T12 = toeplitz(r12(1000:1009));
T31 = toeplitz(r31(1000:1009));
T32 = toeplitz(r32(1000:1009));
% commentaire: sur la diagonale de la matrice de Toeplitz, on observe bien
% la variance du BBG.

%% Q6 
% En théorie, les valeurs propres de la matrice de Toeplitz sont positives
[~, valeurs_propres] = eig(T11);
valeurs_propres = nonzeros(valeurs_propres);
% En pratique, on retrouve des valeurs propres positives
%% Q7
X4 = 3+randn(1000,1)*2;

%% Q8
% La fonction d'autocorrélation de CE processus est:
% Rxx(m) = 4*delta(m)+ 3*3
r41 = xcorr(X4, 'unbiased');
r42 = xcorr(X4, 'biased');
figure, subplot(211), plot(axe_r, r41)
subplot(212), plot(axe_r, r42)
% commentaire:
%   - pour l estimateur biaisé: En zero, le résultat est
%   identique au précédent car valant mu^2 + sigma*2. la moyenne non nulle
%   implique une constante dans la somme qui vaut la moyenne au carré 
%   (9 ici). La somme se faisant sur N-m éléments, le résultat tend vers
%   droite affine de pente -mu^2/N, (N=1000 ici) et d'ordonné à l'origine
%   mu^2. 
%   - pour l estimateur non-biaisé. Resultat identique à ceux de la
%   question 3 mais avec un offset de mu^2

%% Q9
estimateur = 'unbiased';
r = 2000;
N = 1000;
mu = 3;
sigma = 2;

X = mu + randn(r,N)*sigma;
for i =1:r
    R(i,:) = xcorr(X(i,:), estimateur);
end
R_mean = mean(R);
figure, plot(R_mean)

% commentaire: Plus r augmente et moins on a de fluctuations. Donc pour
% l'estimateur non biaisé on se rapproche d'un dirac en zéros d'amplitude
% mu^2 + sigma^2 et un offset valant mu^2 ailleurs. Pour l'estimateur
% biaisé, plus r augmente, moins il y a de fluctuations et donc plus la
% pente se rapproche de la droite affine définie précédemment

%% Q10
clear, close all, clc
U = randn(1,1000);
X = filter([1 2 1], 1, U);
r = xcorr(X, 'biased');
DSP = abs(fft(r));

figure, plot(r)
figure, plot(DSP)

% commentaire: théoriquement la Rxx(t) = 16*mu^2
% + sigma^2*(6*delta(t)+4*delta(t+1) + 4*delta(t-1) + delta(t+2)
% + delta(t-2))
% On retrouve ce résultat en traçant r (=Rxx) et en regardant les valeurs
% en 0, 1 et 2. 
%% Chaines de Markov
close all, clear, clc
%% Q1
% u0 = [0.3 0.3 0.4];
u0 = [0.9 0.05 0.05];
P = [0.6 0.2 0.2; 0.15 0.7 0.15; 0.05 0.05 0.9];
% commentaire :
%   - homogène: tous les élements de la chaine sont de mêmes
%   caractéristiques (tous les noeuds correspondent à un opérateur)
%   - irréductible : on peut atteindre n'importe quel état à partir de
%   n'importe quel état
%   - apériodique : aucune périodicité dans la chaine de Markov

%% Q2

u2 = u0*P*P;
u3 = u2*P;
u4 = u3*P;

%% Q3
u_5 = ch_markov( u0, P, 5 );
u_10 = ch_markov( u0, P, 10 );

% commentaire: Le système à tendance à converger car la chaine est
% apériodique, homogène et irréductible

%methodes proposées : - Le régime permanent [p1 p2 p3] est caractérisé
%par [p1 p2 p3] = [p1 p2 p3]*P. Ce qui donne lieu à un système d'équation.
%Celui-ci est irrésolvable car il reste un degré de liberté. La seconde
%condition a vérifier est : p1 + p2 + p3 = 1
% En résolvant l'équation on trouve p3 = 12/19, p1 = 3/19 et p2 = 4/19
%                     - La seconde méthode consiste à faire une boucle while
%qui calcule pour chaque année la nouvelle matrice P(n). Le programme
%s'arrête lorsque P(n+1) = P(n)

n = 0;
P_n1 = zeros(3,1);
P_n2 = ones(3,1);

proba_n(1,:) =u0;
while (sum(P_n1==P_n2)~=3)
    n = n+1;
    P_n1 = ch_markov( u0, P, n );
    P_n2 = ch_markov( u0, P, n+1 );
    proba_n(n+1,:) = P_n1;
end

figure,
plot((0:30), proba_n(1:31,1), 'b')
hold on, plot((0:30), proba_n(1:31,2),'r')
hold on, plot((0:30), proba_n(1:31,3), 'g')
legend('O1', 'O2', 'O3')

% commentaire: quelque soit la répartition initiale, on tend vers le même
% état stationnaire