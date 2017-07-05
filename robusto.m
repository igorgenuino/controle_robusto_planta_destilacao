%  Controle Robusto Multivariavel - Barreiras

tic

%% Definiçao das matrizes da planta do modelo da coluna de destilação
Ap=[-0.005131 0 0 0 0;
0 -0.07366 0 0 0;
0 0 -0.1829 0 0;
0 0 0 -0.462 0.9895;
0 0 0 -0.9895 -0.462];

Bp=[-0.629 0.624;
    0.055 -0.172;
    0.03 -0.108;
    -0.186 -0.139;
    -1.23 -0.056];

Cp=[-0.7223 -0.517 0.3386 -0.1633 0.1121;
    -0.8913 0.4728 0.9876 0.8425 0.2186];

%% Parâmetros do compensador para o desempenho e a estabilidade
alfa_r = 0.05;
wr = 0.5;

Passos_wL = 5;   %  0.5 rad/s   (Frequência baixa)
Passos_wH = 100; %  10   rad/s   (Frequência alta)
w_var = 0.1:0.1:500;

%% Simplificação da notação
A = Ap;
B = Bp;
C = Cp;

%% Loop de variações no parametro A33 e varredura pelas frequências estipuladas
k = 1;
EM = zeros();
for w = w_var
    E = zeros(1,20);
    for m = 1:20
       
        A1=[-0.005131 0 0 0 0;
            0 -0.07366 0 0 0;
            0 0 -0.16-0.002*m 0 0;
            0 0 0 -0.462 0.9895;
            0 0 0 -0.9895 -0.462];

        ROT  = 400 / (-w^2 + 20*1j*w + 625) * eye(2);
        GN   = C / (1j*w*eye(5) - A) * B;
        G1   = C / (1j*w*eye(5) - A1) * B;
        GR   = ROT*G1;
        ERR  = (GN - GR) / GN;
        ERR1 = svd(ERR);  % Cálculo do valores singulares
        E(m) = max(ERR1); % Construção do vetor de valores singulares máximos para cada variação dentro de uma frequência
    end
    
    % Construção do vetor de valores singulares máximos para todas as frequências
    EM(k) = max(E);
    k = k+1;

end

%% Construção das barreiras de estabilidade e desempenho
BARR_ESTAB =  1 ./ EM;
BARR_DES = 20 ./ (1 - EM);

BARR_GERAL = NaN(size(w_var));
BARR_GERAL( 1:Passos_wL ) = 20 * log10(BARR_DES(1:Passos_wL));
BARR_GERAL(Passos_wH:end) = 20 * log10(BARR_ESTAB(Passos_wH:end));

% Plotando as barreiras
semilogx(w_var, BARR_GERAL,'LineWidth',2.3);
grid

%% Definição das matrizes expandidas do sistema

% Ap = 5x5
% Bp = 5x2
% Cp = 2x5
% Dp = 2x2

A2 = [zeros(2,2)  zeros(2,5)
               B           A];
B2 = [    eye(2)
      zeros(5,2)];
C2 = [zeros(2,2)  C];
D2 = zeros(2,2);

LL = -inv(C / A * B);
LH = -A \ B * LL;

L = [LL
     LH];

% Escolha do valor de mi
mi = 0.0027;

% Plotagem dos valores singulares para o sistema com matrizes expandidas

SYS = ss(A2, L, 1/sqrt(mi)*C2, D2);

hold on
sigma(SYS,[0.1,500])
grid

%% Resolução da Equação Algébrica de Ricatti para o filtro de Kalman

EAR = are(A2', 1/mi*(C2'*C2), L*L');

% Cálculo da matriz de ganhos
H = 1/mi * EAR * C2';

% Plotagem dos valores singulares da malha objetivo
SYS2 = ss(A2, H, C2, D2);

hold on
sigma(SYS2, '--', {0.1, 500});

%% Cálculo dos pólos e zeros de Gn
SYS3 = ss(A2, B2, C2, D2);
ZE = tzero(SYS3);
PO = pole(SYS3);

%% Escolha do valor de rho
rho = 1e-10;    % Geralmente não precisa passar de 1e-10

% Resolução da Equação Algébrica de Ricatti
EAR2 = are(A2, 1/rho*(B2*B2'), C2'*C2);

% Calculo da matriz de ganhos de realimentação de estados
G = 1/rho * B2' * EAR2;

%% Plotagem dos valores de GnK (verificação da proximidade com a malha objetivo)
w2_var = 0.1:0.1:500;
k = 1;

for w2 = w2_var
    COMP =  G / (1j*w2*eye(7) - A2 + B2*G + H*C2) * H;
    PLAN = C2 / (1j*w2*eye(7) - A2) * B2;
    SIST  = PLAN * COMP;
    VS_GnK   = svd(SIST);
    VS_MAX_GnK(k) = max(VS_GnK);
    VS_MIN_GnK(k) = min(VS_GnK);
    k = k + 1;
end

hold on
grid
semilogx(w2_var, 20*log10(BARR_ESTAB), 'k--','LineWidth',2);
semilogx(w2_var, 20*log10(VS_MAX_GnK), 'g')
semilogx(w2_var, 20*log10(VS_MIN_GnK), 'g-')

%% Plotar os valores máximos do sistema com compensador em malha fechada e comparar com a barreira de robustez da estabilidade
k = 1;
for w2 = w2_var
    COMP =  G / (1j*w2*eye(7) - A2 + B2*G + H*C2) * H;
    PLAN = C2 / (1j*w2*eye(7) - A2) * B2;
    SIST  = (eye(2) + PLAN * COMP) \ PLAN * COMP;
    VS_Cn   = svd(SIST);
    VS_MAX_Cn(k) = max(VS_Cn);
    k = k + 1;
end

hold on
grid
semilogx(w2_var, 20*log10(VS_MAX_Cn), 'r')

%% Definição das matrizes no espaço de estados do sistema com compensador em malha fechada
AF = [(A2 - B2*G)         B2*G
       zeros(7,7)  (A2 - H*C2)];

BF = [zeros(7,2)
               H];

CF = [C2  zeros(2,7)];

DF = zeros(2,2);

SYSF = ss(AF, BF, CF, DF);

%% Inclusão de integradores
t = 0:0.01:20;
for i = 1:2001
    r(1,i) = 1;
    r(2,i) = 1;
end

%% Plotar a resposta ao degrau
figure
lsim(SYSF, r, t);

toc;
