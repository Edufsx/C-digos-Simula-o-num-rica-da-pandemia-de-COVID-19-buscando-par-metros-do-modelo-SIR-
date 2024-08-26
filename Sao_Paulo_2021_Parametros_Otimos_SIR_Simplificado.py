import numpy as np
from matplotlib import pyplot as plt

def S_dot(t, S, I, R, alpha, beta):
    return -alpha*I

def I_dot(t, S, I, R, alpha, beta):
    return (alpha - beta)*I

def R_dot(t, S, I, R, alpha, beta):
    return beta*I
  
def rk4_passo_3por3(t, S_dot, I_dot, R_dot, S, I, R, alpha, beta, h):
    
    # PVI :: S' = S_dot , I' = I_dot e R' = R_dot
    # S_dot, I_dot e R_dot sao as três funcoes que definem o sistema
    # t eh o tempo no instante anterior
    # S, I, R sao as solucoes no instante de tempo anterior
    # h eh o passo de integracao
    
    k1_S = S_dot(t,       S,            I,            R,            alpha, beta)
    k1_I = I_dot(t,       S,            I,            R,            alpha, beta)
    k1_R = R_dot(t,       S,            I,            R,            alpha, beta)
    
    k2_S = S_dot(t + h/2, S + h/2*k1_S, I + h/2*k1_I, R + h/2*k1_R, alpha, beta)
    k2_I = I_dot(t + h/2, S + h/2*k1_S, I + h/2*k1_I, R + h/2*k1_R, alpha, beta)  
    k2_R = R_dot(t + h/2, S + h/2*k1_S, I + h/2*k1_I, R + h/2*k1_R, alpha, beta)     
    
    k3_S = S_dot(t + h/2, S + h/2*k2_S, I + h/2*k2_I, R + h/2*k2_R, alpha, beta)
    k3_I = I_dot(t + h/2, S + h/2*k2_S, I + h/2*k2_I, R + h/2*k2_R, alpha, beta)  
    k3_R = R_dot(t + h/2, S + h/2*k2_S, I + h/2*k2_I, R + h/2*k2_R, alpha, beta)     

    k4_S = S_dot(t + h,   S + h*k3_S,   I + h*k3_I,   R + h*k3_R,   alpha, beta)      
    k4_I = I_dot(t + h,   S + h*k3_S,   I + h*k3_I,   R + h*k3_R,   alpha, beta)   
    k4_R = R_dot(t + h,   S + h*k3_S,   I + h*k3_I,   R + h*k3_R,   alpha, beta)     

    S_novo = S + h*(k1_S + 2*k2_S + 2*k3_S + k4_S)/6    
    I_novo = I + h*(k1_I + 2*k2_I + 2*k3_I + k4_I)/6 
    R_novo = R + h*(k1_R + 2*k2_R + 2*k3_R + k4_R)/6     

    return S_novo, I_novo, R_novo


def rk4_sistemas3por3(S_dot, I_dot, R_dot, condicao_inicial, a, b, alpha, beta, h):
    
    # S_dot, I_dot e R_dot sao as tres funcoes que definem o sistema
    # condicao inicial (x0, y0, z0) para t = a
    # intevalo de integracao [a, b]
    # passo de integracao h
    
    # Passo 0 : Iniciar o vetor de tempo e os vetores para S, I e R
    t_vetor = [a]
    S0, I0, R0 = condicao_inicial
    S_vetor = [S0]
    I_vetor = [I0]
    R_vetor = [R0]    
    
    # Passo 1 :: Calcular N (numero de repeticoes)
    N = int((b-a)/h)

    # Passo 2 :: Executar N vezes o rk4_passo_3por3
    for i in range(N):
        t_novo = t_vetor[i] + h
        S_novo, I_novo, R_novo = rk4_passo_3por3(t_vetor[i], S_dot, I_dot, R_dot, \
                                                 S_vetor[i], I_vetor[i], R_vetor[i],\
                                                 alpha, beta, h)
        t_vetor.append(t_novo)
        S_vetor.append(S_novo)
        I_vetor.append(I_novo)
        R_vetor.append(R_novo)
        
    return t_vetor, S_vetor, I_vetor, R_vetor

#Sao Paulo Casos busca Parâmetros Otimos 2021 SIR ARTIGO

# Dados Reais 
dado_real = np.loadtxt('sao_paulo_casos_mm7d_2021_sub.txt')

alpha=0.541
beta=0.522

# Integracao Numerica
h = 10**(-3)
a = 0
b = 365


# Condições Iniciais
S0 = 11466825
I0 = 1117
R0 = 401718

condicao_inicial = [S0, I0, R0]


# Busca I0 Otimo
I0 = 1121
linha = 0
erro_vec=np.zeros([1122,2])
for k in range(0,1122):
    print(linha)
    I0=I0-1
    condicao_inicial = [S0, I0, R0]
    t_vector, S_vector, I_vector, R_vector = rk4_sistemas3por3(S_dot, I_dot, R_dot, \
                                                   condicao_inicial, a, b, \
                                                   alpha, beta, h)
    # Pontos Discretos que Interessam 0, 999, 2*999,...
    len_I = len(I_vector)
    I_diario = []
    for i in range(0, 365):
        t_discreto = i*999
        I_discreto = I_vector[t_discreto]
        I_diario.append(I_discreto)
    
    erro = 0
    for p in range(365):
        erro = erro + (I_diario[p] - dado_real[p])**2
    
    resultado = np.array([I0, erro])               
    erro_vec[k, :] = resultado    
    linha = linha + 1

indice_valor_minimo=np.argmin(erro_vec[:,1])
print(erro_vec[indice_valor_minimo,:])

"""
#Plotagem
t_vector, S_vector, I_vector, R_vector = rk4_sistemas3por3(S_dot, I_dot, R_dot, \
                                               condicao_inicial, a, b, \
                                               alpha, beta, h)

# Pontos Discretos que Interessam 0, 999, 2*999, ...
len_I = len(I_vector)
I_diario = []
for i in range(0, 365):
    t_discreto = i * int(len_I / 365)
    I_discreto = I_vector[t_discreto]
    I_diario.append(I_discreto)        
  
# Infectados Simulados x Infectados Oficiais
plt.figure()
plt.plot([0,I_vector[-1]], [0,I_vector[-1]], label='y = x', linestyle='--')
plt.plot(dado_real, I_diario)
plt.xlabel('Número de casos Oficiais Confirmados em São Paulo 2021')
plt.ylabel('Número de casos Simulados Confirmados em São Paulo 2021')

# Dados Simulados e Dados Oficiais x Tempo
plt.figure()
plt.plot(t_vector, I_vector)
plt.plot(dado_real)
plt.xlabel('Dias - 2021')
plt.ylabel('Número de casos Confirmados em São Paulo')
"""

"""
# Salvar dados simulados em um arquivo texto
with open('sao_paulo_casos_Simulados_2021.txt', 'w') as file:
     file.write(' '.join(map(str, I_diario))) 
"""  

"""
# Busca de Parametros alpha e beta com menor erro quadrático acumulado
alpha_vec = np.linspace(0.3, 0.7, 401)
beta_vec = np.linspace(0.3, 0.7, 401)

linha = 0
matriz_simulacao = np.zeros([len(alpha_vec)*len(beta_vec), 3])
for m in alpha_vec:
    for n in beta_vec:
        if m > n:
           print(linha)
           t_vec, S_vec, I_vec, R_vec = rk4_sistemas3por3(S_dot, I_dot, R_dot, \
                                                          condicao_inicial, a, b, \
                                                          m, n, h)
           # Pontos Discretos que Interessam 0, 999, 2*999,...
           len_I = len(I_vec)
           I_diario = []
           for i in range(0, 365):
               t_discreto = i*999
               I_discreto = I_vec[t_discreto]
               I_diario.append(I_discreto)
               # 
           # Calculo do Erro...
           erro = 0
           for k in range(365):
               erro = erro + (I_diario[k] - dado_real[k])**2
           resultado = np.array([m, n, erro])           
           matriz_simulacao[linha, :] = resultado
           linha = linha + 1
           
matriz_recortada= matriz_simulacao[0:80200,:]
indice_valor_minimo=np.argmin(matriz_recortada[:,2])

print(matriz_simulacao[indice_valor_minimo,:])
"""