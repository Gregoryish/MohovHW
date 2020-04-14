#!/usr/bin/env python
# coding: utf-8

# In[6]:


import numpy as np


# <b>Состав Газа и его использование для нахождения физических характеристик </b>

# In[7]:


#нормальный молярный объём в нормальных условиях
normal_volume = 22.4 
#
standart_volume = 24.05


# In[8]:


# массовая доля компонента в смеси
G = np.arange(1,3)
# объёмная (молярная) доля компонента в смеси
y = np.arange(2,4)
# молекулярная масса i-го компонента
M = np.arange(3,5)


# In[9]:


#средняя молекулярная масса газа 

def _mean_mol_massa (yG, M, how, si):
    if how == 'y':
        if si == '%':
            return _mean_mol_massa_volume_1(yG, M)
        else:
            return _mean_mol_massa_volume_2(yG, M)
    if how == 'G':
        if si == '%':
            return _mean_mol_massa_mas_1 (yG, M)
        else:
            return _mean_mol_massa_mas_2 (yG, M)

#если объёмная доля в процентах
def _mean_mol_massa_volume_1(y, M):
    y_m_multiply = np.multiply(y, M).sum()
    return round(y_m_multiply/100,3)

#если объёмная доля в долях единиц
def _mean_mol_massa_volume_2(y, M):
    y_m_multiply = np.multiply(y, M).sum()
    return round(y_m_multiply,3)

#если массовая доля в процентах
def _mean_mol_massa_mas_1 (G, M):
    g_m_multiply = np.multiply(G, 1/M).sum()
    return round(100/g_m_multiply,3)

#если массовая доля в долях единиц
def _mean_mol_massa_mas_2 (G, M):
    g_m_multiply = np.multiply(G, 1/M).sum()
    return round(1/g_m_multiply,3)


# In[10]:


# средняя плотность газа по вычисленной средней молекулярной массе

# в нормальных условиях
def _mean_density_gas_norm(mean_mol_massa):
    return round(mean_mol_massa/22.414,3)

# в стандартных условиях
def _mean_density_gas_std(mean_mol_massa):
    return round(mean_mol_massa/24.05,3)

#относительная плотность газа по воздуху
def _relative_density (mean_mol_massa, how=None):
    if how == 'norm':
        mean_density = mean_density_gas_norm(mean_mol_massa)
        return round(mean_density/1.293, 3)
    elif how == 'std':
        mean_density = mean_density_gas_std(mean_mol_massa)
        return round(mean_density/1.205, 3)
    else:
        return round(mean_mol_massa/28.98, 3)


# <b>Уравнения состояние и их использование для расчёта физических свойств газов </b>

# In[11]:


p0 = 101.325 * 1000 # , Па
T0 = 273.15 # , К
V_m = 22.414 # , м3/кмоль
R_m = round(p0*V_m/T0,3) # Дж/ (кмоль*К) универсальная газовая постоянная
def _R(mean_mol_massa, R_m=8.314*10**3):
    return round(R_m/mean_mol_massa, 3)


# In[12]:


#приведенное давление
def _pressure_priv(p, p_psev_crit):
    return round(p/p_psev_crit,3)

#приведенная температура
def _temperature_priv(T, T_psev_crit):
    return round(T/T_psev_crit,3)

# псевдокритические параметры газа
def _T_psev_crit(y, T_crit):
    y_T_crit_multiply = np.multiply(y, T_crit)
    y_T_crit_multiply = y_T_crit_multiply.sum()
    return y_T_crit_multiply

def _P_psev_crit(y, P_crit):
    y_P_crit_multiply = np.multiply(y, P_crit)
    y_P_crit_multiply = y_P_crit_multiply.sum()
    return y_P_crit_multiply

# приведенное давление по аппроксимационным формулам П.Д. Ляпкова 
def _Lyapkov_pressure_priv (p, relative_density_hydrocarbon):
    pressire_priv = p / (10**5*(46.9 - 2.06*relative_density_hydrocarbon**2))
    return round(pressire_priv, 3)
# приведенная температура по аппроксимационным формулам П.Д. Ляпкова
def _Lyapkov_temperature_priv (T, relative_density_hydrocarbon):
    #temperature_priv = T / (97 + 172*relative_density_hydrocarbon**2)
    temperature_priv = T / (97 + 172*relative_density_hydrocarbon)
    
    return round(temperature_priv, 3)

# относительная по воздуху плотность смеси УВ газов кроме азота
# на входе относительная плотность всего воздуха вычисляемая, молярная доля азота при стандартных условиях
def _relative_density_hydrocarbon(relative_density, y_a, p_a=0.970):
    relative_density_hydrocarbon = (relative_density - p_a*y_a)/(1-y_a)
    return round (relative_density_hydrocarbon, 3)


# In[13]:


# функцтя коэффициента сверхжимаемости
# y_a - объёменая доля азота
# y_oil - объёмная доля углеводородных компонентов
# z_o, z_a - к сверхжимаемости УВ, азота
def _z_supercompress(z_o, y_o, z_a, y_a):
    z_supercompress = z_o*y_o + z_a*y_a
    return z_supercompress 


# In[1]:


# для p = 0-20 МПа, T = 273-355 K

def _z_y(pressure_priv, temperature_priv):
    if temperature_priv<1.05:
        temperature_priv = 1.05
    if (0<=pressure_priv<=5) and (1.17<=temperature_priv):
        z_y = 1 - pressure_priv*(0.18/(temperature_priv-0.73)-0.135)                + 0.016*pressure_priv**3.45/temperature_priv**6.1

    if (0<=pressure_priv<=1.45) and (1.05<=temperature_priv<1.17):
        z_y = 1 - 0.23*pressure_priv - (1.88 - 1.6*temperature_priv)*pressure_priv**2
    
    if (1.45<=pressure_priv<=5) and (1.05<=temperature_priv<=1.17):
        z_y = 0.13*pressure_priv + (6.05*temperature_priv - 6.26)*temperature_priv/pressure_priv**2  
    return z_y


# In[15]:


# коэффициент сфержимаемости азота от текущих p, T

def _z_a(p, T):
    z_a = 1 + 0.564*10**-10*(T-273)**3.71*(p/10**6)**(14.7/(T-273)**0.5)
    return z_a


# In[16]:



def _density_gas_p_T(p,T,z_supercompress, density_gas_norm, p_0=101325, T_0=273):
    density_gas_p_T = density_gas_norm*p*T_0/(z_supercompress*p_0*T)
    return density_gas_p_T
def _Volume_gas_p_T (p,T,z_supercompress, V_0, p_0=101325, T_0=273):
    density_gas_p_T = V_0*z_supercompress*p_0*T/(p*T_0)
    return density_gas_p_T


# In[17]:


# y_percent = np.array([35.5, 23.9, 19.4, 2.5, 6.7, 1.8, 1.7, 1.1, 0.5, 6.9])
# y_relant = y_percent/100
# M = np.array([16.043, 30.070, 44.097, 58.124, 58.124, 72.151, 72.151, 88.178, 44.011, 28.016])
# G_relat = np.array([0.176, 0.222, 0.264, 0.045, 0.120, 0.040, 0.038, 0.029, 0.007, 0.059])
# G_percent = G_relat*100


# In[30]:


# Функция температуры потока от P_current текущего давления
def T_function_p (p_current, p_reservoir, p_flowline, T_reservoir, T_flowline):
    T_function_p = T_flowline + ((T_reservoir-T_flowline)*(p_current-p_flowline))/( p_reservoir - p_flowline)
    return round(T_function_p,3)


# In[ ]:





# In[31]:


# #well 5PM
# # скважина 2Ф(D): 
# liquid_rate = 40 #м3/сут
# B = 0 #% обводнённость 0%
# D_cas_in = 0.133 # м внутренний диаметр обсадной колонны
# d_tub_in = 0.0503 # м внутренний диаметр НКТ
# p_flowline = 0.9 # МПа давление на устье
# p_wf = 11.435 # МПа забойное давление
# H_well = 1605 # м глубина скважины
# H_tubing = 1205 # м глубина спуска труб

# gas_saturation = 56 #м3/м3 газосодержание пластовой нефти
# pressure_bubble_point_form = 9 # МПа , давление насыщения МПа
# density_oil_without_gas = 860 # кг/м3 плотность нефти дегазированной
# density_oil_reservoir = 800 #кг/м3 плотность нефти в пластовых условиях
# T_reservoir = 313 # К

# y_a, y_c1 = 0.08, 0.4 # мольная доля в процентах азота, метана
# density_gas = 1.45 # кг/м3 плотность газа при ОСР 

# #из плотности газа можно найти молекулярную массу, небходимые параметры для расчёта св-в газа при p,T

# #параметры недостающие
# T_flowline = 293
# pressure_reservoir = H_well * 1000*9.81/10**6 +2 #пластовое давление расчётное

# PI = liquid_rate/(pressure_reservoir - p_wf)/10 # м3/атм


# In[32]:


# p_current_list = np.arange(p_flowline, p_wf,0.5)
# T_current_list = np.array([T_function_p  (p_i, pressure_reservoir, p_flowline, T_reservoir, T_flowline) for p_i in p_current_list])


# In[33]:


# y_a = 0.08
# relative_density = 1.45*24.05/28.98
# p = 0.9*10**6
# T =T_flowline = 293


# In[34]:



# relative_density_hydrocarbon = _relative_density_hydrocarbon(relative_density, y_a, p_a=0.970)

# # приведенное давление по аппроксимационным формулам П.Д. Ляпкова 
# _Lyapkov_pressure_priv (p, relative_density_hydrocarbon)

# # приведенная температура по аппроксимационным формулам П.Д. Ляпкова
# _Lyapkov_temperature_priv (T, relative_density_hydrocarbon)


# In[36]:


# priv={'pressure_priv':[], 'temperature_priv':[]}
# for i in range(len(p_current_list)):
#     p = p_current_list[i]*10**6
#     T = T_current_list[i]
#     pressure_priv = _Lyapkov_pressure_priv (p, relative_density_hydrocarbon)
#     temperature_priv = _Lyapkov_temperature_priv (T, relative_density_hydrocarbon)
#     priv['pressure_priv']+=[pressure_priv]
#     priv['temperature_priv'] +=[temperature_priv]


# In[40]:


# priv


# In[43]:



# for i in range(len(priv['temperature_priv'])):
#     z=0
#     temperature_priv = priv['temperature_priv'][i]
#     pressure_priv = priv['pressure_priv'][i]
#     if temperature_priv<1.05:
#         temperature_priv = 1.05
#     if (0<=pressure_priv<=3.8) and (1.17<=temperature_priv<2):
#         z+=1
#     if (0<=pressure_priv<=1.45) and (1.05<=temperature_priv<1.17):
#         z+=1
#     if (1.45<=pressure_priv<=4) and (1.05<=temperature_priv<=1.17):
#         z_y = 0.13*pressure_priv + (6.05*temperature_priv - 6.26)*temperature_priv/pressure_priv**2  
#         z+=1
#     print(i, z)


# In[ ]:




