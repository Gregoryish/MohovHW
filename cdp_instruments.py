#!/usr/bin/env python
# coding: utf-8

# In[1]:


import gron
import Oil
import Gas

import numpy as np


# In[2]:


#функция под конкретные исходные данные, для получения значений объем газа, плотность газа при p,T
def get_gas_parametrs (p,t, q_gas_std, density_gas_std, relative_density_gas, y_a, y_c):
    """
    p - МПа,
    t - К,
    q_gas_std - суммарный расход газа в стандартных условиях выделившегося при p, t
    relative_density_gas - относительная плотность газа
    y_a, y_c - содержание азота, ув части в д.ед.
    
    volume_gas_p_t, density_gas_p_t - расход газа при [p,t], плотность газа при [p,t]
    """
    
    p=p*10**6
    relative_density_hydrocarbon = Gas._relative_density_hydrocarbon(relative_density_gas, y_a)
    p_priv = Gas._Lyapkov_pressure_priv(p, relative_density_hydrocarbon)
    t_priv = Gas._Lyapkov_temperature_priv(t, relative_density_hydrocarbon)
    
    z_c = Gas._z_y(p_priv,t_priv)
    z_a = Gas._z_a(p, t)
    z_supercompress = Gas._z_supercompress(z_c, y_c, z_a, y_a)
    density_gas_p_t = Gas._density_gas_p_T(p,t,z_supercompress,density_gas_std,T_0=294)
    volume_gas_p_t = Gas._Volume_gas_p_T(p,t,z_supercompress,q_gas_std,T_0=294)
    
    return {'q_gas_p_t':round(volume_gas_p_t,2),'density_gas_p_t':round(density_gas_p_t,3)}


# In[3]:



def get_oil_parametrs(p, t, p_bbp, t_reservoir, gas_saturation, y_a, y_c1, density_oil_without_gas, relative_density_gas):
    
    #2
    pressure_bbp_T = Oil._pressure_bbp_T (p_bbp, t_reservoir, t, gas_saturation, y_a, y_c1)

    #3
    if p < pressure_bbp_T:
        R = Oil._R(p, pressure_bbp_T)
        m = Oil._m (t, density_oil_without_gas, relative_density_gas)
        D = Oil._D (t, density_oil_without_gas, relative_density_gas)
        volume_separate_gas = Oil._volume_separate_gas(gas_saturation, R,m,D)
    else:
        volume_separate_gas = 0


    #4
    if volume_separate_gas == 0:
        volume_dissolved_gas = gas_saturation
    else:
        volume_dissolved_gas = Oil._volume_dissolved_gas(gas_saturation, m, volume_separate_gas)
        if volume_dissolved_gas>gas_saturation:
            volume_dissolved_gas = gas_saturation

    #5
    a = Oil._a(t)
    u = Oil._u(density_oil_without_gas, gas_saturation)

    if volume_separate_gas == 0:
        relative_density_separate_gas = 0
    else: 
        relative_density_separate_gas = Oil._relative_density_separate_gas(relative_density_gas, a, u, R)

    #6
    m = Oil._m (t, density_oil_without_gas, relative_density_gas)
    relative_density_dissolved_gas = Oil._relative_density_dissolved_gas (gas_saturation, a, m, relative_density_gas, relative_density_separate_gas, volume_separate_gas, volume_dissolved_gas)

    #7
    alpha_n = Oil._alpha_n(density_oil_without_gas)
    lambda_T = Oil._lambda_T(density_oil_without_gas, relative_density_dissolved_gas, a, volume_dissolved_gas)
    b_oil = Oil._b_oil(p ,t ,density_oil_without_gas, volume_dissolved_gas, lambda_T, m, alpha_n)



    #8
    density_oil_with_gas = Oil._density_oil_with_gas(
        density_oil_without_gas,
        relative_density_dissolved_gas,
        volume_dissolved_gas,
        a,
        m,
        b_oil)

    #9 
    if t<293:
        t_m = 293.1
    else:
        t_m = t
    viscosity_without_gas, a_viscosity, b_viscosity = Oil._Dunyshkin_viscosity_without_gas (density_oil_without_gas), Oil._a_viscosity(t_m), Oil._b_viscosity (density_oil_without_gas, t_m)
    
    viscosity_without_gas_T = Oil._viscosity_without_gas_T(density_oil_without_gas, t_m, a_viscosity, b_viscosity, viscosity_without_gas)

    #10 
    volume_dissolved_prived = Oil._volume_dissolved_prived(volume_dissolved_gas, density_oil_without_gas, alpha_n)
    A_visc_dissolved = Oil._A_visc_dissolved(volume_dissolved_prived)
    B_visc_dissolved = Oil._B_visc_dissolved (volume_dissolved_prived)
    viscosity_oil_dissolved_gas = round(Oil._viscosity_dissolved_gas(
        A_visc_dissolved,
        B_visc_dissolved,
        viscosity_without_gas_T),4)
    

    #11
    sigma_oil_gas = Oil._sigma_oil_gas(p,t)
    return {'viscosity_oil_dissolved_gas':viscosity_oil_dissolved_gas, 
            'density_oil_with_gas':density_oil_with_gas, 
            'b_oil':b_oil, 
            'relative_density_separate_gas':relative_density_separate_gas, 
            'volume_separate_gas':volume_separate_gas}


# In[4]:


def get_all_parametrs(p, t, p_bbp, t_reservoir, gas_saturation, y_a, y_c1, y_c, density_oil_without_gas, density_gas_std, relative_density_gas, density_water, B, q_liquid, D_tubing, a_gas = 0):

    #расчёт свойств нефти
    oil_parametrs = get_oil_parametrs(p, t, p_bbp, t_reservoir, gas_saturation, y_a, y_c1, density_oil_without_gas, relative_density_gas)

    viscosity_liquid_p_T = oil_parametrs['viscosity_oil_dissolved_gas']
    
    density_liquid_p_T = oil_parametrs['density_oil_with_gas']*(1-B) + density_water*B
    q_liquid_p_T = q_liquid*(1-B)*oil_parametrs['b_oil']+q_liquid*B
    #суммарный объём газа выделившегося в ст. условиях, м3
    q_gas_std = q_liquid*(1-B)*oil_parametrs['volume_separate_gas'] + a_gas*q_liquid
    #относительная плотность газа выделвшиегося при p,t
    relative_density_gas_p_t = oil_parametrs['relative_density_separate_gas']

    #расчёт параметров газа
    gas_parametrs = get_gas_parametrs (p,t, q_gas_std,density_gas_std,relative_density_gas_p_t, y_a, y_c)
    q_gas_p_T, density_gas_p_T = gas_parametrs['q_gas_p_t'], gas_parametrs['density_gas_p_t']

    gron_parametrs = gron.Gron_parametrs(q_liquid_p_T, q_gas_p_T, density_gas_p_T, density_liquid_p_T, viscosity_liquid_p_T, D_tubing).return_dict()

    return oil_parametrs, gas_parametrs, gron_parametrs


# In[5]:


def func_list_temperature (t_reservoir, t_head_tubing, p_reservoir, p_head_tubing, p_list, method=0):
    t_list=list()
    for p in p_list:
        t = t_head_tubing+((t_reservoir - t_head_tubing)*(p - p_head_tubing))/(p_reservoir - p_head_tubing)
        t_list.append(round(t,2))
    if method == 0:
        delta = (t_reservoir - t_list[-1])
    else:
        delta = 25
    return np.array(t_list) + delta


# In[6]:


def get_crd(q_liquid, 
            PI,  
            h_well, h_tubing,
            p_bbp, p_reservoir,
            t_reservoir, t_head_tubing, p_head_tubing,
            gas_saturation, 
            y_a, y_c1, y_c, 
            density_oil_without_gas, 
            density_gas_std,
            relative_density_gas, 
            density_water, 
            B, 
            D_tubing_in, D_casing_in, p_step = 0.1):
    """
    Расчёт КРД от забоя к устью
    Внутри функции расчитывается забойное давление через уравнение притока для q_liquid 
    
    return 
    p_list - список давлений от забоя до устья, 
    list_H - список глубин соотвествующих списку давлений, 
    list_all_parametrs - список состоящий из всех основных расчётов на каждой точке
    """
    
    #создание листов p_list, t_list
    p_wf = round(p_reservoir - q_liquid/PI ,1)
    p_list = np.arange(p_wf, 0.2, -p_step)
    p_list = np.round(p_list,3)
    t_list = func_list_temperature(t_reservoir, t_head_tubing, p_reservoir, p_head_tubing, p_list)
    t_list = np.round(t_list,3)


    #получение значений list_H 
    list_all_parametrs = []
    list_H = []
    list_H.append(h_well)
    h = h_well
    for p, t in zip(p_list, t_list):
        if h<h_tubing:
            D_tubing = D_tubing_in
        else:
            D_tubing = D_casing_in
        if len(list_all_parametrs) == 0:
            all_parametrs = get_all_parametrs(p, t, p_bbp, t_reservoir, gas_saturation, y_a, y_c1, y_c, density_oil_without_gas,density_gas_std, relative_density_gas, density_water, B, q_liquid, D_tubing)
            list_all_parametrs.append(all_parametrs)

        else:
            all_parametrs = get_all_parametrs(p, t, p_bbp, t_reservoir, gas_saturation, y_a, y_c1, y_c, density_oil_without_gas,density_gas_std, relative_density_gas, density_water, B, q_liquid, D_tubing)
            list_all_parametrs.append(all_parametrs)
            mean_gradient = (list_all_parametrs[-2][2]['gradient'] + list_all_parametrs[-1][2]['gradient'])/2
            dh = round((1/ mean_gradient)*p_step,1)
            h = h - dh
            list_H.append(round(h,1))

    list_H = np.array(list_H)
    
    return p_list, list_H, list_all_parametrs


# In[7]:


def get_crd_from_head(q_liquid,
            h_well, h_tubing,
            p_bbp, p_reservoir,
            t_reservoir, t_head_tubing, p_head_tubing,
            gas_saturation, 
            y_a, y_c1, y_c, 
            density_oil_without_gas,
            density_gas_std,
            relative_density_gas, 
            density_water, 
            B, 
            D_tubing_in, D_casing_in, p_step = 0.1, a_gas = 0):
    """
    Расчёт КРД от устья к забою
    Внутри функции расчитывается забойное давление через уравнение притока для q_liquid 
    a_gas - удельный расход газа в ст. условиях закачиваемый в скважину, для расчёта газлифтной установки
    
    return 
    p_list - список давлений от забоя до устья, 
    list_H - список глубин соотвествующих списку давлений, 
    list_all_parametrs - список состоящий из всех основных расчётов на каждой точке
    
    внимание используется подстройка температуры для схождения с тестовым исследованием
    причина - приближенный расчёт температуры дает ошибку, поэтому используется подстройка
    """
        
    #создание листов p_list, t_list
    p_list = np.arange(p_head_tubing, p_reservoir, p_step)
    p_list = np.round(p_list,3)
    t_list = func_list_temperature(t_reservoir, t_head_tubing, p_reservoir, p_head_tubing, p_list, method = 1)
    t_list = np.round(t_list,3)

    
    #получение значений list_H 
    list_all_parametrs = []
    list_H = []
    h = 0
    list_H.append(h)
    
    for p, t in zip(p_list, t_list):
        
        if True:
            D_tubing = D_tubing_in

        if len(list_all_parametrs) == 0:
            all_parametrs = get_all_parametrs(p, t, p_bbp, t_reservoir, gas_saturation, y_a, y_c1, y_c, density_oil_without_gas,density_gas_std, relative_density_gas, density_water, B, q_liquid, D_tubing, a_gas)
            list_all_parametrs.append(all_parametrs)

        else:
            all_parametrs = get_all_parametrs(p, t, p_bbp, t_reservoir, gas_saturation, y_a, y_c1, y_c, density_oil_without_gas,density_gas_std,relative_density_gas, density_water, B, q_liquid, D_tubing, a_gas)
            list_all_parametrs.append(all_parametrs)
            mean_gradient = (list_all_parametrs[-2][2]['gradient'] + list_all_parametrs[-1][2]['gradient'])/2
            dh = round((1/ mean_gradient)*p_step,1)
            h = h + dh
            list_H.append(round(h,1))

    list_H = np.array(list_H)
    
    return p_list, list_H, list_all_parametrs


# In[ ]:


import pandas as pd 
def get_data_frame(p_list, list_H, list_all_parametrs):
    
    
    new_list = list()
    for i, b in enumerate(list_all_parametrs):
        b = {**b[0], **b[1], **b[2]}
        new_list.append(b)
    df = pd.DataFrame(new_list)
    df['p'] = p_list
    df['h'] = list_H
    return df

# df.to_excel('123.xlsx')


# _____________________________________
# Исходные свойства <b>нефти, газа</b> информация по месторождению

# In[8]:


# #газосодрежание пластовой нефти, газовый фактор, м3/м3
# gas_saturation = 56 
# #давление насыщения, МПа
# p_bbp = 9 
# #плотность дегазированной нефти, кг/м3
# density_oil_without_gas = 860
# #плотность нефти в пластовых условиях, кг/м3
# density_oil_reservoir = 800

# #плотность воды, кг/м3
# density_water = 1000

# #плотность газа стандартные условия, кг/м3
# density_gas_std = 1.45
# #содержание в газе метана, %
# y_c1 = 0.4
# #содержание в газе азота, %
# y_a = 0.08
# #содержание в газе УВ части, %
# y_c = 1-y_a-0.1
# #пластовая температура, К
# t_reservoir = 40+273


# 1. Исходные даные о тех.режиме скважины
# 2. Данные исследования КРД

# Исследование - 5РМ
# скважина 2Ф(D) 

# In[9]:


# #дебит по жидкости м3/сут
# q_liquid = 40
# #обводнённость, д.ед
# B = 0
# #диаметр экспл. колонны внутренний, м
# D_casing_in = 0.133
# #диаметр нкт колонны внутренний, м
# D_tubing_in = 0.0503
# #давление на устье в НКТ, МПа
# p_head_tubing = 0.9
# #Глубина скважины, м
# h_well = 1605
# #Глубина спуска труб, м
# h_tubing = 1205

#дебит по жидкости м3/сут
q_liquid = 40
#обводнённость, д.ед
B = 0
#диаметр экспл. колонны внутренний, м
D_casing_in = 0.133
#диаметр нкт колонны внутренний, м
D_tubing_in = 0.0503
#давление на устье в НКТ, МПа
p_head_tubing = 0.4
#Глубина скважины, м
h_well = 1670
#Глубина спуска труб, м
h_tubing = 1099


# дополнительные параметры

# In[10]:


# # #давление на забое , МПа
# # p_wf = 11.435
# # #пластовое давление, МПа
# # p_reservoir = round((h_well*9.81*1000/10**6 + 2),2)
# # #температура на устье, К
# # t_head_tubing = 274
# # # К продуктивности, м3/МПа
# # PI = round(q_liquid/(p_reservoir - p_wf),3)
# # #молекулярная масса газа
# # mean_mol_mass = density_gas_std*24.05
# # #относительная плотность газа, кг/м3
# # relative_density_gas = Gas._relative_density(mean_mol_mass)
# #давление на забое , МПа

# p_wf = 11.72
# #пластовое давление, МПа
# p_reservoir = round((h_well*9.81*1000/10**6 + 2),2)
# #температура на устье, К
# t_head_tubing = 274
# # К продуктивности, м3/МПа
# PI = round(q_liquid/(p_reservoir - p_wf),3)
# #молекулярная масса газа
# mean_mol_mass = density_gas_std*24.05
# #относительная плотность газа, кг/м3
# relative_density_gas = Gas._relative_density(mean_mol_mass)


# In[11]:


# # #данные исследований
# # h_test = [0, 105, 305, 505, 705, 905, 1105, 1305, 1405, 1505, 1605]
# # p_test = [0.9, 1.12, 1.83, 2.957, 4.355, 5.785, 7.3, 8.953, 9.863, 10.176, 11.435]

# #данные исследований
# h_test = [0, 263, 463, 663, 863, 1063, 1263, 1463, 1663]
# p_test = [0.4, 1.28, 2.29, 3.82, 5.11, 6.75, 8.23, 10.0, 11.72]


# _____________________________________________________________
# создание расчёта КРД

# In[14]:


# q_liquid = 40
# a_gas = 80


# In[ ]:





# In[15]:


# p_list, list_H, list_all_parametrs = get_crd_from_head(
#             q_liquid, 
#             h_well, h_tubing,
#             p_bbp, p_reservoir,
#             t_reservoir, t_head_tubing, p_head_tubing,
#             gas_saturation, 
#             y_a, y_c1, y_c, 
#             density_oil_without_gas, 
#             density_gas_std,
#             relative_density_gas, 
#             density_water, 
#             B, 
#             D_tubing_in, D_casing_in, a_gas=a_gas)


# In[ ]:





# In[ ]:





# In[16]:


# import plotly.graph_objects as go

# # Create traces
# fig = go.Figure()
# fig.add_trace(go.Scatter(x=p_test, y=h_test,
#                     mode='lines+markers',
#                     name='КРД-Исследование'))

# fig.add_trace(go.Scatter(x=p_list, y=list_H,
#                     mode='lines+markers',
#                     name='КРД-Метод Грона'))
# fig.update_layout(
#     width = 800,
#     height = 800,
#     title = "КРД",
#     yaxis = dict(
#       range = (h_well+100,0)
#     ),
#     xaxis = dict(
#       range = (0,p_reservoir)))
# fig.show()


# In[ ]:





# In[ ]:





# In[ ]:




