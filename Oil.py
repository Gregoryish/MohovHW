#!/usr/bin/env python
# coding: utf-8

# Иходные данные
# 
#     1. density_oil_without_gas - плотность дегазированной нефти ( p_0 = 0.1 *10^6 МПа, Tст = 293 К), кг/м^3
#     2. viscosity_without_gas - вязкость нефти в стандартных условиях , мПа*с
#     3. gas_saturation  - газонасыщенность (газосодержание) пластовой нефти, т.е. отношение объёма газа, растворённого в нефти, к массе сепарированной нефти м3/т (объём газа приведен к нормальным условиям)
#     4. relative_density_gas (relative_density) - относительная плотность газа по воздуху
#     5. T_formation - пластовая температура, К
#     6. pressure_formation - пластовое давление, МПа
#     7. pressure_bubble_point_form - давление насыщения пластовой нефти при пластовой температуре, МПа
#     8. y_a, y_c1 - молярная доли азота и метана в попутном газе однократного разгазирования нефти до (0.1 МПа, 293 К)
#     
#     

# In[4]:


import numpy as np
import math


# In[2]:


#1 определяем термодинамические условия разгазирования p, T
#TODO

#2 равновесное давление насыщения при T<=Tпл

def _pressure_bbp_T (pressure_bubble_point_form, T_formation, T_current, gas_saturation, y_a, y_c1):
    pressure_bbp_T = pressure_bubble_point_form - (T_formation - T_current)/(9.157 + 701.8/(gas_saturation*(y_c1 - 0.8*y_a)))
    return round(pressure_bbp_T, 3)


# In[8]:


# 3
def _R(p, pressure_bbp_T):
    R = (1+math.log10(p))/(1+math.log10(pressure_bbp_T)) -1
    return round(R, 4)

def _m (T, density_oil_without_gas, relative_density):
    m = 1+0.029*(T-293)*(density_oil_without_gas*relative_density*10**-3 - 0.7966)
    return round(m, 4)

def _D (T, density_oil_without_gas, relative_density):
    D = (10**-3)*density_oil_without_gas*relative_density*(4.5 - 0.00305*(T-293))-4.785
    return round(D, 4)

# приведённый к нормальным условиям удельный объём выделившегося газа

def _volume_separate_gas (gas_saturation, R, m, D,):
    volume_separate_gas = gas_saturation*R*m*(D*(1+R)-1)
    return round(volume_separate_gas, 4)

    


# In[14]:


#4 рассчёт остаточной газонасыщенности нефти (удельный объём растворенного раза) в процессе её разгазирвоания
def _volume_dissolved_gas (gas_saturation, m, volume_separate_gas):
    volume_dissolved_gas = gas_saturation*m - volume_separate_gas
    return round(volume_dissolved_gas, 2)


# In[18]:


#5 относительная плотность выделившегося газа (p, T)
def _relative_density_separate_gas (relative_density, a, u, R):
    relative_density_separate_gas = a*(relative_density - 0.0036*(1+R)*(105.7 + u*R))
    return round(relative_density_separate_gas, 4)

def _a (T):
    a = 1 + 0.0054*(T-293)
    return a

def _u(density_oil_without_gas, gas_saturation):
    u = 10**-3*density_oil_without_gas*gas_saturation - 186
    return u


# In[19]:


#6 находим относительную плотность растворённого раза, остающегося в нефти при данных условиях её разгазирования (p, T)
def _relative_density_dissolved_gas (gas_saturation, a, m, relative_density_gas, relative_density_separate_gas, volume_separate_gas, volume_dissolved_gas):
    relative_density_dissolved_gas = gas_saturation*(a*m*relative_density_gas - relative_density_separate_gas*volume_separate_gas/gas_saturation)/volume_dissolved_gas
    return round (relative_density_dissolved_gas , 4)


# In[39]:


#7 рассчитываем объёмный коэффициент, 
def _b_oil(p ,T ,density_oil_without_gas, volume_dissolved_gas, lambda_T, m, alpha_n):
    b_oil = 1 + 1.0733*10**-3*density_oil_without_gas*volume_dissolved_gas*lambda_T/m        +alpha_n*(T-293) - 6.5*10**-4*p
    return round(b_oil,3)

# предварительно определив удельное приращение объёма нефти за счёт единичного изменения газонасыщенности lambda_T
def _lambda_T(density_oil_without_gas, relative_density_dissolved_gas, a, volume_dissolved_gas):
    lambda_T = 10**-3*(4.3-3.54*10**-3*density_oil_without_gas + 1.0337*relative_density_dissolved_gas/a                      +5.581*10**-6*density_oil_without_gas*(1-1.61*10**-6*density_oil_without_gas*volume_dissolved_gas)                       * volume_dissolved_gas)
    return round(lambda_T, 6)
# температурный коэффициент объёмного расширения дегазированной нефти при стандартном давлении
def _alpha_n(density_oil_without_gas):
    if 780<=density_oil_without_gas<=860:
        alpha_n = 10**-3*(3.083-2.638*10**-3*density_oil_without_gas)
    if 860<=density_oil_without_gas<=960:
        alpha_n = 10**-3*(2.513-1.975*10**-3*density_oil_without_gas)
    return round(alpha_n, 8)
        


# In[25]:


#8 определяем плотность газонасыщенной нефти

def _density_oil_with_gas (density_oil_without_gas, relative_density_dissolved_gas, volume_dissolved_gas, a, m, b_oil,):
    density_oil_with_gas = density_oil_without_gas*(1+1.293*10**-3*relative_density_dissolved_gas*volume_dissolved_gas/(a*m))/b_oil
    return round(density_oil_with_gas, 3)


# In[32]:


#9 определяем вязкость дегазированной нефти при p_0 , T

def _viscosity_without_gas_T ( density_without_gas, T, a, b,viscosity_without_gas=0):
    if viscosity_without_gas == 0:
        viscosity_without_gas = _Dunyshkin_viscosity_without_gas(density_without_gas)
    viscosity_without_gas = viscosity_without_gas*(T-293)**a*math.exp(1)**(b*(293-T))
    return round(viscosity_without_gas,3)
def _Dunyshkin_viscosity_without_gas (density_without_gas):
    if 845<density_without_gas<924:
        viscosity_without_gas = ((0.658*density_without_gas**2)/(886*10**3-density_without_gas**2))**2
    if 780<density_without_gas<=845:
        viscosity_without_gas = ((0.456*density_without_gas**2)/(833*10**3-density_without_gas**2))**2
    return round (viscosity_without_gas, 3)

def _a_viscosity (T):
    a = 10**(-0.0175*(293-T)-2.58)
    return a

def _b_viscosity (density_without_gas, T, viscosity_without_gas=0):
    if viscosity_without_gas == 0:
        viscosity_without_gas = _Dunyshkin_viscosity_without_gas(density_without_gas)
    b_viscosity = (8*10**-5*density_without_gas-0.047)*viscosity_without_gas**(0.13+0.002*(T-293))
    return round(b_viscosity,4)


# In[41]:


#10 определяем вязкость газонасыщенной нефти

def _viscosity_dissolved_gas(A_visc_dissolved, B_visc_dissolved, viscosity_without_gas_T):
    viscosity_dissolved_gas = A_visc_dissolved*viscosity_without_gas_T**B_visc_dissolved
    return viscosity_dissolved_gas

def _A_visc_dissolved (volume_dissolved_prived):
    A_visc_dissolved = 1 + 0.0129*volume_dissolved_prived - 0.0364*volume_dissolved_prived**0.85
    return round(A_visc_dissolved,4)

def _B_visc_dissolved (volume_dissolved_prived):
    B_visc_dissolved = 1 + 0.0017*volume_dissolved_prived - 0.0228*volume_dissolved_prived**0.667
    return round(B_visc_dissolved, 4)

#приведенный объём газа растворенного в нефти к стандартным условиям
def _volume_dissolved_prived (volume_dissolved_gas, density_oil_without_gas, alpha_n):
    volume_dissolved_prived = 1.055*10**-3*(1+5*alpha_n)*volume_dissolved_gas*density_oil_without_gas
    return round(volume_dissolved_prived,3)


# In[38]:


#11 рассчёт повернхностного натяжения

def _sigma_oil_gas(p,T):
    sigma_oil_gas = (1/10)**(1.58+0.05*p) - 72*10**-6*(T-305)
    return round(sigma_oil_gas, 4)


# In[ ]:





# In[ ]:




