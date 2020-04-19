#!/usr/bin/env python
# coding: utf-8

# In[7]:


import numpy as np
import math

import Oil
import Gas


# In[8]:


10**-3


# In[9]:


np.exp(-1)


# In[10]:


def coef_1 (relative_viscosity, D_tubing):
    """
    relative_viscosity - относительная вязкость жидкости, ед.
    D_tubing - внутренний диаметр трубы НКТ
    """
    D_0 = 0.015
    c_1 = 2.2361*np.exp(0.049*relative_viscosity)/(1+1.1002*np.exp(0.049*relative_viscosity))
    c_1 = c_1 - 8.17*(10**-3)*(relative_viscosity**(0.6))*(D_tubing/D_0 - 1)
    return round(c_1, 4)


# In[11]:


def coef_2 (relative_viscosity, D_tubing):
    """
    relative_viscosity - относительная вязкость жидкости, ед.
    D_tubing - внутренний диаметр трубы НКТ
    """
    
    D_0 = 0.015
    
    if relative_viscosity <=40:
        c_2 = (1+0.1082*np.exp(0.049*relative_viscosity))/(1+1.1002*np.exp(0.049*relative_viscosity))
        c_2 = c_2 - (0.1006 - 2.52*(10**-3)*(relative_viscosity-1))*(D_tubing/D_0 - 1)
    
    if relative_viscosity>40:
        c_2 = (1+0.1082*np.exp(0.049*relative_viscosity))/(1+1.1002*np.exp(0.049*relative_viscosity))
    return round(c_2, 4)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


def func_fi_gas(beta_gas, fruda_fluid, c_1, c_2):
    """
    beta_gas - объемное газосодержание,
    fruda_fluid - число фруда, 
    c_1 - коеф с_1, 
    c_2 - коеф с_2
    """
    fi = beta_gas/(c_1 + c_2*fruda_fluid**(-0.5))
    return fi


# In[ ]:


def func_velocity_fluid(q_liquid_p_T, q_gas_p_T, D_tubing):
    """
    q_liquid_p_T - расход жидкости при [p,T], м3/сут 
    q_gas_p_T - расход газа при  [p,T], м3/cут
    D_tubing - диаметр трубы, м
    """
    q_liquid_p_T = q_liquid_p_T/86400
    q_gas_p_T = q_gas_p_T/86400
    velocity_fluid = (q_liquid_p_T + q_gas_p_T)/((np.pi*D_tubing**2)/4)
    return velocity_fluid

def func_fruda_fluid (velocity_fluid, D_tubing):
    """
    velocity_fluid - скорость ГЖС, м/с, 
    D_tubing - диаметр трубы, м
    
    """    
    g = 9.81
    fruda_fluid = velocity_fluid**2/(g*D_tubing)
    return fruda_fluid

def func_beta_gas (q_gas_p_T, q_liquid_p_T):
    """
    q_gas_p_T, q_liquid_p_T - расход газа, жидкости при p,T условиях
    
    """
    
    beta_gas = q_gas_p_T/(q_gas_p_T+q_liquid_p_T)
    return beta_gas


# In[ ]:


def func_relative_viscosity(viscosity_liquid_p_T):
    """
    viscosity_liquid_p_T - вязкость жидкости при [p, T], мПа*с
    
    relative_density - на выходе относительная вязкость жидкости
    """
    # 1*10**-3 вязкость воды, Па*с
    relative_viscosity = viscosity_liquid_p_T*10**-3 / (1*10**-3)
    
    return relative_viscosity


# In[7]:


#градиент потерь на трение

def func_gradient_friction (labmda_friction, velocity_fluid, density_fluid_p_T, D_tubing):
    """
    labmda_friction - гидравлич. коэф потерь на трение
    velocity_fluid - скорость смеси ГЖС, м/c
    density_fluid_p_T - плотность ГЖС, кг/м3
    D_tubing - внутренний диаметр НКТ, м
    
    gradient_friction - градиент потерь на трение, МПа/м
    """
    gradient_friction = labmda_friction*density_fluid_p_T*velocity_fluid**2*(10**-6)/(2*D_tubing)
    return gradient_friction


# In[8]:


def func_reinolds_liquid (velocity_fluid, density_liquid_p_T, viscosity_liquid_p_T, D_tubing):
    """
    velocity_fluid - скорость ГЖС, м/с 
    viscosity_liquid_p_T - вязкость жидкости [p,T], мПа*с
    D_tubing - диаметр НКТ, м
    
    
    reinolds_liquid - критерий Рейнольдса потока жидкости, движущегося со скоростью, равной смеси (w ж = w см)
    """
    
    reinolds_liquid = velocity_fluid*D_tubing*density_liquid_p_T/(viscosity_liquid_p_T*10**(-3))
    
    return reinolds_liquid


# In[9]:


def func_labmda_friction(reinolds_liquid, D_tubing, eps = 0.0015*10**-3):
    """
    reinolds_liquid - число рейнольдса для жидкости со скоростью смеси, 
    eps = 0.0015*10**-3
    
    lambda_friction - коэф. гидравлических потерь
    """
    lambda_friction = 0.067*(158/reinolds_liquid+2*(eps/D_tubing))**(1/5)
    
    return round(lambda_friction,4)


# In[10]:


def func_gradient_static (density_fluid_p_T, cos_a = 1):
    """
    density_fluid_p_T - плотность ГЖС при [p,T], кг/м3
    cos_a - косинус а, где а - угол отклонения траектории скважины от вертикали
    
    gradient_static - градиент потерь на гидростатику, МПа/м
    """
    gradient_static = density_fluid_p_T*9.81*10**-6*cos_a
    
    return gradient_static    


# In[11]:


def func_gradient(gradient_static, gradient_friction):
    return gradient_static + gradient_friction


# In[12]:


def func_density_fluid (density_liquid_p_T, density_gas_p_T, fi_gas):
    """
    density_liquid_p_T - плотность жидкости [p, T], 
    density_gas_p_T - плотность газа [p, T], 
    fi_gas - истинное газосодержание
    """
    density_fluid = (1-fi_gas)*density_liquid_p_T+fi_gas*density_gas_p_T
    
    return round(density_fluid,1)
  


# In[13]:


# def get_gron_parametrs(q_liquid_p_T, q_gas_p_T, density_gas_p_T, density_liquid_p_T, viscosity_liquid_p_T, D_tubing):
#     relative_viscosity = func_relative_viscosity(viscosity_liquid_p_T)
#     c_1 = coef_1(relative_viscosity, D_tubing)
#     c_2 = coef_2(relative_viscosity, D_tubing)
#     velocity_fluid = func_velocity_fluid(q_liquid_p_T, q_gas_p_T, D_tubing)
#     beta_gas = func_beta_gas (q_gas_p_T, q_liquid_p_T)
#     fruda_fluid = func_fruda_fluid(velocity_fluid,  D_tubing)
#     fi_gas = func_fi_gas (beta_gas, fruda_fluid, c_1, c_2)
#     density_fluid_p_T = func_density_fluid (density_liquid_p_T, density_gas_p_T, fi_gas)
#     reinolds_liquid = func_reinolds_liquid (velocity_fluid, density_liquid_p_T, viscosity_liquid_p_T, D_tubing)
#     labmda_friction = func_labmda_friction(reinolds_liquid, D_tubing, eps = 0.0015*10**-3)
#     gradient_friction = func_gradient_friction (labmda_friction, velocity_fluid, density_fluid_p_T, D_tubing)
#     gradient_static = func_gradient_static (density_fluid_p_T, cos_a = 1)
#     gradient = func_gradient(gradient_static, gradient_friction)
    
#     return 


# In[13]:


class Gron_parametrs:
    def __init__(self, q_liquid_p_T, q_gas_p_T, density_gas_p_T, density_liquid_p_T, viscosity_liquid_p_T, D_tubing):
        self.relative_viscosity = func_relative_viscosity(viscosity_liquid_p_T)
        self.c_1 = coef_1(self.relative_viscosity, D_tubing)
        self.c_2 = coef_2(self.relative_viscosity, D_tubing)
        self.velocity_fluid = func_velocity_fluid(q_liquid_p_T, q_gas_p_T, D_tubing)
        self.beta_gas = func_beta_gas (q_gas_p_T, q_liquid_p_T)
        self.fruda_fluid = func_fruda_fluid(self.velocity_fluid,  D_tubing)
        self.fi_gas = abs(func_fi_gas (self.beta_gas, self.fruda_fluid, self.c_1, self.c_2))
        self.density_fluid_p_T = func_density_fluid (density_liquid_p_T, density_gas_p_T, self.fi_gas)
        self.reinolds_liquid = func_reinolds_liquid (self.velocity_fluid, density_liquid_p_T, viscosity_liquid_p_T, D_tubing)
        self.labmda_friction = func_labmda_friction(self.reinolds_liquid, D_tubing)
        self.gradient_friction = func_gradient_friction (self.labmda_friction, self.velocity_fluid, self.density_fluid_p_T, D_tubing)
        self.gradient_static = func_gradient_static (self.density_fluid_p_T)
        self.gradient = func_gradient(self.gradient_static, self.gradient_friction)
    def return_dict(self):
        return self.__dict__


# In[ ]:


# a = Gron_parametrs(q_liquid_p_T, q_gas_p_T, density_gas_p_T, density_liquid_p_T, viscosity_liquid_p_T, D_tubing)
# b = Gron_parametrs(q_liquid_p_T, q_gas_p_T, density_gas_p_T, density_liquid_p_T, viscosity_liquid_p_T, D_tubing)


# In[ ]:


# q_liquid_p_T, q_gas_p_T, density_gas_p_T, density_liquid_p_T, viscosity_liquid_p_T, D_tubing = 30,20, 28.2, 941.2, 112.6, 0.0635


# In[ ]:


# relative_viscosity = func_relative_viscosity(viscosity_liquid_p_T)
# c_1 = coef_1(relative_viscosity, D_tubing)
# c_2 = coef_2(relative_viscosity, D_tubing)
# velocity_fluid = func_velocity_fluid(q_liquid_p_T, q_gas_p_T, D_tubing)
# beta_gas = func_beta_gas (q_gas_p_T, q_liquid_p_T)
# fruda_fluid = func_fruda_fluid(velocity_fluid,  D_tubing)
# fi_gas = func_fi_gas (beta_gas, fruda_fluid, c_1, c_2)
# density_fluid_p_T = func_density_fluid (density_liquid_p_T, density_gas_p_T, fi_gas)
# reinolds_liquid = func_reinolds_liquid (velocity_fluid, density_liquid_p_T, viscosity_liquid_p_T, D_tubing)
# labmda_friction = func_labmda_friction(reinolds_liquid, D_tubing, eps = 0.0015*10**-3)
# gradient_friction = func_gradient_friction (labmda_friction, velocity_fluid, density_fluid_p_T, D_tubing)
# gradient_static = func_gradient_static (density_fluid_p_T, cos_a = 1)
# gradient = func_gradient(gradient_static, gradient_friction)


# In[ ]:




