
# coding: utf-8

# ### Fonctions de résolution explicite

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
#paramètres
a=0.03
b=0.001
c=0.05
d=0.0002    

ordre=1
ti=0
te=750
step=1


v_ini = np.array([200, 50])

#équations Lokta-Verra
def derivee_v(v, t):
    v=np.array([v[0]*(a-b*v[1]),v[1]*(-c+d*v[0])])
    return v

#méthodes de résolution 
def Euler(start, end, step, v_ini, derivee):
    interval = end - start                     # Intervalle
    num_points = int(interval / step) #+1      # Nombre d'éléments
    t = np.linspace(start, end, num_points)    # Tableau temps t

    v = np.empty((2, num_points))
    v[:, 0] = v_ini

    for i in range(num_points - 1):
        v[:, i + 1] = v[:, i] + step * derivee(v[:, i], t[i])

    return t, v 
def rk2(start, end, step, v_ini, derivee, ordre):
    interval = end - start                     # Intervalle
    num_points = int(interval / step)#+1       # Nombre d'éléments
    t = np.linspace(start, end, num_points)    # Tableau temps t

    v = np.empty((2, num_points))
    v[:,0]=v_ini

    for i in range(num_points - 1):
        d1 = derivee(v[:, i], t[i])
        v1 = v[:, i] + step / 2 * d1
        d2 = derivee(v1, t[i] + step / 2)
        v[:, i + 1] = v1 + step / 2 * d2

    return t, v


def rk4(start, end, step, v_ini, derivee, ordre):
    interval = end - start                     # Intervalle
    num_points = int(interval / step)#+1      # Nombre d'éléments
    t = np.linspace(start, end, num_points)    # Tableau temps t

    v = np.empty((2, num_points))
    v[:,0] = v_ini

    for i in range(num_points - 1):
        d1 = derivee(v[:, i], t[i])
        d2 = derivee(v[:, i] + step / 2 * d1, t[i] + step / 2)
        d3 = derivee(v[:, i] + step / 2 * d2, t[i] + step / 2)
        d4 = derivee(v[:, i] + step * d3, t[i] + step)
        v[:, i + 1] = v[:, i] + step / 6 * (d1 + 2 * d2 + 2 * d3 + d4)

    return t, v


# # Méthodes explicites

# In[2]:


t, v = Euler(ti,te,step, v_ini, derivee_v)
t, Drk2 = rk2(ti, te, step, v_ini, derivee_v, ordre)
t, Drk4 = rk4(ti, te,step, v_ini, derivee_v, ordre)
#plot des proies pour les méthodes explicites
plt.plot(t,v[0],label = 'Euler',color='red') 
plt.plot(t, Drk2[0, :],label = 'RK2',color='orange')
plt.plot(t, Drk4[0, :],label = 'RK4',color='black')
plt.xlabel('temps')
plt.ylabel('population proie')
plt.title("Evolution de la population des proies dans le temps")
plt.axis("equal") # Pour avoir des axes isométriques
plt.grid()
plt.legend()
plt.show()


# In[3]:


t, v = Euler(ti,te,step, v_ini, derivee_v)
t, Drk2 = rk2(ti, te, step, v_ini, derivee_v, ordre)
t, Drk4 = rk4(ti, te,step, v_ini, derivee_v, ordre)

plt.plot(t,v[1],label = 'Euler',color='red')
plt.plot(t, Drk2[1, :],label = 'RK2',color='orange')
plt.plot(t, Drk4[1, :],label = 'RK4',color='black')
plt.xlabel('temps')
plt.ylabel('population prédateur')
plt.title("Evolution de la population des prédateurs dans le temps")
plt.grid()
plt.legend()
plt.show()


# ## Evolution proie predateur méthode rk4

# In[4]:


t, Drk4 = rk4(ti, te,step, v_ini, derivee_v, ordre)
plt.plot(t, Drk4[0, :],label = 'proie',color='green')
plt.plot(t, Drk4[1, :],label = 'prédateur',color='black')
plt.xlabel('temps')
plt.ylabel('population prédateur')
plt.title("Evolution de la population des prédateurs dans le temps")
plt.axis("equal") # Pour avoir des axes isométriques
plt.grid()
plt.legend()
plt.show()


# # Méthode simplectique : méthode d'Euler semi-implicite.

# In[5]:


ti = 0
te = 750 
step = 1

alpha=0.03
beta=0.001
gamma=0.05
delta =0.0002


v_ini = np.array([200, 50])

def Lokta(x,y,t):
    
    dx = x*(alpha-beta*y)
    dy = y*(delta*x-gamma)
    return dx,dy

def Euler_simplectique(start, end, step,v_ini, derivee):
   
    t = np.arange(start,end,step)
    
    x = np.empty(t.size)
    y = np.empty(t.size)

    x[0] = v_ini[0]
    y[0] = v_ini[1]

    for i in range(t.size-1):
        y[i+1] = y[i] + step * derivee(x[i],y[i],t[i])[1]
        x[i+1] = x[i] + step * derivee(x[i],y[i+1],t[i])[0] # semi implicite : on utilise y_i+1 et non y_i pour calculer x_i+1.
    return t, x, y


# ### Proie en fonction de prédateur

# In[7]:


v_ini = np.array([[200, 50],[250,62.5],[300,75],[350,87.5],[400,100],[450,112.5],[500,125]])
for i in range(int(v_ini.size/2)):
    t, Drk4 = rk4(ti, te, step, v_ini[i], derivee_v, ordre)
    plt.plot(Drk4[0,:],Drk4[1,:])
    plt.scatter(v_ini[i][0],v_ini[i][1]) # condition initiale
plt.title('population des proies en fonction de celle des predateurs pour différentes conditions initiales.')
plt.xlabel('population des proies')
plt.ylabel('population des predateurs')
plt.show()


# # Etude sur les pas de temps

# ## Méthode Euler

# In[8]:


steps =(0.5,1,2,3)
v_ini = np.array([200, 50])

for i in range(len(steps)):
    t, v = Euler(ti,te,steps[i], v_ini, derivee_v)
    plt.plot(t,v[0])
plt.title("Evolution population proie méthode d'Euler pour différents pas de temps")
plt.show()


# # Etude sur les méthode d'intégration : comparaison des maximums

# In[9]:


step = 1
time = 750
t, v = Euler(ti,te,step, v_ini, derivee_v)
t, Drk2 = rk2(ti, te, step, v_ini, derivee_v, ordre)
t, Drk4 = rk4(ti, te,step, v_ini, derivee_v, ordre)
t,x1,y1 = Euler_simplectique(ti,te,step,v_ini,Lokta)

plt.plot(t,v[0],label = 'Euler (RK1)',color='red')
plt.plot(t, Drk2[0, :],label = 'RK2',color='orange')
plt.plot(t, Drk4[0, :],label = 'RK4',color='black')
plt.plot(t,x1,label ='simplectique ordre 1',color = 'blue')
plt.xlabel('temps [anneés]')
plt.ylabel('populations des proies')
plt.title("Evolution de la population des proies dans le temps pour différentes méthodes de résolution")
#plt.axis("equal") # Pour avoir des axes isométriques
plt.grid()
plt.legend()
legend = plt.legend(loc="upper left", edgecolor="black")
legend.get_frame().set_alpha(None)
legend.get_frame().set_facecolor((0, 0, 1, 0.1))
plt.show()


# In[16]:


def maximum_proie(start, end, step, v_ini, derivee,ordre,derivee_simplectique):
    t, v = Euler(start, end, step, v_ini, derivee)
    t, Drk2 = rk2(start, end, step, v_ini, derivee, ordre)
    t, Drk4 = rk4(start, end, step, v_ini, derivee, ordre)
    t,x1,y1 = Euler_simplectique(start, end, step, v_ini,derivee_simplectique)
    
    max_euler = ()
    max_rk2 = ()
    max_rk4 = ()
    max_simplectique = ()
    
    pop_proie = v[0]
    # si la valeur 
    for i in range(pop_proie.size-2):
        if pop_proie[i+1] >= pop_proie[i] and pop_proie[i+1] >= pop_proie[i+2]: 
            max_euler = np.append(max_euler,pop_proie[i+1])
    pop_proie = Drk2[0,:]
    for i in range(pop_proie.size-2):
        if pop_proie[i+1] >= pop_proie[i] and pop_proie[i+1] >= pop_proie[i+2]:
            max_rk2 = np.append(max_rk2,pop_proie[i+1])
    pop_proie = Drk4[0,:]
    for i in range(pop_proie.size-2):
        if pop_proie[i+1] >= pop_proie[i] and pop_proie[i+1] >= pop_proie[i+2]:
            max_rk4 = np.append(max_rk4,pop_proie[i+1])
    pop_proie = x1
    for i in range(pop_proie.size-2):
        if pop_proie[i+1] >= pop_proie[i] and pop_proie[i+1] >= pop_proie[i+2]:
            max_simplectique = np.append(max_simplectique,pop_proie[i+1])
            
    maximums = np.array([max_euler,max_rk2,max_rk4,max_simplectique]) 
    
    return maximums


# In[18]:


def maximum_pred(start, end, step, v_ini, derivee,ordre,derivee_simplectique):
    t, v = Euler(start, end, step, v_ini, derivee)
    t, Drk2 = rk2(start, end, step, v_ini, derivee, ordre)
    t, Drk4 = rk4(start, end, step, v_ini, derivee, ordre)
    t,x1,y1 = Euler_simplectique(start, end, step, v_ini,derivee_simplectique)
    
    max_euler = ()
    max_rk2 = ()
    max_rk4 = ()
    max_simplectique = ()
    
    pop_pred = v[1]
    for i in range(pop_pred.size-2):
        if pop_pred[i+1] >= pop_pred[i] and pop_pred[i+1] >= pop_pred[i+2]:
            max_euler = np.append(max_euler,pop_pred[i+1])
    pop_pred = Drk2[1,:]
    for i in range(pop_pred.size-2):
        if pop_pred[i+1] >= pop_pred[i] and pop_pred[i+1] >= pop_pred[i+2]:
            max_rk2 = np.append(max_rk2,pop_pred[i+1])
    pop_pred = Drk4[1,:]
    for i in range(pop_pred.size-2):
        if pop_pred[i+1] >= pop_pred[i] and pop_pred[i+1] >= pop_pred[i+2]:
            max_rk4 = np.append(max_rk4,pop_pred[i+1])
    pop_pred = y1
    for i in range(pop_pred.size-2):
        if pop_pred[i+1] >= pop_pred[i] and pop_pred[i+1] >= pop_pred[i+2]:
            max_simplectique = np.append(max_simplectique,pop_pred[i+1])
            
    maximums = np.array([max_euler,max_rk2,max_rk4,max_simplectique]) 
    
    return maximums


# In[19]:


def max_plot(maxs):
    
    pts = np.arange(1,maxs[0].size+1,1)
    plt.plot(pts,maxs[0],label='Euler') #Euler
    plt.scatter(pts,maxs[0])

    pts = np.arange(1,maxs[1].size+1,1)
    plt.plot(pts,maxs[1],label='rk2') #rk2
    plt.scatter(pts,maxs[1])

    pts = np.arange(1,maxs[2].size+1,1)
    plt.plot(pts,maxs[2],label='rk4') #rk4
    plt.scatter(pts,maxs[2])

    pts = np.arange(1,maxs[3].size+1,1)
    plt.plot(pts,maxs[3],label='simplectique') #simplectique
    plt.scatter(pts,maxs[3])
    


# In[22]:


te = 750
maxs_pred = maximum_pred(ti,te,step, v_ini, derivee_v,ordre,Lokta) # euler, rk2, rk4, simplectique
max_plot(maxs_pred)
plt.legend()
plt.ylabel('maximums des populations des prédateurs ')
plt.xlabel("'pic' de population, numéro du maximum pour un temps t [années] ")
plt.title('évolution des maximums des populations des prédateurs pour différentes méthodes de résolution')
plt.show()


# In[ ]:




