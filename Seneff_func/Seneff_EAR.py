def Seneff (input,sr,largo):

  import numpy as np
  import math
  import copy

  #ETAPA DE RECTIFICACIÓN 

  #Constantes del modelo de Seneff

  def rectificacion_seneff(input):
    
    A=10
    B=65
    G=2.35

    #Función de rectificación Seneff

    y=input  
    y_rect=np.piecewise(y, [y<=0, y>0],[lambda y:np.exp(A*B*y), lambda y:1+A*np.arctan(B*y)])

    #Ganancia usada en el modelo G=2.35

    y_rect_amp=y_rect*G  

    return y_rect_amp

  #MODELO DE MEMBRANA 

  def membrana_Seneff(y_rect_amp,sr):

    #constantes del modelo

    tau1=15*0.001
    tau2=120*0.001
    ua=((1/tau1)-(1/tau2))
    ub=(1/tau2)
    start_Cn = 0
    
    s=y_rect_amp
    largo=len(s)
    c= np.zeros(largo)
    c[0]=0
    T=1/sr
    k=ua+ub
    
    """
    #método diferencial

    #for i in range(1,largo):
      #if c[i-1]<s[i]:
        #c[i]=((T/(1+k*T))*(ua*s[i]+(c[i-1]/T)))
      #else:
        #c[i]=(1/(1+(ub*T)))*c[i-1]

     """ 

    #Método de Euler

    for i in range(1,largo):
      if c[i-1]<s[i-1]:
        c[i]=c[i-1]+(T*(-k*c[i-1]+ua*s[i-1]))
      else:
        c[i]=c[i-1]+(T*(-ub*c[i-1]))
      
    return c,ua

  #ETAPA DE FILTRO PASA BAJO (pérdida de sincronía)

  def LPF(c,sr,largo):

    #Constantes del modelo

    tauLP=0.04*0.001
    nLP=sr*tauLP
    alpha=math.exp(-1/nLP)
    
    import copy

    y_fil= np.zeros(largo)
    y_fil[0]=0

    for j in range (4):
      for i in range (1,largo):
        y_fil[i]=(1-alpha)*c[i]+y_fil[i-1]*alpha
      c = copy.deepcopy(y_fil) 

    return y_fil

    # ETAPA DE COMPRESIÓN

  def AGC(y_fil,sr,largo):

    #Constantes del modelo 
    
    tauAGC=3*0.001
    nAGC=sr*tauAGC
    alphaAGC=math.exp(-1/nAGC)
    KAGC=0.002
    
    #Output de la etapa de adaptación filtrados con el filtro de primer orden 

    y_fil_AGC= np.zeros(largo)
    y_fil_AGC[0]=0.23071276
    
    for i in range (1,largo):

      y_fil_AGC[i]=(1-alphaAGC)*y_fil[i]+y_fil_AGC[i-1]*alphaAGC

    #Modelo AGC 

    y_final=np.zeros(largo)
    div=1+KAGC*y_fil_AGC
    y_final=y_fil/div
    
    return y_final


  #Se invoca la función rectificación 

  y_rect_amp=rectificacion_seneff(input)

  #Se invoca la función membrana (método de Euler)

  c,ua=membrana_Seneff(y_rect_amp,sr)
  y_rect_amp=y_rect_amp.reshape(largo,)
    
  flow=np.subtract(y_rect_amp,c)

  flow_ua=np.multiply(ua,flow)
  flow_final=np.clip(flow_ua,0,None)

  #Se invoca la función LPF

  y_fil=LPF(flow_final,sr,largo)

  #Se invoca a función de compresión (AGC)

  y_final=AGC(y_fil,sr,largo)

  return y_final

