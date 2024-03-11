# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 14:14:45 2022
@author: sirik-baranov

Програма для дослідження шліфів двофазних іррегулярних евтектик типу китайське
письмо методом січних ліній.

- Сканування чорно-білого эвтектичного мотиву на підготовленому знимку системою 
паралельних січних ліній під заданим кутом розвороту.

- Створення графічної візуалізації системи січний ліній на знімку.

- Створення масиву довжин чорних, білих або бінарних (біле поле - чорне поле),
доменів, що були відсічені усіма січними лініями.

- Статистична обробка отриманих довжин доменів. Побудування гістограм, 
визначення характерних розмірів домену (середній, найчастіший, максимальний),
розрахунок площин фаз та періметру границі між ними. 

- Створення усередненої профілограми. для вивчення просторового розподілу 
розмірів елементів євтектичного мотиву. Створюються масиви доменів окремо для 
кожної січної лінії. Далі кожна січна лінія розбивається на однакові інтервали 
(канали) таким чином, щоб вони збігалися з відповідними каналами інших ліній.
Довжина домену, що потрапляє у канал, сумується з доменами, що потрапляють у 
відповідні канали інших ліній, та відбувається усереднення.
 



"""
import sys
import os
import matplotlib.image as mpimg
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import math
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d
from scipy import interpolate 

def Hyst (df,aa,hyst_edge_min,hyst_edge_max,Nbins):
    
    'Статистична обробка довжин доменів'
      
    def Half (A,B):
        B1=pd.Series(gaussian_filter1d(B, sigma=1))
        x_points = np.array(A)
        y_points = np.array(B1)
        tck = interpolate.splrep(x_points, y_points)
        j=min(A)
        a=[]
        b=[]
        p=0.5
        while j<max(A):
            a+=[j]
            b+=[ interpolate.splev(j, tck)]
            j=j+p 
        m=b.index(max(b)) 
        mx=a[m]  
        b_=list(filter(lambda y: y  > max(b)/2, b))
        half=len(b_)*p
        return (mx,half)
    'нахождение площади фазы'
    s=sum(df['L'])*aa
    'нахождение L среднего'
    mid=round(sum(df['L'])/len(df['L']),2)
    'нахождение периметра границы'
    P=len(df['L'])*2*aa
    'нахождение сегмента максимальной длины'
    max_len = round(max(df['L']), 2)
   
    delta= hyst_edge_max/Nbins 
    hh=df['L'].tolist()
    h= np.array(hh)
    hyst=np.histogram(h, bins=Nbins, range=(hyst_edge_min,
                                            hyst_edge_max) ) 
    df=pd.DataFrame(hyst)
    l=len(df.columns)-1
    df1=df.drop(df.columns[l], axis=1)
    a=df1.iloc[0]*100/len(h)#sum(df1.iloc[0])
    b=df1.iloc[1]
    df2=pd.DataFrame({'L':b,'Prob':a})
    df2['Prob2']=((b+delta/2)*aa*df1.iloc[0]*100)/s
    mx=Half(df2['L'],df2['Prob'])[0]
    half=Half(df2['L'],df2['Prob'])[1]     
    mx2=Half(df2['L'],df2['Prob2'])[0]
    half2=Half(df2['L'],df2['Prob2'])[1] 
    return(df2,mx,half,mx2,half2,s,mid,P,max_len )

def angle (df,angle,b,K):
    '''
    Побудування масиву пікселів вздовж січної лінії c тангенсом кута нахилу k
    та висотою над початком координат b
    '''
    h=len(df)
    l=len(df.columns)
    n=0
    I=[]
    while n<h:
        I+=[h-n-1]
        n=n+1
    df=df.reindex(I) 
    h=len(df)-1
    a=math.radians(angle)
    k=math.tan(a)
    Y = lambda x: k*x+b
    X = lambda y:(y-b)/k
    
    #*******************************************
    #1
    #*******************************************
    
    if angle < 0:
        if b>=h:
            x0=int(X(h))
            y0=h
            if Y(l)>=0:  
                lteor=math.hypot(h-Y(l),l-X(h))
            else:
                lteor=math.hypot(h,X(0)-X(h))
        else:
            x0=0
            y0=int(b)
            if Y(l)>=0: 
                lteor=math.hypot(b-Y(l),l)
            else:
                lteor=math.hypot(b,X(0))
                
        x=x0
        y=y0
        m=int(Y(x))
        T=[]  
        while y >= 0 :
           try: 
            if y == m:
                T+=[df.iloc[y,x]]
                x=x+1
                m=int(Y(x))
            else:
                y=y-1
                T+=[df.iloc[y,x]]
                m=int(Y(x))
           except IndexError : break 
       
    #*******************************************
    #2
    #******************************************* 
      
    elif angle >0:
        if b<0:
             x0=int(X(0))
             y0=0
             if Y(l)<=h:
                 lteor=math.hypot(Y(l),l-X(0))
             else:
                 lteor=math.hypot(h,X(h)-X(0))
                
        else:
            x0=0
            y0=int(b)
            if Y(l)<=h:
                 lteor=math.hypot(Y(l)-b,l)
            else:       
                 lteor=math.hypot(h-b,X(h))
        x=x0
        y=y0
        m=int(Y(x))
        T=[]      
        while y <= h :
           try: 
            if y == m:
                T+=[df.iloc[y,x]]
                x=x+1
                m=int(Y(x))
            else:
                y=y+1
                T+=[df.iloc[y,x]]
                m=int(Y(x))
           except IndexError : break 
       
    #*******************************************
    #3
    #*******************************************  
    
    else:
            x0=0
            y0=int(b)
            lteor=l
            T=df.iloc[y0].tolist()   
           
    K2=lteor/len(T)  
    
        
    df=pd.DataFrame({'I':T})
    df['x']=df.index
    df['x']=df['x']*K
    df['x']=df['x']*K2
    df=df.round(2)
    L=lteor*K
    Lsec=len(T)*K*K2
    print('^^^^^^^^^^^')
    print('l: ',l)
    print('h: ',h)
    print('lteor, pix: ',round(lteor,2))
    print('lteor, mkm: ',round(L,2))
    print('Lsec,mkm: ',round(Lsec,2))
    print('max x, mkm: ',max(df['x']))
    print('K2: ',round(K2,3))
    print('K: ',K)
    print('^^^^^^^^^^^')    
    return(df,K2,round(Lsec,2))
def get_leg (name):
    
    a=os.path.splitext(name)[0]
    pos=a.find('K')
    nam=a[:pos-1]
    return (nam) 
def Drow (name,title,labelx,labely,data,picname,clip,ΔA,variant,method,limY,limY0):
    
    'Побудування усередненої профілограми'
    
    h=max(data['channel'])
    if method==1:
        l=max(data['AverⅠ'])
    else: 
        l=max(data['AverⅡ'])
    if variant == 1:
        fig = plt.figure(figsize=(5,2))
    else:
        fig = plt.figure(figsize=(5,4))
    ax=fig.add_subplot()
    plt.rcParams.update({'font.size': 12})
    if order > 400:
        plt.ticklabel_format(style='sci', axis='x', 
                         scilimits=(0,0),useMathText=True) 
    else:
        plt.ticklabel_format(style='plain',  axis='x',
                             scilimits=(0,0)) 
    if variant == 1:    
        ax.set_ylim(0, h*1.4)
    else:
        ax.set_ylim(limY0, limY)
        
    ax.set_xlabel(labelx,fontsize=16) 
    ax.set_ylabel(labely,fontsize=16)
    ax.grid(True)
    jjj=['-','dashed',':']
    x=data['channel']
    if method==1:    
        y=data['AverⅠ']
    else:y=data['AverⅡ']     
    m='□'+str(int(max(data['channel'])))+'mkm ' 
    ax.set_title(title,fontsize=12) 
    ax.step(x,y,
            linewidth = 1.1,linestyle =jjj[0],
            color=(0, 0, 0, 1),label=get_leg (name),
            where='pre')
    (markers, stemlines, baseline)=ax.stem(x,y,
    basefmt='none', markerfmt = 'none')       
    plt.setp(stemlines, linestyle="-", 
    color="black",linewidth=0.5 )
     
    ax.legend()
    plt.savefig('профили/'+get_leg(name)+m+ΔA+picname,
                bbox_inches='tight',dpi=300) 
    data.to_excel('профили/'+get_leg(name)+m+ΔA+picname+'.xlsx')
    print(clip) 
    plt.show()
    
def Drow1 (name,title,labelx,labely,data,picname,clip,ΔA,variant):
    
    'Побудування загального профілю по всіх січних лініях'
    
    h=max(data['L'])
    l=max(data['x'])
    if variant == 1:
        fig = plt.figure(figsize=(5,2))
    else:
        fig = plt.figure(figsize=(5,4))
    ax=fig.add_subplot()
    plt.rcParams.update({'font.size': 12})
    if order > 50:
        plt.ticklabel_format(style='sci', axis='x', 
                         scilimits=(0,0),useMathText=True) 
    else:
        plt.ticklabel_format(style='plain',  axis='x',
                             scilimits=(0,0)) 
    if variant == 1:    
        ax.set_ylim(0, h*1.4)
    else:
        ax.set_ylim(0, 200)
    ax.set_xlabel(labelx,fontsize=16) 
    ax.set_ylabel(labely,fontsize=16)
    ax.grid(True)
    level=math.fsum(data['L'])/len(data['L'])
    strukt=pd.DataFrame({'x':[0,l],'y':[level,level]})
    jjj=['-','dashed',':']
    x=data['x']
    y=data['L']
    m='□'+str(int(max(data['x'])))+'mkm ' 
    ax.set_title(title+', '+get_leg (name),fontsize=12) 
    #ax.step(x,y,
    #        linewidth = 1.1,linestyle =jjj[0],
    #        color=(0, 0, 0, 1),
    #        where='mid')
    (markers, stemlines, baseline)=ax.stem(x,y,
    basefmt='none', markerfmt = 'none')       
    plt.setp(stemlines, linestyle="-", 
    color="black",linewidth=0.8 )
    if ph=='d':
        ax.plot(strukt['x'],strukt['y'],color='b',linestyle="--",
            linewidth=1,label='main eutectic distance level',)
    else:
        ax.plot(strukt['x'],strukt['y'],color='b',linestyle="--",
            linewidth=1,label='middle hord level',)    
    ax.legend()
    plt.savefig('профили/'+get_leg(name)+m+ΔA+picname,
                bbox_inches='tight',dpi=300) 
    data.to_excel('профили/'+get_leg(name)+m+ΔA+picname+'.xlsx')
    print(clip) 
    plt.show()
    return(m)
    
def Drow2 (name,i,A,b,Δb,a,k,signΔ,t,ΔA): 
    
    'Створення графічної візуалізації системи січний ліній на знімку'
    
    def irsis(a):
       b=50
       c=a/b
       C=math.modf(c)[0]
       if C==0:
           A=int(a)
       else:
           i=0 
           while C!=0:
               a=a+1
               c=a/b
               C=math.modf(c)[0]
               i=i+1
               if i>1:
                   break
           A=int(a)
       return(A)    
    name_rez='__'+name 
    fig = plt.figure()
    ax=fig.add_subplot()
    img = mpimg.imread('sci_pic_in/'+name) 
    lum_img = img[:, :, 0]
    df1 = pd.DataFrame(lum_img)
    ax.imshow(df1, cmap='gray' )
    ax.axis('off')
    l=len(df1.columns)
    h=len(df1)
    n=0
    ax.set_ylim(0, h)
    while n<i:
        
        x = np.linspace(0, l, l)
        y = lambda x: k*x+b+Δb*n*signΔ
        ax.set_ylim(0, h)
        if n==0:
            co='r'
            al=1
        elif n>0 and order < 3:
            co='orangered'
            al=1
        else:
            co='r'
            al=0.1
        ax.plot(x, y(x),color=co,alpha=al)
        n=n+1  
    i=50
    n=1
    L=n*i/K
    while int(L)<int(l*0.25):
        L=L+n*i/K
        n=n+1
    m =  str(irsis(int(L*K)))+'mkm'
    rect = Rectangle((0.65*l, 0.83*h), L ,0.15*L,
                     facecolor='white',
                     edgecolor='black')
    ax.add_patch(rect) 
    plt.text(0.68*l, 0.84*h,m,color='black')
    if t!= '#':
        plt.savefig('скан_шоу/'+t+ΔA+name_rez,
                    bbox_inches='tight',
                    pad_inches = 0)
    plt.show()
    
def Drow3 (title,labelx,labely,data,col,path,clip):
    
    'Побудування гістограм розподілу доменів по їх довжинам'
    
    fig, ax = plt.subplots( )
    plt.rcParams.update({'font.size': 14})
    ax.set_title(title,fontsize=16)
    plt.ticklabel_format(style='plain', axis='y', 
                         scilimits=(0,0))
    #ax.set_yscale("log", base=10)
    ax.set_xlabel(labelx,fontsize=16) 
    ax.set_ylabel(labely,fontsize=16)
    ax.grid(True)
    fff=['^','o','x']
    jjj=['-','dashed',':']        
    x=data['L']
    y=data[col]
    ax.step(x,y,marker = fff[0],markersize = 4,
                linewidth = 1,linestyle =jjj[0],
                label = get_leg (name),color=(0, 0, 0, 1),
                where='mid')   
    (markers, stemlines, baseline)=ax.stem(x,y,
    basefmt='none', markerfmt = 'none')
    plt.setp(stemlines, linestyle="-", 
    color="grey",linewidth=0.5 )
    ax.legend()
    plt.savefig(path,bbox_inches='tight',dpi=300) 
    print(clip) 
    plt.show()    
         
def do_w (xin,K,df,K2):
    
    'Створення масиву білих доменів на січної лінії'
    
    datapix=df['I'].tolist()
    countwhite = 0
    i=0
    datarezwhite = []
    nrez = []
    while i<len(datapix)-1: 
               i=i+1  
               num = float(datapix[i])             
               if num >= 0.8:
                    countwhite += 1 
               elif num < 0.8 and countwhite != 0:
                    datarezwhite+=[round(countwhite*K*K2,2)]
                    nrez+=[round((i - countwhite/2 - 0.5)
                                 *K*K2,2)]
                    countwhite = 0
    df1=pd.DataFrame({'x':nrez,'L':datarezwhite})  
    df1['x']=df1['x']+xin
    xout=max(df1['x'])+1   
    return (df1,xout)

def do_b (xin,K,df,K2):
    
    'Створення масиву чорних доменів на січної лінії'
    
    datapix=df['I'].tolist()
    countblack = 0
    i=0
    datarezblack = []
    nrez = []
    while i<len(datapix)-1: 
               i=i+1  
               num = float(datapix[i])             
               if num < 0.8:
                    countblack += 1 
               elif num >= 0.8 and countblack != 0:
                    datarezblack+=[round(countblack*K*K2,2)]
                    nrez+=[round((i - countblack/2 - 0.5)
                                 *K*K2,2)]
                    countblack = 0
    df1=pd.DataFrame({'x':nrez,'L':datarezblack})  
    df1['x']=df1['x']+xin
    xout=max(df1['x'])+1      
    return (df1,xout)

def do_d (xin,K,df,K2):
    
    'Створення масиву бінарних доменів на січної лінії'
    
    def P(n):
        if n<0.8:
            i=0
        else:
            i=1
        return(i)    
    
    datapix=df['I'].tolist()
    i=0
    istart=0
    while i<len(datapix)-1:
        num=datapix[i]
        if P(num) == 0:
            istart=i
            break
        i=i+1
    
    count = 0
    i=istart
    j=istart+1
    f=0
    datarez = []
    nrez = []
    while j<len(datapix)-1: 
               i=i+1  
               j=j+1
               num = datapix[i]
               num2= datapix[j]
               if P(num) == P(num2) :
                   count += 1
               elif  P(num) != P(num2) and f<1:  
                   count += 1
                   f=f+1           
               elif f>=1:          
                    datarez+=[round(count*K*K2,2)]
                    nrez+=[round((i - count/2 - 0.5)
                                 *K*K2,2)]
                    count = 0
                    f=0        
    df1=pd.DataFrame({'x':nrez,'L':datarez})  
    df1['x']=df1['x']+xin
    try:
        xout=max(df1['x'])+1
    except ValueError:
        xout =xin+1
    return (df1,xout)

def getdf1 (df,A,b,K,Δb,signΔ,lim,ph,N):
    
    'Утворення спільного масиву з масивів даних окремих січних ліній'
    
    def get_channels (df,lsec,N,m):
        n=0
        lchannel=round(lsec/N,0)
        Φ_x=[]
        Φ_L=[]
        while n<len(df):
            i=0
            x=df['x'][n]
            L=df['L'][n]
            while i<N:
                λ=int(i*lchannel)
                if  λ >= x:
                    Φ_x+=[λ]
                    Φ_L+=[L]
                    break
                i=i+1
            n=n+1    
        df = pd.DataFrame({'channel':Φ_x,'L':Φ_L})  
        df=df.sort_values(by=['L'])
        df=df.drop_duplicates(subset=['channel'],keep='last')
        df=df.sort_values(by=['channel']) 
        Ndrops=max(df.index)-len(df)
        
        n=0
        Nzero=0
        Φ_x=[]
        Φ_L=[]
        while n<N:
            λ=int(n*lchannel)
            if λ in df['channel'].tolist() :
                Φ_x+=[λ]
                Φ_L+=[float(df.loc[df['channel'] == λ]['L'])]
                n=n+1
            else: 
                Φ_x+=[λ]
                Φ_L+=[0]
                Nzero=Nzero+1
                n=n+1              
        df = pd.DataFrame({'channel':Φ_x,'L':Φ_L}) 
        df[m]=df['L']
        del(df['L'])
        return df,Ndrops,Nzero
    def find_0(s):
        i=0
        I=0
        while i<len(s):
            if int(s[i])!=0:
                I+=1
            i+=1
        return I    
                
    xin=0
    order=1
    T=angle(df,A,b,K)
    df0=T[0]
    K2=T[1]
    lsec=T[2]
    if ph=='w':
        T=do_w(xin,K,df0,K2)
    elif ph=='b': 
        T=do_b(xin,K,df0,K2)
    elif ph=='d':
        T=do_d(xin,K,df0,K2)
    df1=T[0]
    xin=T[1]
    T=get_channels(df1,lsec,N,0)
    df_ch=T[0]
    Ndrops=[T[1]]
    Nzero=[T[2]]
    NL=[len(df1)]
    Nd=pd.DataFrame({'Ndrops':[Ndrops],'Nzero':[Nzero]})
    Nd.to_excel('Ndrops.xlsx')
    while max(df1['x']) < lim:          
            T=angle(df,A,b+order*Δb*signΔ,K)
            df01=T[0]
            lsec=T[2]
            if ph=='w':
                T=do_w(xin,K,df01,K2)
                T1=do_w(0,K,df01,K2)
            elif ph=='b':
                T=do_b(xin,K,df01,K2)
                T1=do_b(0,K,df01,K2)
            elif ph=='d':
                T=do_d(xin,K,df01,K2)
                T1=do_d(0,K,df01,K2)
            df02=T[0]
            df03=T1[0]
            xin=T[1]
            
            
            T=get_channels(df03,lsec,N,order)
            df_ch_=T[0]
            df_ch[order]=df_ch_[order]
            Ndrops+=[T[1]]
            Nzero+=[T[2]]
            NL+=[len(df03)]
            df1=pd.concat([df1,df02],ignore_index=True)
            order=order+1
            st=str(order)+' '+'##### '+str(A)
            print(st)
    df_zero=pd.DataFrame({'l, mkm':int(lsec),
                          'δl, mkm':int(lsec/N),
                          'M':order,
                          'N':N,
                          'Nm':[int(sum(NL)/len(NL))],
                          'Ndrops':[int(sum(Ndrops)/len(Ndrops))],
                          'Nzero':[int(sum(Nzero)/len(Nzero))],
                          'Method':['ⅠorⅡ']}).T
    df_zero['item']=df_zero.index
    tt=['lenth of secant line','channel width',
         'number of profilograms','number of channels',
         'average number of domains per secant',
         'average number of missing domains per secant',
         'average number of empty measurements per channel',
         'average method']
    
    
    df_zero.index=tt
    df_zero['comment']=df_zero.index
    df_zero.index=df_zero['item']
    del(df_zero['item'])
    df_zero.rename(columns={0: 'value'}, inplace=True)
    
   
   
            
    df_ch1=df_ch.copy() 
    del(df_ch1['channel'])
    i=0
    I=[]
    while i<len(df_ch1):
        I+=[find_0(df_ch1.iloc[i])]
        i=i+1
    df_ch['nonzero'] = I
    df_ch['AverⅠ']=round(df_ch1.sum(axis=1)/df_ch['nonzero'],2)
    df_ch['AverⅡ']=round(df_ch1.sum(axis=1)/len(df_ch1.columns),2)
    D=pd.DataFrame()
    D['channel']=df_ch['channel']
    D['AverⅠ']=df_ch['AverⅠ']
    D['AverⅡ']=df_ch['AverⅡ']
    D.to_excel('ch.xlsx')
           
    return(df1,order,D,df_zero)  
def get_K (name):
    
    a=os.path.splitext(name)[0]
    pos=a.find('K')
    K=float(a[pos+1:])*0.01
    return (K)
     

mylist = os.listdir('sci_pic_in.') 
name=mylist[0]
K=get_K(name) 
'''
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

В І К Н О   К Е Р У В А Н Н Я   П Р О Г Р А М О Ю
введення робочих параметрів

''' 

A=89.9 #цей параметр задає кут нахилу січних ліній над віссю абсцис в градусах

b=-25000
    # задає висоту першої січної лінії над початком координат в мкм

a=2         # задає відстань (крок) між паралельними січними лініями під час  
            # заповнення ними знімка, мкм

signΔ=-1    # задає напрямок заповнення знимку лініями, від першої лінії

lim=40000   # задає максимальну кількість доменів, відсічених лініями
            # від цього параметру також залежить кількість січних ліній
            
ph = 'd'    # визначає, які домени обчислюються          
            # (білі 'w', чорні 'b', або бінарні 'd')  
            
hyst_edge_min = 0   # завдає лівий край гістограми

hyst_edge_max = 130 # завдає правий край гістограми

Nbins=80            # завдає число стовбчиків гістограми 



N=10                   # завдає число каналів профілограми

limY0, limY = (5,40)   # діапазон побудови профілограми по осі y, в мкм


key='##-'   # цей параметр звадає режим роботи програми
# якщо значення key='#' - програма працює в режимі вибору області сканування
# значення key='' - статистична обробка (пбудова гістограм)
# значення key='@' - пбудова усереднених профілограм

# при значеннях key='##-' та key='##+' виконується серія сканувань 
# з автоматичною послідовною  зміною кута нахилу A
'''
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
'''
method=1
ΔA='_r'+str(A)+'°'
aa=a/K
p=math.radians(A)
k=math.tan(p)
if A!=0:
    Δb = aa/math.sin(math.pi/2 - math.radians(A))
else:
    Δb = aa    
img = mpimg.imread('sci_pic_in./'+name) 
lum_img = img[:, :, 0]
df=pd.DataFrame(lum_img)
order = 1
xin=0
if key =='##-':
    
    i=-3
    I=[]
    J=[]
    M=[]
    try:
     while A>=-89:
       p=math.radians(A)
       k=math.tan(p) 
       Δb = aa/math.sin(math.pi/2 - math.radians(A))
       T=getdf1 (df,A,b,K,Δb,signΔ,lim,ph,N)
       df1=T[0]
       order=T[1]
       d=Hyst(df1,a,hyst_edge_min,hyst_edge_max,Nbins)
       I=I+[d[6]]
       M=M+[d[8]]
       J=J+[A]
       df1.to_excel('df1./A '+str(A)+'.xlsx')
       mm='#'
       Drow2(name,order,A,b,Δb,a,k,signΔ,mm,ΔA)
       A=A+i
       
       
       
    except ValueError:
        df=pd.DataFrame({'угол':J,'Lсредн':I,'Lмакс':M})  
        df.to_excel('find_L_min.xlsx')
        df1.to_excel('df1./A '+str(A)+'.xlsx')
        sys.stdout.write('\a')
        sys.stdout.flush()   
        
       
    df=pd.DataFrame({'угол':J,'Lсредн':I,'Lмакс':M})  
    df.to_excel('find_L_min.xlsx')
    sys.stdout.write('\a')
    sys.stdout.flush()   
if key =='##+':
    
    i=-3
    I=[]
    J=[]
    M=[]
    
    while A>=0:
       p=math.radians(A)
       k=math.tan(p) 
       Δb = aa/math.sin(math.pi/2 - math.radians(A))
       T=getdf1 (df,A,b,K,Δb,signΔ,lim,ph,N)
       df1=T[0]
       order=T[1]
       d=Hyst(df1,a,hyst_edge_min,hyst_edge_max,Nbins)
       I=I+[d[6]]
       M=M+[d[8]]
       J=J+[A]
       df1.to_excel('df1./A '+str(A)+'.xlsx')
       A=int(A+i)
       mm='#'
       Drow2(name,order,A,b,Δb,a,k,signΔ,mm,ΔA)
       
    df=pd.DataFrame({'угол':J,'Lсредн':I,'Lмакс':M})  
    df.to_excel('find_L_min.xlsx')
    sys.stdout.write('\a')
    sys.stdout.flush()        
elif key == '#':
    T=getdf1 (df,A,b,K,Δb,signΔ,lim,ph,N)
    df1=T[0]
    order=T[1]
    d=Hyst(df1,a,hyst_edge_min,hyst_edge_max,Nbins)
    data=d[0]
    L=d[6]
    Lmax=d[8]
    mm='#'
    Drow2(name,order,A,b,Δb,a,k,signΔ,mm,ΔA)
    print('середня міжфазна відстань')    
    print('L   ',L, 'mkm')
    print('максимальна міжфазна відстань') 
    print('Lmax',Lmax,'mkm')
elif key=='@':
    T=getdf1 (df,A,b,K,Δb,signΔ,lim,ph,N)
    df1=T[0]
    order=T[1]
    D=T[2]
    df=T[3]
    if method==1:
        df.loc['Method','value']='Ⅰ'
        df.loc['Method','comment']='average method:   Λ = ΣLᵐ / (M - Nzero)'
        
    elif  method==2: 
        df.loc['Method','value']='Ⅱ'
        df.loc['Method','comment']='average method:   Λ = ΣLᵐ / M '
    df.to_excel('AVER_DATA.xlsx')    
    labelx='channel, mkm'
    labely='Λ, mkm'
    if ph=='w':
        title='white hords profile'
        picname='w_profil'
    elif ph=='b':
        title='black hords profile'
        picname='b_profil'
    elif ph=='d':
        title='average profilogram'
        picname='d_profil'        
    clip='@@@@@@@@@@'
    picname1 = picname + '_var1.png'
    picname2 = picname + '_var2.png'
    picname3 = picname + '_AVER.png'
    Drow(name,title,labelx,labely,D,picname3,clip,ΔA,2,method,limY,limY0) 
else:    
    T=getdf1 (df,A,b,K,Δb,signΔ,lim,ph,N)
    df1=T[0]
    order=T[1]
    D=T[2]        
    labelx='x, mkm'
    labely='L, mkm'
    if ph=='w':
        title='white hords profile'
        picname='w_profil'
    elif ph=='b':
        title='black hords profile'
        picname='b_profil'
    elif ph=='d':
        title='eutectic distance profile'
        picname='d_profil'        
    clip='@@@@@@@@@@'
    picname1 = picname + '_var1.png'
    picname2 = picname + '_var2.png'
    picname3 = picname + '_AVER.png'
    #Drow1(name,title,labelx,labely,df1,picname3,clip,ΔA,2)
    mm=Drow1(name,title,labelx,labely,df1,picname1,clip,ΔA,1)
    #Drow1(name,title,labelx,labely,df1,picname2,clip,ΔA,2)
    Drow2(name,order,A,b,Δb,a,k,signΔ,mm,ΔA)    
    if ph=='w':
        title='White phase hord'
        path='гист граф1/'+get_leg (name)+ΔA+mm+'_w.png'
        path_ex='гистограммы_в_excel/'+get_leg (name)+mm+ΔA+'_w.xlsx'
        path_ex2='таблицы_в_excel/'+get_leg (name)+mm+ΔA+'_w.xlsx'
    elif ph=='b':
        title='Black phase hord'
        path='гист граф1/'+get_leg (name)+ΔA+mm+'_b.png'
        path_ex='гистограммы_в_excel/'+get_leg (name)+mm+ΔA+'_b.xlsx'
        path_ex2='таблицы_в_excel/'+get_leg (name)+mm+ΔA+'_b.xlsx'
    elif ph=='d':
        title=''
        path='гист граф1/'+get_leg (name)+ΔA+mm+'_d.png'
        path_ex='гистограммы_в_excel/'+get_leg (name)+mm+ΔA+'_d.xlsx'
        path_ex2='таблицы_в_excel/'+get_leg (name)+mm+ΔA+'_d.xlsx'
    labelx='L, mkm'
    labely='F$_{L}$, %'
    d=Hyst(df1,a,hyst_edge_min,hyst_edge_max,Nbins)
    data=d[0]
    L=d[6]
    col='Prob'
    clip='@.@.@.@.@.'
    Drow3(title,labelx,labely,data,col,path,clip) 
    data.to_excel(path_ex)
    file=get_leg (name)
    Lmid=d[6]
    Lmax=d[8]
    df=pd.DataFrame({'file':[file],'Lср':[Lmid],'Lмакс':[Lmax]})
    df.to_excel(path_ex2)
    if ph=='w':
        title='White phase hord'
        path='гист граф2/'+get_leg (name)+ΔA+mm+'_w_pr.png'
        path_ex='гистограммы_в_excel_2/'+get_leg (name)+mm+ΔA+'_w_pr.xlsx'
    elif ph=='b':
        title='Black phase hord'
        path='гист граф2/'+get_leg (name)+ΔA+mm+'_b_pr.png'
        path_ex='гистограммы_в_excel_2/'+ get_leg (name)+mm+ΔA+'_b_pr.xlsx'
    elif ph=='d':
        title=''
        path='гист граф2/'+get_leg (name)+ΔA+mm+'_d_pr.png'
        path_ex='гистограммы_в_excel_2/'+ get_leg (name)+mm+ΔA+'_d_pr.xlsx'    
    labelx='L, mkm'
    labely='F$_{S}$, %'
    data=Hyst(df1,a,hyst_edge_min,hyst_edge_max,Nbins)[0]
    col='Prob2'
    clip='@.@.@.@.@.'
    Drow3(title,labelx,labely,data,col,path,clip) 
    data.to_excel(path_ex)
    
    
    
    

