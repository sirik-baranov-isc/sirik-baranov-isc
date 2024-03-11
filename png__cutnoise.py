"""
реконструювання двох або трифазної структури на шліфі кристала
по знімку поля зору мікроскопа   

"""
import matplotlib.image as mpimg
import pandas as pd
import matplotlib.pyplot as plt
import os
import winsound

def cutnoise6 (s):  # фільтрація 6-піксельного шуму
    s=s.round(1).tolist()
    [i,j,k,l,a,b,c,n]=[0,1,2,3,4,5,6,7]
    while n<len(s):
        if s[n]==0 and s[i]==0 :
            [s[j],s[k],s[l],s[a],s[b],s[c]]=[0,0,0,0,0,0]
            [i,j,k,l,a,b,c,n]=[i+1,j+1,k+1,l+1,a+1,b+1,c+1,n+1]
        elif s[n]==1 and s[i]==1 : 
            [s[j],s[k],s[l],s[a],s[b],s[c]]=[1,1,1,1,1,1]
            [i,j,k,l,a,b,c,n]=[i+1,j+1,k+1,l+1,a+1,b+1,c+1,n+1]
        elif s[n]==0.1 and s[i]==0.1 : 
            [s[j],s[k],s[l],s[a],s[b],s[c]]=[0.1,0.1,0.1,0.1,0.1,0.1]
            [i,j,k,l,a,b,c,n]=[i+1,j+1,k+1,l+1,a+1,b+1,c+1,n+1]
        else :
           [i,j,k,l,a,b,c,n]=[i+1,j+1,k+1,l+1,a+1,b+1,c+1,n+1]
    s=pd.Series(s)
    print('6')
    return s

def cutnoise5 (s):   # фільтрація 5-піксельного шуму
    s=s.round(1).tolist()
    [i,j,k,l,a,b,n]=[0,1,2,3,4,5,6]
    while n<len(s):
        if s[n]==0 and s[i]==0 :
            [s[j],s[k],s[l],s[a],s[b]]=[0,0,0,0,0]
            [i,j,k,l,a,b,n]=[i+1,j+1,k+1,l+1,a+1,b+1,n+1]    
        elif s[n]==1 and s[i]==1 : 
            [s[j],s[k],s[l],s[a],s[b]]=[1,1,1,1,1]
            [i,j,k,l,a,b,n]=[i+1,j+1,k+1,l+1,a+1,b+1,n+1]
        elif s[n]==0.1 and s[i]==0.1 : 
            [s[j],s[k],s[l],s[a],s[b]]=[0.1,0.1,0.1,0.1,0.1]
            [i,j,k,l,a,b,n]=[i+1,j+1,k+1,l+1,a+1,b+1,n+1]   
        else :
           [i,j,k,l,a,b,n]=[i+1,j+1,k+1,l+1,a+1,b+1,n+1]
    s=pd.Series(s)
    print('5')
    return s        
def cutnoise4 (s):   # фільтрація 4-піксельного шуму
    s=s.round(1).tolist()
    [i,j,k,l,a,n]=[0,1,2,3,4,5]
    while n<len(s):
        if s[n]==0 and s[i]==0 :
            [s[j],s[k],s[l],s[a]]=[0,0,0,0]
            [i,j,k,l,a,n]=[i+1,j+1,k+1,l+1,a+1,n+1]       
        elif s[n]==1 and s[i]==1 : 
            [s[j],s[k],s[l],s[a]]=[1,1,1,1]
            [i,j,k,l,a,n]=[i+1,j+1,k+1,l+1,a+1,n+1]
        elif s[n]==0.1 and s[i]==0.1 : 
            [s[j],s[k],s[l],s[a]]=[0.1,0.1,0.1,0.1]
            [i,j,k,l,a,n]=[i+1,j+1,k+1,l+1,a+1,n+1]     
        else :
           [i,j,k,l,a,n]=[i+1,j+1,k+1,l+1,a+1,n+1]
    s=pd.Series(s)
    print('4')
    return s               

def cutnoise3 (s):   # фільтрація 3-піксельного шуму
    s=s.round(1).tolist()
    [i,j,k,l,n]=[0,1,2,3,4]
    while n<len(s):
        if s[n]==0 and s[i]==0 :
            [s[j],s[k],s[l]]=[0,0,0]
            [i,j,k,l,n]=[i+1,j+1,k+1,l+1,n+1]      
        elif s[n]==1 and s[i]==1 : 
            [s[j],s[k],s[l]]=[1,1,1]
            [i,j,k,l,n]=[i+1,j+1,k+1,l+1,n+1]
        elif s[n]==0.1 and s[i]==0.1 : 
            [s[j],s[k],s[l]]=[0.1,0.1,0.1]
            [i,j,k,l,n]=[i+1,j+1,k+1,l+1,n+1]      
        else :
            [i,j,k,l,n]=[i+1,j+1,k+1,l+1,n+1]
    s=pd.Series(s)
    print('3')
    return s

def cutnoise2 (s):   # фільтрація 2-піксельного шуму
    s=s.round(1).tolist()
    [i,j,k,n]=[0,1,2,3]
    while n<len(s):
        if s[n]==0 and s[i]==0 :
            [s[j],s[k]]=[0,0]
            [i,j,k,n]=[i+1,j+1,k+1,n+1]  
        elif s[n]==1 and s[i]==1 : 
            [s[j],s[k]]=[1,1]
            [i,j,k,n]=[i+1,j+1,k+1,n+1]
        elif s[n]==0.1 and s[i]==0.1 : 
            [s[j],s[k]]=[0.1,0.1]
            [i,j,k,n]=[i+1,j+1,k+1,n+1] 
        else :
           [i,j,k,n]=[i+1,j+1,k+1,n+1]
    s=pd.Series(s)
    print('2')    
    return s

def cutnoise1 (s):   # фільтрація 1-піксельного шуму
    s=s.round(1).tolist()
    [i,j,n]=[0,1,2]
    while n<len(s):
        if s[n]==0 and s[i]==0 :
            s[j]=0
            [i,j,n]=[i+1,j+1,n+1]   
        elif s[n]==1 and s[i]==1 : 
            [s[j]]=[1]
            [i,j,n]=[i+1,j+1,n+1]
        elif s[n]==0.1 and s[i]==0.1 : 
            [s[j]]=[0.1]
            [i,j,n]=[i+1,j+1,n+1]   
        else :
           [i,j,n]=[i+1,j+1,n+1]
    s=pd.Series(s)
    print('1')    
    return s  
                      
def do (name,RGB1,RGB2):  # основний алгоритм
    
    img = mpimg.imread('sci_pic_rez_bin/'+name) 
    
    img = img[:-250,:,:]
    
    img[:,:,0][(img[:,:,0]<RGB1[0])] = 0 
    img[:,:,1][(img[:,:,1]<RGB1[1])] = 0 
    img[:,:,2][(img[:,:,2]<RGB1[2])] = 0 
    
    img[:,:,0][(img[:,:,0]>RGB2[0])] = 1 
    img[:,:,1][(img[:,:,1]>RGB2[1])] = 1 
    img[:,:,2][(img[:,:,2]>RGB2[2])] = 1  
    
    img[:,:,0][(img[:,:,0]<=RGB2[0])&(img[:,:,0]>=RGB1[0])] = 0 
    img[:,:,1][(img[:,:,1]<=RGB2[0])&(img[:,:,1]>=RGB1[0])] = 0.1 
    img[:,:,2][(img[:,:,2]<=RGB2[0])&(img[:,:,2]>=RGB1[0])] = 1 
    
    df=pd.DataFrame(img[:,:,1])
    
    # за допомогою рещітки можна вилучати зайві функції фільтрації шуму, 
    # з ціллю отримання найкращого результату
    
    #df=df.apply(cutnoise6,axis=1)
    #df=df.apply(cutnoise6,axis=0)
    
    #df=df.apply(cutnoise5,axis=1)
    #df=df.apply(cutnoise5,axis=0)
    
    #df=df.apply(cutnoise4,axis=1)
    #df=df.apply(cutnoise4,axis=0)
    
    #df=df.apply(cutnoise3,axis=1)
    #df=df.apply(cutnoise3,axis=0)
    
    df=df.apply(cutnoise2,axis=1)
    df=df.apply(cutnoise2,axis=0)
    
    df=df.apply(cutnoise1,axis=1)
    df=df.apply(cutnoise1,axis=0)
    
    #df=df.apply(cutnoise6,axis=0)
    #df=df.apply(cutnoise5,axis=0)
    #df=df.apply(cutnoise4,axis=0)
    #df=df.apply(cutnoise3,axis=1)
    
    #df=df.apply(cutnoise2,axis=1)
    #df=df.apply(cutnoise2,axis=0)
    
    #df=df.apply(cutnoise1,axis=1)
    #df=df.apply(cutnoise1,axis=0)
    
    #df=df.apply(cutnoise6,axis=1)
    
    img[:,:,1]=df.to_numpy()
    img[:,:,0][(img[:,:,1]==0)] = 0
    img[:,:,2][(img[:,:,1]==0)] = 0
    
    img[:,:,0][(img[:,:,1]==1)] = 1
    img[:,:,2][(img[:,:,1]==1)] = 1
    
    img[:,:,0][(img[:,:,1]==0.1)] = 0
    img[:,:,2][(img[:,:,1]==0.1)] = 1
    
    
    
    plt.imshow(img )
    plt.axis('off')
    plt.savefig('sci_pic_rez/'+name,bbox_inches='tight',pad_inches = 0,dpi=600) 
    
mylist = os.listdir('sci_pic_rez_bin.') 

s1=pd.Series((1,1,1))
s2=s1.copy()

s1=s1*0.49 
s2=s2*0.49
# якщо параметр s1 = s2, зображення буде оброблене як двофазна структура

RGB1=s1.tolist()
RGB2=s2.tolist()


i=0
while i<len(mylist):
   
   do(mylist[i],RGB1,RGB2)
   i=i+1
winsound.Beep(1000,2000)  