"""
вирівнювання яскравості фону зображення png

"""

import matplotlib.image as mpimg
import pandas
from scipy.ndimage.filters import gaussian_filter1d
import matplotlib.pyplot as plt
import os


def gauss (serie):
    ysmoothed = gaussian_filter1d(serie, sigma=1000)
    s = pandas.Series(ysmoothed)
    return s 


def Drow (df1,df2):
    '''
    за допомогою цієї функції можлива побудова графіків вихідного, згладженого 
    та скоригованого розподілів яскравості пікселів у рядку зображення
    '''
    h=int(len(df1)/2)
    fig = plt.figure()
    ax=fig.add_subplot()
    ax.set_ylim(0, 1.0)
    ax.set_xlabel('distance, pix',fontsize=16) 
    ax.set_ylabel('I, a.u.',fontsize=16)
    x=df1.columns.tolist()
   
    y1=df1.iloc[h]
    y2=df2.iloc[h]
    y3=(y1-y2)+df1.iloc[h,0]/2+0.15
    ax.plot(x,y1)
    ax.plot(x,y2)
    ax.plot(x,y3)


def do (name):  
    '''
   основний алгоритм програми
    '''
    name_rez='rez_'+name+'.png' 
    img = mpimg.imread('sci_pic_input/'+name) 
    lum_img = img[:, :, 0]
    df = pandas.DataFrame(lum_img)
    h=int(len(df.index)*0.88) 
    h0=len(df.index)
    df_=df.iloc[h:h0]
    df=df.iloc[0:h]
    df1=df.apply(gauss,axis=0)
    df=df+df.iloc[0,0]*0.5-df1  
    df=pandas.concat([df, df_])
    img_rez = df.to_numpy()
    
    plt.axis('off')
    plt.imshow(img_rez, cmap='gray' )
    plt.savefig('sci_pic_rez_bin/'+name_rez,
                bbox_inches='tight',pad_inches = 0,dpi=600) 
    plt.show()
    
key=1 # якщо параметр key=0 - відбувається побудова графіків розподілів,
      # інше значення - режим обробки зображень

if key!=0:
    mylist = os.listdir('sci_pic_input.') 
    for filename in mylist:
       do(filename)
else:
    name=os.listdir('sci_pic_input.') [0]
    img = mpimg.imread('sci_pic_input/'+name) 
    lum_img = img[:, :, 0]
    df = pandas.DataFrame(lum_img)
    df1=df.apply(gauss,axis=0)
    Drow(df,df1)
       
  