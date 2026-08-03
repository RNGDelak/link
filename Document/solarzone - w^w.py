import pygame
from fractions import *
# from LBMS2 import *
# from math import gcd

def ubi(a,n):
    if a==[1]:return list(range(n))
    j=-a[::-1].index(a[-1]-1)-1;return a[:j]+a[j:-1]*n

# def ubi(a,n):
#   #print(a)
#   if a==[1]:return[0,n+1]
#   def p(n):
#     if a[n]==0:return 0
#     r=n
#     while a[r]>=a[n]:r-=1
#     return r
#   r=p(-1)
#   l=len(a)
#   if a[r]==a[-1]-1:return a[:r]+a[r:-1]*n
#   t=[a[i]-a[p(i)]for i in range(l)]
#   q=[]
#   for i in range(l):
#     if t[i]<t[-1]:q.append(i)
#   k=[l-1]
#   while k[-1]!=0:k.append(p(k[-1]))
#   r=max(set(q)&set(k))
#   d=a[-1]-a[r]-1
#   o=a[:r]
#   for i in range(n):o+=[j+d*i for j in a[r:-1]]
#   return o

def getString(x):
    if x<0:return 'negative'
    if x==0:return'*'
    else:
        a=''
        while x!=1:
            if x.numerator%4==1:
                a='0'+a
                x=Fraction((x.numerator+1)//2,x.denominator//2)
            else:
                a='1'+a
                x=Fraction((x.numerator-1)//2,x.denominator//2)
        return a

def strBMS(s):
    if s=='*':return[]
    x=[]
    y=limit#[[],[[[],[]]]]
    d=[]
    for i in s:
        t=0
        if i=='0':y=d.copy()
        else:x=d.copy()
        if y[-1]==0:return'Error'
        while list(ubi(y,t))<=x:t+=1
        d=list(ubi(y,t))
    return d
def strBMS2(s):
    if s=='*':return[]
    x=[]
    y=limit#[[],[[[],[]]]]
    d=[]
    for i in s:
        t=0
        if i=='0':y=d.copy()
        else:x=d.copy()
        if not y:return d[:-1]
        if y[-1]==0:return d[:-1]
        while list(ubi(y,t))<=x:t+=1
        d=list(ubi(y,t))
    return d

def getOrdinal(x):return strBMS(getString(x/2+Fraction(1,2))[1:])
def getOrdinal2(x):return strBMS2(getString(x/2+Fraction(1,2))[1:])

def gz(x):return Fraction(2**(x//64)*e[x%64],1024) if x>0 else 1

pygame.init()

width=pygame.display.Info().current_w
#width=int('1'+'0'*(len(bin(width))-3),2)
height=pygame.display.Info().current_h
canvas = pygame.display.set_mode((width,height))
pygame.display.set_caption("Ordinal number line - ε_0")
limit = [0,1,2]
exit = False

e=[1024,1035,1046,1058,1069,1081,1093,1105,1117,1129,1141,1154,1166,1179,1192,1205,1218,1231,1244,1258,1272,1286,1300,1314,1328,1342,1357,1372,1387,1402,1417,1433,1448,1464,1480,1496,1512,1529,1545,1562,1579,1596,1614,1631,1649,1667,1685,1704,1722,1741,1760,1779,1798,1818,1838,1858,1878,1898,1919,1940,1961,1983,2004,2026]

pos=width//2
z=0
d=0

def showtext(canvas,s,m,x,y):
    font = pygame.font.Font(r'C:\Windows\Fonts\Arial.ttf', m)
    text=font.render(s,True,(255,255,255),(0,0,0))
    canvas.blit(text,(x,y))

def toString(s): # from waffle
 if not s:return'0'
 o,c,r,l="",1,0,None
 for i in range(0,len(s)+1):
  if s[0]==(s[0 if i+2>len(s)else i+1]):
   S=s[r+1:i+1];t="1"if i==r else"ω"if i-1==r else("ω^%s"if sum(v==s[r+1]for v in S)<2 or toString(S).isdigit()else"ω^(%s)")%toString(S)
   if t==l and i<len(s):c+=1
   else:
    if l:o+=" + "+(l if c<2 else str(c)if l=="1"else l+(""if l=="ω"else"*")+str(c))
    l=t;c=1
   r=i+1
 return o[3:]
# def toString(s):
#     if not s:return 'Zero'
#     return str(s)[1:-1].replace(' ','')

C=pygame.time.Clock()
while not exit:
    canvas.fill((0,0,0))
    pygame.draw.line(canvas,(0,0,255),(width//2,0), (width//2,height-100))
    t=z
    for event in pygame.event.get():
        #print(event)
        if event.type == pygame.QUIT:
            exit = True
        if event.type == pygame.MOUSEBUTTONDOWN:d=1
        if event.type == pygame.MOUSEBUTTONUP:d=0
        if d and event.type==pygame.MOUSEMOTION:pos+=event.rel[0];z-=event.rel[1]
        if event.type==pygame.MOUSEWHEEL:z+=event.y*4
        if event.type==pygame.KEYDOWN:
            if event.unicode=='q':exit=True
        if z<0:z=0
    zoom=gz(z)
    ozoom=gz(t)
    pos=(pos*zoom)//ozoom
    rpos=Fraction(pos,width)
    left=Fraction(pos-width//2,zoom*width)
    center=Fraction(pos,zoom*width)
    right=Fraction(pos+width//2,zoom*width)
    ep=Fraction(1,2**(len(bin(int(zoom*width)))-2))
    epx=Fraction(1,zoom*width)
    q=[None for i in range(2*(width//2))]
    for i in range(0,len(bin(int(zoom*width)))-1):
        ep1=Fraction(1,2**i)
        for j in range(left//ep1,right//ep1+1):
            k=(j*ep1)//epx-left//epx
            if 0<=k<len(q):
                if q[k]==None:q[k]=ep1*j
    #d0=q[0]-right
    #q=[ep*((i-d0)//ep) for i in q]
    for i in range(len(q)):
        if not 0<=q[i]<=1:continue
        v=getOrdinal(q[i])
        if v=='Error':continue
        t=max(len(getString(q[i]/2+Fraction(1,2)))-1,0)
        pygame.draw.line(canvas,(255,255,255),(i,height-100),(i,max(int(height-100-Fraction(height,2**t)*zoom),0)))
    if q[width//2]<0:g='0'
    elif q[width//2]>=1:g='ε_0'
    else:g=toString(getOrdinal2(q[width//2]))
    showtext(canvas,g,15,10,height-95)
    showtext(canvas,f'{int(zoom):,}'+'x',12,10,height-75)
    showtext(canvas,str(f'{C.get_fps():.1f} fps'),12,200,height-75)
    showtext(canvas,str(q[width//2]),12,10,height-60)
    showtext(canvas,str(float(q[width//2])),12,200,height-60)
    showtext(canvas,str(center),12,10,height-45)
    showtext(canvas,str(float(center)),12,200,height-45)
    showtext(canvas,'diff. '+str('{:.1f}'.format(1000000*float(zoom*(q[width//2]-center))))+' μpx',12,10,height-30)
    C.tick(10**100)
    pygame.display.update()
