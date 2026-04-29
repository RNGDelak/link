from functools import total_ordering # psi v4.2
import sys
veblen_cutoff='max'#'X**X**w'
klammersymbolen=True
@total_ordering
class Ordinal:
   def __init__(self,I_subscript=0,subscript=0,arg=1,copies=1,addend=0):
       self.Isub=I_subscript
       self.sub=subscript
       self.arg=arg
       self.copies=copies
       self.addend=addend
   def copy(self):
       return Ordinal(self.Isub, self.sub, self.arg, self.copies, self.addend)
   def __repr__(self):
       return str(self)
   def __str__(self):
       term=''
       def k(ord):
           if ord < w+1:
               return str(ord)
           if ord.copies==1 and ord.addend==0 and ((ord<W and (W*divW(ord.arg,1)==ord.arg)) or (ord>=W and omega(ord.Isub,ord.sub+1)*div(ord.arg,omega(ord.Isub,ord.sub+1))==ord.arg)):
               return str(ord)
           return f'{{{str(ord)}}}'
       if self.arg==0:
           if self.Isub==0 and self.sub<w:
               if isinstance(self.sub, int):
                   #term='Ω'+''.join(['₀₁₂₃₄₅₆₇₈₉'[int(i)] for i in str(self.sub)])*(self.sub!=1) #ψ_α(0) = Ω_α
                   term='Ω'+('_'+str(self.sub))*(self.sub!=1)
           else:
               if self.sub==0:
                   if isinstance(self.Isub, int):
                       #term='I'+''.join(['₀₁₂₃₄₅₆₇₈₉'[int(i)] for i in str(self.Isub)])*(self.Isub!=1)  # ψ_{α,0}(0) = I_α
                       term='I'+('_'+str(self.Isub))*(self.Isub!=1)
                   else:
                       term='I_'+k(self.Isub)
               else:
                   if self.sub>=omega(self.Isub+1,0):
                       if terms(self.sub)[-1]>=omega(self.Isub+1,0):
                           if terms(self.sub)[-1]<omega(self.Isub+1,0)**2:
                               x=terms(divW(self.sub,omega(self.Isub+1,0))).copy()
                               y=0
                               while x[-1]<omega(self.Isub+1,0):
                                   y=x[-1]+y;x=x[:-1]
                                   if len(x)==0:break
                               q=omega(self.Isub,omega(self.Isub+1,0)*sum(x))+y
                               if self.sub<I*w:term='Λ'+('_'+str(q-1))*(q!=2)#term='Λ'+''.join(['₀₁₂₃₄₅₆₇₈₉'[int(i)] for i in str(divW(self.sub,I))])*(divW(self.sub,I)!=1)
                               else:term=f'Λ_{k(q)}'
                           else:
                               term=f'Ω_{{{str(self.Isub)+","+str(self.sub)}}}'
                       else:
                           x=terms(self.sub).copy()
                           y=0
                           while x[-1]<omega(self.Isub+1,0):y=x[-1]+y;x=x[:-1]
                           q=omega(self.Isub,sum(x))+y
                           term=f'Ω_{k(q)}'
                   else:
                       term = f'Ω_{k(omega(self.Isub,0)+self.sub)}'
       else:
           X=omega(self.Isub,self.sub+1)
           if veblen_cutoff.lower()!='max':p=eval(veblen_cutoff)
           else:p=psi3(self.Isub,self.sub+1,omega(self.Isub,self.sub+2))
           if self.sub==0 and self.Isub==0 and self.arg<w:term='ω'+('^'+str(self.arg))*(self.arg!=1)#term='ω'+''.join(['⁰¹²³⁴⁵⁶⁷⁸⁹'[int(i)] for i in str(self.arg)])*(self.arg!=1)
           elif terms(self.arg)[-1]>=p or self>W and terms(self.arg)[-1]<X:
               x=terms(self.arg)
               y=0
               while x[-1]<omega(self.Isub,self.sub+1):
                   y=x[-1]+y;x=x[:-1]
                   if len(x)==0:break
               if terms(self.arg)[-1]>=omega(self.Isub,self.sub+1) or w**psi3(self.Isub,self.sub,sum(x))!=psi3(self.Isub,self.sub,sum(x)):
                   if isinstance(self.sub, int) and self.Isub==0:
                       #term= 'ψ' +''.join(['₀₁₂₃₄₅₆₇₈₉'[int(i)] for i in str(self.sub)]) + '(' + str(self.arg) + ')'
                       term= 'ψ'+('_'+str(self.sub))*(self.sub>0)+'(' + str(self.arg) + ')'
                   else:
                       if self.Isub==0:
                           term= 'ψ_' +k(self.sub)+ '(' + str(self.arg) + ')'
                       else:
                           term = 'ψ_' + '{' + str(self.Isub)+','+str(self.sub)+ '}(' + str(self.arg) + ')'
               else:
                   q=log(self)
                   q1=div(q,omega(self.Isub,self.sub))
                   if q1>w:
                       if q1.arg>0:p=terms(q1.arg)[-1]>=omega(q1.Isub,q1.sub)
                       else:p=1
                       if q1.addend==0 and q1.copies==1 and((q1.Isub+q1.sub>0 and p)or(q1.sub==0 and q1.Isub==0)):q1=str(q1)
                       else:q1=f'({q1})'
                   else:
                       q1=str(q1)
                   t=terms(self)
                   l=[]
                   while t[0]>=omega(self.Isub,self.sub)**div(q,omega(self.Isub,self.sub)):
                       l.append(t[0])
                       t=t[1:]
                       if not t:break
                   q2=div(sum(l),omega(self.Isub,self.sub)**div(q,omega(self.Isub,self.sub)))
                   q3=q2
                   if q2<w+1:q2=str(q2)
                   elif (q2.addend>0 and q2<W)or(q2.addend>0 and terms(q2)[-1]<omega(q2.Isub,q2.sub)**q3):q2=f'({q2})'
                   else:q2=f'{q2}'
                   if psi3(self.Isub, self.sub, 0)**div(q,omega(self.Isub,self.sub))==div(q,omega(self.Isub,self.sub)):term=f'{q1}{f"·{q2}"*(q2!="1")}'
                   else:term=f'{psi3(self.Isub, self.sub, 0)}{f"^{q1}"*(q1!="1")}{f"·{q2}"*(q2!="1")}'
                   if sum(t)!=self:
                       if t:term+=f'+{sum(t)}'
                       return term
                   else:
                       if isinstance(self.sub, int) and self.Isub==0:
                           #term= 'ψ' +''.join(['₀₁₂₃₄₅₆₇₈₉'[int(i)] for i in str(self.sub)]) + '(' + str(self.arg) + ')'
                           term= 'ψ'+('_'+str(self.sub))*(self.sub>0)+'(' + str(self.arg) + ')'
                       else:
                           if self.Isub==0:
                               term= 'ψ_' +k(self.sub)+ '(' + str(self.arg) + ')'
                           else:
                               term = 'ψ_' + '{' + str(self.Isub)+','+str(self.sub)+ '}(' + str(self.arg) + ')'
           else:
               a1,a2=phi_inv(self.arg,X,omega(self.Isub,self.sub))
               if a1==X:a1=-1
               j=k(a2)
               if a1==0 and len(terms(a2))<2:term=f'ω^{a2}'
               elif a1==0:term=f'ω^({a2})'
               elif a1<=3 and a2<w:term=' εζηΓ'[a1]+'_'+str(a2)#''.join(['₀₁₂₃₄₅₆₇₈₉'[int(i)] for i in str(a2)])
               elif a1<=3:term=' εζηΓ'[a1]+f'_{j}'
               elif a1<X:
                   if a1<w:term=f'φ_{a1}({a2})'#+''.join(['₀₁₂₃₄₅₆₇₈₉'[int(i)] for i in str(a1)])+f'({a2})'
                   else:term=f'φ_{k(a1)}({a2})'
               else:
                   def f(a,h,str,sp):
                       if a<h:return f'@{str(a)}'
                       return f'@({array(a,h,str,f,sp)})'
                   term=f'φ({array(X*a1+a2,X,str,f,",")})'
       if self.copies!=1:
           term+=f'·{self.copies}'
       if self.addend!=0:
           term+=f'+{self.addend}'
       return term
   def as_tuple(self):
       return (self.Isub ,self.sub, self.arg, self.copies, self.addend)
   def __eq__(self,other):
       if isinstance(other,Ordinal):
           return self.as_tuple() == other.as_tuple()
       return False
   def __lt__(self,other):
       if isinstance(other,Ordinal):
           return self.as_tuple()<other.as_tuple()
       return False
   def __add__(self,other):
       if isinstance(other,Ordinal):
           if (self.Isub,self.sub,self.arg)<(other.Isub,other.sub,other.arg):
               return other
           else:
               if (self.Isub,self.sub,self.arg)==(other.Isub,other.sub,other.arg):
                   k=self.copy()
                   k.copies+=other.copies
                   if Ordinal(self.Isub,self.sub,self.arg,self.copies)<Ordinal(other.Isub,other.sub,other.arg,other.copies):k.addend+=other.addend
                   else:k.addend=other.addend
                   return k
               else:
                   k=self.copy()
                   k.addend+=other
                   return k
       else:
           k=self.copy()
           k.addend+=other
           return k
   def __radd__(self,other):
       return self
   def __sub__(self,other):
       return unadd(other,self)
   def __rsub__(self, other):
       return unadd(self,other)
   def __mul__(self, other):
       if other==0: return 0
       if other==1: return self
       if self.addend==0:
           if other<w: q=self.copy();q.copies*=other;return q
           if other.addend==0:return w**(log(self)+log(other))*other.copies
       if type(other)==int:return Ordinal(self.Isub,self.sub,self.arg,self.copies)*other+self.addend
       if other.addend>0:return self*Ordinal(other.Isub,other.sub,other.arg,other.copies)+self*other.addend
       return Ordinal(self.Isub,self.sub,self.arg,self.copies)*other
   def __rmul__(self, other):
       return Ordinal(self.Isub,self.sub,self.arg,self.copies)+other*self.addend if other!=0 else 0
   def __pow__(self, other):
       if self==w:
           if other<e0:return psi(0,other)
           return fromlist(exp(tolist(other)))
       if self.addend==0 and self.copies==1:
           return w**(log(self)*other)
       if other<w:
           if other==0: return 1
           if other==1: return self
           if other%2==0: return self**(other//2)*self**(other//2)
           return self**(other//2)*self**(other//2)*self
       if other.addend==0:
           return terms(self)[0]**other
       return self**(w**log(other)*other.copies)*self**other.addend
   def __rpow__(self, other):
       if other==1:return 1
       if other==0:return 0
       if self.addend!=0 or self.copies!=1:
           if self.addend==0:
               return other**(w**log(self)*(self.copies-1))*other**w**log(self)
           return other**(w**log(self)*self.copies)*other**self.addend
       if self>=w**w:
           return w**self
       if self==w:return w
       return w**w**(log(self)-1)
   def __hash__(self):
       return hash(self.as_tuple())
   def __getitem__(self,n):
       if cof(self)==1:return sum(terms(self)[:-1])
       if len(terms(self))>1:return sum(terms(self)[:-1])+terms(self)[-1][n]
       if cof(self)==self:
           return n
       if self.arg==0:
           if self.sub==0:return omega(self.Isub[n],0)
           if cof(self.sub)<=self:return omega(self.Isub,self.sub[n])
           if cof(self.sub).sub==0:
               def g(n):return omega(sum(terms(cof(self.sub).Isub)[:-1]),0)if n==0 else omega(sum(terms(cof(self.sub).Isub)[:-1]),self.sub[g(n-1)])
               return omega(self.Isub,self.sub[g(n)])
           def g(n):
               if cof(self.sub).sub!=0:
                   if n==0:return 0
                   return psi3(cof(self.sub).Isub,sum(terms(cof(self.sub).sub)[:-1]),self.sub[g(n-1)])
               else:return omega(sum(terms(cof(self.sub).Isub)[:-1]),0)if n==0 else omega(sum(terms(cof(self.sub).Isub)[:-1]),self.sub[g(n-1)])
           return omega(self.Isub,self.sub[g(n)])
       if cof(self.arg)==1:return psi3(self.Isub,self.sub,sum(terms(self.arg)[:-1]))*n
       if cof(self.arg)<=self:return psi3(self.Isub,self.sub,self.arg[n])
       def g(n):
            if n==0:return 0
            return psi3(cof(self.arg).Isub,sum(terms(cof(self.arg).sub)[:-1]),self.arg[g(n-1)])
       return psi3(self.Isub,self.sub,self.arg[g(n)])# if terms(self.arg)[-1]>=W**2 else g(n)

def cof(a):
   if a==0:return 0
   if a<w: return 1
   if len(terms(a))!=1: return cof(terms(a)[-1])
   if a.arg==0:
       #if omega(a.Isub+1,0)<a.sub<omega(a.Isub+2,0):return a
       if a.sub==0:
           if a.Isub==0:return w
           return a if cof(a.Isub)==1 else cof(a.Isub)
       elif cof(a.sub)==1:return a
       else:return cof(a.sub) if cof(a.sub)<a else w
   if cof(a.arg)in[1,w]:return w
   return cof(a.arg) if cof(a.arg)<a else w
def std(a):return fromlist(SF(tolist(a)))
def omega(n,k):
   if (n,k)==(0,0):
       return 1
   return Ordinal(n,k,0)
def psi(n,k):
   if (n,k)!=(0,0):
       return Ordinal(0,n,k)
   return 1
def psi3(n,k,r):
   if (n,k,r)!=(0,0,0):
       return Ordinal(n,k,r)
   return 1
def divW(k,q=1):
   # k/(W^q)
   if k<W**q:
       return 0
   if q==0:
       return k
   return w**(unadd(W*q,log(k)))+divW(sum(terms(k)[1:]),q)
def div(k,q):
   if k<q:
       return 0
   return w**(unadd(log(q),log(k)))+div(sum(terms(k)[1:]),q)
def csub1(k):
   if k>=w:
       return k
   return k-1
def is_successor(self):
   if type(self)==int:return self>0
   return is_successor(self.addend)
def terms(ord):
   if ord==0:
       return []
   if ord<w:
       if ord>100:
           raise Exception
       return [1]*int(ord)
   return [Ordinal(ord.Isub,ord.sub,ord.arg)]*ord.copies+terms(ord.addend)
def tolist(ord):
   if ord<=0: return []
   if isinstance(ord,int): return [[],[],[],tolist(ord-1)]
   return [tolist(terms(ord)[0].Isub),tolist(terms(ord)[0].sub),tolist(terms(ord)[0].arg),tolist(sum(terms(ord)[1:]))]
def fromlist(Q):
   if Q==[]:
       return 0
   return psi3(fromlist(Q[0]),fromlist(Q[1]),fromlist(Q[2]))+fromlist(Q[3])
def add(a,b): # will be converted soon - 275 - 368 will be gone
   if a==[]:
       return b
   if b==[]:
       return a
   if [a[0],a[1],a[2]]<[b[0],b[1],b[2]]:
       return b
   return [a[0],a[1],a[2],add(a[3],b)]
def _unadd(a,b):
   if a==[]:
       return b
   if b==[]:
       return []
   if [a[0],a[1],a[2]]<[b[0],b[1],b[2]]:
       return b
   if [a[0],a[1],a[2]]>[b[0],b[1],b[2]]:
       return []
   return _unadd(a[3],b[3])
def unadd(a,b):
   return fromlist(_unadd(tolist(a),tolist(b)))
def mul(a,b):
   if b==[]:
       return []
   if a==[]:
       return []
   [x,y,z,w]=b
   if [x,y,z,[]]==[[],[],[],[]]:
       return add(a,mul(a,w))
   return add(exp(add(_log(a),_log(b))),mul(a,w))
def exp(a):
   if a<[[],[[],[],[],[]],[],[]]:
       return SF([[],[],a,[]])
   else:
       return SF([a[0],a[1],_unadd([a[0],a[1],[],[]],a),[]])
def _log(a):
    if a==[]:return []
    [x,y,z,w]=a
    p,q=split(z,[x,add(y,[[],[],[],[]]),[],[]])
    if x==y==p==[]:return q
    return [x,y,p,q]
def split(a,l):
   # split n into terms >=l and terms <l
   if a==[]:
       return [],[]
   [x,y,z,w]=a
   if [x,y,z,[]]<l:
       return [],a
   p,q=split(w,l)
   return [x,y,z,p],q
def SFpsi(a,b,p,c,q):
   # assuming [n,a,b,[]], p, and q are standard,
   # return standard form of [n,add(a,p),add(b,q),[]]
   one=[[],[],[],[]]
   if p==[]:
       if q<[a,b,[],[]]:
           return [a,b,add(c,q),[]]
       [x,y,z,w]=q
       nc=add(c,[x,y,z,[]])
       m,f=findmaxsub([q[0],q[1],q[2],[]],[a,b,[],[]])
       if m>nc:
           t,_=split(m,[x,add(y,one),[],[]])
           nc=add(t,_unadd(f(t),[x,y,z,[]]))
       return SFpsi(a,b,p,nc,w)
   # not [n+1,[],[],[]] since stuff like psi_Wfp
   if p<[a,[],[],[]]:
       return SFpsi(a,add(b,p),[],c,q)
   [x,y,z,w]=p
   nb=add(b,[x,y,z,[]])
   m,f=findmaxsub([p[0],p[1],p[2],[]],[add(a,one),[],[],[]])
   if m>nb:
       t,_=split(m,[x,add(y,one),[],[]])
       nb=add(t,_unadd(f(t),[x,y,z,[]]))
   return SFpsi(a,nb,w,c,q)
def SF(a):
   if a==[]:
       return []
   [x,y,z,w]=a
   return add(SFpsi(SF(x),[],SF(y),[],SF(z)),SF(w))
def findmaxsub(a,b):
   # return the largest subterm in n, ignoring terms below a
   # and n function that places something at that position,
   #  deleting everything afterward
   #  eg psi(W_2+psi_1(W_3)+3)->lambda x:W_2+psi_1(x)
   if a<b or a==[]:
       return [],None
   ret=a
   retf=lambda x:x
   for i in [0,1,2,3]:
       v,f=findmaxsub(a[i],b)
       if v>ret:
           ret=v
           # friggin python closure bullshit
           retf=lambda x,f=f,i=i:a[:i]+[f(x)]+[[]]*(3-i)
   return ret,retf
def log(ord):
   if ord<w:
       return 0
   if ord<e0:
       return ord.arg
   return fromlist(_log(tolist(ord)))
def logW(ord):
   return divW(fromlist(_log(tolist(ord)))) if ord>=W else 0
def subterms2(a,j=Ordinal(0,1,0)):
    if a==0:return []
    if a<j:return [(a,j)]
    x=terms(a)
    m=[] # [exp,coef]
    for i in x:
        c=div(log(i),j)
        if m==[]:m+=[[div(log(i),j),div(i,j**div(log(i),j))]]
        elif m[-1][0]!=c:m+=[[div(log(i),j),div(i,j**div(log(i),j))]]
        else:m[-1][1]+=div(i,j**div(log(i),j))
    r=[]
    m1=[]
    for i in m:
        m1+=[i]
        t=sum([j**q[0]*q[1]for q in m1[:-1]])
        r+=[(p,t+W**q) for p,q in subterms2(i[0],j)]
        if i[1]>1 or i[0]==0:r+=[(i[1],sum([j**q[0]*q[1]for q in m1[:-1]])+j**(i[0]+1))]
    return r
def phi(*args):
    if args==():return 1
    #if len(args)==1:return psi(0,args[0])
    X=list(args[:-1])
    Y=[]
    b=args[-1]
    for i in range(len(X)):
        Y.append(W**i*X[::-1][i])
    a=sum(Y[::-1])
    if b>=W:
        c=div(b,W)
        b-=W*c
        a+=c
    X=subterms2(a) if a>0 else [(0,W)]
    c,d=X[-1]
    l=W**d*div(c.arg if c>=w else 0,W**d)
    Y=X[:-1] or [(0,0)]
    if psi(0,l)==c and c>=w and max(Y)[0]<c:
        u=0
        t,v=max(X)
    else:
        t,v=max(X)
        if t<G0:u=-1 # ideal: t==1, this is more efficient
        elif t>=phi(v,0):u,l=1,W**v*div(t.arg,W**v)
        else:u=-1
    if b<w:p=0
    else:
        if t>=w:r=l
        else:r=0
        r=psi(0,r+W**(a+1))
        p=W**(a+1)*div(b.arg,W**(a+1))*(b>=r)
        if p==1:p=0
    b=b-(psi(0,p)if p>0 else 0)
    if p==0 and u and a:b=1+b
    return psi(0,(p or l)+W**a*b)
def vartheta(*args):
    if args==():return 1
    #if len(args)==1:return psi(0,args[0])
    X=list(args[:-1])
    Y=[]
    b=args[-1]
    for i in range(len(X)):
        Y.append(W**i*X[::-1][i])
    a=sum(Y[::-1])
    if b>=W:
        c=div(b,W)
        b-=W*c
        a+=c
    X=subterms2(a) if a>0 else [(0,W)]
    c,d=X[-1]
    l=W**d*div(c.arg if c>=w else 0,W**d)
    Y=X[:-1] or [(0,0)]
    if psi(0,l)==c and c>=w and max(Y)[0]<c:
        u=0
        t,v=max(X)
    else:
        t,v=max(X)
        if t<G0:u=-1 # ideal: t==1, this is more efficient
        elif t>=phi(v,0):u,l=1,W**v*div(t.arg,W**v)
        else:u=-1
    if b<w:p=0
    else:
        if t>=w:r=l
        else:r=0
        r=psi(0,r+W**(a+1))
        p=W**(a+1)*div(b.arg,W**(a+1))*(b>=r)
        if p==1:p=0
    b=b-(psi(0,p)if p>0 else 0)
    if a or p>0:b=1+b
    return psi(0,(p or l)+W**a*b)

def prod(A):
    x=1
    for i in A:x*=i
    return x
def check(a,b,j): # a = P(W), b = beta, j = W
    if a<j:
        if a<b:return 1
        if a==b:return 2
        if a>b:return 3
    #if b==0:return 1
    x=terms(a)
    m=[] # [exp,coef]
    for i in x:
        c=div(log(i),j)
        if m==[]:m+=[[div(log(i),j),div(i,j**div(log(i),j))]]
        elif m[-1][0]!=c:m+=[[div(log(i),j),div(i,j**div(log(i),j))]]
        else:m[-1][1]+=div(i,j**div(log(i),j))
    s=[]
    for i in m:s+=i
    if max([check(i,b,j) for i in s])==1:return 1
    s=[0]+s
    if max([check(i,b,j) for i in s[:-2]])==1 and ((s[-2]==b and s[-1]==1) or (check(s[-2],b,j)==1 and s[-1]==b)):return 2
    return 3
def phi_inv(a,h,s): # credit to regtar
    if s==1:s=w
    if a<h:return 0,a
    if terms(a)[-1]<h:
        return 0,psi3(s.Isub,s.sub,h*div(a,h))+(a-h*div(a,h))
    x=terms(a)
    y=[]
    for i in x:
        m=div(log(i),h)
        if y==[]:y+=[[div(log(i),h),div(i,h**div(log(i),h))]]
        elif y[-1][0]!=m:y+=[[div(log(i),h),div(i,h**div(log(i),h))]]
        else:y[-1][1]+=div(i,h**div(log(i),h))
    X=y[-1][1]
    b=psi3(s.Isub,s.sub,sum([h**i*j for i,j in y[:-1]]))if len(y)>1 else 0
    if s==w:s=0
    if b==0 and a>W:b=s
    l=div(log(h**y[-1][0]*y[-1][1]),h)
    if check(l,b,h)==1:X=b+X # case 1
    elif check(l,b,h)==2:X=1+X # case 2
    X-=1
    return y[-1][0],X

def multerms(a):
    if a==0:return []
    if a<w:return [a]
    return [Ordinal(a.Isub,a.sub,a.arg,a.copies,0)]+multerms(a.addend)
def flatten(x):
    while list in [type(i)for i in x]:
        m=[]
        for i in x:
            if type(i)==list:m+=i
            else:m+=[i]
        x=m[::]
    return x
def subterms(a,O):
    h=multerms(a)
    if a==0:return []
    if a<O:return [0,a]
    if O**a==a and a.arg!=0:
        q=subterms(a.arg,omega(0,multerms(a.arg)[-1].sub))
    else:
        O1=omega(0,O.sub+1)
        p,q1=sum([i for i in h if i>=O1]),sum([i for i in h if i<O1])
        q=[]
        while q1!=0:
            e1=div(log(q1),O)
            e2=div(q1,O**e1)
            q+=subterms(e1,O)+[e2]
            q1-=O**e1*e2
    while max(flatten(q))>=O if len(q) else False:
        q=[subterms(i,omega(0,multerms(i)[-1].sub)) if i>=O else i for i in q]
    return subterms(a.sub,O)+flatten(q)
def NF(x,v=0):
    if type(x)==int:return x
    y=[NF(i) for i in terms(x.arg)]
    t=0
    if y!=[]:
        while y[0]>=omega(0,x.sub+2):
            t+=y[0]
            y=y[1:]
            if not y:break
    if t>=w:
        z=[i.sub>=x.sub+w and terms(i.sub)[-1]==1 or div(x.sub,w)==div(i.sub,w) and i.sub-div(i.sub,w)>x.sub-div(x.sub,w) for i in terms(t)]
        if sum(z):
            s=len(z)-z[::-1].index(True)-1
            t=psi(sum(terms(terms(t)[s].sub)[:-1]),t)
    q=psi(x.sub,t+sum(y))
    if q>=w:
        a,b=q.sub,q.arg
        q=psi(NF(a,v),NF(b,v))
    return q+NF(sum(terms(x)[1:]),v)
def theta_inv(a): # credit to regtar
    s=w
    h=W
    x=terms(a)
    y=[]
    for i in x:
        m=div(log(i),h)
        if y==[]:y+=[[div(log(i),h),div(i,h**div(log(i),h))]]
        elif y[-1][0]!=m:y+=[[div(log(i),h),div(i,h**div(log(i),h))]]
        else:y[-1][1]+=div(i,h**div(log(i),h))
    X=y[-1][1]
    b=psi3(s.Isub,s.sub,sum([h**i*j for i,j in y[:-1]]))if len(y)>1 else 0
    if s==w:s=0
    if b==0 and a>W:b=s
    l=div(log(h**y[-1][0]*y[-1][1]),h)
    if check(l,b,h)==1:X=b+(X-1) # case 1
    X-=1
    return y[-1][0],X
def array(x,j,str,sub,sp):
    if(x<j):return str(x)
    a=x
    s=omega(x.Isub,x.sub)
    x=terms(x)
    y=[]
    for i in x:
        m=div(log(i),j)
        if y==[]:y+=[[div(log(i),j),div(i,s**div(log(i),j))]]
        elif y[-1][0]!=m:y+=[[div(log(i),j),div(i,s**div(log(i),j))]]
        else:y[-1][1]+=div(i,s**div(log(i),j))
    if a<j**w:
        z=[0]*(div(log(sum(x)),j)+1)
        for i in y:
            z[-i[0]-1]=i[1]
        return ','.join([str(i)for i in z])
    elif klammersymbolen:
        return sp.join([f'{str(i)}{sub(k,j,str,sp)}'for k,i in y])
    else:
        q=[i[0]for i in y]
        d=[q[i]-q[i+1]for i in range(len(y)-1)]+[q[-1]]
        o=''
        for i in range(len(y)):
            o+=str(y[i][1])
            for k in terms(d[i])[::-1]:
                o+=f'{{{array(log(k),j,str,sub,sp)}}}0'
            if i<len(y)-1:o=o[:-1]
        return o.replace('{0}',',')
def aslatex(self):
    def sup(ord):return f'^{{{aslatex(ord)}}}'
    def sub(ord):return f'_{{{aslatex(ord)}}}'
    if type(self)!=Ordinal:
        return str(self)
    term=''
    if self.sub==0 and self.Isub==0:
        if terms(self.arg)[-1]>=W**(W+1):
           term=f'\\psi_0({aslatex(self.arg)})'
        else:
           a1,a2=phi_inv(self)
           if a1==0 and a2==1:term='\\omega'
           elif a1==0:term=f'\\omega{sup(a2)}'
           elif a1<=3:term=[None,'\\varepsilon','\\zeta','\\eta'][a1]+f'{sub(a2)}'
           elif a1==terms(self)[0]:
               x=terms(self.arg)
               y=0
               while x[-1]<W**(W+1):
                   y=divW(x[-1],W)+y
                   x=x[:-1]
                   if not x:break
               y=psi(0,sum(x))+y-2
               term=f'\\Gamma{sub(y)}'
           else:
               term=f'\\varphi{sub(a1)}({aslatex(a2)})'
    else:
        if self.arg==0:
           if self.sub==0:
                term=f'I{sub(self.Isub)*(self.Isub>1)}'
           else:
               if self.sub>=omega(self.Isub+1,0):
                   if terms(self.sub)[-1]>=omega(self.Isub+1,0):
                       if terms(self.sub)[-1]<omega(self.Isub+1,0)**2:
                           x=terms(divW(self.sub,omega(self.Isub+1,0))).copy()
                           y=0
                           while x[-1]<omega(self.Isub+1,0):
                               y=x[-1]+y;x=x[:-1]
                               if len(x)==0:break
                           q=omega(self.Isub,omega(self.Isub+1,0)*sum(x))+y-1
                           term=f'\\Lambda{sub(q)*(q>1)}'
                       else:
                           term=f'\\Omega_{{{aslatex(self.Isub)},{aslatex(self.sub)}}}'
                   else:
                       x=terms(self.sub).copy()
                       y=0
                       while x[-1]<omega(self.Isub+1,0):y=x[-1]+y;x=x[:-1]
                       q=omega(self.Isub,sum(x))+y
                       term=f'\\Omega{sub(q)}'
               else:
                   term = f'\\Omega{sub(omega(self.Isub,0)+self.sub-1)*((self.Isub,self.sub)!=(0,1))}'
        else:
            x=terms(self.arg)
            y=0
            while x[-1]<omega(self.Isub,self.sub+1):
                y=x[-1]+y;x=x[:-1]
                if len(x)==0:break
            if terms(self.arg)[-1]>=omega(self.Isub,self.sub+1) or w**psi3(self.Isub,self.sub,sum(x))!=psi3(self.Isub,self.sub,sum(x)):
                if self.Isub==0:
                    term=f'\\psi{sub(self.sub)}({aslatex(self.arg)})'
                else:
                    term = f'\\psi_{{{aslatex(self.Isub)},{aslatex(self.sub)}}}({aslatex(self.arg)})'
            else:
                q = psi3(self.Isub,self.sub,sum(x))+y
                q1=div(q,omega(self.Isub,self.sub))
                q2=w**(q-omega(self.Isub,self.sub)*q1)
                if psi3(self.Isub,self.sub,0)**div(q,omega(self.Isub,self.sub))==div(q,omega(self.Isub,self.sub)):term = f'{aslatex(q1)}{f"{chr(92)}cdot {aslatex(q2)}"*(q2!=1)}'
                else:term = f'{{{aslatex(psi3(self.Isub, self.sub, 0))}}}{f"{sup(q1)}"*(q1!=1)}{f"{chr(92)}cdot {aslatex(q2)}"*(q2!=1)}'
    if self.copies != 1:
        term += f'\\cdot {self.copies}'
    if self.addend != 0:
        term += f'+{aslatex(self.addend)}'
    return term
def psi_str(ord):
   def k(ord):
       if ord < psi(0,1):
           return str(ord)
       if ord.copies==1 and ord.addend==0:
           return psi_str(ord)
       return f'{{{psi_str(ord)}}}'
   if type(ord)!=Ordinal:
       return str(ord)
   term=''
   if ord.Isub==0:
       if ord.sub==0:
           term+=f'ψ({psi_str(ord.arg)})'
       else:
           term += f'ψ_{k(ord.sub)}({psi_str(ord.arg)})'
   else:
       term += f'ψ_{k(ord.Isub)}({psi_str(ord.sub)},{psi_str(ord.arg)})'
   if ord<psi(0,2):
       term='ω'
   term = ((term+'+')*ord.copies)[:-1]
   if ord.addend != 0:
       term += f'+{psi_str(ord.addend)}'
   return term
def ascii(self):
    def sup(t):return f'^({t})'
    def sub(t):return f'({t})'
    if type(self)!=Ordinal:
        return str(self)
    term = ''
    if self.sub==0 and self.Isub==0:
        if terms(self.arg)[-1]>=psi(1,omega(0,2)):
           term=f'psi({ascii(self.arg)})'
        else:
           a1,a2=phi_inv(self.arg,W,w)
           if a1==W:a1=-1
           if a1==0 and a2==1:term='w'
           elif a1==0:term=f'w{sup(ascii(a2))}'
           elif a1<=3:term=' eznG'[a1]+f'{sub(ascii(a2))}'
           elif a1<W:term=f'phi({ascii(a1)},{ascii(a2)})'
           else:
               def f(a,h,str,sp):
                   if a<h:return f'@{str(a)}'
                   return f'@({array(a,h,ascii,f,sp)})'
               term=f'phi({array(W*a1+a2,W,ascii,f,",")})'
    else:
        if self.arg==0:
           if self.sub==0:
                term=f'I{sub(ascii(self.Isub))*(self.Isub>1)}'
           else:
               if self.Isub==0:
                   term=f'W{sub(ascii(self.sub))*(self.sub>1)}'
               else:
                   term='W'+sub(f'{ascii(self.Isub)},{ascii(self.sub)}')
        else:
            x=terms(self.arg)
            y=0
            while x[-1]<omega(self.Isub,self.sub+1):
                y=x[-1]+y;x=x[:-1]
                if len(x)==0:break
            if terms(self.arg)[-1]>=omega(self.Isub,self.sub+1) or w**psi3(self.Isub,self.sub,sum(x))!=psi3(self.Isub,self.sub,sum(x)):
                if self.Isub==0:
                    term=f'psi({ascii(self.sub)},{ascii(self.arg)})'
                else:
                    term = f'psi({ascii(self.Isub)},{ascii(self.sub)},{ascii(self.arg)})'
            else:
                q=log(self)
                q1=div(q,omega(self.Isub,self.sub))
                if q1>w:
                    if q1.arg>0:p=terms(q1.arg)[-1]>=omega(q1.Isub,q1.sub)
                    else:p=1
                    if q1.addend==0 and q1.copies==1 and((q1.Isub+q1.sub>0 and p)or(q1.sub==0 and q1.Isub==0)):q1=ascii(q1)
                    else:q1=f'({ascii(q1)})'
                else:
                    q1=ascii(q1)
                t=terms(self)
                l=[]
                while t[0]>=omega(self.Isub,self.sub)**div(q,omega(self.Isub,self.sub)):
                    l.append(t[0])
                    t=t[1:]
                    if not t:break
                q2=div(sum(l),omega(self.Isub,self.sub)**div(q,omega(self.Isub,self.sub)))
                if q2<w+1:q2=ascii(q2)
                elif q2.copies>1 or(q2.addend>0 and q2<W)or(q2.addend>0 and terms(q2)[-1]<omega(q2.Isub,q2.sub)):q2=f'({ascii(q2)})'
                else:q2=ascii(q2)
                if psi3(self.Isub, self.sub, 0)**div(q,omega(self.Isub,self.sub))==div(q,omega(self.Isub,self.sub)):term=f'{q1}{f"*{q2}"*(q2!="1")}'
                else:term=f'{ascii(psi3(self.Isub, self.sub, 0))}{f"^{q1}"*(q1!="1")}{f"*{q2}"*(q2!="1")}'
                if t:term+=f'+{ascii(sum(t))}'
                return term
    if self.copies != 1:
        term += f'*{self.copies}'
    if self.addend != 0:
        term += f'+{ascii(self.addend)}'
    return term
def F(self):
    def sup(ord):
        if ord<w+1:return f'^{F(ord)}'
        if ord.arg>0:p=terms(ord.arg)[-1]>=omega(ord.Isub,ord.sub)
        else:p=1
        if ord.addend==0 and ord.copies==1 and((ord.Isub+ord.sub>0 and p)or(ord.sub==0 and ord.Isub==0)):return f'^{F(ord)}'
        else:return f'^({F(ord)})'
    def sub(ord):
       if ord<w+1:return f'_{F(ord)}'
       if ord.copies==1 and ord.addend==0 and ((ord<W and (W*divW(ord.arg,1)==ord.arg)) or ord>=W and (ord.arg==0 or ord.arg>=omega(ord.Isub,ord.sub+1))):return f'_{F(ord)}'
       return f'_{{{F(ord)}}}'

    if type(self)!=Ordinal:
        return str(self)
    term=''
    if self.sub==0 and self.Isub==0:
        if terms(self.arg)[-1]>=W**W**w:
           term=f'ψ_0({F(self.arg)})'
        else:
           a1,a2=phi_inv(self)
           if a1==W:a1=-1
           if a1==0 and a2==1:term='ω'
           elif a1==0:term=f'ω{sup(a2)}'
           elif a1<=3:term=' εζηΓ'[a1]+f'{sub(a2)}'
           elif a1<W:term=f'φ{sub(a1)}({F(a2)})'
           else:term=f'φ({",".join([F(i)for i in array(a1)])},{F(a2)})'
    else:
        if self.arg==0:
           if self.sub==0:
                term=f'I{sub(self.Isub)*(self.Isub>1)}'
           else:
               if self.sub>=omega(self.Isub+1,0):
                   if terms(self.sub)[-1]>=omega(self.Isub+1,0):
                       if terms(self.sub)[-1]<omega(self.Isub+1,0)**2:
                           x=terms(divW(self.sub,omega(self.Isub+1,0))).copy()
                           y=0
                           while x[-1]<omega(self.Isub+1,0):
                               y=x[-1]+y;x=x[:-1]
                               if len(x)==0:break
                           q=omega(self.Isub,omega(self.Isub+1,0)*sum(x))+y-1
                           term=f'Λ{sub(q)*(q>1)}'
                       else:
                           term=f'Ω_{{{F(self.Isub)},{F(self.sub)}}}'
                   else:
                       x=terms(self.sub).copy()
                       y=0
                       while x[-1]<omega(self.Isub+1,0):y=x[-1]+y;x=x[:-1]
                       q=omega(self.Isub,sum(x))+y
                       term=f'Ω{sub(q)}'
               else:
                   term = f'Ω{sub(omega(self.Isub,0)+self.sub-1)*((self.Isub,self.sub)!=(0,1))}'
        else:
            x=terms(self.arg)
            y=0
            while x[-1]<omega(self.Isub,self.sub+1):
                y=x[-1]+y;x=x[:-1]
                if len(x)==0:break
            if terms(self.arg)[-1]>=omega(self.Isub,self.sub+1) or w**psi3(self.Isub,self.sub,sum(x))!=psi3(self.Isub,self.sub,sum(x)):
                if self.Isub==0:
                    term=f'ψ{sub(self.sub)}({F(self.arg)})'
                else:
                    term = f'ψ_{{{F(self.Isub)},{F(self.sub)}}}({F(self.arg)})'
            else:
                q=log(self)
                q1=div(q,omega(self.Isub,self.sub))
                if q1>w:
                    if q1.arg>0:p=terms(q1.arg)[-1]>=omega(q1.Isub,q1.sub)
                    else:p=1
                    if q1.addend==0 and q1.copies==1 and((q1.Isub+q1.sub>0 and p)or(q1.sub==0 and q1.Isub==0)):q1=F(q1)
                    else:q1=f'({F(q1)})'
                else:
                    q1=F(q1)
                t=terms(self)
                l=[]
                while t[0]>=omega(self.Isub,self.sub)**div(q,omega(self.Isub,self.sub)):
                    l.append(t[0])
                    t=t[1:]
                    if not t:break
                q2=div(sum(l),omega(self.Isub,self.sub)**div(q,omega(self.Isub,self.sub)))
                if q2<w+1:q2=F(q2)
                elif(q2.addend>0 and q2<W)or(q2.addend>0 and terms(q2)[-1]<omega(q2.Isub,q2.sub)):q2=f'({F(q2)})'
                else:q2=F(q2)
                if psi3(self.Isub, self.sub, 0)**div(q,omega(self.Isub,self.sub))==div(q,omega(self.Isub,self.sub)):term=f'{q1}{f"*{q2}"*(q2!="1")}'
                else:term=f'{F(psi3(self.Isub, self.sub, 0))}{f"^{q1}"*(q1!="1")}{f"*{q2}"*(q2!="1")}'
                if t:term+=f'+{F(sum(t))}'
                return term
    if self.copies != 1:
        term += f'*{self.copies}'
    if self.addend != 0:
        term += f'+{F(self.addend)}'
    return term
def lvs(y):
    t=[0]
    for i in y:
        if i in '[}':t.append(t[-1]+1)
        if i in ']{':t.append(t[-1]-1)
    return[max(t),abs(min(t))]
def sup(a):
    return '['*(lvs(a)[1]+1)+a+']'*(lvs(a)[1]+1)
def sub(a):
    return '{'*(lvs(a)[0]+1)+a+'}'*(lvs(a)[0]+1)
def gamestuffs(self):
    if type(self)!=Ordinal:
        return str(self)
    term = ''
    X=omega(self.Isub,self.sub+1)
    if veblen_cutoff.lower()!='max':p=eval(veblen_cutoff)
    else:p=psi3(self.Isub,self.sub+1,omega(self.Isub,self.sub+2))

    if (not(terms(self.arg)[-1]>=p or self>W and terms(self.arg)[-1]<X))if self.arg>0 else (self.sub==self.Isub==0):
           a1,a2=phi_inv(self.arg,X,omega(self.Isub,self.sub))
           if a1==W:a1=-1
           if a1==0 and a2==1:term='ω'
           elif a1==0:term=f'ω{sup(gamestuffs(a2))}'
           elif a1<=3:term=' εζηΓ'[a1]+f'{sub(gamestuffs(a2))}'
           elif a1<W:term=f'φ{sub(gamestuffs(a1))}({gamestuffs(a2)})'
           else:
               def f(a,h,str,sp):
                   if a<h:return '@'+str(a)
                   return '@('+array(a,h,str,f,sp)+')'
               term=f'φ({array(W*a1+a2,W,gamestuffs,f,",")})'
    else:
        if self.arg==0:
           if self.sub==0:
                term=f'I{sub(gamestuffs(self.Isub))*(self.Isub>1)}'
           else:
               if self.sub>=omega(self.Isub+1,0):
                   if terms(self.sub)[-1]>=omega(self.Isub+1,0):
                       if terms(self.sub)[-1]<omega(self.Isub+1,0)**2:
                           x=terms(divW(self.sub,omega(self.Isub+1,0))).copy()
                           y=0
                           while x[-1]<omega(self.Isub+1,0):
                               y=x[-1]+y;x=x[:-1]
                               if len(x)==0:break
                           q=omega(self.Isub,omega(self.Isub+1,0)*sum(x))+y-1
                           term=f'Λ{sub(gamestuffs(q))*(q>1)}'
                       else:
                           term=f'Ω{sub(gamestuffs(self.Isub))}{{,}}{sub(gamestuffs(self.sub))}'
                   else:
                       x=terms(self.sub).copy()
                       y=0
                       while x[-1]<omega(self.Isub+1,0):y=x[-1]+y;x=x[:-1]
                       q=omega(self.Isub,sum(x))+y
                       term=f'Ω{sub(gamestuffs(q))}'
               else:
                   term = f'Ω{sub(gamestuffs(omega(self.Isub,0)+self.sub-1))*((self.Isub,self.sub)!=(0,1))}'
        else:
            x=terms(self.arg)
            y=0
            while x[-1]<omega(self.Isub,self.sub+1):
                y=x[-1]+y;x=x[:-1]
                if len(x)==0:break
            if terms(self.arg)[-1]>=omega(self.Isub,self.sub+1) or w**psi3(self.Isub,self.sub,sum(x))!=psi3(self.Isub,self.sub,sum(x)):
                if self.Isub==0:
                    term=f'ψ{sub(gamestuffs(self.sub))*(self.sub>0)}({gamestuffs(self.arg)})'
                else:
                    term = f'ψ{sub(gamestuffs(self.Isub))}{{,}}{sub(gamestuffs(self.sub))}({gamestuffs(self.arg)})'
            else:
                q=log(self)
                q1=gamestuffs(div(q,omega(self.Isub,self.sub)))
                t=terms(self)
                l=[]
                while t[0]>=omega(self.Isub,self.sub)**div(q,omega(self.Isub,self.sub)):
                    l.append(t[0])
                    t=t[1:]
                    if not t:break
                q2=div(sum(l),omega(self.Isub,self.sub)**div(q,omega(self.Isub,self.sub)))
                if q2<w+1:q2=gamestuffs(q2)
                elif q2.copies>1 or(q2.addend>0 and q2<W)or(q2.addend>0 and terms(q2)[-1]<omega(q2.Isub,q2.sub)):q2=f'·{"("*(q2.addend>0)+gamestuffs(q2)+")"*(q2.addend>0)}'
                else:q2=gamestuffs(q2)
                if psi3(self.Isub, self.sub, 0)**div(q,omega(self.Isub,self.sub))==div(q,omega(self.Isub,self.sub)):term=f'{q1}{f"{q2}"*(q2!="1")}'
                else:term=f'{gamestuffs(omega(self.Isub,self.sub))}{f"{sup(q1)}"*(q1!="1")}{f"{q2}"*(q2!="1")}'
                if t:term+=f'+{gamestuffs(sum(t))}'
                return term
    if self.copies != 1:
        term += f'{self.copies}'
    if self.addend != 0:
        term += f'+{gamestuffs(self.addend)}'
    return term

def HTML(self):
    def sup(a):
        return f'<sup>{a}</sup>'
    def sub(a):
        return f'<sub>{a}</sub>'
    if type(self)!=Ordinal:
        return str(self)
    term = ''
    X=omega(self.Isub,self.sub+1)
    p=eval(veblen_cutoff)if veblen_cutoff.lower()!='max'else psi3(self.Isub,self.sub+1,omega(self.Isub,self.sub+2))
    if self<W or self>=W and self.arg>0:
        if X<=terms(self.arg)[-1]<p or self<W:
            if terms(self.arg)[-1]>=p:
               term=f'ψ({HTML(self.arg)})'
            else:
               a1,a2=phi_inv(self.arg,X,omega(self.Isub,self.sub))
               if a1==X:a1=-1
               if a1==0 and a2==1:term='ω'
               elif a1==0:term=f'ω{sup(HTML(a2))}'
               elif a1<=3:term=' εζηΓ'[a1]+f'{sub(HTML(a2))}'
               elif a1<W:term=f'φ{sub(HTML(a1))}({HTML(a2)})'
               else:
                   def f(a,h,str,sp):
                       if a<h:return'<sub>'+str(a)+'</sub>'
                       return sub('('+array(a,h,str,f,sp)+')')
                   term=f'φ({array(X*a1+a2,X,HTML,f,",")})'
    if self>=W:
        if self.arg==0:
           if self.sub==0:
                term=f'I{sub(HTML(self.Isub))*(self.Isub>1)}'
           else:
               if self.sub>=omega(self.Isub+1,0):
                   if terms(self.sub)[-1]>=omega(self.Isub+1,0):
                       if terms(self.sub)[-1]<omega(self.Isub+1,0)**2:
                           x=terms(divW(self.sub,omega(self.Isub+1,0))).copy()
                           y=0
                           while x[-1]<omega(self.Isub+1,0):
                               y=x[-1]+y;x=x[:-1]
                               if len(x)==0:break
                           q=omega(self.Isub,omega(self.Isub+1,0)*sum(x))+y-1
                           term=f'Λ{sub(HTML(q))*(q>1)}'
                       else:
                           term=f'Ω{sub(HTML(self.Isub))}<sub>,</sub>{sub(HTML(self.sub))}'
                   else:
                       x=terms(self.sub).copy()
                       y=0
                       while x[-1]<omega(self.Isub+1,0):y=x[-1]+y;x=x[:-1]
                       q=omega(self.Isub,sum(x))+y
                       term=f'Ω{sub(HTML(q))}'
               else:
                   term = f'Ω{sub(HTML(omega(self.Isub,0)+self.sub-1))*((self.Isub,self.sub)!=(0,1))}'
        elif terms(self.arg)[-1]<X:
            x=terms(self.arg)
            y=0
            while x[-1]<omega(self.Isub,self.sub+1):
                y=x[-1]+y;x=x[:-1]
                if len(x)==0:break
            if terms(self.arg)[-1]>=omega(self.Isub,self.sub+1) or w**psi3(self.Isub,self.sub,sum(x))!=psi3(self.Isub,self.sub,sum(x)):
                if self.Isub==0:
                    term=f'ψ{sub(HTML(self.sub))}({HTML(self.arg)})'
                else:
                    term = f'ψ{sub(HTML(self.Isub))}<sub>,</sub>{sub(HTML(self.sub))}({HTML(self.arg)})'
            else:
                q=log(self)
                q1=HTML(div(q,omega(self.Isub,self.sub)))
                t=terms(self)
                l=[]
                while t[0]>=omega(self.Isub,self.sub)**div(q,omega(self.Isub,self.sub)):
                    l.append(t[0])
                    t=t[1:]
                    if not t:break
                q2=div(sum(l),omega(self.Isub,self.sub)**div(q,omega(self.Isub,self.sub)))
                if q2<w+1:q2=HTML(q2)
                elif q2.copies>1 or(q2.addend>0 and q2<W)or(q2.addend>0 and terms(q2)[-1]<omega(q2.Isub,q2.sub)):q2=f'·{"("*(q2.addend>0)+gamestuffs(q2)+")"*(q2.addend>0)}'
                else:q2=HTML(q2)
                if psi3(self.Isub, self.sub, 0)**div(q,omega(self.Isub,self.sub))==div(q,omega(self.Isub,self.sub)):term=f'{q1}{f"{q2}"*(q2!="1")}'
                else:term=f'{HTML(omega(self.Isub,self.sub))}{f"{sup(q1)}"*(q1!="1")}{f"{q2}"*(q2!="1")}'
                if t:term+=f'+{HTML(sum(t))}'
                return term
    if self.copies != 1:
        term += f'{self.copies}'
    if self.addend != 0:
        term += f'+{HTML(self.addend)}'
    return term
def evil(ord):
    if ord<w:
        return str(ord)
    term=''
    if log(ord)==1:
        term+='ω'
    elif ord<W:
        term+=f'ψ({evil(ord.arg)})'
    else:
        term+=f'ψ{sub(evil(ord.sub))}({evil(ord.arg)})'
    if ord.copies>1: term+=f'{ord.copies}'
    if ord.addend>0: term+=f'+{evil(ord.addend)}'
    return term
def eviler(ord):
    if ord==0: return '0'
    if ord==1: return '00ψ'
    x=terms(ord)
    if len(x)==1: return f'{eviler(ord.arg)}{eviler(ord.sub)}ψ'
    return f'{eviler(x[0])}{eviler(sum(x[1:]))}+'
def eviler2(ord):
    if ord==0: return '0'
    if ord==1: return 'ψ{0}(0)'
    x=terms(ord)
    if len(x)==1: return f'ψ{{{eviler2(ord.sub)}}}({eviler2(ord.arg)})'
    return f'{eviler2(x[0])}+{eviler2(sum(x[1:]))}'
def hydra(ord):
    if ord==0:
        return ''
    if ord==1:
        return '()'
    if len(terms(ord))>1:
        return ''.join([hydra(i) for i in terms(ord)])
    if ord<W:
        return f'({hydra(ord.arg)})'
    return f'({ascii(ord.sub)}{hydra(ord.arg)})'
def hydra2(ord):
    if ord==0:
        return ''
    if ord==1:
        return '()'
    if len(terms(ord))>1:
        return ''.join([hydra2(i) for i in terms(ord)])
    if ord<W:
        return f'({hydra2(ord.arg)})'
    if ord<I:
        return f'({hydra2(ord.sub)},{hydra2(ord.arg)})'
    return f'({hydra2(ord.Isub)},{hydra2(ord.sub)},{hydra2(ord.arg)})'
def hydraAlt(ord):
    if ord==0:
        return ''
    if ord==1:
        return ('0+'*ord)[:-1]
    if len(terms(ord)) > 1:
        return '+'.join([hydraAlt(i) for i in terms(ord)])
    if ord.arg==0: return ascii(ord.sub)
    return f'{ascii(ord.sub)}^{"("*bool(len(terms(ord.arg))-1)}{hydraAlt(ord.arg)}{")"*bool(len(terms(ord.arg))-1)}'
def hydraAlt2(ord):
    if ord==0:
        return '0'
    if ord==1:
        return ('(0)^0+'*ord)[:-1]
    if len(terms(ord)) > 1:
        return '+'.join([hydraAlt2(i) for i in terms(ord)])
    if ord.arg==0: return '('+hydraAlt2(ord.sub)+')^(0)'
    return f'({hydraAlt2(ord.sub)})^({hydraAlt2(ord.arg)})'
def hprss(ord):
    if ord<w:
        if ord>100:
            raise Exception
        return [0]*int(ord)
    return [ord.sub]+[1+ord.sub+i for i in hprss(ord.arg)]+hprss(sum(terms(ord)[1:]))
def render(y):
    # [ ] superscript
    # { } subscript
    q=[]
    t=[0]
    for i in y:
        if i in '[}':
            t.append(t[-1]+1)
        if i in ']{':
            t.append(t[-1]-1)
    for i in range(max(t)+abs(min(t))+1):
        q.append([' ']*len(y.replace('[','').replace(']','').replace('{','').replace('}',''))) # weird python list shit
    x=abs(min(t))+1
    j=0
    for i in range(len(y)):
        if y[i] not in '[]{}':
            q[-x][j]=y[i]
            j+=1
        else:
            if y[i] in '[}':
                x+=1
            else:
                x-=1
    return '\n'.join([''.join(i) for i in q])
def prettyprint(ord):
    return render(gamestuffs(ord))
def X(ord):
    a=logW(ord)-W
    b=divW(ord,W+a)
    if ord-W**(W+a)*b>=(W**(W+a)*b)[psi(0,W**(W+a)*b)]:return W**(W+a)*b
    return 0
def repr(ord):
    if type(ord)==int:return str(ord)
    x=f'psi({repr(ord.sub)},{repr(ord.arg)})'if ord>=w**2 else 'w'
    if ord.copies>1:x+=f'*{ord.copies}'
    if ord.addend!=0:x+=f'+{repr(ord.addend)}'
    return x
def voidstyle(self):
    if not isinstance(self,Ordinal):return str(self)
    term = ''

    def k(ord):
        if ord < w + 1:
            return voidstyle(ord)
        if ord.copies == 1 and ord.addend == 0 and ((ord < W and (W * divW(ord.arg, 1) == ord.arg)) or (
                ord >= W and omega(ord.Isub, ord.sub + 1) * div(ord.arg, omega(ord.Isub, ord.sub + 1)) == ord.arg)):
            return voidstyle(ord)
        return f'({voidstyle(ord)})'

    if self.arg == 0:
        term= f'Ω_{k(self.sub)}'if self.sub>1 else'Ω'
    else:
        X = omega(self.Isub, self.sub + 1)
        if veblen_cutoff.lower() != 'max':
            p = eval(veblen_cutoff)
        else:
            p = psi3(self.Isub, self.sub + 1, omega(self.Isub, self.sub + 2))
        if self.sub == 0 and self.Isub == 0 and self.arg < w:
            term = f'ω^{self.arg}'if self.arg>1 else'ω'
        elif terms(self.arg)[-1] >= p or self > W and terms(self.arg)[-1] < X:
            x = terms(self.arg)
            y = 0
            while x[-1] < omega(self.Isub, self.sub + 1):
                y = x[-1] + y
                x = x[:-1]
                if len(x) == 0: break
            if terms(self.arg)[-1] >= omega(self.Isub, self.sub + 1) or w ** psi3(self.Isub, self.sub, sum(x)) != psi3(
                    self.Isub, self.sub, sum(x)):
                term = f'ψ_{voidstyle(omega(0,self.sub+1))}'+'(' + voidstyle(self.arg) + ')'
            else:
                q = log(self)
                q1 = div(q, omega(self.Isub, self.sub))
                if q1 > w:
                    if q1.arg > 0:
                        p = terms(q1.arg)[-1] >= omega(q1.Isub, q1.sub)
                    else:
                        p = 1
                    if q1.addend == 0 and q1.copies == 1 and (
                            (q1.Isub + q1.sub > 0 and p) or (q1.sub == 0 and q1.Isub == 0)):
                        q1 = voidstyle(q1)
                    else:
                        q1 = f'({voidstyle(q1)})'
                else:
                    q1 = voidstyle(q1)
                t = terms(self)
                l = []
                while t[0] >= omega(self.Isub, self.sub) ** div(q, omega(self.Isub, self.sub)):
                    l.append(t[0])
                    t = t[1:]
                    if not t: break
                q2 = div(sum(l), omega(self.Isub, self.sub) ** div(q, omega(self.Isub, self.sub)))
                f=0
                if q2 < w + 1:
                    q2 = str(q2)
                elif (q2.addend > 0 and q2 < W) or (q2.addend > 0 and terms(q2)[-1] < omega(q2.Isub, q2.sub)):
                    q2 = f'({voidstyle(q2)})'
                    f=1
                else:
                    q2 = f'{voidstyle(q2)}'
                if psi3(self.Isub, self.sub, 0) ** div(q, omega(self.Isub, self.sub)) == div(q, omega(self.Isub,
                                                                                                      self.sub)):
                    term = f'{q1}{("*"*(q1[-1]in"0123456789ωΩ")+q2 or f) * (q2 != "1")}'
                else:
                    term = f'{voidstyle(psi3(self.Isub, self.sub, 0))}{f"^{q1}" * (q1 != "1")}'
                    term+=(("*"*(term[-1]in"0123456789ωΩ" or f))+q2) * (q2 != "1")
                if sum(t) != self:
                    if t: term += f'+{voidstyle(sum(t))}'
                    return term
                else:
                    term = f'ψ_{voidstyle(omega(0,self.sub+1))}'+'(' + voidstyle(self.arg) + ')'
        else:
            a1, a2 = phi_inv(self.arg, X, omega(self.Isub, self.sub))
            if a1 == X: a1 = -1
            j = k(a2)
            if a1 == 0 and len(terms(a2)) < 2:
                term = f'ω^{voidstyle(a2)}'
            elif a1 == 0:
                term = f'ω^({voidstyle(a2)})'
            elif a1 <= 3:
                term = ' εζηΓ'[a1] + f'_{k(j)}'
            elif a1 < X:
                term = f'φ({voidstyle(a1)},{voidstyle(a2)})'
            else:
                def f(a, h, str, sp):
                    if a < h: return f'@{str(a)}'
                    return f'@({array(a, h, str, f, sp)})'

                term = f'φ({array(X * a1 + a2, X, voidstyle, f, ",")})'
    if self.copies != 1:
        term += '*'*(voidstyle(psi(self.sub,self.arg))[-1]in'0123456789ω'and self.arg!=1)+str(self.copies)
    if self.addend != 0:
        term += f'+{voidstyle(self.addend)}'
    return term.replace('ψ_Ω(','ψ(')
def as_theta(ord):
    if ord<w:return str(ord)
    if ord==w:return 'ω'
    h=omega(0,ord.sub)
    if ord.sub>0 and multerms(ord.arg)[-1]<omega(0,ord.sub+1) if ord.arg>0 else True:
        x=terms(ord)
        y=[]
        for i in x:
            m=div(log(i),h)
            if y==[]:y+=[[div(log(i),h),div(i,h**div(log(i),h))]]
            elif y[-1][0]!=m:y+=[[div(log(i),h),div(i,h**div(log(i),h))]]
            else:y[-1][1]+=div(i,h**div(log(i),h))
        v=''
        for i in y:
            q1,q2=i
            if q1==0:v+=f'+{as_theta(q2)}';continue
            v+='+Ω'+f'_{ord.sub}'*(ord.sub>1)
            if q1>w:
               if q1.arg>0:p=terms(q1.arg)[-1]>=omega(q1.Isub,q1.sub)
               else:p=1
               if q1.addend==0 and q1.copies==1 and((q1.Isub+q1.sub>0 and p)or(q1.sub==0 and q1.Isub==0)):q1=as_theta(q1)
               else:q1=f'({as_theta(q1)})'
            else:
               q1=as_theta(q1)
            q3=q2
            if q2<w+1:q2=as_theta(q2)
            elif (q2.addend>0 and q2<h)or(q2.addend>0 and terms(q2)[-1]<omega(q2.Isub,q2.sub)**q3):q2=f'({as_theta(q2)})'
            else:q2=f'{as_theta(q2)}'
            if q1!='1':v+=f'^{q1}'
            if q2!='1':v+=f'·{q2}'[(q1=='1' and q2[0]!='(' and h==W):]
        return v[1:]
    m=theta_inv(ord.arg)
    m=W*m[0]+m[1]
    v=f'ϑ({as_theta(m)})'
    if ord.copies>1:v+=f'{ord.copies}'
    if ord.addend>0:v+=f'+{as_theta(ord.addend)}'
    return v
w=Ordinal()
W=omega(0,1)
W_2=omega(0,2)
I=omega(1,0)
L=omega(0,I)
e0=psi(0,W)
z0=psi(0,psi(1,W))
G0=psi(0,psi(1,psi(1,W)))

#
# def f(a, h, str, sp):
#     if a < h: return sub(str(a))
#     return sub('(' + array(a, h, str, f, sp) + ')')
# print(array(W**w,W,HTML,f,','))
# print(theta_inv(W+4,W,w))
# print(psi(0,W**W**2+W**(W*psi(0,W**W**2))))
# print(check(W**(W*psi(0,W**W**2)),psi(0,W**W**2),W))
# print(psi(0,W**W+W**G0))
# print(psi(0,W**(W+1)+W**W))
# print(subterms(omega(0,omega(0,omega(0,2))),W))
# print(as_theta(psi(1,W_2**2)))

# print('ω^(φ(φ(1@(1,0))@(φ(1@(1,0)),φ(1@(1,0))),φ(1@(1,0))@1,φ(1@(1,0))@0)+φ(1@(1,0)))'.replace('φ(1@(1,0))','LVO'))
# print(theta_inv(psi(1,psi(2,psi(3,0)))*2,0))
# print(theta_inv(psi(1,psi(2,0))+1))
# print(phi(W**phi(W**W),phi(W**W)).arg)
# print(phi(G0+1,phi(G0+2,0)).arg)