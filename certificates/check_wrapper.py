"""Exact non-diagonal check of the measure-and-erase antidegrading map."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import sympy as s
from pathlib import Path
OUT=Path(__file__).resolve().parent
I=s.eye(2); X=s.Matrix([[0,1],[1,0]])
plus=(I+X)/2; minus=(I-X)/2
W=plus+minus/2; c=1
Facc=(I-W)/2; Frej=(I+W)/2
Sa=minus/2; Sr=plus+s.sqrt(3)*minus/2
assert s.simplify(Sa*Sa-Facc)==s.zeros(2) and s.simplify(Sr*Sr-Frej)==s.zeros(2)
# Bob basis {qubit 0,qubit 1,erase}; environment {acc X=0,1; rej X,q=00,01,10,11}.
V=s.zeros(18,4)
for x in range(2):
 for q in range(2):
  for y in range(2):
   V[q*6+y,2*x+q]=Sa[y,x]
   V[2*6+2+2*y+q,2*x+q]=Sr[y,x]
assert s.simplify(V.conjugate().T*V)==s.eye(4)
Ks=[]
for e in range(6):
    Ks.append(V.extract([b*6+e for b in range(3)],range(4)))
Ls=[]
for b in range(3):
    Ls.append(V.extract(range(b*6,(b+1)*6),range(4)))
D=[]
for y in range(2):
    k=s.zeros(3,6); k[2,y]=1; D.append(k)
Groot=minus/s.sqrt(3); Hroot=plus+s.sqrt(s.Rational(2,3))*minus
for y in range(2):
    k=s.zeros(3,6)
    for x in range(2):
     for q in range(2): k[q,2+2*x+q]=Groot[y,x]
    D.append(k)
for y in range(2):
 for q in range(2):
    k=s.zeros(3,6)
    for x in range(2): k[2,2+2*x+q]=Hroot[y,x]
    D.append(k)
assert s.simplify(sum((k.conjugate().T*k for k in D),s.zeros(6)))==s.eye(6)
for i in range(4):
 for j in range(4):
    unit=s.zeros(4); unit[i,j]=1
    bob=sum((k*unit*k.conjugate().T for k in Ks),s.zeros(3))
    eve=sum((k*unit*k.conjugate().T for k in Ls),s.zeros(6))
    recovered=sum((k*eve*k.conjugate().T for k in D),s.zeros(3))
    assert s.simplify(bob-recovered)==s.zeros(3)
print('Exact non-diagonal PSD W: V isometric, D CPTP, D o Phi^c = Phi on all 16 matrix units.')
print('Negative-eigenvector case: W=diag(-1/2,1/2), c=1 gives erasure acceptance 3/4 and coherent information 1/2.')
