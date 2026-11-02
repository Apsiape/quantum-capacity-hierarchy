import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[k]='1'
import itertools as it,json,time,sys,hashlib
from pathlib import Path
from fractions import Fraction as F
OUT=Path(__file__).resolve().parent
DATA=OUT
contexts=[tuple(range(3*r,3*r+3)) for r in range(3)]+[tuple(range(c,9,3)) for c in range(3)]
V=[(q,b) for q in range(6) for b in it.product((0,1),repeat=3) if sum(b)%2==int(q==5)]
labels=[dict(zip(contexts[q],b)) for q,b in V]
def edge(u,v):
 return u!=v and (V[u][0]==V[v][0] or any(labels[u][x]!=labels[v][x] for x in set(labels[u])&set(labels[v])))
A=[[edge(u,v) for v in range(24)] for u in range(24)]
def finite():
 import sympy as s
 import networkx as nx
 G=nx.Graph(); G.add_nodes_from(range(24)); G.add_edges_from((u,v) for u in range(24) for v in range(u) if A[u][v])
 cliques=list(nx.find_cliques(G)); assert max(map(len,cliques))==4
 ind=list(nx.find_cliques(nx.complement(G))); assert max(map(len,ind))==5
 print('Fresh graph: 24 vertices, 108 edges, regular degree 9, omega=4, alpha=5.',flush=True)
 cert=json.loads((DATA/'dimension-dual.json').read_text()); core=[V.index((q,tuple(b))) for q,b in cert['core_vertices']]
 assert len(core)==len(set(core))==10
 T=(1<<10)-1; acc={}; rhs=F(0)
 def add(row,w):
  for z,c in row.items(): acc[z]=acc.get(z,F(0))+w*c
 for item in cert['inequalities']:
  w=F(item['weight']); assert w>=0; S=item['S']; assert 0<=S<=T
  if item['type']=='mono':
   i=item['index']; assert 0<=i<10 and not (S>>i)&1
   row={S:-1,S|1<<i:1}
  else:
   assert item['type']=='submod'; i,j=item['index']; assert 0<=i<10 and 0<=j<10 and i!=j and not (S>>i)&1 and not (S>>j)&1
   row={S:-1,S|1<<i|1<<j:-1,S|1<<i:1,S|1<<j:1}
  assert {int(k):-v for k,v in item['row'].items()}==row
  add(row,w)
 for item in cert['equalities']:
  i=item['index']; w=F(item['weight']); typ=item['type']
  if typ=='empty': row={0:1}; target=0
  elif typ=='singleton': row={1<<i:1}; target=1
  else:
   assert typ=='independence'; S=sum(1<<j for j in range(10) if j!=i and not A[core[i]][core[j]])
   row={S|1<<i:1,S:-1}; target=1
  assert row=={int(k):v for k,v in item['row'].items()} and target==item['rhs']; add(row,w); rhs+=w*target
 assert {z:v for z,v in acc.items() if v}=={T:F(1)} and rhs==6
 B=G.subgraph(core); isolated=[v for v in core if B.degree(v)==0]; H=B.subgraph([v for v in core if B.degree(v)])
 mobius=nx.cycle_graph(8); mobius.add_edges_from((i,i+4) for i in range(4)); assert len(isolated)==2 and nx.is_isomorphic(H,mobius)
 print('Exact all-field dual: 33 nonnegative inequality rows + 10 equality rows, r(core)-6; induced core 2K1 + Mobius ladder.',flush=True)
 X=s.Matrix([[0,1],[1,0]]); Y=s.Matrix([[0,-s.I],[s.I,0]]); Z=s.diag(1,-1); I=s.eye(2); I4=s.eye(4)
 kron=s.kronecker_product
 obs=[kron(X,I),kron(I,X),kron(X,X),kron(I,Z),kron(Z,I),kron(Z,Z),kron(X,Z),kron(Z,X),kron(Y,Y)]
 P=[]
 for q,bits in V:
  p=I4
  for cell,bit in zip(contexts[q],bits): p=p*(I4+(-1)**bit*obs[cell])/2
  assert p==p.conjugate().T and p*p==p and s.trace(p)==1; P.append(p)
 for q in range(6): assert sum((P[u] for u in range(24) if V[u][0]==q),s.zeros(4))==I4
 for u,v in G.edges: assert P[u]*P[v]==s.zeros(4)
 print('Exact rank-one Pauli projections give six-message projective code; six context cliques bound theta and fractional clique cover by 6; omega=4 gives fractional-cover lower bound 6.',flush=True)
 ms=json.loads((DATA/'lower29.json').read_text()); words=[tuple(w) for w in ms['code']]
 assert len(words)==len(set(words))==29
 assert all(0<=u<24 and 0<=v<24 for u,v in words)
 assert all(not ((a==c or A[a][c]) and (b==d or A[b][d])) for (a,b),(c,d) in it.combinations(words,2))
 print('Lower29: all 406 pairs verified from fresh graph.',flush=True)
 return {'graph_edges':G.number_of_edges(),'alpha':5,'omega':4,'dimension_bound':6,'quantum_code':6,'lower_square':29}
def symmetry():
 lookup={v:i for i,v in enumerate(V)}; qgroup=set(); mappings=[]
 for rp in it.permutations(range(3)):
  for cp in it.permutations(range(3)):
   for tr in (0,1):
    cellmap={3*r+c:3*cp[c]+rp[r] if tr else 3*rp[r]+cp[c] for r in range(3) for c in range(3)}
    qp=tuple(next(i for i,Q in enumerate(contexts) if set(Q)=={cellmap[c] for c in C}) for C in contexts)
    for flips in it.product((0,1),repeat=9):
     if any(sum(flips[c] for c in contexts[qp[q]])%2 != ((q==5)^(qp[q]==5)) for q in range(6)): continue
     p=[]
     for q,bits in V:
      assignment={cellmap[c]:b^flips[cellmap[c]] for c,b in zip(contexts[q],bits)}
      p.append(lookup[(qp[q],tuple(assignment[c] for c in contexts[qp[q]]))])
     assert len(set(p))==24 and all(A[u][v]==A[p[u]][p[v]] for u in range(24) for v in range(24)); mappings.append(p); qgroup.add(qp)
 assert len({tuple(p) for p in mappings})==1152 and len(qgroup)==72
 kernel=[p for p in mappings if all(V[p[u]][0]==V[u][0] for u in range(24))]
 assert len(kernel)==16 and all(len({p[4*q] for p in kernel})==4 for q in range(6))
 todo=set(it.permutations(range(6))); reps=[]
 while todo:
  sig=min(todo); orbit=set()
  for a in qgroup:
   ainv=[a.index(q) for q in range(6)]
   for b in qgroup: orbit.add(tuple(b[sig[ainv[q]]] for q in range(6)))
  assert orbit<=set(it.permutations(range(6))); todo-=orbit; reps.append((sig,len(orbit)))
 assert [n for _,n in reps]==[72,648]
 print('All 1152 explicit automorphisms checked; vertex transitive; 16 kernel maps answer-transitive; 720 empty patterns covered by double cosets 72+648.',flush=True)
 return reps
def sat(which):
 from pysat.formula import CNF
 from pysat.card import CardEnc,EncType
 from pysat.solvers import Solver
 reps=symmetry(); sig,orb=reps[which]; clauses=[]
 for x,y in it.combinations(range(576),2):
  a,b=divmod(x,24); c,d=divmod(y,24)
  if (a==c or A[a][c]) and (b==d or A[b][d]): clauses.append([-x-1,-y-1])
 assert len(clauses)==28512
 for x in range(576):
  a,b=divmod(x,24)
  if V[b][0]==sig[V[a][0]]: clauses.append([-x-1])
 top=576
 for coord in (0,1):
  for q in range(6):
   lits=[x+1 for x in range(576) if V[divmod(x,24)[coord]][0]==q]
   enc=CardEnc.equals(lits,bound=5,top_id=top,encoding=EncType.seqcounter); top=enc.nv; clauses.extend(enc.clauses)
 assert sig[0]!=1; clauses.append([5])
 tag=''.join(map(str,sig)); path=OUT/f'independent-upper30-{tag}.cnf'
 with path.open('w') as f:
  f.write(f'p cnf {top} {len(clauses)}\n')
  for c in clauses: f.write(' '.join(map(str,c))+' 0\n')
 start=time.monotonic()
 with Solver(name='lingeling',bootstrap_with=clauses,with_proof=True) as solver:
  result=solver.solve(); assert result is False
  proof=solver.get_proof(); (OUT/f'independent-upper30-{tag}.drat').write_text('\n'.join(proof)+'\n')
  info={'solver':'Lingeling','result':'UNSAT','seconds':time.monotonic()-start,'proof_lines':len(proof),'variables':top,'clauses':len(clauses),'orbit':orb,'cnf_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
  (OUT/f'independent-upper30-{tag}.json').write_text(json.dumps(info,indent=2)); print(info,flush=True)
  os._exit(0)
if __name__=='__main__':
 if sys.argv[1]=='finite': result=finite(); result['symmetry']=symmetry(); (OUT/'independent-finite.json').write_text(json.dumps(result,indent=2))
 else: sat(int(sys.argv[1]))
 sys.stdout.flush(); os._exit(0)
