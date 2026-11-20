"""Find a short obstruction to a regular action on class 24."""
import json,subprocess,tempfile
from pathlib import Path
from compare import HERE,read_graph
import igraph as ig
def main():
    rows,_=read_graph(HERE/'classes/class_24.g6')
    g=ig.Graph(n=len(rows),edges=[(i,j) for j,r in enumerate(rows) for i in range(j) if (r>>i)&1])
    gens=g.automorphism_group(sh='fl');perms=','.join('PermList('+str([x+1 for x in p])+')' for p in gens)
    with tempfile.TemporaryDirectory() as d:
        path=Path(d)/'short.g'
        path.write_text('G:=Group(['+perms+']);;\n'+'''H:=Stabilizer(G,1);;
t:=First(Elements(H),x->not IsOne(x));;
Print("group ",Size(G),"  point stabilizer ",Size(H),"  derived subgroup ",Size(DerivedSubgroup(G)),"  stabilizing involution in derived group ",t in DerivedSubgroup(G),"\\n");
x:=First(Elements(G),g->g^2=t);;
if x<>fail then
  Print("square witness order ",Order(x),"  basis images ",List([1,82,163,325,649,2,4,10,28],v->v^x-1),"\\n");
fi;
Print("DONE\\n");QUIT;''')
        out=subprocess.run(['gap','-q','-b',str(path)],check=True,capture_output=True,text=True)
        assert 'DONE' in out.stdout,out.stderr
        print(out.stdout,end='')
if __name__=='__main__': main()
