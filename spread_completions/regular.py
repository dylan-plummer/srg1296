"""Test all possible regular subgroups of the transitive mixed-spread graph."""
from pathlib import Path
import argparse,json,subprocess,tempfile
import numpy as np
import igraph as ig
from compare import HERE, read_graph

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('class_number',type=int,choices=(18,24))
    n=ap.parse_args().class_number
    p=HERE/'classes'/f'class_{n:02}.g6'
    rows,_=read_graph(p)
    edges=[(i,j) for j,r in enumerate(rows) for i in range(j) if (r>>i)&1]
    g=ig.Graph(n=len(rows),edges=edges)
    generators=g.automorphism_group(sh='fl')
    assert g.count_automorphisms(sh='fl')=={18:5184,24:2592}[n]
    with tempfile.TemporaryDirectory() as d:
        inp=Path(d)/'test.g'
        perms=',\n'.join('PermList('+str([x+1 for x in gen])+')' for gen in generators)
        neighbours=[i+1 for i in range(1296) if (rows[0]>>i)&1]
        header='G:=Group(['+perms+']);;\nN0:='+str(neighbours)+';;\n'
        header+='outputDir:='+json.dumps(str(HERE))+';;\n'
        inp.write_text(header+(HERE/'regular.g').read_text())
        run=subprocess.run(['gap','-q','-b',str(inp)],text=True,capture_output=True,check=True)
        if 'FINISHED' not in run.stdout or 'Error' in run.stderr:
            raise RuntimeError(run.stderr+'\n'+run.stdout)
        print(run.stdout,end='')
if __name__=='__main__':main()
