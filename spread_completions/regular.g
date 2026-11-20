# G acts on the original graph6 labels, augmented by one for GAP.
iso:=IsomorphismPcGroup(G);;
P:=Image(iso);;
if not Size(P) in [2592,5184] then Error("wrong full group order"); fi;
front:=[P];;
candidates:=[];;
while Length(front)>0 do
  H:=Remove(front);
  for K in MaximalSubgroupClassReps(H) do
    if Size(K)=1296 then
      Add(candidates,K);
    elif Size(K)=2592 then
      Add(front,K);
    fi;
  od;
od;
seen:=[];; types:=[];;
for K in candidates do
  if not ForAny(seen,L->IsConjugate(P,L,K)) then
    Add(seen,K);
    S:=PreImage(iso,K);
    if IsTransitive(S,[1..1296]) then
      id:=IdGroup(K);
      if not id in types then
        model:=SmallGroup(id[1],id[2]);;
        transport:=IsomorphismGroups(model,K);;
        if transport=fail then Error("no isomorphism to library group"); fi;
        D:=Filtered(Elements(model), d->1^PreImagesRepresentative(iso,Image(transport,d)) in N0);;
        exps:=List(D,d->ExponentsOfPcElement(Pcgs(model),d));;
        if Length(D)<>595 then Error("incorrect connection set"); fi;
        file:=Concatenation(outputDir,"/pds_",String(id[2]),".json");;
        PrintTo(file,exps,"\n");
        AddSet(types,id);
      fi;
      Print("regular ",id,"  abelian=",IsAbelian(K),"  orbit=",Length(Orbit(S,1)),"  stabilizer=",Size(Stabilizer(S,1)),"\n");
    fi;
  fi;
od;
Print("FINISHED candidate subgroup classes ",Length(seen),"  regular types ",types,"\n");
QUIT;
