# Run with the working directory containing the exported exponent lists.
VerifyOnePDS := function(id)
  local H, codes, D, elements, counts, x, y, pos, i, expected;
  H:=SmallGroup(1296,id);;
  codes:=EvalString(StringFile(Concatenation("pds_",String(id),".json")));;
  D:=Set(List(codes,e->PcElementByExponents(Pcgs(H),e)));;
  if Length(D)<>595 or One(H) in D or Set(List(D,Inverse))<>D then
    Error("invalid regular PDS input");
  fi;
  elements:=AsSSortedList(H);;
  counts:=List(elements,x->0);;
  for x in D do
    for y in D do
      pos:=PositionSorted(elements,x/y);
      counts[pos]:=counts[pos]+1;
    od;
  od;
  for i in [1..Length(elements)] do
    if IsOne(elements[i]) then expected:=595;
    elif elements[i] in D then expected:=274;
    else expected:=272; fi;
    if counts[i]<>expected then Error("incorrect difference multiplicity"); fi;
  od;
  Print("VERIFIED ",IdGroup(H)," inverse-closed, identity-free, k=595, all ",Length(D)^2," differences; multiplicities ",Collected(counts),"\n");
end;;
for id in [2797,2848] do VerifyOnePDS(id); od;
QUIT;
