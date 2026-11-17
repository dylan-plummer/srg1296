// Exact two-eigenvalue audit of unions of saved connection-set pieces.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>
int main(){
 int n,m;std::cin>>n>>m;if(!std::cin||n>24)return 2;
 std::vector<std::vector<int>> parts(n,std::vector<int>(m));
 for(auto& p:parts)for(int& x:p)std::cin>>x;
 std::vector<int> eigen(m);uint32_t mask=0;
 for(uint32_t t=1;t<(1u<<n);t++){
  int j=__builtin_ctz(t);mask^=1u<<j;int sign=(mask>>j&1)?1:-1;
  for(int x=0;x<m;x++)eigen[x]+=sign*parts[j][x];
  int k=eigen[0];if(k<1||k>647)continue;
  int a=eigen[1],b=a;bool ok=true;
  for(int x=2;x<m;x++)if(eigen[x]!=a){if(b==a)b=eigen[x];else if(eigen[x]!=b){ok=false;break;}}
  if(!ok||a==b)continue;int r=std::max(a,b),s=std::min(a,b);
  int mu=k+r*s,lambda=mu+r+s;
  if(mu<1||lambda<0||r>=k||k*(k-lambda-1)!=(1295-k)*mu)continue;
  std::cout<<mask<<' '<<k<<' '<<lambda<<' '<<mu<<'\n';
 }
}
