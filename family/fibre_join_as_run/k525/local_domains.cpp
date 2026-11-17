// Complete binary moment-code domains by a 20+20 meet in the middle.
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <chrono>
using namespace std;
struct Point{array<int,4> x;array<int,5> hit;array<int,10> moment;int diagonal;};
struct State{int weight=0;array<int,5> hit{};array<int,10> moment{};};
uint64_t key(const State& s){
    uint64_t k=0,p=1;for(int x:s.moment){k+=p*x;p*=3;}
    k|=uint64_t(s.weight)<<16;for(int j=0;j<5;j++)k|=uint64_t(s.hit[j])<<(21+5*j);
    return k;
}
void change(State& s,const Point& p,int sign){
    s.weight+=sign;for(int j=0;j<5;j++)s.hit[j]+=sign*p.hit[j];
    for(int j=0;j<10;j++){s.moment[j]+=sign*p.moment[j];if(s.moment[j]<0)s.moment[j]+=3;if(s.moment[j]>=3)s.moment[j]-=3;}
}
int main(int argc,char**argv){
    if(argc!=3)return 2;ifstream in(argv[1]);int n;if(!(in>>n)||n!=40)return 3;
    array<Point,40> p;
    for(auto& q:p){for(int& a:q.x)in>>a;for(int& a:q.hit)in>>a;in>>q.diagonal;
        int k=0;for(int i=0;i<4;i++)for(int j=i;j<4;j++)q.moment[k++]=q.x[i]*q.x[j]%3;}
    if(!in)return 4;
    array<uint64_t,40> planes{};
    for(int i=0;i<40;i++)for(int j=0;j<40;j++){
        int dot=0;for(int k=0;k<4;k++)dot+=p[i].x[k]*p[j].x[k];
        if(dot%3==0)planes[i]|=uint64_t(1)<<j;
    }
    auto started=chrono::steady_clock::now();unordered_multimap<uint64_t,uint32_t> left;left.max_load_factor(.8);left.reserve(1<<20);
    State s;uint32_t mask=0;
    for(uint32_t t=0;t<(1u<<20);t++){
        if(t){int j=__builtin_ctz(t);mask^=1u<<j;change(s,p[j],mask>>j&1?1:-1);}
        left.emplace(key(s),mask);
    }
    array<uint64_t,6> raw{},kept{};ofstream out(argv[2]);State r;mask=0;
    for(uint32_t t=0;t<(1u<<20);t++){
        if(t){int j=__builtin_ctz(t);mask^=1u<<j;change(r,p[20+j],mask>>j&1?1:-1);}
        for(int type=0;type<6;type++){
            const bool large=type==5;const int weight=large?18:13,shift=large?6:4;
            State want;want.weight=weight-r.weight;if(want.weight<0||want.weight>20)continue;
            bool valid=true;
            for(int j=0;j<5;j++){want.hit[j]=(large?6:(j==type?10:4))-r.hit[j];if(want.hit[j]<0||want.hit[j]>20)valid=false;}
            if(!valid)continue;for(int j=0;j<10;j++)want.moment[j]=r.moment[j]?3-r.moment[j]:0;
            auto range=left.equal_range(key(want));
            for(auto it=range.first;it!=range.second;++it){
                raw[type]++;uint64_t whole=(uint64_t(mask)<<20)|it->second;bool good=true;
                for(int j=0;j<40;j++){
                    int value=__builtin_popcountll(whole&planes[j])-shift;
                    if(value%3){cerr<<"moment error\n";return 5;}
                    value/=3;int bound=p[j].diagonal==-6?1:2;if(value < -bound||value>bound){good=false;break;}
                }
                if(good){kept[type]++;out<<(large?36:27)<<'\t'<<(large?-1:type)<<'\t'<<hex<<whole<<dec<<'\n';}
            }
        }
    }
    out.close();cout<<"COMPLETE left_subsets="<<(1<<20)<<" right_subsets="<<(1<<20)<<" seconds="<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<'\n';
    for(int i=0;i<6;i++)cout<<"type="<<i<<" raw="<<raw[i]<<" kept="<<kept[i]<<'\n';
}
