// Complete local projective moment-code catalogs, without fixed matching channels.
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <unordered_map>
using namespace std;
struct Point { array<int,4> x; array<int,10> moment; };
struct State { int weight=0; array<int,10> moment{}; };
uint32_t key(const State& s) {
    uint32_t v=0,p=1;for(int a:s.moment){v+=p*a;p*=3;}
    return v|(uint32_t(s.weight)<<16);
}
void change(State& s,const Point& p,int sign) {
    s.weight+=sign;
    for(int j=0;j<10;j++)s.moment[j]=(s.moment[j]+sign*p.moment[j]+3)%3;
}
int main(int argc,char** argv) {
    if(argc!=3)return 2;
    ifstream in(argv[1]);array<Point,40> points;array<array<int,40>,3> bounds;
    for(auto& p:points){for(int& x:p.x)in>>x;int k=0;
        for(int i=0;i<4;i++)for(int j=i;j<4;j++)p.moment[k++]=p.x[i]*p.x[j]%3;}
    for(auto& b:bounds)for(int& x:b)in>>x;
    if(!in)return 3;
    array<uint64_t,40> planes{};
    for(int i=0;i<40;i++)for(int j=0;j<40;j++) {
        int dot=0;for(int k=0;k<4;k++)dot+=points[i].x[k]*points[j].x[k];
        if(dot%3==0)planes[i]|=uint64_t(1)<<j;
    }
    auto started=chrono::steady_clock::now();unordered_multimap<uint32_t,uint32_t> left;
    left.max_load_factor(.8);left.reserve(1<<20);
    State state;uint32_t mask=0;
    for(uint32_t t=0;t<(1u<<20);t++) {
        if(t){int j=__builtin_ctz(t);mask^=1u<<j;change(state,points[j],(mask>>j&1)?1:-1);}
        left.emplace(key(state),mask);
    }
    array<ofstream,3> output;for(int c=0;c<3;c++)output[c].open(string(argv[2])+"_"+to_string(c)+".tsv");
    array<uint64_t,2> raw{};array<array<uint64_t,2>,3> kept{};
    State right;mask=0;
    for(uint32_t t=0;t<(1u<<20);t++) {
        if(t){int j=__builtin_ctz(t);mask^=1u<<j;change(right,points[20+j],(mask>>j&1)?1:-1);}
        for(int type=0;type<2;type++) {
            int weight=type?18:13,shift=type?6:4;
            State want;want.weight=weight-right.weight;if(want.weight<0||want.weight>20)continue;
            for(int j=0;j<10;j++)want.moment[j]=(3-right.moment[j])%3;
            auto range=left.equal_range(key(want));
            for(auto it=range.first;it!=range.second;++it) {
                raw[type]++;uint64_t whole=(uint64_t(mask)<<20)|it->second;
                array<bool,3> good{true,true,true};
                for(int q=0;q<40;q++) {
                    int value=__builtin_popcountll(whole&planes[q])-shift;
                    if(value%3)return 4;value/=3;
                    for(int c=0;c<3;c++)if(value < -bounds[c][q] || value>bounds[c][q])good[c]=false;
                    if(!good[0]&&!good[1]&&!good[2])break;
                }
                for(int c=0;c<3;c++)if(good[c]) {
                    kept[c][type]++;output[c]<<(type?36:27)<<'\t'<<hex<<whole<<dec<<'\n';
                }
            }
        }
    }
    cout<<"COMPLETE seconds="<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<'\n';
    for(int t=0;t<2;t++)cout<<"raw size="<<(t?36:27)<<" count="<<raw[t]<<'\n';
    for(int c=0;c<3;c++)for(int t=0;t<2;t++)cout<<"case="<<c<<" size="<<(t?36:27)<<" kept="<<kept[c][t]<<'\n';
}
