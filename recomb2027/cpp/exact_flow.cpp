// Exact-integer Dinic reference backend. This is standard maximum flow, not a
// new flow algorithm or a claim of whole-slide runtime. Arbitrary precision
// avoids silently rounding the rational compiler's scaled capacities.
#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <vector>
using Integer=boost::multiprecision::cpp_int;
struct Edge{int to,reverse;Integer residual;};
struct Original{int from,to,index;Integer capacity;};
struct Network {
 std::vector<std::vector<Edge>> adj;
 std::vector<int> level,next;
 explicit Network(int n):adj(n),level(n),next(n){}
 int add(int u,int v,Integer c){
  if(u==v)return -1;
  int i=static_cast<int>(adj[u].size()),j=static_cast<int>(adj[v].size());
  adj[u].push_back({v,j,c});adj[v].push_back({u,i,0});return i;
 }
 bool bfs(int source,int sink){
  std::fill(level.begin(),level.end(),-1);level[source]=0;std::queue<int> q;q.push(source);
  while(!q.empty()){int u=q.front();q.pop();for(auto const&e:adj[u])if(e.residual>0&&level[e.to]<0){level[e.to]=level[u]+1;q.push(e.to);}}
  return level[sink]>=0;
 }
 Integer send(int u,int sink,Integer amount){
  if(u==sink)return amount;
  for(int &i=next[u];i<static_cast<int>(adj[u].size());++i){
   Edge &e=adj[u][i];if(e.residual<=0||level[e.to]!=level[u]+1)continue;
   Integer cap=amount<e.residual?amount:e.residual;
   Integer f=send(e.to,sink,cap);
   if(f>0){e.residual-=f;adj[e.to][e.reverse].residual+=f;return f;}
  }
  return 0;
 }
 Integer flow(int source,int sink,Integer upper){
  Integer total=0;
  while(bfs(source,sink)){std::fill(next.begin(),next.end(),0);while(true){Integer f=send(source,sink,upper);if(f==0)break;total+=f;}}
  return total;
 }
};
int main(int argc,char**argv){
 try{
  if(argc!=3)throw std::runtime_error("usage: exact_flow input.graph output.flow");
  std::ifstream in(argv[1]);std::ofstream out(argv[2]);if(!in||!out)throw std::runtime_error("cannot open input/output");
  int n,m,s,t;if(!(in>>n>>m>>s>>t)||n<2||n>200000||m<0||m>10000000||s<0||s>=n||t<0||t>=n||s==t)throw std::runtime_error("invalid graph header");
  Network net(n);std::vector<Original> original;original.reserve(m);Integer upper=1;
  for(int k=0;k<m;++k){int u,v;Integer c;if(!(in>>u>>v>>c)||u<0||v<0||u>=n||v>=n||c<0)throw std::runtime_error("invalid edge");original.push_back({u,v,net.add(u,v,c),c});upper+=c;}
  Integer value=net.flow(s,t,upper);std::vector<int> side;for(int i=0;i<n;++i)if(net.level[i]>=0)side.push_back(i);
  out<<value<<'\n'<<side.size();for(int i:side)out<<' '<<i;out<<'\n';
  for(auto const&e:original)out<<(e.index<0?Integer(0):e.capacity-net.adj[e.from][e.index].residual)<<'\n';
  if(!out)throw std::runtime_error("write failed");
 }catch(std::exception const&e){std::cerr<<e.what()<<'\n';return 2;}
 return 0;
}
