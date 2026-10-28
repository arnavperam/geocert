// Streaming marginal bounds for independent unit-RNA candidate assignments.
// This accelerates accumulation AFTER candidate sets exist. It does not solve
// coherent geometry, extract candidates from images, or certify anatomical truth.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
struct Totals { std::uint64_t source=0, removable=0, incoming=0; };
int main(int argc,char** argv) {
    try {
        if(argc!=4) throw std::runtime_error("usage: count_bounds input.tsv output.tsv RNA_budget_or_all");
        std::ifstream in(argv[1]); std::ofstream out(argv[2]);
        if(!in||!out) throw std::runtime_error("cannot open input/output");
        std::uint64_t K=std::string(argv[3])=="all"?std::numeric_limits<std::uint64_t>::max():std::stoull(argv[3]);
        std::map<std::pair<std::int64_t,std::uint32_t>,Totals> table;
        std::uint32_t gene;std::int64_t source,cell;std::size_t count;std::uint64_t records=0;
        while(in>>gene>>source>>count) {
            if(count==0||count>10000000)throw std::runtime_error("invalid candidate count");
            std::vector<std::int64_t> candidates; candidates.reserve(count);
            for(std::size_t j=0;j<count;++j){if(!(in>>cell))throw std::runtime_error("truncated record");candidates.push_back(cell);}
            std::sort(candidates.begin(),candidates.end());candidates.erase(std::unique(candidates.begin(),candidates.end()),candidates.end());
            if(!std::binary_search(candidates.begin(),candidates.end(),source))throw std::runtime_error("source not a candidate");
            auto &original=table[{source,gene}];++original.source;
            if(candidates.size()>1)++original.removable;
            for(auto c:candidates)if(c!=source)++table[{c,gene}].incoming;
            ++records;
        }
        if(!in.eof())throw std::runtime_error("invalid input token");
        out<<"cell\tgene\tsource\tlower\tupper\n";
        for(auto const &[key,v]:table)out<<key.first<<'\t'<<key.second<<'\t'<<v.source<<'\t'<<(v.source-std::min(K,v.removable))<<'\t'<<(v.source+std::min(K,v.incoming))<<'\n';
        if(!out)throw std::runtime_error("output write failed");
        std::cerr<<"records="<<records<<" entries="<<table.size()<<"\n";
    } catch(std::exception const& e){std::cerr<<e.what()<<"\n";return 2;}
    return 0;
}
