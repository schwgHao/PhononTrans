#include <iostream>
#include <vector>
#include <numeric>

using namespace std;
int main()
{
     const size_t nscz = 3;
     const size_t na = 6;
     vector<vector<int> > res(2, vector<int>(na));

     if(nscz == 3){
        iota(std::begin(res[0]), std::end(res[0]), na);
        iota(std::begin(res[1]), std::end(res[1]), 2*na);
     }
     else if(nscz == 1){
        res[0].resize(na/2);
        res[1].resize(na/2);
        iota(std::begin(res[0]), std::begin(res[0]) + na/2, 0);
        iota(std::begin(res[1]), std::begin(res[1]) + na/2, na/2);
     }
     else{
        std::cerr << "Wrong number of super cells." << std::endl;
        abort();
     }
     for(auto i : res){
         for(auto j : i)cout << j << " ";
         cout << endl;
     }
}

