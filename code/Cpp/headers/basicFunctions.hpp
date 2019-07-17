#ifndef GUARD_basicFunctions_hpp
#define GUARD_basicFunctions_hpp

#include <vector>
#include <string>
#include <cmath>

//#define EXP_VALUES 1000
//extern double exponential[EXP_VALUES];

void createDir(std::string dir_name);
template <class T> T  mean(const std::vector<T>& vec);
template <class T> T  stdev(const std::vector<T>& vec);
//void computeExponentials(double mu);


/* Template functions definition */

template <class T> T mean(const std::vector<T>& vec) {
    return accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
}

template <class T> T  stdev(const std::vector<T>& vec) {
    std::vector<double>::size_type size = vec.size();
    if(size < 2) return 0;
    double vec_mean = mean(vec);
    double s = 0;
    for(double elem : vec) {
        s += (elem - vec_mean)*(elem - vec_mean);
    }
    s = ( 1./(size-1) ) * sqrt(s);
    return s;
}

#endif