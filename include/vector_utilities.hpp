#pragma once
#include <ostream>
#include <vector>

namespace vectorext {
    //coefficientvector-like multiplication
    template<typename T>
    std::vector<T> operator*(const std::vector<T>& a,const std::vector<T>& b){
        std::vector<T> resultvector(a.size()+b.size()-1);
        for(std::size_t i=0; i<a.size();++i)
            for(std::size_t j=0; j<b.size();++j)
                resultvector[i+j]+=a[i]*b[j];
        return resultvector;
    }

    template<class T>
    void MultiliedByScalar(std::vector<T>& a,const T& num){
        for (std::size_t i=0; i<a.size();++i) a[i]*=num;
    }

    template<typename T>
    std::ostream& operator<<(std::ostream& os, const std::vector<T>& v){
        os << "[";
        for (size_t i = 0; i < v.size(); ++i) {
            os << v[i];
            if (i + 1 < v.size()) os << ", ";
        }
        os << "]";
        return os;
    }
}
