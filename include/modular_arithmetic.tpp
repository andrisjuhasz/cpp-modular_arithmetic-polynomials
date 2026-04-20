template<int P>
IntMod<P>::IntMod(): nMod_(0){}

template<int P>
IntMod<P>::IntMod(int n){
    nMod_ = n % P;
    if (nMod_ < 0)
        nMod_ += P;
}

template<int P>
bool IntMod<P>::operator==(const IntMod<P>& other) const{
    return nMod_==other.nMod_;
}

template<int P>
IntMod<P>& IntMod<P>::operator+=(const IntMod<P>& other){
    nMod_+=other.nMod_;
    nMod_%=P;
    return *this;
}

template<int P>
IntMod<P>& IntMod<P>::operator++(){
    nMod_++;
    nMod_%=P;
    return *this;
}

template <int P>
IntMod<P> IntMod<P>::operator+(const IntMod<P>& other) const{
    IntMod<P> temp=*this;
    temp+=other;
    return temp;
}

template<int P>
IntMod<P>& IntMod<P>::operator-=(const IntMod<P>& other){
    nMod_=nMod_-other.nMod_+P;
    nMod_%=P;
    return *this;
}

template <int P>
IntMod<P> IntMod<P>::operator-(const IntMod<P>& other) const{
    int result = nMod_ - other.nMod_;
    if (result < 0)
        result += P;
    return IntMod(result);
}

template <int P>
IntMod<P> IntMod<P>::operator-() const{
    return nMod_ == 0 ? *this : IntMod(P - nMod_);
}

template <int P>
IntMod<P>& IntMod<P>::operator*=(const IntMod<P>& other){
    nMod_*=other.nMod_;
    nMod_%=P;
    return *this;
}

template <int P>
IntMod<P> IntMod<P>::operator*(const IntMod<P>& other) const{
    IntMod<P> temp=*this;
    temp*=other;
    return temp;
}

template <int P>
IntMod<P> IntMod<P>::operator*(int num) const{
    return IntMod<P>{(nMod_*num)%P};
}

template<int P>
IntMod<P> IntMod<P>::FastExponentiation(int exp) const{
    int result = 1;
    int base = nMod_;
    while (exp > 0) {
        if (exp & 1)
            result = (result * base) % P;
        base = (base * base) % P;
        exp >>= 1;
    }
    return IntMod<P>(result);
}

//Fermat's Little Theorem
template<int P>
IntMod<P> IntMod<P>::Inverse() const{
    return FastExponentiation(P-2);
}

