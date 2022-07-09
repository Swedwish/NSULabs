#pragma once
#include <tuple>
#include <functional>
#include <iostream>
using namespace std;

template<const int index, typename ...Args>
class iterate {
public:
    static void next(tuple<Args...> &Tuple, ostream &stream) {
        iterate<index-1, Args...>::next(Tuple, stream);
        stream<< get<index>(Tuple)<<" ";
    }
};
template<typename Ch, typename Tr, typename ... Args>
auto &operator<<(std::basic_ostream<Ch, Tr> &os, std::tuple<Args...> &t)
{
    int const length = tuple_size<tuple<Args...>>::value;
    iterate<length-1, Args...>::next(t,os);
    return os;
}
template<typename ...Args>
class iterate<0, Args...> {
public:
    static void next(tuple<Args...> &Tuple, ostream& stream) {
        stream<<get<0>(Tuple)<<" ";
    }
};