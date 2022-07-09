#ifndef UNTITLED_CSVPARSER_H
#define UNTITLED_CSVPARSER_H
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iosfwd>
using namespace std;
template<class Arg_One, class... others>
Arg_One& return_type(int num = 0) {
    if (num == 0) {
        return new Arg_One;
    }
    return return_type<others...>(num-1);
}
template<class ...Args>
class CSVparser;

template <class ...Args>
class Iterator{
private:
    ifstream &IT_stream;
    int is_end = 0;
    tuple<Args...> data;
public:
    friend ostream& operator<< (ostream& off, Iterator<Args...>& it);

    Iterator(int is_end, ifstream& off, int fk) : is_end(is_end), IT_stream(off){};

    Iterator(int is_end, ifstream& str):is_end(is_end),IT_stream(str){
        data = mk_tuple<Args...>(IT_stream);
    };

    template<class Arg1>
    tuple<Arg1> mk_tuple(ifstream& stream){
        char c;
        string str;
        stream.get(c);
        while (c != '\n' && c!=  EOF){
            str.push_back(c);
            c = EOF;
            stream.get(c);
        };
        istringstream ss(str);
        Arg1 a;
        ss>>a;
        tuple<Arg1> k = make_tuple(a);
        if (c == EOF) is_end = 2;
        return k;
    }

    template<typename first, typename second, typename...Arguments>
    tuple<first,second, Arguments...> mk_tuple(ifstream& stream){
        char c;
        string str;
        stream.get(c);
        while (c != ','){
            str.push_back(c);
            stream.get(c);
        };
        istringstream ss(str);
        first a;
        getline(ss,a);
        tuple<first> k = make_tuple(a);
        tuple<first,second, Arguments...> f = tuple_cat(k, mk_tuple < second, Arguments... > (stream));
        return f;
    }

    void operator++(){
        if (this->is_end == 2)
            this->is_end = 1;
        else
            this->data = mk_tuple<Args...>(IT_stream);
    }

    tuple<Args...> operator *(){
        return data;
    }

    bool operator!=(Iterator<Args...>& f){
        return (this->is_end != f.is_end);
    }
};
template <class ...Args>
ostream& operator<< (ostream& off, Iterator<Args...>& it){
    off<< it->data;
    return off;
}

template<class ...Args>
class CSVparser {
private:
    ifstream& m_stream;
    tuple<Args...> data;
public:
    CSVparser(ifstream& stream, int skip = 0): m_stream(stream){
        string str;
        for (int i = 0; i < skip; i++)
            getline(stream,str);
    }

    CSVparser(ifstream& stream, tuple<Args...> a): m_stream(stream), data(a) {};

    Iterator<Args...>& end(){
        auto a = new Iterator<Args...>(1,m_stream,0);
        return *a;
    }

    Iterator<Args...>& begin(){
        auto a = new Iterator<Args...>(0,m_stream);
        return *a;
    }
};

#endif //UNTITLED_CSVPARSER_H
