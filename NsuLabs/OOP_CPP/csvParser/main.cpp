#pragma once
#include <tuple>
#include <functional>
#include <iostream>
#include <iosfwd>
#include <vector>
#include <fstream>
#include "CSVparser.cpp"

using namespace std;
#include "print_tuple.h"

int main(){
//    tuple<int, string, string> a = {3,"hello", "world"};
//    cout<<a;
//    string b = "hello";
//    auto c = make_tuple(b);
//    tuple<int> d = {3};
//    tuple<string,int> e = tuple_cat(c,d);
//    cout<<endl<<e;

    ifstream file("test.csv");
    CSVparser<string, string,int > parser(file, 0 /*skip first lines count*/);
    for (tuple<string,string,int> rs : parser){
        cout<<rs<<endl;
    }

}
