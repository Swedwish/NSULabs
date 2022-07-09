#include "Sort.h"

bool cmp(const string& a,const string& b){
    return a<b;
}

bool Sort::exec(vector<string> *input, vector<string> *output) {
    sort((*output).begin(), (*output).end(), greater<>());
    output->clear();
    *output = *input;
    return true;
}