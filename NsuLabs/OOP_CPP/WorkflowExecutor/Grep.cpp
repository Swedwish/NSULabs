#include "Grep.h"

Grep::Grep(string word) {
    args.push_back(word);
}

bool Grep::exec(vector<string> *input, vector<string> *output) {
    for (auto &i : *input){
        if(i.find(args[0])!=string::npos){
            output->push_back(i);
        }
    }
    return true;
}

