//
// Created by Admin on 05.12.2021.
//

#include "Replace.h"

Replace::Replace(const string& word1, const string& word2) {
    args.push_back(word1);
    args.push_back(word2);
}

bool Replace::exec(vector<string> *input, vector<string> *output) {
    int index;
    output->clear();
    if (args[1].find(args[0])!= string::npos)
        Validator::badFile("In \"replace\" first string can't be substring of the second.");
    for (auto& i: *input){
        while ((index = i.find(args[0]))!=string::npos){
            i.replace(index,args[0].length(),args[1]);
        }
        output->push_back(i);
    }
    return true;
}