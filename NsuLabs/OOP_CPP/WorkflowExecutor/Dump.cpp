//
// Created by Admin on 06.12.2021.
//

#include "Dump.h"

Dump::Dump(string filename) {
    args.push_back(filename);
}

bool Dump::exec(vector<string> *input, vector<string> *output) {
    ofstream stream(args[0]);
    for (auto& i : *input){
        stream<<i<<endl;
        output->push_back(i);
    }
    stream.close();
    return true;
}
