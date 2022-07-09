#include "ReadFile.h"

bool ReadFile::exec(vector<string> *input, vector<string> *output) {
    string line;
    ifstream stream(args[0]);
    while (getline(stream,line)){
        output->push_back(line);
    }
    stream.close();
    return true;
}

ReadFile::ReadFile(string input) {
    args.push_back(input);
}