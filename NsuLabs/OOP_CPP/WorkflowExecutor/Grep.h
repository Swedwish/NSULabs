#ifndef LAB2_GREP_H
#define LAB2_GREP_H
#include "Worker.h"

class Grep: public Worker{
    private:
    vector <string> args;
public:
    Grep(string word);
    bool exec(vector<string> *input, vector<string> *output) override;
};


#endif //LAB2_GREP_H
