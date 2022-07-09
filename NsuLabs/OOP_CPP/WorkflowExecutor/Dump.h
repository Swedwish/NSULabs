#ifndef LAB2_DUMP_H
#define LAB2_DUMP_H
#include "Worker.h"

class Dump: public Worker {
    private:
    vector <string> args;
public:
    Dump(string filename);
    bool exec(vector<string> *input, vector<string> *output) override;
};


#endif //LAB2_DUMP_H
