#ifndef LAB2_SORT_H
#define LAB2_SORT_H
#include "Worker.h"
#include <algorithm>

class Sort: public Worker{
private:
    vector <string> args;
public:
    bool exec(vector<string> *input, vector<string> *output) override;
};


#endif //LAB2_SORT_H
