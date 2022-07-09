#ifndef LAB2_VALIDATOR_H
#define LAB2_VALIDATOR_H
#pragma once
#include <string>
#include <sys/stat.h>
#include <stdexcept>
using namespace std;

class Validator
{
public:
    bool file_Exists(const string& file);
    void badFile(const string& message);
};

#endif //LAB2_VALIDATOR_H
