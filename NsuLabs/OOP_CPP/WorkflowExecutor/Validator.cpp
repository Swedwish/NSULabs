#include "Validator.h"

bool Validator::file_Exists(const string& file)
{
    struct stat buffer{};
    return (stat (file.c_str(), &buffer) == 0);
}
void Validator::badFile(const string& message)
{
    throw invalid_argument("Bad program file! " + message);
}