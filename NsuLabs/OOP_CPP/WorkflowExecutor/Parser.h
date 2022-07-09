#ifndef LAB2_PARSER_H
#define LAB2_PARSER_H
#include "ReadFile.cpp"
#include "WriteFile.cpp"
#include "Sort.cpp"
#include "Grep.cpp"
#include "Replace.cpp"
#include "Dump.cpp"
#include "unordered_map"

class Parser {
public:
    ifstream& stream;
    Parser(ifstream&);
    unordered_map<int,Worker*> commands;
    Worker* make_block(const string& arg0);
    void Parse_desc();
    void Exec_queue();
};


#endif