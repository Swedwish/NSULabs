#include "Parser.h"


Worker *Parser::make_block(const string &arg0) {
    vector<string> args;
    string a;
    args.push_back(arg0);
    if (args[0] == "readfile"){
        stream>>a;
        args.push_back(a);
        auto* temp = new ReadFile(args[1]);
        Worker* prod = temp;
        return prod;
    }
    else if (args[0] == "writefile"){
        stream>>a;
        args.push_back(a);
        auto *temp = new WriteFile(args[1]);
        Worker* prod = temp;
        return prod;
    }
    else if (args[0] == "grep"){
        stream>>a;
        args.push_back(a);
        auto* temp = new Grep(args[1]);
        Worker* prod = temp;
        return prod;
    }
    else if (args[0] == "sort"){
        auto* temp = new Sort;
        Worker* prod = temp;
        return prod;
    }
    else if (args[0] == "replace"){
        stream>>a;
        args.push_back(a);
        stream>>a;
        args.push_back(a);
        auto *temp = new Replace(args[1],args[2]);
        Worker* prod = temp;
        return prod;
    }
    else if (args[0] == "dump"){
        stream>>a;
        args.push_back(a);
        auto *temp = new Dump(args[1]);
        Worker* prod = temp;
        return prod;
    }
    else{
        Validator::badFile("Wrong function name. Doesn't match to existing functions.");
    }
}

void Parser::Parse_desc() {
    string str;
    stream>>str;
    if (str != "desc"){
        Validator::badFile("No \"desc\", wrong description format");
    }
    int index;
    while (true){
        stream >> str;
        if (str == "csed") {
            break;
        }
        index = stoi(str);
        stream >> str;
        if (str != "=") {
            Validator::badFile("Wrong description format.");
        }
        stream >> str;
        commands[index] = make_block(str);
    }
}

Parser::Parser(ifstream& stream1):stream(stream1) {

}

void Parser::Exec_queue() {
    string str;
    auto *input = new vector<string>, *output = new vector<string>;
    while (stream>>str){
        if (str == "->")
            continue;
        commands[stoi(str)]->exec(input, output);
        input->clear();
        input = output;
        output = new vector<string>;
    }
    output->clear();
}