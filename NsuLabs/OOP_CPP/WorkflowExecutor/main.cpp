#include "Parser.cpp"
#include <iostream>

int main() {
    ifstream stream("workflow.txt");
    if (!Validator::file_Exists("workflow.txt")){
        throw invalid_argument("workflow.txt does not exist");
    }
    Parser parser(stream);
    parser.Parse_desc();
    parser.Exec_queue();

    return 0;
}
