package ru.nsu.ccfit.saburov.task2.materials;

import java.util.Map;
import java.util.Stack;
import java.util.Vector;

public class Context {
    Stack<Double> Data;
    Map<String, Double> definitions;
    Vector<String> args;

    public Context(Stack<Double> Data,
                   Map<String, Double> definitions,
                   Vector<String> args) {
        this.Data = Data;
        this.definitions = definitions;
        this.args = args;
    }
    public void argsAdd(String s){
        args.add(s);
    }
    public void argsClear(){
        args.clear();
    }
    public boolean dataEmpty(){
        return Data.empty();
    }
    public Double dataPeek()
    {
        return Data.peek();
    }
    public void argsRemove(int index){
        args.remove(index);
    }

    public void dataClear() {
        Data.clear();
    }

    public void argsSet(int i, String s) {
        args.set(i,s);
    }
}
