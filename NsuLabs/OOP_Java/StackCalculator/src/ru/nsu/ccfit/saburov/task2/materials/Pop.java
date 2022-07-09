package ru.nsu.ccfit.saburov.task2.materials;

import java.io.InputStreamReader;

public class Pop implements Operation {
    public Pop(){}

    @Override
    public void makeOperation(Context context, InputStreamReader reader) {
        if (!context.Data.empty()){
            context.Data.pop();
        }
        else{
            System.out.println("Stack is empty when attempting to pop.");
        }

    }
}
