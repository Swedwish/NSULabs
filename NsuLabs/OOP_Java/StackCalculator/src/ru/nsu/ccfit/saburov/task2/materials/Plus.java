package ru.nsu.ccfit.saburov.task2.materials;

import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.Objects;

public class Plus implements Operation{

    public Plus() {
    }

    @Override
    public void makeOperation(Context context, InputStreamReader reader) {

        double a,b;
        if (context.Data.empty()){
            System.out.println("Stack is empty when attempting to make a \"+\" operation");
        }
        else{
            a = context.Data.pop();
            if (context.Data.empty()) {
                System.out.println("Stack has only one value (" + a + ") when attempting to make a \"-\" operation");
                context.Data.push(a);
                return;
            }
            b = context.Data.pop();
            context.Data.push(a + b);
        }
    }
}