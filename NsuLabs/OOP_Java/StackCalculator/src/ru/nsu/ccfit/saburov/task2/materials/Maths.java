package ru.nsu.ccfit.saburov.task2.materials;

import java.io.InputStreamReader;
import java.util.Objects;

public class Maths implements Operation {

    public Maths() {
    }

    @Override
    public void makeOperation(Context context, InputStreamReader reader) {

        double a,b;
        if (context.Data.empty()){
            System.out.println("Stack is empty when attempting to make an \""+context.args.get(0) +"\" operation");
        }
        else{
            a = context.Data.pop();
            if (context.Data.empty()) {
                System.out.println("Stack has only one value (" + a + ") when attempting to make an \"" + context.args.get(0) + "\" operation");
                context.Data.push(a);
            }
            b = context.Data.pop();
            if (Objects.equals(context.args.get(0), "+")) {
                context.Data.push(a + b);
            } else if (Objects.equals(context.args.get(0), "-")) {
                context.Data.push(b - a);
            } else if (Objects.equals(context.args.get(0), "*")) {
                context.Data.push(a * b);
            } else if (Objects.equals(context.args.get(0), "/")) {
                context.Data.push(b / a);
            }
        }
    }
}
