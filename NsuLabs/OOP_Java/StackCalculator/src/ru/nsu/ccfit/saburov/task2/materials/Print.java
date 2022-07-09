package ru.nsu.ccfit.saburov.task2.materials;

import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;

public class Print implements Operation{

    public Print() {
    }

    @Override
    public void makeOperation(Context context, InputStreamReader reader) {
        if (!context.Data.empty()) {
            System.out.println(context.Data.peek());
            try {
                OutputStreamWriter s = new OutputStreamWriter(new FileOutputStream("src/test/ru/nsu/ccfit/saburov/task2/materials/output-log.txt", true));
                s.write(context.Data.peek()+"\n");
                s.close();
            }
            catch (Exception e){
                System.out.println(e);
            }
        }
        else{
            System.out.println("Stack is empty when attempting to print.");
        }
    }
}
