package ru.nsu.ccfit.saburov.task2.materials;

import org.apache.log4j.Logger;

import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Stack;
import java.util.Vector;
import static ru.nsu.ccfit.saburov.task2.materials.Factory.readWord;

public class Calculator {
    private final Context context = new Context(new Stack<>(),new HashMap<>(),new Vector<>());
    private static final Logger logger = Logger.getLogger(Calculator.class);
    private final InputStreamReader reader;

    public Calculator(InputStreamReader reader) {
        this.reader = reader;
    }
    public void calculate() {
        try {
            String word;
            Operation worker;
            Factory factory= new Factory();

            while (!(word = readWord(reader)).equals("")) {
                context.args.clear();
                worker = factory.makeClass(word);
                if (worker == null){
                    continue;
                }
                logger.info("Executing "+ word+"...");
                worker.makeOperation(context, reader);
            }
        } catch (Exception e) {
            System.out.println("Exception in Factory." + e);
        }
    }
}
