package ru.nsu.ccfit.saburov.task2.materials;

import java.io.InputStreamReader;
import java.util.Map;
import java.util.Vector;

import static java.lang.Character.isDigit;
import static ru.nsu.ccfit.saburov.task2.materials.Factory.readWord;

public class Push implements Operation {
    public Push(){}

    static double getValue(Vector<String> args, Map<String, Double> definitions){
        int dotCounter = 0;
        if(definitions.containsKey(args.get(0)))
            return definitions.get(args.get(0));
        for (int i = 0; i<args.get(0).length(); i++){
            if (args.get(0).toCharArray()[i] == '.'){
                dotCounter++;
            }
            if(!isDigit(args.get(0).toCharArray()[i]) && !(i==0 && (args.get(0).toCharArray()[i] == '-' || args.get(0).toCharArray()[i] == '+')) && !(args.get(0).toCharArray()[i] == '.') || dotCounter>1){
                System.out.println("Failed to recognise param of the function. Expected defined variable or number. Got \""+args.get(0)+ "\" instead.");
                //System.exit(5);
                return Double.NaN;
            }
        }
        return Double.parseDouble(args.get(0));
    }
    @Override
    public void makeOperation(Context context, InputStreamReader reader) {
        context.args.add(readWord(reader));
        double a = getValue(context.args,context.definitions);
        if (!Double.isNaN(a))
            context.Data.push(a);
    }
}
