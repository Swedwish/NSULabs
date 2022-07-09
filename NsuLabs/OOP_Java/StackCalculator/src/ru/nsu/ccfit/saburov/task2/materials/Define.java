package ru.nsu.ccfit.saburov.task2.materials;

import java.io.InputStreamReader;
import java.util.Map;
import java.util.Vector;

import static java.lang.Character.isDigit;
import static ru.nsu.ccfit.saburov.task2.materials.Factory.readWord;

public class Define implements Operation{

    public Define() {
    }

    static double getValue(Vector<String> args, Map<String, Double> definitions){
        int dotCounter = 0;
        if(definitions.containsKey(args.get(1)))
            return definitions.get(args.get(1));
        for (int i = 0; i<args.get(1).length(); i++){
            if (args.get(1).toCharArray()[i] == '.'){
                dotCounter++;
            }
            if(!isDigit(args.get(1).toCharArray()[i]) && !(i==0 && (args.get(1).toCharArray()[i] == '-' || args.get(1).toCharArray()[i] == '+')) && !(args.get(1).toCharArray()[i] == '.') || dotCounter>1){
                System.out.println("Failed to recognise param of the function. Expected defined variable or number. Got \""+args.get(1)+ "\" instead.");
                //System.exit(5);
                return Double.NaN;
            }
        }
        return Double.parseDouble(args.get(1));
    }

    @Override
    public void makeOperation(Context context, InputStreamReader reader) {
        context.args.add(readWord(reader));
        context.args.add(readWord(reader));
        double a = getValue(context.args,context.definitions);
        if (!Double.isNaN(a))
            context.definitions.put(context.args.get(0),a);
    }


}
