package ru.nsu.ccfit.saburov.task2.materials;

import org.apache.log4j.Logger;

import java.io.*;

public class Main {
    private static final Logger logger = Logger.getLogger(Main.class);
    public static void main(String[] args) {
        logger.info("Program start.");
        try {
            InputStreamReader reader;
            if (args.length!=0) {
                reader = new InputStreamReader(new FileInputStream(args[0]));
            } else reader = new InputStreamReader(System.in);
            var Calculator = new Calculator(reader);
            Calculator.calculate();
            logger.info("Program end.");
        }
        catch (Exception e){
            System.out.println("Exception in main." + e);
        }
    }
}
