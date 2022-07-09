package ru.nsu.ccfit.saburov.task2.materials;

import java.io.InputStreamReader;
import java.util.*;

public class Factory {

    static public String readWord(InputStreamReader i) {
        try {
            return new Scanner(i).next();
        }
        catch(Exception e){
            System.out.println("IOException." + e);
            //System.exit(5);
            return "";
        }
    }

    public Operation makeClass(String word) {
        try {
            Properties pps = new Properties();
            pps.load(Factory.class.getResourceAsStream("prop.properties"));
            if(pps.containsKey(word)) {
                return (Operation) Class.forName(pps.getProperty(word)).getDeclaredConstructor().newInstance();
            }
            else if(pps.containsKey(String.valueOf(word.charAt(0)))){
                return (Operation) Class.forName(pps.getProperty(String.valueOf(word.charAt(0)))).getDeclaredConstructor().newInstance();
            }
            else{
                System.out.println("Can't recognise word \"" + word + "\".\n");
                return null;
            }
        } catch (Exception e) {
            System.out.println("Exception in Factory." + e);
        }
        return null;
    }
}
