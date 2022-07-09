package com.company;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

public class Main {

    public static void main(String[] args) {
        Map<String,Integer> dict = new HashMap<>();
        Reader reader = null;
        try
        {
            reader = new InputStreamReader(new FileInputStream("C:\\Users\\katri\\IdeaProjects\\lab0\\src\\com\\company\\input.txt"));
            int letter;
            letter = reader.read();
            var builder = new StringBuilder();
            int wordCounter = 0;
            while (true){
                if (letter == ' ' || letter == -1){
                    if (dict.containsKey(builder.toString())){
                        dict.put(builder.toString(),dict.get(builder.toString())+1);
                    }
                    else{
                        dict.put(builder.toString(),1);
                    }
                    wordCounter+=1;
                    builder.delete(0,builder.length());
                }
                else {
                    builder.append((char) letter);
                }
                if (letter == -1){
                    break;
                }
                letter = reader.read();
            }
            var output = new FileWriter("C:\\Users\\katri\\IdeaProjects\\lab0\\src\\com\\company\\output.csv",false);
            for (var a: dict.keySet()) {
                output.write(a + '.' + dict.get(a) + '.' + (float)dict.get(a)/ (float)wordCounter*100 + "%\n");
            }
            output.close();
        }
        catch (IOException e)
        {
            System.err.println("Error while reading file: " + e.getLocalizedMessage());
        }
        finally
        {
            if (null != reader)
            {
                try
                {
                    reader.close();
                }
                catch (IOException e)
                {
                    e.printStackTrace(System.err);
                }
            }
        }

    }
}
