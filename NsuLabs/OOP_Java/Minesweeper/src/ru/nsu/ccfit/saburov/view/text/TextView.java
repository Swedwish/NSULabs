package ru.nsu.ccfit.saburov.view.text;

import ru.nsu.ccfit.saburov.model.MineField;
import ru.nsu.ccfit.saburov.model.PairInt;
import ru.nsu.ccfit.saburov.view.View;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.Scanner;

public class TextView implements View {
    static void printDashLine(int size){
        for(int i = 0; i <= size*4; i++){
            System.out.print("\u2014");
        }
        System.out.println();
    }

    MineField mineField;

     public TextView(MineField mineField){
        this.mineField = mineField;
    }

    @Override
    public String[] makeMove() {
        String[] ans = new String[3];
        System.out.print("Write command:");
        Scanner scanner = new Scanner(System.in);
        String command = scanner.next();
        if (command.equals("help")) {
            help();
            String a = new Scanner(System.in).nextLine();
        } else if (scanner.hasNextInt()) {
            if (command.equals("g")) {
                String firstCoordinate = scanner.next();
                if (scanner.hasNextInt()) {
                    String secondCoordinate = scanner.next();
                    ans[0] = "g";
                    ans[1] = firstCoordinate;
                    ans[2] = secondCoordinate;
                    return ans;
                } else wrong();
            } else if (command.equals("f")) {
                String firstCoordinate = scanner.next();
                if (scanner.hasNextInt()) {
                    String secondCoordinate = scanner.next();
                    ans[0] = "f";
                    ans[1] = firstCoordinate;
                    ans[2] = secondCoordinate;
                    return ans;
                } else wrong();
            } else if (Character.isDigit(command.charAt(0))) {
                if (scanner.hasNextInt()) {
                    String secondCoordinate = scanner.next();
                    ans[0] = "g";
                    ans[1] = command;
                    ans[2] = secondCoordinate;
                    return ans;
                } else wrong();
            } else {
                wrong();
            }
        }
        ans[0] = "0";
        return ans;
    }

     public void showTextField(MineField mineField){
        int size = mineField.getSize();
        printDashLine(size);
        for (int i = 0; i < size; i++){
            System.out.print("| ");
            for(int j = 0; j< size; j++){
                System.out.print(mineField.getCell(new PairInt(i,j)));
                System.out.print(" | ");
            }
            System.out.println();
            printDashLine(size);
        }
    }

    public void updateGameTextField(){
        int size = mineField.getSize();
        printDashLine(size);
        for (int i = 0; i < size; i++){
            System.out.print("| ");
            for(int j = 0; j< size; j++){
                if(mineField.isOpened(new PairInt(i,j))==1){
                    System.out.print(mineField.getCell(new PairInt(i,j)));
                }
                else if(mineField.isOpened(new PairInt(i,j))==2){
                    System.out.print("\u1F32");
                }
                else{
                    System.out.print("#");
                }
                System.out.print(" | ");
            }
            System.out.println();
            printDashLine(size);
        }
    }

    public void victory(){
        System.out.println("You win!!!!!!!");
        showTextField(mineField);
        System.exit(0);
    }

    public void loss(){
        System.out.println("Game over (9(((((");
        showTextField(mineField);
        System.exit(0);
    }

    public void greetings(){
        System.out.println("""
                Welcome to the Minesweeper!
                Game starts now!
                (Type "help" if you feel lost)""");
    }

    public void printMenu() {
        System.out.println("""
                ----!!!MENU!!!----
                New Game
                High Scores
                About
                Exit""");
    }

    public void about() {
        System.out.println("""
                Made in 24 hours by a 2nd grade student by name V.S. (Visual studio).
                """);
        System.out.println("Type anything to go back to menu.");
    }

    public void printHighScores() {
        try {
            Scanner scanner = new Scanner(new FileInputStream("src/ru/nsu/ccfit/saburov/resoruses/HighScores.txt"));
            for (int i = 0; i < 11; i++) {
                System.out.println(scanner.nextLine());
            }
            System.out.println("Type anything to go back to menu.");
        }
        catch (FileNotFoundException e){
            System.out.println("Whoops. High Scores file not found. Not a big deal tho, u can still play.");
        }
    }

    public void exit() {
        System.out.println("Thanks for playing (or not playing if you weren't really playing)!");
    }

    public void tryAgain() {
        System.out.println("Incorrect input, try again.");
    }

    public void help(){
        System.out.println("""
                If you don't know the rules, go google it, idk.
                g + coords (e.g. g 3 4) or just coords (e.g. 3 4) to open the field.
                f + coords (e.g. f 3 4) to mark the place where you think the mine could be.
                That's it.""");
        System.out.println("Type anything to go back to game.");
    }

    public void wrong() {
        System.out.println("Wrong, try again!");
    }
}