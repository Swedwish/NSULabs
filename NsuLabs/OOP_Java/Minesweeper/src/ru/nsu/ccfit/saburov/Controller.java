package ru.nsu.ccfit.saburov;

import ru.nsu.ccfit.saburov.model.MineField;
import ru.nsu.ccfit.saburov.model.PairInt;
import ru.nsu.ccfit.saburov.view.View;
import ru.nsu.ccfit.saburov.view.gui.GuiView;
import ru.nsu.ccfit.saburov.view.text.TextView;

import java.util.Objects;
import java.util.Scanner;




public class Controller {
    View view;
    MineField mineField;
    static void openNearby(MineField mineField, PairInt coordinates) {
        mineField.opened[coordinates.getFirst()][coordinates.getSecond()] = 1;
        for (int i = coordinates.getFirst() - 1; i <= coordinates.getFirst() + 1; i++) {
            for (int j = coordinates.getSecond() - 1; j <= coordinates.getSecond() + 1; j++) {
                try {
                    if(j == coordinates.getSecond() || i == coordinates.getFirst()) {
                        if (mineField.field[i][j] == '0' && mineField.opened[i][j] == 0) {
                            openNearby(mineField, new PairInt(i, j));
                        }
                        mineField.opened[i][j] = 1;
                    }
                } catch (ArrayIndexOutOfBoundsException ignored) {
                }
            }
        }
    }

    /*void menu(){
        view.printMenu();
        Scanner scanner = new Scanner(System.in);
        String option = scanner.nextLine().toLowerCase(Locale.ROOT);
        if (Objects.equals(option, "new game")){
            play();
        }
        else if(Objects.equals(option, "high scores")){
            view.printHighScores();
            scanner.nextLine();
            menu();
        }
        else if(Objects.equals(option, "about")){
            view.about();
            scanner.nextLine();
            menu();
        }
        else if (Objects.equals(option, "exit")){
            view.exit();
            System.exit(0);
        }
        else{
            view.tryAgain();
            menu();
        }
    }*/

    void play(int mode){
        Scanner in = new Scanner(System.in);
        System.out.println("Choose field size");
        int fieldSize = in.nextInt();
        if (fieldSize > 33){
            System.out.println("Nope, try again.");
            System.exit(0);
        }
        System.out.println("Choose mine count");
        int mineCount = in.nextInt();

        mineField = new MineField(fieldSize,mineCount);
        mineField.fillField();
        if (mode == 'g') {
            view = new GuiView(mineField);
        }
        else if (mode == 't'){
            view = new TextView(mineField);
        }
        String[] task;
        view.greetings();
        view.updateGameTextField();
        while (true) {
            task = view.makeMove();
            if (Objects.equals(task[0], "g")){
                go(task[1],task[2]);
            }
            else if (Objects.equals(task[0], "f")){
                flag(task[1],task[2]);
            }
        }
    }



    void flag(String secondCoordinate, String firstCoordinate){
        if (mineField.opened[Integer.parseInt(firstCoordinate)-1][Integer.parseInt(secondCoordinate)-1] == 0) {
            mineField.opened[Integer.parseInt(firstCoordinate)-1][Integer.parseInt(secondCoordinate)-1] = 2;
        }
        else if (mineField.opened[Integer.parseInt(firstCoordinate)-1][Integer.parseInt(secondCoordinate)-1] == 1){
            view.wrong();
        }
        else{
            mineField.opened[Integer.parseInt(firstCoordinate)-1][Integer.parseInt(secondCoordinate)-1] = 0;
        }
        view.updateGameTextField();
    }

    void go(String secondCoordinate, String firstCoordinate) {
        PairInt coordinate = new PairInt(Integer.parseInt(firstCoordinate)-1, Integer.parseInt(secondCoordinate)-1);
        if (coordinate.getSecond()<0 || coordinate.getSecond()>= mineField.getSize() || coordinate.getFirst()<0 || coordinate.getFirst()>= mineField.getSize()){
            view.wrong();
            return;
        }
        mineField.opened[coordinate.getFirst()][coordinate.getSecond()] = 1;
        if (mineField.field[coordinate.getFirst()][coordinate.getSecond()] == '*') {
            view.loss();
        } else if (mineField.field[coordinate.getFirst()][coordinate.getSecond()] == '0') {
            openNearby(mineField, coordinate);
        }
        view.updateGameTextField();
        int flag = 0;
        for (int i = 0; i < mineField.getSize(); i++) {
            for (int j = 0; j < mineField.getSize(); j++) {
                if (mineField.opened[i][j] != 1){
                    flag++;
                }
            }
        }
        if (flag == mineField.mineCount){
            view.victory();
        }
    }
}
