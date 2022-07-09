package ru.nsu.ccfit.saburov.view.gui;

import ru.nsu.ccfit.saburov.model.MineField;

import javax.swing.*;

public class MyFrame extends JFrame {
    MineField mineField;
    MyButton[] buttons;
    MyFrame(MineField mineField, GuiView view) {
        this.mineField = mineField;
        int size = mineField.getSize();
        buttons = new MyButton[size*size];
        this.setLayout(null);
        for(int i = 0; i < size; i++){
            for(int j = 0; j < size; j++){
                MyButton button = new MyButton(view);
                button.settX(j);
                button.settY(i);
                buttons[i*size+j] = button;
                button.addMouseListener(button);
                button.setBounds(30*j,30*i,30,30);
                this.add(button);
            }
        }

        this.setSize(size*30+15, size*30+39);
        this.setVisible(true);
        this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        this.setResizable(false);


        ImageIcon image = new ImageIcon("src\\ru\\nsu\\ccfit\\saburov\\resoruses\\mine.jpg");
        this.setIconImage(image.getImage());
        this.setTitle("Minesweeper");
    }

    MyButton[] getButtons(){
        return buttons;
    }


}
