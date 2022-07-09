package ru.nsu.ccfit.saburov.view.gui;

import ru.nsu.ccfit.saburov.model.MineField;
import ru.nsu.ccfit.saburov.model.PairInt;
import ru.nsu.ccfit.saburov.view.View;

import javax.swing.*;
import javax.swing.border.Border;
import java.awt.*;

public class GuiView implements View {

    public int getWaitFlag() {
        return waitFlag;
    }

    public void setWaitFlag(int waitFlag) {
        this.waitFlag = waitFlag;
    }

    public int getX() {
        return x;
    }

    public void setX(int x) {
        this.x = x;
    }

    public int getY() {
        return y;
    }

    public void setY(int y) {
        this.y = y;
    }

    volatile int waitFlag;
    int x;
    int y;

    public int getIsLmb() {
        return isLmb;
    }

    public void setIsLmb(int isLmb) {
        this.isLmb = isLmb;
    }

    int isLmb;
    MineField mineField;
    MyFrame frame;
    MyButton[] buttons;

    public GuiView(MineField mineField){
        this.mineField = mineField;
        this.frame = new MyFrame(mineField, this);
        buttons = frame.getButtons();
    }

    @Override
    public void updateGameTextField() {
        int size = mineField.getSize();
        for (int i = 0; i < size; i++){
            for(int j = 0; j < size; j++){
                if (mineField.isOpened(new PairInt(i,j)) == 1){
                    if(mineField.getCell(new PairInt(i,j))=='0'){
                        buttons[i*size+j].setEnabled(false);
                        Border border = BorderFactory.createLineBorder(new Color(0x000000),1);
                        JLabel label = new JLabel();
                        label.setBorder(border);
                        label.setBounds(j*30,i*30,30,30);
                        frame.add(label);
                    }
                    else if (mineField.getCell(new PairInt(i,j))=='*') {
                        buttons[i*size+j].setEnabled(false);
                        buttons[i*size+j].setVisible(false);
                        Border border = BorderFactory.createLineBorder(new Color(0x000000),1);
                        ImageIcon image = new ImageIcon("src\\ru\\nsu\\ccfit\\saburov\\resoruses\\mine2.jpg");
                        JLabel label = new JLabel();
                        label.setBorder(border);
                        label.setIcon(image);
                        label.setBounds(j*30,i*30,30,30);
                        frame.add(label);
                    }
                    else {
                        Border border = BorderFactory.createLineBorder(new Color(0x000000),1);
                        buttons[i * size + j].setVisible(false);
                        JLabel label = new JLabel();
                        label.setText(Character.toString(mineField.getCell(new PairInt(i,j))));
                        label.setBounds(j*30,i*30,30,30);
                        label.setHorizontalTextPosition(JLabel.CENTER);
                        label.setFont(new Font("Comic Sans", Font.BOLD,20));
                        label.setBorder(border);
                        frame.add(label);
                    }
                }
                else if (mineField.isOpened(new PairInt(i,j)) == 2){
                    Border border = BorderFactory.createLineBorder(new Color(0x000000),1);
                    ImageIcon image = new ImageIcon("src\\ru\\nsu\\ccfit\\saburov\\resoruses\\flag2.jpg");
                    buttons[i * size + j].setIcon(image);
                }
                else if (mineField.isOpened(new PairInt(i,j)) == 0){
                    buttons[i * size + j].setIcon(null);
                    buttons[i * size + j].setVisible(true);
                }
            }
        }
    }

    @Override
    public void help() {

    }

    @Override
    public void wrong() {
        JFrame wrong = new JFrame("pls");
        wrong.setResizable(false);
        wrong.setSize(200,100);
        wrong.setVisible(true);
        wrong.setLayout(new FlowLayout());
        wrong.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);

        JLabel label = new JLabel("Dont do dis agan!");
        wrong.add(label);
    }


    @Override
    public void loss() {
        JFrame loss = new JFrame("Sadge");
        loss.setResizable(false);
        loss.setSize(200,100);
        loss.setVisible(true);
        loss.setLayout(new FlowLayout());
        loss.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        JLabel label = new JLabel("You lost!");
        loss.add(label);
        for (int i = 0; i < mineField.getSize();i++){
            for(int j = 0; j < mineField.getSize(); j++) {
                mineField.setOpened(i,j,1);
            }
        }
    }

    @Override
    public void victory() {
        JFrame victory = new JFrame("GJ!");
        victory.setResizable(false);
        victory.setSize(200,100);
        victory.setVisible(true);
        victory.setLayout(new FlowLayout());
        victory.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        JLabel label = new JLabel("You win!");
        victory.add(label);
    }

    @Override
    public void greetings() {

    }

    @Override
    public String[] makeMove() {
        String[] ans = new String[3];
        while (waitFlag == 0) Thread.onSpinWait();
        waitFlag = 0;
        if (isLmb == 1){
            ans[0] = "g";
            ans[1] = Integer.toString(x);
            ans[2] = Integer.toString(y);
        }
        else{
            ans[0] = "f";
            ans[1] = Integer.toString(x);
            ans[2] = Integer.toString(y);
        }
        return ans;
    }
}
