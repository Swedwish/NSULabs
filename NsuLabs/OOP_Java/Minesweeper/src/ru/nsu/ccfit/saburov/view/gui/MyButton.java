package ru.nsu.ccfit.saburov.view.gui;

import javax.swing.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

public class MyButton extends JButton implements MouseListener {
    public int gettX() {
        return x;
    }


    public void settX(int x) {
        this.x = x;
    }

    public int gettY() {
        return y;
    }

    public void settY(int y) {
        this.y = y;
    }

    int x;
    int y;
    GuiView view;

    MyButton(GuiView view){
        this.view = view;
    }
    @Override
    public void mouseClicked(MouseEvent e) {
        if (SwingUtilities.isLeftMouseButton(e)){
            view.setIsLmb(1);
            view.setX(x + 1);
            view.setY(y + 1);
            view.setWaitFlag(1);
        }
        else if (SwingUtilities.isRightMouseButton(e)){
            //System.out.println(x+" "+y);
            view.setIsLmb(0);
            view.setX(x + 1);
            view.setY(y + 1);
            view.setWaitFlag(1);
        }
    }

    @Override
    public void mousePressed(MouseEvent e) {

    }

    @Override
    public void mouseReleased(MouseEvent e) {

    }

    @Override
    public void mouseEntered(MouseEvent e) {

    }

    @Override
    public void mouseExited(MouseEvent e) {

    }
}
